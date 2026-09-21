#!/usr/bin/env python3
"""
ED_Fig_age_model.py
====================
Extended Data Figure: HS4 age-depth chronology.

Two-panel figure:
  (A) Full age-depth model with 31 tie points (30 U-Th ages and the growth surface; 2σ error bars),
      interpolated 1-yr resolution curve, growth rate inset, and
      8.2 ka zone highlight.
  (B) 8.2 ka zone detail showing age reversals in raw U-Th dates
      and the monotonic interpolation.

Inputs:
  - HS4_age_depth.csv      (31 tie points: 30 U-Th ages and the growth surface; depth, age, 2σ error)
  - age_model.csv           (interpolated 1-yr resolution: depth, age_yBP)

Outputs:
  ED_Fig_age_model.png (600 dpi)
  ED_Fig_age_model.pdf (vector)

Usage:
  python ED_Fig_age_model.py

  To override input paths:
    python ED_Fig_age_model.py --tiepoints path/to/HS4_age_depth.csv \
                                --age_model path/to/age_model.csv
"""

import argparse
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', 'manuscript_figures'))
from ngeo_style import (apply_style, style_ax, panel_label, DOUBLE_COL,
                        COL_BG_900, COL_BG_600, COL_BG_400, COL_BG_300, COL_BG_200,
                        COL_TEAL_800, COL_TEAL_600, COL_NI, COL_CO,
                        COL_BROWN_500, COL_DORANGE, COL_RED_900, COL_GREEN_800)
apply_style()

# ── Arguments ────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description='ED age model figure')
parser.add_argument('--tiepoints', default='HS4_age_depth.csv',
                    help='Path to U-Th tie-point CSV')
parser.add_argument('--age_model', default='age_model.csv',
                    help='Path to interpolated age model CSV')
parser.add_argument('--output', default='ExD_6',
                    help='Output filename stem (no extension)')
args = parser.parse_args()

# Nature Communications house style comes from ngeo_style.apply_style() (above).

# ── Load data ────────────────────────────────────────────────────────
tp = pd.read_csv(args.tiepoints)
tp.columns = [c.strip() for c in tp.columns]
tp_depth = tp.iloc[:, 0].values
tp_age = tp.iloc[:, 1].values
tp_err = tp.iloc[:, 2].values

am = pd.read_csv(args.age_model)
am_depth = am['depth'].values
am_age = am['age_yBP'].values
# the 1-yr age grid runs on past the base of the stalagmite with depth held at its maximum;
# keep the model only up to the first time it reaches the base
_n_end = int(np.argmax(am_depth >= am_depth.max())) + 1
am_depth, am_age = am_depth[:_n_end], am_age[:_n_end]

# ── Identify 8.2 ka zone (dense tie points) ─────────────────────────
zone_82_mask = (tp_depth >= 228) & (tp_depth <= 240)
normal_mask = ~zone_82_mask

# ── Figure: 2 panels ────────────────────────────────────────────────
fig, (ax1, ax2) = plt.subplots(
    1, 2, figsize=(DOUBLE_COL, DOUBLE_COL * 0.4),
    gridspec_kw={'width_ratios': [3, 1.2], 'wspace': 0.35})

# ── Panel A: Full age-depth model ────────────────────────────────────
ax1.plot(am_depth, am_age / 1000, color=COL_NI, lw=0.7, zorder=2,
         label='Interpolated model (1-yr resolution)')

# Normal tie points
ax1.errorbar(tp_depth[normal_mask], tp_age[normal_mask] / 1000,
             yerr=tp_err[normal_mask] / 1000,
             fmt='o', markersize=3.5, color=COL_RED_900, ecolor=COL_RED_900, elinewidth=0.4, capsize=1.2, capthick=0.3, markeredgecolor=COL_BG_900, markeredgewidth=0.25,
             zorder=4, label=f'U\u2013Th ages and growth surface (n = {normal_mask.sum()})')

# 8.2 ka zone tie points — different marker
ax1.errorbar(tp_depth[zone_82_mask], tp_age[zone_82_mask] / 1000,
             yerr=tp_err[zone_82_mask] / 1000,
             fmt='s', markersize=3.5, color=COL_DORANGE, ecolor=COL_DORANGE, elinewidth=0.4, capsize=1.2, capthick=0.3, markeredgecolor=COL_BG_900, markeredgewidth=0.25,
             zorder=4, label=f'8.2 ka zone U\u2013Th ages (n = {zone_82_mask.sum()})')

# 8.2 ka zone shading
ax1.axvspan(228, 240, color=COL_DORANGE, alpha=0.10, lw=0, zorder=0)
ax1.text(226, 10.2, '8.2 ka zone', fontsize=5.5, color=COL_DORANGE,
         ha='right', va='top', style='italic')

ax1.set_xlabel('Depth (cm)')
ax1.set_ylabel('Age (ka BP)')
ax1.set_xlim(-5, 260)
ax1.set_ylim(-0.5, 10.5)
ax1.legend(loc='upper left', bbox_to_anchor=(0.07, 1.0), frameon=False, fontsize=5.5)
style_ax(ax1)
panel_label(ax1, 'a')

# Growth rate inset — from the monotonic interpolated model (age_model.csv)
# Enforce depth monotonicity (remove tiny numerical reversals in 8.2 ka zone)
am_depth_mono = np.maximum.accumulate(am_depth)
dd = np.diff(am_depth_mono)
da = np.diff(am_age)
gr_raw = dd / da * 1000  # cm kyr⁻¹
gr_depth_mid = (am_depth_mono[:-1] + am_depth_mono[1:]) / 2
# Smooth with 50-yr rolling window
from scipy.ndimage import uniform_filter1d
gr_smooth = uniform_filter1d(gr_raw, size=50)

ax1_inset = ax1.inset_axes([0.60, 0.13, 0.24, 0.21])   # clear of the data and short of the 8.2 ka band
ax1_inset.fill_between(gr_depth_mid, gr_smooth, color=COL_NI, alpha=0.2, lw=0)
ax1_inset.plot(gr_depth_mid, gr_smooth, color=COL_NI, lw=0.5)
ax1_inset.set_xlabel('Depth (cm)', fontsize=5)
ax1_inset.set_ylabel('Growth rate\n(cm kyr$^{-1}$)', fontsize=5)
ax1_inset.tick_params(labelsize=5)
ax1_inset.set_xlim(0, 253)
ax1_inset.set_ylim(0, None)
ax1_inset.axvspan(228, 240, color=COL_DORANGE, alpha=0.15, lw=0)

# ── Panel B: 8.2 ka zone zoom ───────────────────────────────────────
z_mask_am = (am_depth >= 225) & (am_depth <= 243)
ax2.plot(am_depth[z_mask_am], am_age[z_mask_am], color=COL_NI, lw=0.9,
         zorder=2)

ax2.errorbar(tp_depth[zone_82_mask], tp_age[zone_82_mask],
             yerr=tp_err[zone_82_mask],
             fmt='s', markersize=3.5, color=COL_DORANGE, ecolor=COL_DORANGE, elinewidth=0.4, capsize=1.5, capthick=0.4, markeredgecolor=COL_BG_900, markeredgewidth=0.25,
             zorder=4)

# Annotate each tie point — labels right-aligned at x=250, dashed connectors
LABEL_X = 250  # fixed x for all labels
z_idx = np.where(zone_82_mask)[0]
# Sort by age for clean vertical stacking
z_sorted = sorted(zip(tp_depth[z_idx], tp_age[z_idx], tp_err[z_idx]),
                  key=lambda t: t[1])

# Keep each label as close to its own age as the text height allows (1-D repel), so the
# connectors stay short and do not cross
MIN_GAP = 62.0   # yr: one 5 pt line at this panel's scale
label_ages = np.array([a for _, a, _ in z_sorted], float)
for _ in range(500):
    moved = False
    for i in range(len(label_ages) - 1):
        gap = label_ages[i + 1] - label_ages[i]
        if gap < MIN_GAP - 1e-6:
            push = (MIN_GAP - gap) / 2
            label_ages[i] -= push; label_ages[i + 1] += push; moved = True
    if not moved:
        break

# Assign labels to those slots so that no two connectors cross: swap neighbours whose
# connectors intersect until none do
def _cross(p1, p2, q1, q2):
    def o(a, b, c): return (b[0]-a[0])*(c[1]-a[1]) - (b[1]-a[1])*(c[0]-a[0])
    return (o(p1, p2, q1) * o(p1, p2, q2) < 0) and (o(q1, q2, p1) * o(q1, q2, p2) < 0)
slots = list(label_ages); order = list(range(len(z_sorted)))
for _ in range(200):
    swapped = False
    for i in range(len(order) - 1):
        (d1, a1, _), (d2, a2, _) = z_sorted[order[i]], z_sorted[order[i + 1]]
        if _cross((d1, a1), (LABEL_X, slots[i]), (d2, a2), (LABEL_X, slots[i + 1])):
            order[i], order[i + 1] = order[i + 1], order[i]; swapped = True
    if not swapped:
        break
z_sorted = [z_sorted[k] for k in order]

for (d, a, e), y_label in zip(z_sorted, label_ages):
    # Dashed connector line from point to label
    ax2.plot([d, LABEL_X], [a, y_label], color=COL_BG_300, lw=0.3,
             clip_on=False, zorder=1)
    # Label text
    ax2.text(LABEL_X + 0.5, y_label, f'{a:.0f} \u00b1 {e:.0f}',
             fontsize=5, color=COL_BG_600, va='center', ha='left',
             clip_on=False)

# Mark age reversals with arrows
z_depths = tp_depth[zone_82_mask]
z_ages = tp_age[zone_82_mask]
for i in range(len(z_depths) - 1):
    if z_ages[i + 1] < z_ages[i]:  # reversal
        ax2.annotate('', xy=(z_depths[i + 1], z_ages[i + 1]),
                     xytext=(z_depths[i], z_ages[i]),
                     arrowprops=dict(arrowstyle='->', color=COL_RED_900, lw=0.5, ls='--'))

ax2.set_xlabel('Depth (cm)')
ax2.set_ylabel('Age (yr BP)')
ax2.set_xlim(226, 248)
ax2.set_ylim(min(label_ages.min(), (tp_age[zone_82_mask] - tp_err[zone_82_mask]).min()) - 60,
             max(label_ages.max(), (tp_age[zone_82_mask] + tp_err[zone_82_mask]).max()) + 60)
panel_label(ax2, 'b')

style_ax(ax2); style_ax(ax1_inset)

# ── Save ─────────────────────────────────────────────────────────────
fig.savefig(f'{args.output}.png', dpi=600, bbox_inches='tight',
            facecolor='white')
fig.savefig(f'{args.output}.pdf', bbox_inches='tight', facecolor='white')
fig.savefig(f'{args.output}.eps', bbox_inches='tight', facecolor='white')
plt.close()
print(f'Saved {args.output}.png (600 dpi) and .pdf')
