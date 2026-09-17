#!/usr/bin/env python3
"""
SuppFig_SF_Holocene.py

Generates the Smart & Friedrich Holocene evolution figure for ED Fig 1.
Single-panel figure:
  250-yr binned scatter in S&F classification space (mean vs CV),
  with IQR error bars, whole-record star, 5.2 ka triangle,
  monitoring diamonds (raw and integration-time corrected).

Inputs:
  - drip_rate_summary.csv  (from a fixed-concentration model run)
  - age_model.csv  (depth -> age mapping; required if summary has depth only)
  - Copy_of_Geochemistry_of_HS4_dripwater__and_pool_water_-For_Adam.xlsx
    (monitoring drip rate data from Climate sheet, column D)
  - chart_data.json (for whole-record S&F summary stats)

Usage:
  python SuppFig_SF_Holocene.py \
      --summary drip_rate_summary.csv \
      --age_model age_model.csv \
      --chart_data chart_data.json \
      --monitoring Copy_of_Geochemistry_of_HS4_dripwater__and_pool_water_-For_Adam.xlsx \
      --output SuppFig_SF_Holocene

Outputs:
  SuppFig_SF_Holocene.png (600 dpi)
  SuppFig_SF_Holocene.pdf (vector)
  SuppFig_SF_Holocene.eps (vector)
"""

import argparse
import os
import numpy as np
import pandas as pd
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', 'manuscript_figures'))
from ngeo_style import (apply_style, style_ax, SINGLE_COL, DOUBLE_COL, COL_BG_900, COL_BG_600, COL_BG_400, COL_BG_300, COL_BG_200, COL_BG_50, COL_TEAL_800, COL_TEAL_600, COL_TEAL_200, COL_BROWN_500, COL_BROWN_200, COL_DORANGE, COL_GREEN_800, COL_RED_900)
apply_style()
from matplotlib.patches import Rectangle
from matplotlib.colors import Normalize
import matplotlib.cm as cm

# ── Parse arguments ──────────────────────────────────────────────────────
print('SuppFig_SF_Holocene v6 — single panel, no time series')
parser = argparse.ArgumentParser(description='S&F Holocene evolution figure')
parser.add_argument('--summary', default='drip_rate_summary.csv',
                    help='Path to drip_rate_summary.csv')
parser.add_argument('--age_model', default=None,
                    help='Path to age_model.csv (depth -> age mapping). '
                         'Required if summary has depth but no age column.')
parser.add_argument('--chart_data', default='chart_data.json',
                    help='Path to chart_data.json')
parser.add_argument('--monitoring', default=None,
                    help='Path to dripwater geochemistry xlsx (optional)')
parser.add_argument('--output', default='SuppFig_SF_Holocene',
                    help='Output filename stem (no extension)')
parser.add_argument('--bin_width', type=int, default=250,
                    help='Bin width in years (default: 250)')
parser.add_argument('--mon_mean', type=float, default=None,
                    help='Override monitoring mean drip rate (drips/min)')
parser.add_argument('--mon_cv_raw', type=float, default=None,
                    help='Override monitoring raw CV')
parser.add_argument('--T_sample', type=float, default=25,
                    help='Integration time per sample in years (default: 25)')
parser.add_argument('--f_s', type=float, default=10.8,
                    help='Monitoring sampling frequency yr-1 (default: 10.8)')
parser.add_argument('--event_half', type=int, default=100,
                    help='Half-width of 5.2 ka event window in years (default: 100)')
args = parser.parse_args()


# ── Helper ───────────────────────────────────────────────────────────────
def find_col(df, candidates, label):
    """Find the first matching column (case-insensitive, stripped)."""
    cols_lower = {c.strip().lower(): c for c in df.columns}
    for cand in candidates:
        if cand.lower() in cols_lower:
            return cols_lower[cand.lower()]
    return None


# ── Load data ────────────────────────────────────────────────────────────
df = pd.read_csv(args.summary)

col_age = find_col(df, ['age', 'age_bp', 'age (bp)', 'age_yr_bp'], 'age')
col_depth = find_col(df, ['depth', 'depth_mm', 'depth (mm)'], 'depth')
col_med = find_col(df, ['pc50', 'p50', 'median', 'drip_rate_median', 'q50'], 'median')
col_p25 = find_col(df, ['pc25', 'p25', 'q25', 'drip_rate_q25'], '25th pctl')
col_p75 = find_col(df, ['pc75', 'p75', 'q75', 'drip_rate_q75'], '75th pctl')

if col_med is None or col_p25 is None or col_p75 is None:
    raise KeyError(f"Cannot find percentile columns. Available: {list(df.columns)}")

# ── Resolve age axis ─────────────────────────────────────────────────────
if col_age is not None:
    ages = df[col_age].values
    print(f'Using age column "{col_age}" directly from summary CSV')
elif col_depth is not None:
    # Depth-only summary -> need age model to convert
    am_path = args.age_model
    if am_path is None:
        # Auto-detect: look for age_model.csv alongside the summary
        summary_dir = os.path.dirname(os.path.abspath(args.summary))
        for fname in ['age_model.csv', 'age_model.json']:
            candidate = os.path.join(summary_dir, fname)
            if os.path.isfile(candidate):
                am_path = candidate
                break
    if am_path is None:
        raise FileNotFoundError(
            f"Summary has depth but no age column, and no age model found.\n"
            f"Available columns: {list(df.columns)}\n"
            f"Provide --age_model path/to/age_model.csv (columns: depth, age_yBP)")

    print(f'Loading age model from {am_path}')
    if am_path.endswith('.json'):
        with open(am_path) as f:
            am = json.load(f)
        am_depth = np.array(am['depth'])
        am_age = np.array(am['age_median'])
    else:
        am_df = pd.read_csv(am_path)
        am_depth_col = find_col(am_df, ['depth', 'depth_mm'], 'age model depth')
        am_age_col = find_col(am_df, ['age_ybp', 'age', 'age_bp', 'age_median',
                                       'cal_bp', 'calbp'], 'age model age')
        if am_depth_col is None or am_age_col is None:
            raise KeyError(f"Cannot parse age model. Columns: {list(am_df.columns)}")
        am_depth = am_df[am_depth_col].values
        am_age = am_df[am_age_col].values

    depths = df[col_depth].values
    ages = np.interp(depths, am_depth, am_age)
    print(f'Interpolated {len(ages)} depth values -> age range '
          f'{ages.min():.0f}-{ages.max():.0f} yr BP')
else:
    raise KeyError(f"Cannot find age or depth column. Available: {list(df.columns)}")

med = df[col_med].values
p25 = df[col_p25].values
p75 = df[col_p75].values

with open(args.chart_data) as f:
    cd = json.load(f)
sf = cd['sf']

# ── 250-year bins ────────────────────────────────────────────────────────
bin_width = args.bin_width
n_boot = 1000
rng = np.random.default_rng(42)
bins = []
for lo in np.arange(ages.min(), ages.max(), bin_width):
    hi = lo + bin_width
    mask = (ages >= lo) & (ages < hi)
    v50 = med[mask]
    v25 = p25[mask]
    v75 = p75[mask]
    fin = np.isfinite(v50) & np.isfinite(v25) & np.isfinite(v75)
    v50, v25, v75 = v50[fin], v25[fin], v75[fin]
    if len(v50) > 5:
        mu = np.mean(v50)
        cv = np.std(v50) / mu if mu > 0 else 0
        q25_mean = np.percentile(v50, 25)
        q75_mean = np.percentile(v50, 75)
        # Bootstrap CV: resample within bin, compute CV each time
        boot_cvs = np.empty(n_boot)
        for b in range(n_boot):
            sample = rng.choice(v50, size=len(v50), replace=True)
            bmu = np.mean(sample)
            boot_cvs[b] = np.std(sample) / bmu if bmu > 0 else 0
        cv_lo = np.percentile(boot_cvs, 25)
        cv_hi = np.percentile(boot_cvs, 75)
        bins.append({
            'age_mid': (lo + hi) / 2,
            'mean': mu,
            'cv': cv,
            'mean_lo': q25_mean,
            'mean_hi': q75_mean,
            'cv_lo': cv_lo,
            'cv_hi': cv_hi,
        })
bdf = pd.DataFrame(bins)
bdf['yerr_lo'] = np.maximum(0, bdf['mean'] - bdf['mean_lo'])
bdf['yerr_hi'] = np.maximum(0, bdf['mean_hi'] - bdf['mean'])
bdf['xerr_lo'] = np.maximum(0, bdf['cv'] - bdf['cv_lo'])
bdf['xerr_hi'] = np.maximum(0, bdf['cv_hi'] - bdf['cv'])

print(f'Bins: {len(bdf)}, age range {bdf["age_mid"].min():.0f}-{bdf["age_mid"].max():.0f}')

# ── Monitoring data ──────────────────────────────────────────────────────
if args.mon_mean is not None and args.mon_cv_raw is not None:
    mon_mean = args.mon_mean
    mon_cv_raw = args.mon_cv_raw
    mon_sd = mon_cv_raw * mon_mean
elif args.monitoring:
    import openpyxl
    wb = openpyxl.load_workbook(args.monitoring, data_only=True)
    ws = wb['Climate']
    dr_vals = []
    for row in range(5, ws.max_row + 1):
        v = ws.cell(row=row, column=4).value  # DR column D (HS4)
        if v is not None:
            try:
                fv = float(v)
                if np.isfinite(fv) and fv > 0:
                    dr_vals.append(fv)
            except (ValueError, TypeError):
                pass
    dr_vals = np.array(dr_vals)
    mon_mean = np.mean(dr_vals)
    mon_sd = np.std(dr_vals)
    mon_cv_raw = mon_sd / mon_mean
    print(f'Monitoring: n={len(dr_vals)}, mean={mon_mean:.2f}, '
          f'sd={mon_sd:.2f}, CV={mon_cv_raw:.3f}')
else:
    # Fallback defaults from Heshang Cave 2005-2015
    mon_mean = 16.04
    mon_sd = 6.53
    mon_cv_raw = 0.407
    print('Using default Heshang monitoring values')

# Integration-time corrected CV: sigma_eff = sigma / sqrt(T * f_s)
T_sample = args.T_sample
f_s = args.f_s
mon_cv_corr = (mon_sd / np.sqrt(T_sample * f_s)) / mon_mean
print(f'Monitoring CV: raw={mon_cv_raw:.3f}, corrected={mon_cv_corr:.3f} '
      f'(T={T_sample}yr, f_s={f_s:.1f}/yr)')

# 5.2 ka event — use narrow window around event from raw data,
# not the coarse 250-yr bin which dilutes the short-lived event
event_centre = 5200
event_half = args.event_half  # ±N yr window to capture the ~60-yr event
event_mask = (ages >= event_centre - event_half) & (ages <= event_centre + event_half)
event_med = med[event_mask]
event_fin = np.isfinite(event_med) & (event_med > 0)
event_med = event_med[event_fin]
if len(event_med) > 2:
    ev_mean = np.mean(event_med)
    ev_cv = np.std(event_med) / ev_mean if ev_mean > 0 else 0
    # Bootstrap CV and mean IQR for event window
    ev_boot_cv = np.empty(n_boot)
    ev_boot_mean = np.empty(n_boot)
    for b in range(n_boot):
        s = rng.choice(event_med, size=len(event_med), replace=True)
        bmu = np.mean(s)
        ev_boot_mean[b] = bmu
        ev_boot_cv[b] = np.std(s) / bmu if bmu > 0 else 0
    ev_mean_lo, ev_mean_hi = np.percentile(ev_boot_mean, [25, 75])
    ev_cv_lo, ev_cv_hi = np.percentile(ev_boot_cv, [25, 75])
    print(f'5.2 ka event (±{event_half} yr): n={len(event_med)}, '
          f'mean={ev_mean:.1f} [{ev_mean_lo:.1f}-{ev_mean_hi:.1f}], '
          f'CV={ev_cv:.3f} [{ev_cv_lo:.3f}-{ev_cv_hi:.3f}]')
else:
    # Fallback: use the minimum drip rate point in the 4.8-5.6 ka range
    wide_mask = (ages >= 4800) & (ages <= 5600) & np.isfinite(med) & (med > 0)
    wide_med = med[wide_mask]
    ev_mean = np.min(wide_med) if len(wide_med) > 0 else 5.0
    ev_cv = 0.25
    ev_mean_lo, ev_mean_hi = ev_mean, ev_mean
    ev_cv_lo, ev_cv_hi = ev_cv, ev_cv
    print(f'5.2 ka event (fallback min): mean={ev_mean:.1f}')

ev_yerr = [[max(0, ev_mean - ev_mean_lo)], [max(0, ev_mean_hi - ev_mean)]]
ev_xerr = [[max(0, ev_cv - ev_cv_lo)], [max(0, ev_cv_hi - ev_cv)]]

# Also keep the bin index so we can still exclude it from the regular scatter
idx_52 = (bdf['age_mid'] - 5200).abs().idxmin()
row_52_bin = bdf.loc[idx_52]


# ── Figure ───────────────────────────────────────────────────────────────
fig, ax = plt.subplots(figsize=(DOUBLE_COL, DOUBLE_COL * 0.72), facecolor='white')


# ── Panel a: S&F classification space ────────────────────────────────────
# Background zones
zones = [
    (0, 0.5, 0, 5, COL_TEAL_50 if False else COL_BG_50),
    (0.5, 1.2, 0, 5, COL_BG_50),
    (0, 0.5, 5, 35, COL_TEAL_200),
    (0.5, 1.2, 5, 35, COL_BROWN_200),
]
for x0, x1, y0, y1, col in zones:
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0,
                            fc=col, ec='none', alpha=0.35, zorder=0))
ax.axvline(0.5, color=COL_BG_300, lw=0.4, ls='--', zorder=1)
ax.axhline(5, color=COL_BG_300, lw=0.4, ls='--', zorder=1)

# Bold flow regime labels
label_props = dict(ha='center', va='center', fontsize=6.5, color=COL_BG_600, fontstyle='italic', zorder=20)
ax.text(0.15, 2.0, 'Seepage /\npercolation', **label_props)
ax.text(0.85, 2.0, 'Fracture /\nconduit', **label_props)
ax.text(0.85, 11, 'Flood / conduit\noverflow', **label_props)
ax.text(0.25, 6.4, 'Buffered overflow', ha='center', va='bottom', fontsize=6.5, color=COL_BG_600, fontstyle='italic', zorder=20)

# Age colormap — custom: burnt orange (young) → teal (mid) → dark navy (old)
# All saturated; visible on white; no yellow; avoids red/green marker clash
from matplotlib.colors import LinearSegmentedColormap
norm = Normalize(vmin=0, vmax=10000)
_age_colors = [COL_DORANGE, COL_TEAL_600, COL_TEAL_800]  # young → old (ngeo)
cmap = LinearSegmentedColormap.from_list('age_dark', _age_colors)

# Plot all bins except 5.2 ka (special marker)
for i, row in bdf.iterrows():
    if i == idx_52:
        continue
    c = cmap(norm(row['age_mid']))
    ax.errorbar(row['cv'], row['mean'],
                xerr=[[row['xerr_lo']], [row['xerr_hi']]],
                yerr=[[row['yerr_lo']], [row['yerr_hi']]],
                fmt='o', markersize=4.0, color=c, ecolor=c, elinewidth=0.4, capsize=1.2, capthick=0.3, markeredgecolor=COL_BG_900, markeredgewidth=0.25, alpha=0.9, zorder=5)

# 5.2 ka event bin -- special triangle marker (from narrow event window)
ax.errorbar(ev_cv, ev_mean, xerr=ev_xerr, yerr=ev_yerr,
            fmt='^', markersize=6, color=COL_RED_900, ecolor=COL_RED_900, elinewidth=0.4, capsize=1.2, capthick=0.3, markeredgecolor=COL_BG_900, markeredgewidth=0.3, zorder=9, label='5.2 ka event bin')

# Whole-record star
ax.scatter(sf['cv'], sf['mean'], c=COL_RED_900, s=90, marker='*', edgecolors=COL_BG_900, linewidth=0.3, zorder=8,
           label=f'Whole record (\u03bc={sf["mean"]:.1f}, CV={sf["cv"]:.2f})')
ax.errorbar(sf['cv'], sf['mean'],
            xerr=[[sf['cv'] - sf['cv_lo']], [sf['cv_hi'] - sf['cv']]],
            yerr=[[sf['mean'] - sf['mean_lo']], [sf['mean_hi'] - sf['mean']]],
            fmt='none', ecolor=COL_RED_900, elinewidth=0.4, capsize=1.5, capthick=0.3, zorder=7)

# Monitoring diamonds: raw (open) and corrected (filled)
ax.scatter(mon_cv_raw, mon_mean, c='none', s=45, marker='D', edgecolors=COL_GREEN_800, linewidth=0.9, zorder=8,
           label=f'Monitoring raw (CV={mon_cv_raw:.2f})')
ax.scatter(mon_cv_corr, mon_mean, c=COL_GREEN_800, s=45, marker='D', edgecolors=COL_BG_900, linewidth=0.3, zorder=9,
           label=f'Monitoring corrected (CV={mon_cv_corr:.3f})')

# Colorbar
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array([])
cbar = fig.colorbar(sm, ax=ax, shrink=0.55, pad=0.02, aspect=25)
cbar.set_label('Age (yr BP)', fontsize=6.5)
cbar.ax.tick_params(labelsize=5.5, width=0.25, length=2.5)
cbar.ax.invert_yaxis()

ax.set_xlabel('Coefficient of variation (CV = \u03c3/\u03bc)')
ax.set_ylabel('Mean drip rate (drips min$^{-1}$)')
ax.set_xlim(-0.02, 1.1)
ax.set_ylim(0, 35)


# Legend
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, labels, loc='upper right', frameon=False, bbox_to_anchor=(0.98, 0.98), ncol=1, fontsize=5.5, columnspacing=1.2, handletextpad=0.5)
style_ax(ax)
cbar.outline.set_linewidth(0.25)


# ── Save ─────────────────────────────────────────────────────────────────
plt.savefig(f'{args.output}.png', dpi=600, bbox_inches='tight',
            facecolor='white')
plt.savefig(f'{args.output}.pdf', bbox_inches='tight', facecolor='white')
plt.savefig(f'{args.output}.eps', bbox_inches='tight', facecolor='white',
            format='eps')
print(f'Saved {args.output}.png (600 dpi), .pdf, and .eps')
