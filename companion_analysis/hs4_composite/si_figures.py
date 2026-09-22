#!/usr/bin/env python3
"""
si_figures.py -- the two Supplementary Figures built from the photographic composite of HS4
(Supplementary Methods 14.4; Supplementary Figures 21 and 22).

Inputs
  manuscript_figures/external/HS4_composite.jpg (+ .json)   stitched composite, 100 px/cm,
        row 0 = 0 cm on the pencil depth scale, pencil line at column 1650 (build_composite.py)
  extended_data/Ex_Data_6_age_depth/HS4_age_depth.csv       230Th sample depths
  dr_app/HS4_example_inputs/HS4_TE_canonical.csv            primary-run Ni and Co
  manuscript_figures/external/crosslab_replication_2019_corrected.csv
                                                            second-run Ni and Co, Ca-corrected

Outputs (manuscript_figures/output/)
  FigS_hs4_composite.{png,pdf}     Supplementary Figure 21
  FigS_hs4_5p2ka_photo.{png,pdf}   Supplementary Figure 22
"""
import json
from pathlib import Path

import cv2
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe

REPO = Path(__file__).resolve().parents[2]
EXT = REPO / 'manuscript_figures' / 'external'
OUT = REPO / 'manuscript_figures' / 'output'
NI, CO = '#0C57A0', '#D84315'
HALO = [pe.withStroke(linewidth=2.6, foreground='white', alpha=0.85)]
FRAG = (155.1, 157.3)          # loose fragment between the two breaks, on the pencil scale

plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 7})

meta = json.load(open(EXT / 'HS4_composite.json'))
IMG = cv2.cvtColor(cv2.imread(str(EXT / 'HS4_composite.jpg')), cv2.COLOR_BGR2RGB)
H, W = IMG.shape[:2]
PX, XL = meta['px_per_cm'], meta['depth_line_column_px']
EXTENT = [-XL / PX, (W - XL) / PX, H / PX, 0]      # cm from the pencil line; depth (cm)

te = pd.read_csv(REPO / 'dr_app' / 'HS4_example_inputs' / 'HS4_TE_canonical.csv')
te.columns = ['d', 'Ni', 'Co']
rep = pd.read_csv(EXT / 'crosslab_replication_2019_corrected.csv')
uth = pd.read_csv(REPO / 'extended_data' / 'Ex_Data_6_age_depth' / 'HS4_age_depth.csv', encoding='utf-8-sig')
uth.columns = ['d', 'age', 'err']
uth = uth[uth.err > 1]


def image_panel(ax, d0, d1, im, wcm=None):
    """Draw the composite between depths d0 and d1, cropped to the stone."""
    sc = im.shape[0] / H
    r0, r1 = int(d0 * PX * sc), int(d1 * PX * sc)
    cols = np.where((im[r0:r1].min(2) < 250).mean(0) > 0.02)[0]
    c0, c1 = cols.min(), cols.max() + 1
    e = [EXTENT[0] + c0 / (PX * sc), EXTENT[0] + c1 / (PX * sc), d1, d0]
    ax.imshow(im[r0:r1, c0:c1], extent=e, aspect='equal', interpolation='lanczos')
    ax.set_adjustable('box')
    if wcm:
        mid = (e[0] + e[1]) / 2
        ax.set_xlim(mid - wcm / 2, mid + wcm / 2)
    else:
        ax.set_xlim(e[0], e[1])
    ax.set_ylim(d1, d0)
    ax.set_xticks([])
    ax.set_ylabel('Depth (cm)')
    for sp in ax.spines.values():
        sp.set_linewidth(0.5)


def overlay_axis(ax):
    """A transparent axis over the middle half of the image (25-75% of its width)."""
    a = ax.inset_axes([0.25, 0, 0.5, 1])
    a.patch.set_alpha(0)
    a.set_yticks([])
    for sp in ('left', 'right', 'bottom'):
        a.spines[sp].set_visible(False)
    a.xaxis.set_ticks_position('top')
    a.xaxis.set_label_position('top')
    return a


# ── Supplementary Figure 21: the composite ─────────────────────────────────────
small = cv2.resize(IMG, None, fx=0.25, fy=0.25, interpolation=cv2.INTER_AREA)
cuts = [(0, 62), (62, 124), (124, 186), (186, 246)]
fig, axs = plt.subplots(1, 4, figsize=(180 / 25.4, 112 / 25.4))
for ax, (a, b) in zip(axs, cuts):
    image_panel(ax, a, b, small, wcm=23)
    u = uth[(uth.d >= a) & (uth.d <= b)]
    ax.plot(np.full(len(u), ax.get_xlim()[0]), u.d, marker='>', ls='none', ms=3.5,
            color='#263238', clip_on=False, zorder=5)
    if a <= FRAG[0] <= b:
        ax.plot([ax.get_xlim()[1]] * 2, FRAG, color='#B71C1C', lw=3, solid_capstyle='butt',
                clip_on=False, zorder=5)
        ax.text(ax.get_xlim()[1] + 0.6, np.mean(FRAG), '5.2 ka', color='#B71C1C', fontsize=6,
                va='center', ha='left', fontweight='bold', clip_on=False)
    if a == 0:
        ax.text(ax.get_xlim()[0] + 0.5, 1.2, '▸ ²³⁰Th samples', fontsize=5.5, color='#263238', va='center')
fig.tight_layout(w_pad=2.6)
for ext in ('png', 'pdf'):
    fig.savefig(OUT / f'FigS_hs4_composite.{ext}', dpi=300, bbox_inches='tight', pad_inches=0.04)
plt.close(fig)

# ── Supplementary Figure 22: the 5.2 ka horizon with Ni and Co ─────────────────
d0, d1 = 148, 164
fig, ax = plt.subplots(figsize=(120 / 25.4, 112 / 25.4))
image_panel(ax, d0, d1, IMG)
s = te[(te.d >= d0) & (te.d <= d1)]
r = rep[(rep.depth_cm >= d0) & (rep.depth_cm <= d1)]
a1, a2 = overlay_axis(ax), overlay_axis(ax)
for a, col, v, mx, lab in [(a1, NI, 'Ni', 12, 'Ni (ppm)'), (a2, CO, 'Co', 2.8, 'Co (ppm)')]:
    a.plot(s[v], s.d, '-o', color=col, alpha=0.8, ms=3, lw=1.2, mec='white', mew=0.4,
           path_effects=HALO, zorder=3, label=f'{v}, primary run')
    a.plot(r[v + '_corrected'], r.depth_cm, 'D', color=col, alpha=0.95, ms=4.2, mfc='white',
           mew=1.1, path_effects=HALO, zorder=4, label=f'{v}, second run (Ca-corrected)')
    a.set_xlim(0, mx)
    a.set_ylim(d1, d0)
    a.set_xticks(np.linspace(0, mx, 3))
    a.tick_params(colors=col, labelsize=6, length=2, width=0.5)
    a.set_xlabel(lab, color=col, labelpad=2)
    a.spines['top'].set_color(col)
    a.spines['top'].set_linewidth(0.6)
a2.spines['top'].set_position(('outward', 24))
h1, l1 = a1.get_legend_handles_labels()
h2, l2 = a2.get_legend_handles_labels()
fig.legend(h1 + h2, l1 + l2, loc='lower center', ncol=2, fontsize=6, frameon=False)
fig.tight_layout(rect=[0, 0.06, 1, 1])
for ext in ('png', 'pdf'):
    fig.savefig(OUT / f'FigS_hs4_5p2ka_photo.{ext}', dpi=300, bbox_inches='tight', pad_inches=0.04)
plt.close(fig)
print('saved FigS_hs4_composite and FigS_hs4_5p2ka_photo')
