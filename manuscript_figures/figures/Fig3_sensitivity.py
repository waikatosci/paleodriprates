"""
Fig 3 — forward-model OMC dissociation sensitivity: Ni and Co in calcite vs drip rate.

Reads 'Fig3_sensitivity' from HS4_SourceData.xlsx. Concentrations are in ppm in
calcite. Model, palette, layout and ngeo_style unchanged (scale_ni=25).
"""

import os as _os, sys as _sys
# Resolve everything against the repo root so generate_figures.py can run each script
# from one place, and so a figure run from its own directory behaves identically.
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _ROOT)
_os.chdir(_ROOT)
_os.makedirs('output', exist_ok=True)

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from ngeo_style import (apply_style, style_ax, COL_NI, COL_NI_FILL,
                        COL_CO, COL_CO_FILL, COL_KD, SINGLE_COL, MM,
                        COL_TEAL_300, COL_BG_400, COL_BG_600)
apply_style()

# Override with more saturated tiers — still within the palette (as parent)
COL_NI_FILL = COL_TEAL_300
COL_CO_FILL = COL_BG_400
COL_KD      = COL_BG_600

# ── THE SINGLE SOURCE ────────────────────────────────────────────────────────
SOURCE   = 'HS4_SourceData.xlsx'
SHEET    = 'Fig3_sensitivity'
OUT_STEM = 'output/Fig3_sensitivity'

def _find_header(sheet, first_col_name):
    """Locate a sheet's header row by name rather than by a hardcoded skiprows.

    Sheets grow notes and provenance blocks above the table as findings accumulate.
    Binding to a fixed row number makes every script hostage to that; binding to the
    column name does not.
    """
    raw = pd.read_excel(SOURCE, sheet_name=sheet, header=None)
    hits = raw.index[raw[0].astype(str).str.strip() == first_col_name]
    if len(hits) == 0:
        raise ValueError(f"header '{first_col_name}' not found in sheet '{sheet}'")
    return int(hits[0])


data = pd.read_excel(SOURCE, sheet_name=SHEET,
                     skiprows=_find_header(SHEET, 'n_F'), usecols='A:H')
data = data.dropna(subset=['n_F'])
print(f"Loaded {SOURCE}[{SHEET}], shape: {data.shape}")

# canonical column names (the CSV's *_ppb are gone for good)
NI_M, NI_L, NI_U = 'Ni_mean_ppm', 'Ni_CI_lower_ppm', 'Ni_CI_upper_ppm'
CO_M, CO_L, CO_U = 'Co_mean_ppm', 'Co_CI_lower_ppm', 'Co_CI_upper_ppm'
DRIP = 'drip_rate_drips_per_min'

data = data[data['n_F'].isin([0.01])]
d_new = np.linspace(0, 30, 1000)
d_log = np.logspace(-2, np.log10(30), 1000)

interpolated = []
for lf in [0.01]:
    df = data[data['n_F'] == lf]
    row = {'n_F': lf, DRIP: d_new}
    for col in [NI_M, NI_L, NI_U, CO_M, CO_L, CO_U]:
        f = interp1d(df[DRIP], df[col], bounds_error=False, fill_value='extrapolate')
        row[col] = f(d_new)
    interpolated.append(pd.DataFrame(row))
data = pd.concat(interpolated, ignore_index=True)

scale_ni = 25
kd_values = [0.01 / 60, 0.1 / 60]
labels_kd = [r'$k_d$ = 0.01 s$^{-1}$', r'$k_d$ = 0.1 s$^{-1}$']
linestyles_kd = [':', '--']


def labile_single(kd, tau, lambda_F=0, lambda_I=0):
    return lambda_F + (1 - lambda_F - lambda_I) * (1 - np.exp(-kd * tau))


fig, ax = plt.subplots(figsize=(SINGLE_COL, 70 * MM))
style_ax(ax)

tau = 60 / np.where(d_new > 0, d_new, np.inf)
for kd, label, ls in zip(kd_values, labels_kd, linestyles_kd):
    ax.plot(d_new, scale_ni * labile_single(kd, tau), color=COL_KD, lw=0.5, ls=ls, label=label)

df = data[data['n_F'] == 0.01]
ax.fill_between(df[DRIP], df[NI_L], df[NI_U], color=COL_NI_FILL, alpha=0.55, lw=0)
ax.plot(df[DRIP], df[NI_M], color=COL_NI, lw=0.9, label=r'Ni ($\lambda_F$, $\lambda_I$ 95% CI)')
ax.fill_between(df[DRIP], df[CO_L], df[CO_U], color=COL_CO_FILL, alpha=0.55, lw=0)
ax.plot(df[DRIP], df[CO_M], color=COL_CO, lw=0.9, label=r'Co ($\lambda_F$, $\lambda_I$ 95% CI)')

ax.set_xlabel(r'Drip rate (drips min$^{-1}$)')
ax.set_ylabel('Metal concentration in calcite (ppm)')
ax.set_xlim(0, 30); ax.set_ylim(0, 30)
ax.set_xticks(np.arange(0, 31, 5)); ax.set_yticks([0, 10, 20, 30])

ax2 = ax.twiny()
ax2.spines['top'].set_visible(True); ax2.spines['top'].set_linewidth(0.25)
ax2.tick_params(axis='x', which='major', width=0.25, length=2.5, direction='in')
ax2.set_xlim(0, 30)
drip_ticks = np.array([0, 5, 10, 15, 20, 25, 30])
tau_t = 60 / np.where(drip_ticks > 0, drip_ticks, np.inf)
ax2.set_xticks(drip_ticks)
ax2.set_xticklabels(['∞'] + [f'{t:.1f}' if t < 100 else f'{int(t)}' for t in tau_t[1:]], fontsize=5.5)
ax2.set_xlabel(r'Residence time, $\tau$ (s)', fontsize=6.5)

ax.legend(loc='upper left', bbox_to_anchor=(0.02, 0.98), framealpha=0.9,
          borderpad=0.3, handlelength=1.5, handletextpad=0.4, labelspacing=0.3)

# Inset
inset = ax.inset_axes([0.55, 0.28, 0.42, 0.45])
style_ax(inset)
for lf in [0.01]:
    df0 = data[data['n_F'] == lf]
    for col, cl, cf in [('Ni', COL_NI, COL_NI_FILL), ('Co', COL_CO, COL_CO_FILL)]:
        m, lo, hi = (f'{col}_mean_ppm', f'{col}_CI_lower_ppm', f'{col}_CI_upper_ppm')
        fm = interp1d(df0[DRIP], df0[m], bounds_error=False, fill_value='extrapolate')
        fl = interp1d(df0[DRIP], df0[lo], bounds_error=False, fill_value='extrapolate')
        fu = interp1d(df0[DRIP], df0[hi], bounds_error=False, fill_value='extrapolate')
        inset.fill_between(d_log, fl(d_log), fu(d_log), color=cf, alpha=0.55, lw=0)
        inset.plot(d_log, fm(d_log), color=cl, lw=0.55)
tau_log = 60 / d_log
for kd, ls in zip(kd_values, linestyles_kd):
    inset.plot(d_log, scale_ni * labile_single(kd, tau_log), color=COL_KD, lw=0.4, ls=ls)
inset.set_xscale('log'); inset.set_xlim(0.01, 30); inset.set_ylim(0, 30)
inset.set_xticks([0.01, 0.1, 1, 10]); inset.set_xticklabels(['0.01', '0.1', '1', '10'], fontsize=4.5)
inset.set_yticks([0, 15, 30]); inset.set_yticklabels(['0', '15', '30'], fontsize=4.5)

inset2 = inset.twiny()
inset2.spines['top'].set_visible(True); inset2.spines['top'].set_linewidth(0.25)
inset2.tick_params(axis='x', which='major', width=0.25, length=2, direction='in')
inset2.set_xscale('log'); inset2.set_xlim(60 / 0.01, 60 / 30)
dti = np.array([0.01, 0.1, 1, 10]); tti = 60 / dti
inset2.set_xticks(tti)
inset2.set_xticklabels([f'{int(t)}' if t >= 1000 else f'{t:.0f}' for t in tti], fontsize=4.5)

# NOTE: EPS flattens transparency (the CI fills) — ship PDF. Kept for parity with parent.
for fmt in ('pdf', 'png', 'eps'):
    fig.savefig(f'{OUT_STEM}.{fmt}', dpi=600, bbox_inches='tight')
    print(f"Saved → {OUT_STEM}.{fmt}")
plt.close(fig); print("Done.")
