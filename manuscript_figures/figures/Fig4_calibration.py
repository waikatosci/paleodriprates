"""
Fig 4 — calibration (three panels).

Reads from HS4_SourceData.xlsx:
  panel a   05b_calibration    P/T vs drip-rate calibration (n=20). PT_h is
            recomputed as P/(365.25*T).
  panel b   Fig4_panelB_kd     Ni/Co dissociation-rate distribution parameters
            Fig4_panelB_lab    laboratory replicates (dashed curves)
  panel c   Fig4_panelC_PT     P/T validation series (Yichang)

Panel b is regenerated from parameters: N(mu, sigma) on linspace(-20, 10, 1001), with
sigma = pi/sqrt(6) = 1.282550 (fixed a priori; SIGMA_CELL below). Panel c shows the
absolute-level check only; its R2/p annotation is not drawn.
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
from scipy.stats import linregress

from ngeo_style import (apply_style, style_ax, MM, DOUBLE_COL,
                        COL_NI, COL_NI_FILL, COL_NI_LAB,
                        COL_CO, COL_CO_FILL, COL_CO_LAB,
                        COL_FIT, COL_SCATTER, COL_OBS, COL_RECON, COL_STAT)
apply_style()

# ══════════════════════════════════════════════════════════════════════
# THE SINGLE SOURCE
# ══════════════════════════════════════════════════════════════════════
SOURCE = 'HS4_SourceData.xlsx'
OUT_STEM = 'output/Fig4_calibration'

# Which sigma panel b draws. 'sigma_PER_CAPTION_AND_SI' = pi/sqrt(6) = 1.2825,
# fixed a priori by the dissociation-rate sampling window (caption and SI Table S1/S3).
SIGMA_CELL = 'sigma_PER_CAPTION_AND_SI'

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


print("Loading ...")

# ── Panel a — one-step calibration, 2004-2023, n=20 ──
_cal = pd.read_excel(SOURCE, sheet_name='05b_calibration',
                     skiprows=_find_header('05b_calibration', 'year'), usecols='A:E')
_cal = _cal.dropna(subset=['year'])
_cal['year'] = _cal['year'].astype(int)
_cal = _cal.set_index('year')
# PT_h recomputed from P and T rather than trusted as a stored column
_cal['PT_h'] = _cal['P_mm'] / (365.25 * _cal['T_degC'])
df_a = pd.DataFrame({'P_T_mean':  _cal['PT_h'].values,
                     'drips_mean': _cal['baseflow_drips_min'].values})
print(f"  Panel A: {len(df_a)} years {_cal.index.min()}-{_cal.index.max()} from "
      f"{SOURCE}[05b_calibration]")

slope, intercept, r_value, p_value, std_err = linregress(df_a['P_T_mean'], df_a['drips_mean'])
r_sq = r_value ** 2
x_fit = np.linspace(df_a['P_T_mean'].min(), df_a['P_T_mean'].max(), 100)
y_fit = slope * x_fit + intercept
print(f"  Panel A: slope={slope:.3f}, R²={r_sq:.3f}, p={p_value:.4f}, n={len(df_a)}")

# Operational orientation (the chain actually used downstream)
_a_h, _b_h, _r_h, _p_h, _se_h = linregress(df_a['drips_mean'], df_a['P_T_mean'])
print(f"  Panel A operational: a_h={_a_h:.5f}±{_se_h:.5f}, b_h={_b_h:+.5f}")
assert abs(_a_h - 0.00439) < 5e-5, "a_h drifted from the SI value 0.00439"
assert abs(_b_h - 0.08313) < 5e-5, "b_h drifted from the SI value +0.08313"

# ── Panel b — HS4 k_d distributions REGENERATED from workbook parameters ──
_par_raw = pd.read_excel(SOURCE, sheet_name='Fig4_panelB_kd', usecols='A:B',
                         header=None, names=['k', 'v'])
_par = {str(k).split('  ')[0].strip(): v
        for k, v in zip(_par_raw['k'], _par_raw['v']) if pd.notna(k) and pd.notna(v)}
mu_ni = float(_par['mu_Ni'])
mu_co = float(_par['mu_Co'])
sigma = float(_par[SIGMA_CELL])
g0, g1, gn = float(_par['grid_min']), float(_par['grid_max']), int(_par['grid_n'])
print(f"  Panel B: mu_Ni={mu_ni}, mu_Co={mu_co}, sigma={sigma} ({SIGMA_CELL}), "
      f"grid=[{g0},{g1}] n={gn}")


def gaussian_pdf(x, mu, sd):
    """The raw (unnormalised-to-grid) normal pdf — exactly what Panel_C_updated.xlsx holds."""
    return np.exp(-0.5 * ((x - mu) / sd) ** 2) / (sd * np.sqrt(2 * np.pi))


_grid = np.linspace(g0, g1, gn)
df_ni = pd.DataFrame({'X(Ni)': _grid, 'Prob(Ni)': gaussian_pdf(_grid, mu_ni, sigma)})
df_co = pd.DataFrame({'X(Co)': _grid, 'Prob(Co)': gaussian_pdf(_grid, mu_co, sigma)})

# ── Panel b lab replicates ──
_lab = pd.read_excel(SOURCE, sheet_name='Fig4_panelB_lab', header=None)
_ni_hdr = _lab.index[_lab[0] == 'X1'][0]
_co_hdr = _lab.index[_lab[0] == 'X1'][1]
df_ni_lab = pd.read_excel(SOURCE, sheet_name='Fig4_panelB_lab', skiprows=_ni_hdr,
                          nrows=_co_hdr - _ni_hdr - 3).dropna(how='all', axis=1)
df_co_lab = pd.read_excel(SOURCE, sheet_name='Fig4_panelB_lab',
                          skiprows=_co_hdr).dropna(how='all', axis=1)
print(f"  Panel B lab: Ni {df_ni_lab.shape}, Co {df_co_lab.shape}")

# ── Panel c — P/T validation ──
df_b = pd.read_excel(SOURCE, sheet_name='Fig4_panelC_PT',
                     skiprows=_find_header('Fig4_panelC_PT', 'Age (P/T Yichang)'))
df_b = df_b.dropna(how='all')
valid_b = df_b['P/T (Moving Avg, Yichang)'].notna() & df_b['P/T (Reconstructed, Yichang)'].notna()
slope_b, intercept_b, r_b, p_b, se_b = linregress(
    df_b.loc[valid_b, 'P/T (Moving Avg, Yichang)'],
    df_b.loc[valid_b, 'P/T (Reconstructed, Yichang)'])
r_sq_b = r_b ** 2
print(f"  Panel C: R²={r_sq_b:.3f}, p={p_b:.2e}")


def calc_mean_2sd(x, p):
    v = x.notna() & p.notna()
    xv, pv = x[v].values, p[v].values
    if len(xv) == 0:
        return np.nan, np.nan
    mu = np.sum(pv * xv) / np.sum(pv)
    sd = np.sqrt(np.sum(pv * (xv - mu) ** 2) / np.sum(pv))
    return mu, 2 * sd


ni_mu, ni_2sd = calc_mean_2sd(df_ni['X(Ni)'], df_ni['Prob(Ni)'])
co_mu, co_2sd = calc_mean_2sd(df_co['X(Co)'], df_co['Prob(Co)'])
print(f"  Ni: ln μ = {ni_mu:.2f}, 2σ = {ni_2sd:.2f}")
print(f"  Co: ln μ = {co_mu:.2f}, 2σ = {co_2sd:.2f}")

# ══════════════════════════════════════════════════════════════════════
# FIGURE — double column (180 mm)
# ══════════════════════════════════════════════════════════════════════
fig, axs = plt.subplots(1, 3, figsize=(180 / 25.4, 55 / 25.4),
                        gridspec_kw={'width_ratios': [1, 1.5, 1], 'wspace': 0.35})
for ax in axs:
    style_ax(ax)

# ───────────────────────────────────────── PANEL a
axs[0].scatter(df_a['P_T_mean'], df_a['drips_mean'], s=15,
               facecolors='none', edgecolors=COL_SCATTER, linewidth=0.5, zorder=3)
axs[0].plot(x_fit, y_fit, color=COL_FIT, lw=0.9, zorder=2)

_yr = 2020
if _yr in _cal.index:
    _i = _cal.index.get_loc(_yr)
    axs[0].scatter(df_a['P_T_mean'].iloc[_i], df_a['drips_mean'].iloc[_i],
                   s=62, facecolors='none', edgecolors=COL_FIT, linewidth=0.6, zorder=4)
    print(f"  Ringed {_yr}: P/T={df_a['P_T_mean'].iloc[_i]:.4f}, "
          f"drip={df_a['drips_mean'].iloc[_i]:.2f}")

axs[0].set_xlabel('P/T')
axs[0].set_ylabel(r'Drip rate (drips min$^{-1}$)')
xm = (df_a['P_T_mean'].max() - df_a['P_T_mean'].min()) * 0.05
ym = (df_a['drips_mean'].max() - df_a['drips_mean'].min()) * 0.05
axs[0].set_xlim(df_a['P_T_mean'].min() - xm, df_a['P_T_mean'].max() + xm)
axs[0].set_ylim(df_a['drips_mean'].min() - ym, df_a['drips_mean'].max() + ym)
axs[0].text(0.95, 0.08, f'R² = {r_sq:.2f}\np = {p_value:.4f}',
            transform=axs[0].transAxes, fontsize=5, ha='right', va='bottom', color=COL_STAT)
axs[0].text(0.04, 0.96, 'a', transform=axs[0].transAxes,
            fontsize=9, fontweight='bold', va='top')

# ───────────────────────────────────────── PANEL b
axs[1].fill_between(df_ni['X(Ni)'], df_ni['Prob(Ni)'], color=COL_NI_FILL, alpha=0.55, zorder=3)
axs[1].plot(df_ni['X(Ni)'], df_ni['Prob(Ni)'], color=COL_NI, lw=0.9, label='HS4 Ni', zorder=4)
axs[1].fill_between(df_co['X(Co)'], df_co['Prob(Co)'], color=COL_CO_FILL, alpha=0.55, zorder=2)
axs[1].plot(df_co['X(Co)'], df_co['Prob(Co)'], color=COL_CO, lw=0.9, label='HS4 Co', zorder=3)


def plot_lab_avg(ax, df_lab, x_prefix, p_prefix, col, label):
    x_cols = [c for c in df_lab.columns if str(c).startswith(x_prefix)]
    p_cols = [c for c in df_lab.columns if str(c).startswith(p_prefix)]
    if not x_cols:
        return
    x_grid = np.linspace(df_lab[x_cols].min().min(), df_lab[x_cols].max().max(), 1000)
    avg_p = np.zeros_like(x_grid)
    for xc, pc in zip(x_cols, p_cols):
        v = df_lab[xc].notna() & df_lab[pc].notna()
        if v.any():
            avg_p += np.interp(x_grid, df_lab[xc][v], df_lab[pc][v], left=0, right=0) / len(x_cols)
    if avg_p.max() > 0:
        target_col = 'Prob(Ni)' if 'Ni' in label else 'Prob(Co)'
        target_df = df_ni if 'Ni' in label else df_co
        avg_p *= target_df[target_col].max() / avg_p.max()
    ax.plot(x_grid, avg_p, color=col, alpha=0.7, ls='--', lw=0.6, label=label)
    ax.fill_between(x_grid, avg_p, color=col, alpha=0.15)


plot_lab_avg(axs[1], df_ni_lab, 'X1', 'P_Ni1', COL_NI_LAB, 'Lab Ni')
plot_lab_avg(axs[1], df_co_lab, 'X1', 'P_Co1', COL_CO_LAB, 'Lab Co')

axs[1].axvline(ni_mu, color=COL_NI, ls=':', lw=0.5, alpha=0.7)
axs[1].axvline(co_mu, color=COL_CO, ls=':', lw=0.5, alpha=0.7)
axs[1].axvline(0, color='#263238', ls='-', lw=0.25)
axs[1].text(0.14, 0.62, f'Co\nln μ = {co_mu:.2f}\n2σ = {co_2sd:.2f}',
            transform=axs[1].transAxes, fontsize=5, va='center', ha='center', color=COL_CO)
axs[1].text(0.86, 0.62, f'Ni\nln μ = {ni_mu:.2f}\n2σ = {ni_2sd:.2f}',
            transform=axs[1].transAxes, fontsize=5, va='center', ha='center', color=COL_NI)
axs[1].set_xlabel(r'ln $k_d$ (s$^{-1}$)')
axs[1].set_ylabel('Probability')
axs[1].set_xlim(-20, 10)
all_y = pd.concat([df_ni['Prob(Ni)'].dropna(), df_co['Prob(Co)'].dropna()])
axs[1].set_ylim(0, all_y.max() * 1.05)
axs[1].legend(framealpha=0.9, borderpad=0.3, handlelength=1.2,
              handletextpad=0.3, labelspacing=0.25)
axs[1].text(0.04, 0.96, 'b', transform=axs[1].transAxes,
            fontsize=9, fontweight='bold', va='top')

# ───────────────────────────────────────── PANEL c
axs[2].plot(df_b['Age (P/T Yichang)'], df_b['P/T (Moving Avg, Yichang)'],
            color=COL_OBS, lw=0.9, label='Observed P/T', zorder=2)
axs[2].scatter(df_b['Age (P/T Yichang)'], df_b['P/T (Moving Avg, Yichang)'],
               s=3, color=COL_OBS, zorder=3)
axs[2].plot(df_b['Age (P/T Yichang)'], df_b['P/T (Reconstructed, Yichang)'],
            color=COL_RECON, lw=0.9, label='Reconstructed P/T (HS4)', zorder=2)
axs[2].scatter(df_b['Age (P/T Yichang)'], df_b['P/T (Reconstructed, Yichang)'],
               s=3, color=COL_RECON, zorder=3)
axs[2].set_xlabel('Year CE')
axs[2].set_ylabel('P/T')
axs[2].set_xlim(1950, 2000)
axs[2].set_ylim(0.165, 0.205)
# The R2/p annotation is not drawn: panel c is a check on the absolute level of the
# reconstruction, not an independent test of temporal skill. r_sq_b / p_b are still
# computed above and printed to stdout for the record.
axs[2].legend(loc='lower left', framealpha=0.9, borderpad=0.3,
              handlelength=1.2, handletextpad=0.3, labelspacing=0.25)
axs[2].text(0.04, 0.96, 'c', transform=axs[2].transAxes,
            fontsize=9, fontweight='bold', va='top')

plt.tight_layout()
for fmt in ('pdf', 'png', 'eps'):
    fig.savefig(f'{OUT_STEM}.{fmt}', dpi=600, bbox_inches='tight')
    print(f"Saved → {OUT_STEM}.{fmt}")
plt.close(fig)
print("Done.")
