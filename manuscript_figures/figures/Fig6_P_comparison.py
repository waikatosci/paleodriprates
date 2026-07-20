"""
Fig 6 — precipitation and model comparison, one-step chain.

Reads from HS4_SourceData.xlsx:
  05b_calibration             P/T vs baseflow calibration
  Fig6_WangT                  Wang et al. temperature reconstruction
  Fig6_P_recon_LR / _dense    one-step P reconstruction (σ = π/√6)
  Fig5_driprate_AP            age-propagated drip-rate IQR
  Fig6_Hu / Fig6_TraCE        comparator reconstructions

  a (2/3 width): MC-propagated local P with gradient IQR envelope + δ¹⁸O-differencing P.
  b (1/3 width): data-model anomaly bar chart (anomalies in %, level-invariant).

Calibration is the single regression P/T_h = a_h·baseflow + b_h on median baseflow +
CY precipitation + HS temperature, 2004-2023, n=20 (a_h=0.00439, b_h=+0.08313). Input
reconstructions are the one-step products (73 Wang-T horizons; drip-grid dense).
HU_MODE and SHOW_MODERN_P are author display choices. INSOL_FILE is not read by Fig 6.
"""

import os as _os, sys as _sys
# Resolve everything against the repo root so generate_figures.py can run each script
# from one place, and so a figure run from its own directory behaves identically.
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _ROOT)
_os.chdir(_ROOT)
_os.makedirs('output', exist_ok=True)

import os, warnings
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from matplotlib.patches import Patch
from scipy.interpolate import interp1d
warnings.filterwarnings('ignore')

# ══════════════════════════════════════════════════════════════════════
# FILE PATHS
# ══════════════════════════════════════════════════════════════════════
AGEPROP_FILE  = 'drip_rate_percentiles.xlsx'       # age-propagated DR (high-res)

# One-step MC P reconstructions (2026-07-15 rebuild)
MC_RECON_LR   = 'p_reconstruction_onestep.csv'        # 73 Wang-T horizons (primary, smooth)
MC_RECON_HR   = 'p_reconstruction_onestep_dense.csv'  # drip-grid res (faint overlay)

SHOW_MODERN_P = False   # dashed reference at modern local P (CY 2004–2023 mean)

# δ¹⁸O differencing curve handling — AUTHOR DECISION PENDING:
#   'raw'    plot Hu as published (Holocene mean ~1407 mm/yr; sits ~40% above
#            the local one-step level — levels NOT comparable)
#   'scaled' multiply Hu by (one-step Holocene mean / Hu Holocene mean) so the
#            curves share a level; SHAPE comparison only, factor printed and
#            must be declared in the caption
#   'off'    omit from panel a
HU_MODE = 'raw'

D18O_P_FILE   = 'Hu_et_al_rainfall_reconstruction.csv'
TRACE_FILE    = 'East_Asia_TraCE21KII_PRECT_mm_yr_LGM_to_present.csv'
INSOL_FILE    = 'insolation_30N_Jun21.csv'         # not plotted, retained for reference

OUTPUT_PDF = 'output/Fig6_P_comparison.pdf'
OUTPUT_PNG = 'Fig6_P_comparison.png'

# ── Styling constants (consistent with Fig 5) ────────────────────
# LR — primary smooth signal
LR_LW      = 0.5
LR_ZORDER  = 4

# HR — high-frequency detail, faint background
HR_LW      = 0.3
HR_ALPHA   = 0.25
HR_ZORDER  = 3

# Gradient envelope
N_GRAD_BANDS = 8      # bands from median to each IQR bound
GRAD_ALPHA_MAX = 0.55  # alpha at median edge
GRAD_ALPHA_MIN = 0.08  # alpha at IQR bound

# ══════════════════════════════════════════════════════════════════════
# CALIBRATION (drip rate → precipitation) — ONE-STEP, local
#   P/T_h = a_h·baseflow + b_h  (median baseflow + CY + HS_Tem, 2004–2023)
#   P     = P/T_h · T(age) · 365.25, T from Wang et al. (2018) RAN15-MAAT
# The Yichang transfer (a_y/b_y) is GONE. Do not reintroduce it.
# ══════════════════════════════════════════════════════════════════════
from scipy.stats import linregress, t as t_dist

SOURCE     = 'HS4_SourceData.xlsx'              # <<< the single source


def _find_header(sheet, first_col_name, occurrence=0):
    """Locate a sheet's header row by name rather than a hardcoded skiprows."""
    raw = pd.read_excel(SOURCE, sheet_name=sheet, header=None)
    hits = raw.index[raw[0].astype(str).str.strip() == first_col_name]
    if len(hits) <= occurrence:
        raise ValueError(f"header '{first_col_name}' #{occurrence} not in '{sheet}'")
    return int(hits[occurrence])


def _sheet(name, first_col, occurrence=0):
    return pd.read_excel(SOURCE, sheet_name=name,
                         skiprows=_find_header(name, first_col, occurrence)
                         ).dropna(how='all')


CAL_FILE   = 'calibration_onestep.csv'          # 20-yr annual table (2026-07-15 rebuild)
T_REC_FILE = 'T_recon_Wang_et_al.xlsx'          # palaeo T for the analytical conversion
T_MODERN   = None                               # filled from the calibration table below

_cal = _sheet('05b_calibration', 'year')
_cal = _cal[_cal['year'].apply(lambda v: isinstance(v, (int, float)))]
_cal = _cal.set_index('year').rename(columns={
    'baseflow_drips_min': 'baseflow', 'P_mm': 'P', 'T_degC': 'T'})
_sl1, _int1, _r1, _p1, _se_sl1 = linregress(_cal['baseflow'], _cal['PT_h'])
_n1 = len(_cal)
_resid1 = (_cal['PT_h'] - (_sl1*_cal['baseflow'] + _int1)).std(ddof=2)
_d_bar = _cal['baseflow'].mean()
_SSd = np.sum((_cal['baseflow'] - _d_bar)**2)
_t1 = t_dist.ppf(0.975, _n1 - 2)
T_MODERN  = _cal['T'].mean()          # ~17.5 °C, local HS_Tem 2004–2023
P_MODERN  = _cal['P'].mean()          # ~1002 mm/yr, local CY 2004–2023

print(f"  One-step calibration (baseflow→P/T_h): a_h={_sl1:.5f}±{_se_sl1:.5f}, "
      f"b_h={_int1:+.5f}, R²={_r1**2:.3f}, p={_p1:.4f}, n={_n1}")
print(f"  Modern local reference: P={P_MODERN:.0f} mm/yr, T={T_MODERN:.2f} °C")
assert abs(_sl1 - 0.00439) < 5e-5, "a_h drifted from SI value 0.00439"

# Palaeotemperature for the analytical (non-MC) conversion path
_T_interp = None
if True:
    _tr = _sheet('Fig6_WangT', 'Age (yr, BP)')
    _tr = _tr[['Age (yr, BP)', 'RAN15-MAAT(°C)']].dropna()
    _T_interp = interp1d(_tr.iloc[:, 0].values, _tr.iloc[:, 1].values,
                         bounds_error=False,
                         fill_value=(_tr.iloc[0, 1], _tr.iloc[-1, 1]))

def _T_at(age_bp):
    """Local temperature at age (yr BP): Wang RAN15-MAAT, else modern mean."""
    if _T_interp is not None:
        return np.asarray(_T_interp(age_bp), float)
    return np.full_like(np.asarray(age_bp, float), T_MODERN)

def drip_to_precip(d, age_bp=None):
    """Point estimate: drip → local P (mm yr⁻¹), one-step chain."""
    d = np.asarray(d, float)
    T = _T_at(age_bp) if age_bp is not None else T_MODERN
    return np.maximum(0, (_sl1 * d + _int1) * T * 365.25)

def drip_to_precip_ci(d, age_bp=None):
    """Return (P, P_lower, P_upper) at 95% CI on the single regression."""
    d = np.asarray(d, float)
    T = _T_at(age_bp) if age_bp is not None else T_MODERN
    pt_h = _sl1 * d + _int1
    ci1 = _t1 * _resid1 * np.sqrt(1/_n1 + (d - _d_bar)**2 / _SSd)
    P = np.maximum(0, pt_h * T * 365.25)
    ci_total = ci1 * T * 365.25
    return P, np.maximum(0, P - ci_total), P + ci_total

# ══════════════════════════════════════════════════════════════════════
# MODEL CONSTANTS
# ══════════════════════════════════════════════════════════════════════
PMIP4_AGE  = 6.0
PMIP4_ANOM = 10.0
PMIP4_ERR  = 7.5

KCM_KA   = np.array([9, 6.5, 3, 0])
KCM_ANOM = np.array([44, 44, 33, 0])

PERIODS = [
    ("Early Hol.\n(9.5–8 ka)",   9500, 8000),
    ("Mid Hol.\n(8–5 ka)",       8000, 5000),
    ("5.2 ka\nevent",            5200, 5100),
    ("Late Hol.\n(5–1 ka)",      5000, 1000),
    ("Pre-ind.\n(1 ka–200 BP)",  1000,  200),
    ("Post-1820",                 130,  -70),
]

# ══════════════════════════════════════════════════════════════════════
# STYLE — Nature Communications house
# ══════════════════════════════════════════════════════════════════════
from ngeo_style import (apply_style, style_ax, draw_events, EVENTS_FIG6,
                        DOUBLE_COL, MM,
                        COL_HS4, COL_HS4_IQR,
                        COL_D18O_P, COL_D18O_FILL,
                        COL_BAR_HS4 as COL_BAR_AP, COL_TRACE, COL_PMIP, COL_KCM)
apply_style()
EVENTS = EVENTS_FIG6

def compute_kcm_anom(age_hi, age_lo):
    f = interp1d(KCM_KA, KCM_ANOM, bounds_error=False,
                 fill_value=(KCM_ANOM[0], KCM_ANOM[-1]))
    return f(np.linspace(age_lo / 1000, age_hi / 1000, 50)).mean()

def compute_trace_anom(trace_df, age_hi, age_lo):
    """Compute TraCE anomaly for a period window."""
    if trace_df is None:
        return np.nan
    lo_ka, hi_ka = age_lo / 1000, age_hi / 1000
    mask = (trace_df['Age_ka_BP'] >= lo_ka) & (trace_df['Age_ka_BP'] <= hi_ka)
    for col in ['Precip_anomaly_percent', 'precip_anomaly_percent',
                'P_anomaly_pct', 'anom_pct']:
        if col in trace_df.columns:
            vals = trace_df.loc[mask, col]
            return vals.mean() if len(vals) > 0 else np.nan
    for col in ['Precip_mm_yr', 'PRECT_mm_yr', 'precip']:
        if col in trace_df.columns:
            vals = trace_df.loc[mask, col]
            all_mean = trace_df[col].mean()
            if all_mean > 0 and len(vals) > 0:
                return (vals.mean() - all_mean) / all_mean * 100
    return np.nan


def gradient_envelope(ax, x, y_med, y_lo, y_hi, color, n_bands=N_GRAD_BANDS,
                      alpha_max=GRAD_ALPHA_MAX, alpha_min=GRAD_ALPHA_MIN,
                      zorder=2):
    """
    Draw a gradient-shaded uncertainty envelope that fades from opaque
    near the median to transparent at the IQR bounds.
    Mimics the PDF heatmap density aesthetic from Fig 5.
    """
    alphas = np.linspace(alpha_max, alpha_min, n_bands)
    # Upper half: median → pc75
    for i in range(n_bands):
        frac_lo = i / n_bands
        frac_hi = (i + 1) / n_bands
        band_lo = y_med + frac_lo * (y_hi - y_med)
        band_hi = y_med + frac_hi * (y_hi - y_med)
        ax.fill_between(x, band_lo, band_hi,
                        color=color, alpha=alphas[i], lw=0, zorder=zorder)
    # Lower half: median → pc25
    for i in range(n_bands):
        frac_lo = i / n_bands
        frac_hi = (i + 1) / n_bands
        band_lo = y_med - frac_hi * (y_med - y_lo)
        band_hi = y_med - frac_lo * (y_med - y_lo)
        ax.fill_between(x, band_lo, band_hi,
                        color=color, alpha=alphas[i], lw=0, zorder=zorder)

# ══════════════════════════════════════════════════════════════════════
# LOAD DATA
# ══════════════════════════════════════════════════════════════════════
print("Loading ...")

# ── MC-propagated P reconstruction — LR (primary) and HR (overlay) ──
def load_mc(path, label):
    df = _sheet('Fig6_P_recon_LR' if label == 'LR' else 'Fig6_P_recon_dense', 'age')
    df = df.sort_values('age')
    df['ka'] = df['age'] / 1000
    print(f"  MC {label}: {len(df)} pts, {df.ka.min():.2f}–{df.ka.max():.2f} ka")
    return df

mc_lr = load_mc(MC_RECON_LR, 'LR')
mc_hr = load_mc(MC_RECON_HR, 'HR')

# Bar chart interpolates the DENSE series when available — the 73-horizon
# primary is too sparse for narrow windows (e.g. the 5.2 ka event, 100 yr).
mc = mc_hr if mc_hr is not None else mc_lr
has_mc = mc is not None

if has_mc:
    p_interp = interp1d(mc['age'].values, mc['P_med'].values,
                        kind='linear', bounds_error=False, fill_value=np.nan)
    mc_hol = mc[(mc.age >= 0) & (mc.age <= 9500)]
    P_hol_mean = mc_hol['P_med'].mean()

# ── Drip rate percentiles (analytical fallback / bar chart supplement) ──
ap = _sheet('Fig5_driprate_AP', 'age_calBP')
ap.columns = [c.strip().lower().replace(' ', '_') for c in ap.columns]
age_col = [c for c in ap.columns if 'age' in c][0]
dr_cols = {
    'med':  [c for c in ap.columns if 'med' in c][0],
    'pc25': [c for c in ap.columns if '25' in c][0],
    'pc75': [c for c in ap.columns if '75' in c][0],
}
# NOTE: drip_rate_percentiles.xlsx is the PUBLISHED posterior. Under the
# one-step chain it must be rescaled by s = 0.8379 (the mu re-anchoring;
# applied before conversion.
_S_REANCHOR = 0.8379
ap['P_med']  = drip_to_precip(ap[dr_cols['med']].values  * _S_REANCHOR, ap[age_col].values)
ap['P_pc25'] = drip_to_precip(ap[dr_cols['pc25']].values * _S_REANCHOR, ap[age_col].values)
ap['P_pc75'] = drip_to_precip(ap[dr_cols['pc75']].values * _S_REANCHOR, ap[age_col].values)
ap['age_bp'] = ap[age_col]
ap['ka']     = ap['age_bp'] / 1000

if not has_mc:
    # Fallback: use analytical conversion
    p_interp = interp1d(ap['age_bp'].values, ap['P_med'].values,
                        kind='linear', bounds_error=False, fill_value=np.nan)
    ap_hol = ap[(ap.age_bp >= 0) & (ap.age_bp <= 9500)]
    P_hol_mean = ap_hol['P_med'].mean()

print(f"  Drip rate percentiles: {len(ap)} pts")
print(f"  Holocene mean P: {P_hol_mean:.0f} mm/yr")

# δ¹⁸O differencing P
has_d18o = HU_MODE != 'off'
if has_d18o:
    d18o_p = _sheet('Fig6_Hu', 'age_ka')
    # Detect age units
    if d18o_p.iloc[:, 0].max() > 100:
        d18o_p['age_ka'] = d18o_p.iloc[:, 0] / 1000
    else:
        d18o_p['age_ka'] = d18o_p.iloc[:, 0]
    print(f"  δ¹⁸O differencing P: {len(d18o_p)} pts (HU_MODE='{HU_MODE}')")
else:
    print(f"  δ¹⁸O differencing P: omitted (HU_MODE='{HU_MODE}' or file not found)")

# TraCE-21K-II
has_trace = True
trace_df = None
if has_trace:
    trace_df = _sheet('Fig6_TraCE', 'Age_ka_BP')
    print(f"  TraCE-21K-II: {len(trace_df)} pts")
else:
    print(f"  TraCE-21K-II: not found — will omit from bars")

# ══════════════════════════════════════════════════════════════════════
# FIGURE — single row: a (2/3) | b (1/3)
# ══════════════════════════════════════════════════════════════════════
print("Rendering ...")
fig = plt.figure(figsize=(180 / 25.4, 70 / 25.4))  # 180 mm wide, compact height

gs = gridspec.GridSpec(1, 2, figure=fig, wspace=0.18,
                       width_ratios=[2, 1],
                       left=0.07, right=0.97, top=0.93, bottom=0.18)

# ─────────────────────────────────────────
# PANEL a — Holocene P
# ─────────────────────────────────────────
ax_a = fig.add_subplot(gs[0, 0])
style_ax(ax_a)
draw_events(ax_a, EVENTS)

# Age-propagated P — dual resolution with gradient envelope
if has_mc:
    # ── HR overlay: faint background showing high-frequency detail ──
    if mc_hr is not None:
        gradient_envelope(ax_a, mc_hr['ka'].values,
                          mc_hr['P_med'].values,
                          mc_hr['P_pc25'].values,
                          mc_hr['P_pc75'].values,
                          color=COL_HS4_IQR,
                          alpha_max=GRAD_ALPHA_MAX * 0.35,
                          alpha_min=GRAD_ALPHA_MIN * 0.35,
                          zorder=HR_ZORDER - 1)
        ax_a.plot(mc_hr['ka'].values, mc_hr['P_med'].values,
                  color=COL_HS4, lw=HR_LW, alpha=HR_ALPHA, zorder=HR_ZORDER)

    # ── LR: primary smooth signal with gradient envelope ──
    gradient_envelope(ax_a, mc_lr['ka'].values,
                      mc_lr['P_med'].values,
                      mc_lr['P_pc25'].values,
                      mc_lr['P_pc75'].values,
                      color=COL_HS4_IQR, zorder=LR_ZORDER - 1)
    ax_a.plot(mc_lr['ka'].values, mc_lr['P_med'].values,
              color=COL_HS4, lw=LR_LW, zorder=LR_ZORDER,
              label='HS4 kinetic proxy')

    # Auto-scale y-axis to data range with padding
    all_lo = mc_lr['P_pc25'].min()
    all_hi = mc_lr['P_pc75'].max()
    if mc_hr is not None:
        all_lo = min(all_lo, mc_hr['P_pc25'].min())
        all_hi = max(all_hi, mc_hr['P_pc75'].max())
    y_pad = (all_hi - all_lo) * 0.12
    y_lo = max(0, all_lo - y_pad)
    y_hi = all_hi + y_pad
else:
    # Analytical fallback (drip rate IQR only, fixed T=16°C)
    gradient_envelope(ax_a, ap['ka'].values,
                      ap['P_med'].values,
                      ap['P_pc25'].values,
                      ap['P_pc75'].values,
                      color=COL_HS4_IQR, zorder=LR_ZORDER - 1)
    ax_a.plot(ap['ka'].values, ap['P_med'].values,
              color=COL_HS4, lw=LR_LW, zorder=LR_ZORDER,
              label='HS4 kinetic proxy')
    y_pad = (ap['P_pc75'].max() - ap['P_pc25'].min()) * 0.12
    y_lo = max(0, ap['P_pc25'].min() - y_pad)
    y_hi = ap['P_pc75'].max() + y_pad

# δ¹⁸O differencing P
if has_d18o:
    if HU_MODE == 'raw':
        print("  !! Hu δ¹⁸O differencing P is plotted at its published level"
              " (~40% above local one-step) — levels are NOT comparable."
              " Set HU_MODE='scaled' for a shape-only comparison, or 'off'.")
    cols = d18o_p.columns
    p_col = [c for c in cols if 'rainfall' in c.lower() and 'shift' not in c.lower()
             and 'neg' not in c.lower() and 'pos' not in c.lower()]
    neg_col = [c for c in cols if 'neg' in c.lower()]
    pos_col = [c for c in cols if 'pos' in c.lower()]
    p_name = p_col[0] if p_col else cols[-1]  # fallback to last column

    hu_factor = 1.0
    hu_label = r'$\delta^{18}$O differencing'
    if HU_MODE == 'scaled' and has_mc:
        _hu_hol = d18o_p[(d18o_p.age_ka >= mc['ka'].min()) & (d18o_p.age_ka <= mc['ka'].max())]
        hu_factor = mc['P_med'].mean() / _hu_hol[p_name].mean()
        hu_label += ' (scaled)'
        print(f"  Hu curve scaled by ×{hu_factor:.3f} to match one-step Holocene"
              f" mean — declare this factor in the caption.")

    if neg_col and pos_col:
        ax_a.fill_between(d18o_p['age_ka'].values,
                          d18o_p[neg_col[0]].values * hu_factor,
                          d18o_p[pos_col[0]].values * hu_factor,
                          color=COL_D18O_FILL, alpha=0.40, zorder=1, lw=0)
    ax_a.plot(d18o_p['age_ka'].values, d18o_p[p_name].values * hu_factor,
              color=COL_D18O_P, lw=0.7, zorder=3.6, label=hu_label)
    # zorder 3.6: above both gradient envelopes (2, 3), below the LR median (4).
    # At zorder 2 the curve is buried under the envelope whenever it shares the
    # reconstruction's level — invisible in 'scaled' mode.
    # Hu participates in the auto y-limits
    y_hi = max(y_hi, (d18o_p[p_name].values * hu_factor).max() * 1.03)
    y_lo = min(y_lo, (d18o_p[p_name].values * hu_factor).min() * 0.97)

# Optional modern local reference (CY 2004–2023 mean, ~1002 mm/yr)
if SHOW_MODERN_P:
    ax_a.axhline(P_MODERN, color='#546E7A', lw=0.4, ls=(0, (4, 3)), zorder=1.5)
    ax_a.text(8.9, P_MODERN, f'modern local ({P_MODERN:.0f})', fontsize=4.5,
              color='#546E7A', va='bottom', ha='left', zorder=6)

ax_a.set_ylabel(r'Local precipitation (mm yr$^{-1}$)', fontsize=6.5)
ax_a.set_ylim(y_lo, y_hi)   # auto-scaled — one-step local P sits ~18% below Yichang-equiv.

# Event labels — position relative to axis limits
for label, key in [('8.2 ka', '8.2 ka'), ('5.2 ka', '5.2 ka'), ('Post-1750', 'Post-1750')]:
    ev = EVENTS[key]
    x = ev['c'] if key != 'Post-1750' else 0.25
    ha = 'center' if key != 'Post-1750' else 'right'
    ax_a.text(x, y_hi - 0.02 * (y_hi - y_lo), label, ha=ha, va='top', fontsize=5,
              color='#90A4AE', fontstyle='italic', zorder=6)
ax_a.set_xlim(9, 0)
ax_a.xaxis.set_major_locator(ticker.MultipleLocator(1))
ax_a.xaxis.set_minor_locator(ticker.MultipleLocator(0.5))
ax_a.set_xlabel('Age (ka BP)', fontsize=6.5)
ax_a.legend(fontsize=5, framealpha=0.9, handlelength=1.5,
            loc='lower left', borderpad=0.3, handletextpad=0.4)
ax_a.text(0.01, 0.96, 'a', transform=ax_a.transAxes,
          fontsize=9, fontweight='bold', va='top')

# ─────────────────────────────────────────
# PANEL b — Data–model anomaly bars
# ─────────────────────────────────────────
ax_b = fig.add_subplot(gs[0, 1])
style_ax(ax_b)

labels = []
ap_anoms = []; trace_anoms = []; pmip4_anoms = []; kcm_anoms = []

for plabel, age_hi, age_lo in PERIODS:
    labels.append(plabel)

    # HS4 age-propagated
    dense = np.arange(age_lo, age_hi + 1, 1)
    dense_P = p_interp(dense)
    valid = np.isfinite(dense_P)
    if valid.sum() > 10:
        ap_anoms.append((np.median(dense_P[valid]) - P_hol_mean) / P_hol_mean * 100)
    else:
        ap_anoms.append(np.nan)

    # TraCE
    trace_anoms.append(compute_trace_anom(trace_df, age_hi, age_lo))

    # PMIP4 (mid-Holocene only)
    pmip4_anoms.append(PMIP4_ANOM if age_lo <= 6000 <= age_hi else np.nan)

    # KCM
    kcm_anoms.append(compute_kcm_anom(age_hi, age_lo))

n_per = len(labels)
bar_w = 0.18
x = np.arange(n_per)
offsets = np.array([-1.5, -0.5, 0.5, 1.5]) * bar_w

data_sets  = [ap_anoms, trace_anoms, pmip4_anoms, kcm_anoms]
colours    = [COL_BAR_AP, COL_TRACE, COL_PMIP, COL_KCM]
bar_labels = ['HS4 age-prop.', 'TraCE-21K-II', 'PMIP4', 'KCM S. China']

for i, (vals, col) in enumerate(zip(data_sets, colours)):
    vv = [v if np.isfinite(v) else 0 for v in vals]
    vm = [np.isfinite(v) for v in vals]
    bars = ax_b.bar(x + offsets[i], vv, bar_w, color=col, alpha=0.90,
                    edgecolor='white', lw=0.35, zorder=3)
    for bar, vis in zip(bars, vm):
        if not vis:
            bar.set_visible(False)

# PMIP4 error bars
for j, pm in enumerate(pmip4_anoms):
    if np.isfinite(pm):
        ax_b.errorbar(x[j] + offsets[2], pm, yerr=PMIP4_ERR,
                      fmt='none', ecolor=COL_PMIP, elinewidth=0.5,
                      capsize=1.5, zorder=5)

legend_patches = [Patch(facecolor=c, edgecolor='white', lw=0.35,
                        alpha=0.90, label=l)
                  for c, l in zip(colours, bar_labels)]

ax_b.axhline(0, color='#263238', lw=0.3, ls='-', zorder=2)
ax_b.set_xticks(x)
ax_b.set_xticklabels(labels, fontsize=4.5, ha='center')
ax_b.set_ylabel('P anomaly (%)', fontsize=6.5)
ax_b.set_ylim(-30, 55)
ax_b.legend(handles=legend_patches, fontsize=4.5, loc='lower left',
            framealpha=0.9, handlelength=0.8, borderpad=0.3,
            handletextpad=0.3, labelspacing=0.3)
ax_b.text(0.02, 0.96, 'b', transform=ax_b.transAxes,
          fontsize=9, fontweight='bold', va='top')

# ══════════════════════════════════════════════════════════════════════
# SAVE
# ══════════════════════════════════════════════════════════════════════
for fmt in ('pdf', 'png', 'eps'):
    path = OUTPUT_PDF.replace('.pdf', f'.{fmt}')
    fig.savefig(path, dpi=600, bbox_inches='tight')
    print(f"  Saved → {path}")
plt.close(fig)
print("Done.")
