"""
Fig 5 — Holocene record, two-panel layout.

Every tabular input comes from HS4_SourceData.xlsx:
  02_chronology          age model (depth -> age) and U-Th tie-points
  03_isotopes            δ¹⁸O, δ¹³C
  04_driprate_posterior  1 cm drip-rate percentiles (Fig 5a thicker lines)
  Fig5_insolation        30°N June insolation
  Fig5_driprate_AP       age-propagated drip-rate IQR
  Fig5_d18O_AP           δ¹⁸O age-propagated IQR (precomputed, 1,223 rows)

Bulk arrays live beside the workbook in external/:
  pdf_heatmap_hr.json, drip_rate_summary_hr.csv (drip-rate posterior, σ = π/√6).

Nature Communications production figure — two-panel layout.
  Panel a: Holocene drip rate PDF heatmap (high-res) + dual-resolution
           median/IQR overlay:
             - HR (native-res, ~576 pts): thin lines showing high-frequency
               epikarst / kinetic noise
             - LR (1 cm interpolated, ~258 pts): thicker lines showing the
               smoother signal backbone
           + age-propagated IQR envelope + 30°N June insolation +
           approximate precipitation axis.
  Panel b: δ¹⁸O (best-estimate chronology) + age-propagated IQR band.

Event zoom panels (c1–c3) moved to extended data figure.
Colour palette from ngeo_style for manuscript consistency.
"""

import os as _os, sys as _sys
# Resolve everything against the repo root so generate_figures.py can run each script
# from one place, and so a figure run from its own directory behaves identically.
_ROOT = _os.path.dirname(_os.path.dirname(_os.path.abspath(__file__)))
_sys.path.insert(0, _ROOT)
_os.chdir(_ROOT)
_os.makedirs('output', exist_ok=True)

import json, os, warnings
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.colors as mcolors
from matplotlib.transforms import blended_transform_factory
from scipy.interpolate import interp1d
warnings.filterwarnings('ignore')

# ══════════════════════════════════════════════════════════════════════
# PATHS
# ══════════════════════════════════════════════════════════════════════
# High-resolution (native BayProX spacing, ~576 pts)
HEATMAP_JSON_HR  = 'external/pdf_heatmap_hr.json'
DRIP_SUMMARY_HR  = 'external/drip_rate_summary_hr.csv'

# Low-resolution (1 cm interpolated, ~258 pts)
DRIP_SUMMARY_LR  = 'external/drip_rate_summary_lr.csv'

SOURCE        = 'HS4_SourceData.xlsx'   # <<< the single source
AGE_DEPTH     = 'HS4_age_depth.csv'      # label only; age model reads from SOURCE[02_chronology]
# ²³⁰Th tie-points come from 02_chronology (with 2σ errors), so the plotted
# tie-points and the age model cannot diverge.

# Output stem
INSOL_FILE    = 'insolation_30N_Jun21.csv'   # 30°N 21 June insolation (ka_BP, Wm2)

OUT_STEM      = 'output/Fig5_record'
DRIP_XLSX     = 'Drip_rate.xlsx'
ISO_SHEET     = '3.Isotopes'
AGEPROP_FILE  = 'drip_rate_percentiles.xlsx'
BASEFLOW_DRIP = 17.0
BASEFLOW_AGE  = -0.060   # ka BP

# δ¹⁸O age-propagated IQR (optional — from realisations or percentiles)
D18O_REAL_FILE    = 'd18O_realisations.csv'
D18O_AGEPROP_FILE = None     # e.g. 'd18O_percentiles.xlsx'

# Precipitation conversion (Supp. Table S2)
PRECIP_A_H = 0.0102;  PRECIP_B_H = -0.0151
PRECIP_A_Y = 0.4393;  PRECIP_B_Y =  0.1074
PRECIP_T_MEAN = 16.0

# ══════════════════════════════════════════════════════════════════════
# STYLE — ngeo_style palette
# ══════════════════════════════════════════════════════════════════════
from ngeo_style import (apply_style, style_ax, make_heatmap_cmap,
                        EVENTS_FIG5 as EVENTS, DOUBLE_COL, MM,
                        COL_MEDIAN, COL_IQR, COL_D18O,
                        COL_AP, COL_BASEFLOW, COL_INSOL,
                        COL_BG_500, COL_BG_600, COL_BG_700)
apply_style()

# EVENTS_FIG5 is defined in yr BP in ngeo_style; rescale to ka for this figure.
EVENTS = {k: {**v, 'c': v['c']/1000.0, 'hw': v.get('hw',200)/1000.0}
          for k, v in EVENTS.items()}

# ── Styling constants ─────────────────────────────────────────────
# AP envelope (age-propagated IQR)
AP_ALPHA    = 0.25
AP_EDGE_LW  = 0.6
AP_EDGE_A   = 0.45

# LR (1 cm) — primary smooth signal (original styling)
LR_MEDIAN_LW = 0.5
LR_IQR_LW   = 0.5
LR_IQR_LS   = '-'
LR_IQR_A    = 1.0

# HR (native-res) — high-frequency detail, faint background
HR_MEDIAN_LW = 0.3
HR_MEDIAN_A  = 0.25
HR_IQR_LW   = 0.25
HR_IQR_LS   = '--'
HR_IQR_A    = 0.2

# Event bands
EVENT_ALPHA  = 0.18
EVENT_ZORDER = 2        # above heatmap (1), below data (3+)

# ══════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════
def to_ka(age_values):
    """Coerce an age series to ka BP.

    The x-axis is ka BP (shared with Figs 6 & 7) and depth2age() already returns
    ka. But several inputs carry their OWN age column read straight from disk
    (d18O_realisations.csv is yr BP: -51..9499), and those do not pass through
    depth2age. Interpolating a ka grid onto a yr-BP axis silently clamps every
    query into the top ~1 cm of the core and yields a flat band, so the units
    are normalised here at every such read.

    Heuristic: HS4 spans ~9.5 ka / ~9500 yr, so any series reaching beyond 100
    is in years.
    """
    a = np.asarray(age_values, dtype=float)
    return a / 1000.0 if np.nanmax(np.abs(a)) > 100 else a


def load_insolation(age_grid_ka):
    """30°N 21 June insolation on age_grid (ka BP), read from INSOL_FILE.

    Replaces the former synthetic_insolation() Berger & Loutre approximation,
    which returned 476->513 W m-2 over 0-9.5 ka (span 36.6) against the
    tabulated 481->523 (span 42.0) — a mean offset of -6.4 and max 10.2 W m-2.
    The tabulated series is authoritative; the approximation is removed rather
    than kept as a fallback so the two can never diverge.
    """
    d = _sheet('Fig5_insolation', 'ka_BP')
    d.columns = [c.strip() for c in d.columns]
    a_col = [c for c in d.columns if 'ka' in c.lower()][0]
    w_col = [c for c in d.columns if 'wm2' in c.lower().replace('/', '').replace(' ', '')][0]
    d = d[[a_col, w_col]].apply(pd.to_numeric, errors='coerce').dropna().sort_values(a_col)
    return np.interp(age_grid_ka, d[a_col].values, d[w_col].values)


def drip_to_precip(d):
    pt_h = PRECIP_A_H * np.asarray(d, float) + PRECIP_B_H
    pt_y = PRECIP_A_Y * pt_h + PRECIP_B_Y
    return np.maximum(0, pt_y * PRECIP_T_MEAN * 365.25)


def precip_to_drip(p):
    pt_y = np.asarray(p, float) / (PRECIP_T_MEAN * 365.25)
    pt_h = (pt_y - PRECIP_B_Y) / PRECIP_A_Y
    return (pt_h - PRECIP_B_H) / PRECIP_A_H


def draw_event_bands(axes, events):
    """Semi-transparent vertical bands for key events."""
    for ev in events.values():
        lo = ev['c'] - ev.get('hw', 0.2)
        hi = ev['c'] + ev.get('hw', 0.2)
        for ax in axes:
            ax.axvspan(lo, hi, color=ev['col'], alpha=EVENT_ALPHA,
                       zorder=EVENT_ZORDER, lw=0)


def add_event_labels(ax, events):
    """Two-tier labelling: age labels inside bands + descriptive arrows."""
    y_top = ax.get_ylim()[1]

    # Tier 1: age labels inside / near bands
    ax.text(events['8.2 ka']['c'], y_top * 0.96, '8.2 ka',
            ha='center', va='top', fontsize=5.5, color='#333',
            fontstyle='italic', zorder=8,
            bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.8))
    ax.text(events['5.2 ka']['c'], y_top * 0.96, '5.2 ka',
            ha='center', va='top', fontsize=5.5, color='#333',
            fontstyle='italic', zorder=8,
            bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.8))
    ax.text(0.25, y_top * 0.96, 'Post-1750',
            ha='right', va='top', fontsize=5.5, color='#333',
            fontstyle='italic', zorder=8,
            bbox=dict(boxstyle='round,pad=0.12', fc='white', ec='none', alpha=0.8))

    # Tier 2: descriptive labels with arrows
    desc = dict(fontsize=5.5, color='#444', ha='center', va='bottom', zorder=8,
                bbox=dict(boxstyle='round,pad=0.12', fc='white',
                          ec='#AAA', lw=0.4, alpha=0.9))
    arrow = dict(arrowstyle='->', color='#555', lw=0.6, shrinkA=2, shrinkB=2)

    ax.annotate('Laurentide ice\nsheet collapse',
                xy=(events['8.2 ka']['c'], y_top * 0.70),
                xytext=(6.8, y_top * 0.90),
                arrowprops=arrow, **desc)
    ax.annotate('End of\nGreen Sahara',
    #            xy=(events['5.2 ka']['c'], y_top * 0.35),
		xy=(5.2, 26),
                xytext=(4.0, y_top * 0.78),
                arrowprops=arrow, **desc)
    ax.annotate('Little\nIce Age',
                xy=(0.3, 20),
                xytext=(0.9, y_top * 0.75),
                arrowprops=arrow, **desc)


def load_d18o_ap_iqr(age_grid):
    """
    Load age-propagated δ¹⁸O IQR from realisations or pre-computed file.
    Returns (pc25, pc75) arrays on age_grid, or (None, None).
    """
    if D18O_AGEPROP_FILE is not None and os.path.exists(D18O_AGEPROP_FILE):
        try:
            df = pd.read_excel(D18O_AGEPROP_FILE)
            df.columns = [c.strip().lower().replace(' ', '_') for c in df.columns]
            a_col = [c for c in df.columns if 'age' in c][0]
            p25 = [c for c in df.columns if 'pc25' in c or '25' in c][0]
            p75 = [c for c in df.columns if 'pc75' in c or '75' in c][0]
            _a = to_ka(df[a_col].values)
            f25 = interp1d(_a, df[p25], bounds_error=False, fill_value=np.nan)
            f75 = interp1d(_a, df[p75], bounds_error=False, fill_value=np.nan)
            print(f"  δ¹⁸O AP IQR from {D18O_AGEPROP_FILE}")
            return f25(age_grid), f75(age_grid)
        except Exception as e:
            print(f"  Warning: {D18O_AGEPROP_FILE} — {e}")

    # PRECOMPUTED: the plotted IQR now lives in the workbook. The 1,000-member
    # ensemble (248 MB) is an INPUT to this, not the plotted quantity; it is cited
    # with its sha256 on 06_external_files. precompute_d18o_iqr.py replicates the
    # chain below exactly and is what produced the sheet.
    try:
        _iqr = _sheet('Fig5_d18O_AP', 'age_ka')
        f25 = interp1d(_iqr['age_ka'].values, _iqr['pc25'].values,
                       bounds_error=False, fill_value=np.nan)
        f75 = interp1d(_iqr['age_ka'].values, _iqr['pc75'].values,
                       bounds_error=False, fill_value=np.nan)
        print(f"  d18O AP IQR from {SOURCE}[Fig5_d18O_AP] ({len(_iqr)} pts, precomputed)")
        return f25(age_grid), f75(age_grid)
    except Exception as e:
        print(f"  Warning: Fig5_d18O_AP - {e}")

    if os.path.exists(D18O_REAL_FILE):
        try:
            d18o_real = pd.read_csv(D18O_REAL_FILE).sort_values('age')
            d18o_real['age'] = to_ka(d18o_real['age'].values)   # yr BP -> ka
            real_cols = [c for c in d18o_real.columns if c.startswith('r')]
            mat = np.full((len(age_grid), len(real_cols)), np.nan)
            for ri, col in enumerate(real_cols):
                f = interp1d(d18o_real['age'].values, d18o_real[col].values,
                             kind='linear', bounds_error=False,
                             fill_value=(d18o_real[col].values[0],
                                         d18o_real[col].values[-1]))
                mat[:, ri] = f(age_grid)
            print(f"  δ¹⁸O AP IQR from {D18O_REAL_FILE} ({len(real_cols)} realisations)")
            return np.nanpercentile(mat, 25, axis=1), np.nanpercentile(mat, 75, axis=1)
        except Exception as e:
            print(f"  Warning: {D18O_REAL_FILE} — {e}")

    return None, None



# ══════════════════════════════════════════════════════════════════════
# WORKBOOK READERS
# ══════════════════════════════════════════════════════════════════════
def _find_header(sheet, first_col_name):
    """Locate a sheet's header row by name, not by a hardcoded skiprows.

    Sheets grow provenance notes above their tables as findings accumulate; binding
    to a row number makes every script hostage to that. (Learned the hard way on
    05_monitoring.)
    """
    raw = pd.read_excel(SOURCE, sheet_name=sheet, header=None)
    hits = raw.index[raw[0].astype(str).str.strip() == first_col_name]
    if len(hits) == 0:
        raise ValueError(f"header '{first_col_name}' not found in sheet '{sheet}'")
    return int(hits[0])


def _sheet(name, first_col):
    return pd.read_excel(SOURCE, sheet_name=name,
                         skiprows=_find_header(name, first_col)).dropna(how='all')


# ══════════════════════════════════════════════════════════════════════
# LOAD DATA
# ══════════════════════════════════════════════════════════════════════
print("Loading ...")

# Age model
ad = _sheet('02_chronology', 'depth_cm')
ad.columns = [c.strip() for c in ad.columns]
d_col = [c for c in ad.columns if 'dist' in c.lower() or 'depth' in c.lower()][0]
a_col = [c for c in ad.columns if 'age' in c.lower() and 'error' not in c.lower()
         and 'err2s' not in c.lower()][0]
ad_d = pd.to_numeric(ad[d_col], errors='coerce').values
ad_a = pd.to_numeric(ad[a_col], errors='coerce').values
# 2σ age error, if the chronology carries it (used for the tie-point display only)
e_cols = [c for c in ad.columns if 'error' in c.lower() or 'err2s' in c.lower()]
ad_e = pd.to_numeric(ad[e_cols[0]], errors='coerce').values if e_cols else None
v = np.isfinite(ad_d) & np.isfinite(ad_a)
# NOTE: the age model uses ALL rows including the collection-date anchor
# (err2s ~ 0). Dropping it shifts the core top by ~238 yr.
_depth2age_yr = interp1d(ad_d[v], ad_a[v], kind='linear',
                         bounds_error=False, fill_value='extrapolate')
# Migrated to ka BP (shared x-axis with Figs 6 & 7): convert at source so every
# downstream age column is in ka. Only the hard-coded year literals below change.
def depth2age(d):
    return _depth2age_yr(d) / 1000.0

# PDF heatmap — use HR for maximum detail in pcolormesh.
# Bulk array: stays a file beside the workbook, not a sheet (~118k cells).
if not os.path.exists(HEATMAP_JSON_HR):
    raise FileNotFoundError(
        f"{HEATMAP_JSON_HR} is required by Fig 5a (the drip-rate PDF heatmap). "
        f"See external/README.md.")
with open(HEATMAP_JSON_HR) as f:
    hm = json.load(f)
V_span = np.array(hm['V_span'])
hm_depths = np.array(hm['ages'])
V_pdf = np.array(hm['V_pdf']).T
hm_ages = depth2age(hm_depths)
ok = np.isfinite(hm_ages) & (hm_depths > 0)
hm_ages, V_pdf = hm_ages[ok], V_pdf[ok]
for i in range(V_pdf.shape[0]):
    mx = V_pdf[i].max()
    if mx > 0:
        V_pdf[i] /= mx
idx = np.argsort(hm_ages)
hm_ages, V_pdf = hm_ages[idx], V_pdf[idx]

ae = np.zeros(len(hm_ages) + 1)
ae[1:-1] = 0.5 * (hm_ages[:-1] + hm_ages[1:])
ae[0]  = hm_ages[0]  - (hm_ages[1]  - hm_ages[0])  / 2
ae[-1] = hm_ages[-1] + (hm_ages[-1] - hm_ages[-2]) / 2
de = np.zeros(len(V_span) + 1)
de[1:-1] = 0.5 * (V_span[:-1] + V_span[1:])
de[0]  = max(0, V_span[0] - (V_span[1] - V_span[0]) / 2)
de[-1] = V_span[-1] + (V_span[-1] - V_span[-2]) / 2
AGE_E, DRIP_E = np.meshgrid(ae, de)

# Drip rate summaries — HR and LR
if not os.path.exists(DRIP_SUMMARY_HR):
    raise FileNotFoundError(
        f"{DRIP_SUMMARY_HR} is required by Fig 5a (native-resolution median/IQR). "
        f"See external/README.md.")
dr_hr = pd.read_csv(DRIP_SUMMARY_HR)
dr_hr.columns = [c.lower() for c in dr_hr.columns]
dr_hr['age'] = depth2age(dr_hr['depth'].values)
dr_hr = dr_hr[np.isfinite(dr_hr['age']) & (dr_hr['pc50'] > 0)].sort_values('age')
print(f"  Drip rate HR: {len(dr_hr)} pts, native resolution")

dr_lr = _sheet('04_driprate_posterior', 'depth')
dr_lr.columns = [c.lower() for c in dr_lr.columns]
dr_lr['age'] = depth2age(dr_lr['depth'].values)
dr_lr = dr_lr[np.isfinite(dr_lr['age']) & (dr_lr['pc50'] > 0)].sort_values('age')
print(f"  Drip rate LR: {len(dr_lr)} pts, 1 cm interpolated")

# Isotope data
iso_data = _sheet('03_isotopes', 'depth_cm')[['depth_cm', 'd18O_permil']].copy()
iso_data.columns = ['depth', 'd18O']
iso_data = iso_data.apply(pd.to_numeric, errors='coerce').dropna()
iso_data['age'] = depth2age(iso_data['depth'].values)
iso_data = iso_data[np.isfinite(iso_data['age'])].sort_values('age')

# Age-propagated IQR — drip rate
has_ap = False
if True:
    try:
        ap = _sheet('Fig5_driprate_AP', 'age_calBP')
        ap.columns = [c.strip().lower().replace(' ', '_') for c in ap.columns]
        age_col_ap = [c for c in ap.columns if 'age' in c][0]
        ap[age_col_ap] = to_ka(ap[age_col_ap].values)           # yr BP -> ka
        p25_col = [c for c in ap.columns if 'pc25' in c or '25' in c][0]
        p75_col = [c for c in ap.columns if 'pc75' in c or '75' in c][0]
        has_ap = True
        print(f"  Drip rate AP IQR: {len(ap)} pts")
    except:
        pass

# Age-propagated IQR — δ¹⁸O
d18o_ap_ages = iso_data['age'].values
d18o_ap_pc25, d18o_ap_pc75 = load_d18o_ap_iqr(d18o_ap_ages)
has_d18o_ap = d18o_ap_pc25 is not None

age_grid = dr_lr['age'].values
insol = load_insolation(age_grid)
print(f"  Insolation: {INSOL_FILE} -> {insol.min():.1f}-{insol.max():.1f} W m-2 on age_grid")


# ══════════════════════════════════════════════════════════════════════
# FIGURE — two panels (v20 layout)
# ══════════════════════════════════════════════════════════════════════
print("Rendering ...")
fig = plt.figure(figsize=(180 / 25.4, 140 / 25.4))
gs = gridspec.GridSpec(2, 1, figure=fig, hspace=0.12,
                       height_ratios=[4.0, 2.5],
                       left=0.16, right=0.82, top=0.96, bottom=0.08)

ax_hero = fig.add_subplot(gs[0])
ax_d18o = fig.add_subplot(gs[1], sharex=ax_hero)

for ax in (ax_hero, ax_d18o):
    style_ax(ax)

draw_event_bands([ax_hero, ax_d18o], EVENTS)

# ══════════════════════════════════════════════════════════════════════
# PANEL a — drip rate hero
# ══════════════════════════════════════════════════════════════════════
hm_cmap = make_heatmap_cmap()
im = ax_hero.pcolormesh(AGE_E, DRIP_E, V_pdf.T, cmap=hm_cmap,
                        vmin=0, vmax=1, shading='flat',
                        rasterized=True, zorder=1)

# ── HR overlay: high-frequency detail (thin, semi-transparent) ──
ax_hero.plot(dr_hr['age'].values, dr_hr['pc25'].values,
             color=COL_IQR, ls=HR_IQR_LS, lw=HR_IQR_LW, alpha=HR_IQR_A, zorder=3)
ax_hero.plot(dr_hr['age'].values, dr_hr['pc75'].values,
             color=COL_IQR, ls=HR_IQR_LS, lw=HR_IQR_LW, alpha=HR_IQR_A, zorder=3)
ax_hero.plot(dr_hr['age'].values, dr_hr['pc50'].values,
             color=COL_MEDIAN, lw=HR_MEDIAN_LW, alpha=HR_MEDIAN_A, zorder=3)

# ── LR overlay: smooth backbone (primary, solid) ──
ax_hero.plot(dr_lr['age'].values, dr_lr['pc50'].values,
             color=COL_MEDIAN, lw=LR_MEDIAN_LW, zorder=4,
             label='Median (best-estimate)')
ax_hero.plot(dr_lr['age'].values, dr_lr['pc25'].values,
             color=COL_IQR, ls=LR_IQR_LS, lw=LR_IQR_LW, alpha=LR_IQR_A, zorder=4)
ax_hero.plot(dr_lr['age'].values, dr_lr['pc75'].values,
             color=COL_IQR, ls=LR_IQR_LS, lw=LR_IQR_LW, alpha=LR_IQR_A, zorder=4,
             label='IQR (best-estimate)')

# AP IQR envelope
if has_ap:
    ax_hero.fill_between(ap[age_col_ap].values,
                         ap[p25_col].values, ap[p75_col].values,
                         color=COL_AP, alpha=AP_ALPHA, lw=0, zorder=5,
                         label='Age-propagated 25th–75th %ile')
    ax_hero.plot(ap[age_col_ap].values, ap[p25_col].values,
                 color=COL_AP, lw=AP_EDGE_LW, alpha=AP_EDGE_A, zorder=5)
    ax_hero.plot(ap[age_col_ap].values, ap[p75_col].values,
                 color=COL_AP, lw=AP_EDGE_LW, alpha=AP_EDGE_A, zorder=5)

# Baseflow anchor
ax_hero.plot(BASEFLOW_AGE, BASEFLOW_DRIP, '*', color=COL_BASEFLOW,
             markersize=7, zorder=7,
             markeredgecolor='white', markeredgewidth=0.3)

ax_hero.set_ylim(0, 45)
ax_hero.yaxis.set_major_locator(ticker.MultipleLocator(10))
ax_hero.yaxis.set_minor_locator(ticker.MultipleLocator(5))
ax_hero.set_ylabel(r'Drip rate (drips min$^{-1}$)', labelpad=8)

# Precipitation axis (secondary left)
ax_precip = ax_hero.secondary_yaxis('left',
                                     functions=(drip_to_precip, precip_to_drip))
ax_precip.set_ylabel(r'$\sim$P (mm yr$^{-1}$)', fontsize=5, color='#777',
                     labelpad=1, fontstyle='italic')
ax_precip.tick_params(axis='y', colors='#777', labelsize=5,
                      width=0.4, length=2, direction='in', pad=1)
ax_precip.spines['left'].set_position(('outward', 30))
ax_precip.spines['left'].set_color('#777')
ax_precip.spines['left'].set_linewidth(0.4)

# Insolation overlay
ax_insol = ax_hero.twinx()
ax_insol.plot(age_grid, insol, color=COL_INSOL, lw=0.8, alpha=0.7,
              zorder=2, label='30°N June insolation')
ax_insol.set_ylabel(r'Insolation (W m$^{-2}$)', color=COL_INSOL,
                    fontsize=6, labelpad=6)
ax_insol.tick_params(axis='y', colors=COL_INSOL, labelsize=5,
                     width=0.4, length=2, direction='in')
ax_insol.spines['right'].set_color(COL_INSOL)
ax_insol.spines['right'].set_linewidth(0.4)
ax_insol.set_ylim(450, 530)

# Combined legend
h1, l1 = ax_hero.get_legend_handles_labels()
h2, l2 = ax_insol.get_legend_handles_labels()
ax_hero.legend(h1 + h2, l1 + l2, loc='lower left', frameon=True,
               framealpha=0.9, fontsize=5.5, ncol=1,
               bbox_to_anchor=(0.01, 0.02))

# Event labels (two-tier)
add_event_labels(ax_hero, EVENTS)

# Panel label
ax_hero.text(0.015, 0.96, 'a', transform=ax_hero.transAxes,
             fontsize=10, fontweight='bold', va='top',
             bbox=dict(boxstyle='round,pad=0.2',
                       fc='white', ec='none', alpha=0.8))
plt.setp(ax_hero.get_xticklabels(), visible=False)

# Colourbar
pos_a = ax_hero.get_position()
cbar_ax = fig.add_axes([pos_a.x1 + 0.08, pos_a.y0, 0.010, pos_a.height])
cbar = fig.colorbar(im, cax=cbar_ax, ticks=[0, 0.5, 1])
cbar.set_label('Normalised PDF', fontsize=5, labelpad=4)
cbar.ax.tick_params(labelsize=4)
cbar.outline.set_linewidth(0.4)

# ══════════════════════════════════════════════════════════════════════
# PANEL b — δ¹⁸O with age-propagated IQR
# ══════════════════════════════════════════════════════════════════════

# AP IQR band
if has_d18o_ap:
    ax_d18o.fill_between(d18o_ap_ages, d18o_ap_pc25, d18o_ap_pc75,
                         color=COL_AP, alpha=AP_ALPHA, lw=0, zorder=2,
                         label='Age-propagated IQR')
    ax_d18o.plot(d18o_ap_ages, d18o_ap_pc25,
                 color=COL_AP, lw=AP_EDGE_LW, alpha=AP_EDGE_A, zorder=2)
    ax_d18o.plot(d18o_ap_ages, d18o_ap_pc75,
                 color=COL_AP, lw=AP_EDGE_LW, alpha=AP_EDGE_A, zorder=2)

# Raw δ¹⁸O
ax_d18o.plot(iso_data['age'].values, iso_data['d18O'].values,
             color=COL_D18O, lw=0.6, zorder=4,
             label='δ¹⁸O (best-estimate)')

ax_d18o.set_ylabel(r'$\delta^{18}$O (‰ VPDB)', labelpad=8)
y_lo = iso_data['d18O'].quantile(0.995) + 0.3
y_hi = iso_data['d18O'].quantile(0.005) - 0.3
ax_d18o.set_ylim(y_lo, y_hi)

ax_d18o.legend(loc='upper right', frameon=True, framealpha=0.9, fontsize=5.5)

ax_d18o.text(0.015, 0.94, 'b', transform=ax_d18o.transAxes,
             fontsize=10, fontweight='bold', va='top',
             bbox=dict(boxstyle='round,pad=0.2',
                       fc='white', ec='none', alpha=0.8))

# ══════════════════════════════════════════════════════════════════════
# SHARED
# ══════════════════════════════════════════════════════════════════════
ax_hero.set_xlim(hm_ages.max(), hm_ages.min())
ax_d18o.set_xlabel('Age (ka BP)', fontsize=7, fontweight='bold')
ax_d18o.xaxis.set_major_locator(ticker.MaxNLocator(12))
ax_d18o.xaxis.set_minor_locator(ticker.AutoMinorLocator(2))

# ══════════════════════════════════════════════════════════════════════
# ²³⁰Th TIE-POINTS — above the panel-a top spine.
# Ported verbatim from Fig. 7 (build_f1_events.py): same marker (v, s=5),
# same greys, same ±2σ errorbars, same axes-fraction offsets, same label.
# NOTE: Fig. 7 works in ka and filters/plots on uth.ka; this figure's
# x-axis is yr BP, so ages are used unconverted (age_bp) and xerr is err2s
# in years. The 0–9500 filter is Fig. 7's 0–9.5 ka filter in yr BP.
# ══════════════════════════════════════════════════════════════════════
if ad_e is not None:
    uth = pd.DataFrame({'age_bp': ad_a, 'err2s': ad_e}).dropna()
    uth['ka']     = uth['age_bp'] / 1000.0
    uth['err_ka'] = uth['err2s']  / 1000.0
    # Fig. 7's filter. err2s > 1 drops the collection-date anchor, which
    # constrains the age model but is not a ²³⁰Th date.
    uth = uth[(uth['ka'] >= 0) & (uth['ka'] <= 9.5) & (uth['err2s'] > 1)]

    trans    = blended_transform_factory(ax_hero.transData, ax_hero.transAxes)
    MARKER_Y = 1.04
    LABEL_Y  = 1.09
    ax_hero.scatter(uth['ka'], np.full(len(uth), MARKER_Y), marker='v', s=5,
                    color=COL_BG_700, zorder=8, clip_on=False, transform=trans)
    ax_hero.errorbar(uth['ka'], np.full(len(uth), MARKER_Y), xerr=uth['err_ka'],
                     fmt='none', ecolor=COL_BG_500, elinewidth=0.4,
                     capsize=1.2, capthick=0.4, zorder=7, clip_on=False,
                     transform=trans)
    # Fig. 7 places the label 0.15 ka inside the left (older) edge; same inset here.
    ax_hero.text(ax_hero.get_xlim()[0] - 0.15, LABEL_Y, '\u00b2\u00b3\u2070Th (\u00b12\u03c3)',
                 fontsize=4.8, color=COL_BG_600, va='center', ha='left',
                 clip_on=False, transform=trans)
    print(f"  \u00b2\u00b3\u2070Th tie-points: {len(uth)} plotted "
          f"(2\u03c3 {uth.err2s.min():.0f}\u2013{uth.err2s.max():.0f} yr)")
else:
    print("  WARNING: no age-error column in AGE_DEPTH \u2014 tie-points skipped")

# ══════════════════════════════════════════════════════════════════════
# SAVE
# ══════════════════════════════════════════════════════════════════════
for fmt in ('pdf', 'png', 'eps'):
    fig.savefig(f'{OUT_STEM}.{fmt}', dpi=600, bbox_inches='tight')
    print(f"  Saved → {OUT_STEM}.{fmt}")
plt.close(fig)
print("Done.")
