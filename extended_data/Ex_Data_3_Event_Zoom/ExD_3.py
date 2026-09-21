"""
ED_Fig_event_zooms.py
=====================
Extended Data Figure — event zoom panels for stalagmite HS4.
Three panels side by side:
  a: 8.2 ka event (7,800–8,600 yr BP)
  b: 5.2 ka event (4,800–5,600 yr BP)
  c: Post-1750 CE (1,680–2,010 CE)

Each panel shows the drip rate PDF heatmap (Gaussian kernel from
native-res percentiles), median, IQR fill, and co-located δ¹⁸O
(right axis, inverted). Panel a includes Liu et al. (2013) high-
resolution δ¹⁸O. Panel c includes historical event markers (i 1788 flood, ii Tambora,
iii 1870 flood, iv 1931 flood) and a
post-1850 decline trend.

Colour palette from ngeo_style for manuscript consistency.
"""
import json, os, warnings
import numpy as np, pandas as pd, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from scipy.interpolate import interp1d
from scipy.stats import linregress
warnings.filterwarnings('ignore')

# ══════════════════════════════════════════════════════════════════════
# PATHS
# ══════════════════════════════════════════════════════════════════════
HEATMAP_JSON  = 'pdf_heatmap.json'
DRIP_SUMMARY  = 'drip_rate_summary.csv'
AGE_DEPTH     = 'HS4_age_depth.csv'
DRIP_XLSX     = 'Drip_rate.xlsx'
import os as _os
if not _os.path.exists(DRIP_XLSX):   # the isotope workbook lives at the repository root
    DRIP_XLSX = _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', 'Drip_rate.xlsx')
ISO_SHEET     = '3.Isotopes'
LIU_FILE      = 'liu2013_hs4_full.csv'
BASEFLOW_DRIP = 16.7   # 2004-2023 annual-baseflow mean

# ══════════════════════════════════════════════════════════════════════
# STYLE
# ══════════════════════════════════════════════════════════════════════
from ngeo_style import (apply_style, style_ax, make_heatmap_cmap,
                        EVENTS_FIG5 as EVENTS, MM,
                        COL_MEDIAN, COL_IQR, COL_D18O_C, COL_LIU,
                        COL_TREND, COL_BASEFLOW, COL_MARKER)
apply_style()

D18O_COMMON = (-10.7, -7.2)
COMMON_YLIM = (0, 45)

# ══════════════════════════════════════════════════════════════════════
# LOAD DATA
# ══════════════════════════════════════════════════════════════════════
print("Loading ...")

# Age model
ad = pd.read_csv(AGE_DEPTH, encoding='utf-8-sig')
ad.columns = [c.strip() for c in ad.columns]
d_col = [c for c in ad.columns if 'dist' in c.lower() or 'depth' in c.lower()][0]
a_col = [c for c in ad.columns if 'age' in c.lower() and 'error' not in c.lower()][0]
ad_d = pd.to_numeric(ad[d_col], errors='coerce').values
ad_a = pd.to_numeric(ad[a_col], errors='coerce').values
v = np.isfinite(ad_d) & np.isfinite(ad_a)
depth2age = interp1d(ad_d[v], ad_a[v], kind='linear',
                     bounds_error=False, fill_value='extrapolate')

# Drip rate summary
dr = pd.read_csv(DRIP_SUMMARY)
dr.columns = [c.lower() for c in dr.columns]
dr['age'] = depth2age(dr['depth'].values)
dr['ce'] = 1950 - dr['age']
dr = dr[np.isfinite(dr['age'])].sort_values('age')

# Isotope data
iso_raw = pd.read_excel(DRIP_XLSX, sheet_name=ISO_SHEET, header=None)
iso_data = iso_raw.iloc[13:, :2].copy()
iso_data.columns = ['depth', 'd18O']
iso_data = iso_data.apply(pd.to_numeric, errors='coerce').dropna()
iso_data['age'] = depth2age(iso_data['depth'].values)
iso_data['ce'] = 1950 - iso_data['age']
iso_data = iso_data[np.isfinite(iso_data['age'])].sort_values('age')

# Liu 2013 high-res δ¹⁸O
liu = pd.read_csv(LIU_FILE).sort_values('age_bp')
print(f"  Liu: {len(liu)} pts, {liu.age_bp.min():.0f}–{liu.age_bp.max():.0f} BP")

# Interpolator for c3 markers
f_dr_ce = interp1d(dr['ce'].values, dr['pc50'].values,
                   bounds_error=False, fill_value=np.nan)


# ══════════════════════════════════════════════════════════════════════
# HELPERS
# ══════════════════════════════════════════════════════════════════════
hm_cmap = make_heatmap_cmap()


def make_event_panel(ax, age_lo, age_hi):
    """Draw drip rate PDF heatmap + median + IQR for an age window."""
    m = (dr.age >= age_lo) & (dr.age <= age_hi)
    sub = dr[m].sort_values('age')
    x_vals = sub.age.values
    y_grid = np.linspace(0, COMMON_YLIM[1], 80)
    Z = np.zeros((80, len(x_vals)))
    for j, (_, row) in enumerate(sub.iterrows()):
        mu = row.pc50
        sigma = max((row.pc75 - row.pc25) / 1.35, 0.5)
        pdf = np.exp(-0.5 * ((y_grid - mu) / sigma) ** 2)
        pdf /= (pdf.max() + 1e-10)
        Z[:, j] = pdf
    ax.pcolormesh(x_vals, y_grid, Z, cmap=hm_cmap, vmin=0, vmax=1,
                  shading='nearest', rasterized=True, zorder=2)
    ax.plot(x_vals, sub.pc50.values, color=COL_MEDIAN, lw=0.8, zorder=4)
    ax.fill_between(x_vals, sub.pc25.values, sub.pc75.values,
                    color=COL_IQR, alpha=0.20, lw=0, zorder=3)
    style_ax(ax)


def add_d18o(ax, age_lo, age_hi, show_label=False):
    """Add δ¹⁸O on a twin right axis (inverted)."""
    m_i = (iso_data.age >= age_lo) & (iso_data.age <= age_hi)
    ax2 = ax.twinx()
    ax2.plot(iso_data.loc[m_i, 'age'].values,
             iso_data.loc[m_i, 'd18O'].values,
             color=COL_D18O_C, lw=0.4, alpha=0.5, zorder=1)
    ax2.set_ylim(*D18O_COMMON)
    ax2.invert_yaxis()
    if show_label:
        ax2.set_ylabel(r'$\delta^{18}$O (‰)', fontsize=5.5, color='#78909C')
        ax2.tick_params(axis='y', labelsize=5, colors='#78909C',
                        width=0.3, length=2, direction='in')
    else:
        ax2.set_yticklabels([])
        ax2.tick_params(axis='y', length=0)
    ax2.spines['right'].set_linewidth(0.3)
    return ax2


# ══════════════════════════════════════════════════════════════════════
# FIGURE — three panels side by side
# ══════════════════════════════════════════════════════════════════════
print("Rendering ...")
fig = plt.figure(figsize=(180 / 25.4, 65 / 25.4))
gs = gridspec.GridSpec(1, 3, figure=fig, wspace=0.10,
                       left=0.08, right=0.92, top=0.92, bottom=0.18)

# ─────────────────────────────────────────
# Panel a: 8.2 ka event
# ─────────────────────────────────────────
ax_a = fig.add_subplot(gs[0, 0])
make_event_panel(ax_a, 7800, 8600)
ax_a.axvspan(8100, 8250, color=EVENTS['8.2 ka']['col'],
             alpha=0.15, zorder=1.5, lw=0)
ax_a_r = add_d18o(ax_a, 7800, 8600)

# Liu et al. high-res δ¹⁸O
m_liu = (liu.age_bp >= 7800) & (liu.age_bp <= 8600)
ax_a_r.plot(liu.loc[m_liu, 'age_bp'].values,
            liu.loc[m_liu, 'd18O'].values,
            color=COL_LIU, lw=0.45, alpha=0.7, zorder=1.5)

ax_a.set_xlim(8600, 7800)
ax_a.set_ylim(*COMMON_YLIM)
ax_a.set_ylabel(r'Drip rate (drips min$^{-1}$)', fontsize=6)
ax_a.set_xlabel('Age (yr BP)', fontsize=6)
ax_a.text(0.04, 0.95, 'a', transform=ax_a.transAxes,
          fontsize=9, fontweight='bold', va='top')

# ─────────────────────────────────────────
# Panel b: 5.2 ka event
# ─────────────────────────────────────────
ax_b = fig.add_subplot(gs[0, 1])
make_event_panel(ax_b, 4800, 5600)
ax_b.axvspan(5000, 5300, color=EVENTS['5.2 ka']['col'],
             alpha=0.15, zorder=1.5, lw=0)
add_d18o(ax_b, 4800, 5600)

ax_b.set_xlim(5600, 4800)
ax_b.set_ylim(*COMMON_YLIM)
ax_b.set_yticklabels([])
ax_b.set_xlabel('Age (yr BP)', fontsize=6)
ax_b.text(0.04, 0.95, 'b', transform=ax_b.transAxes,
          fontsize=9, fontweight='bold', va='top')

# ─────────────────────────────────────────
# Panel c: Post-1750 CE
# ─────────────────────────────────────────
ax_c = fig.add_subplot(gs[0, 2])
m_post = dr['ce'] >= 1680
sub_p = dr[m_post].sort_values('ce')
ce_vals = sub_p.ce.values

# Heatmap
y_grid = np.linspace(0, COMMON_YLIM[1], 80)
Z = np.zeros((80, len(ce_vals)))
for j, (_, row) in enumerate(sub_p.iterrows()):
    mu = row.pc50
    sigma = max((row.pc75 - row.pc25) / 1.35, 0.5)
    pdf = np.exp(-0.5 * ((y_grid - mu) / sigma) ** 2)
    pdf /= (pdf.max() + 1e-10)
    Z[:, j] = pdf
ax_c.pcolormesh(ce_vals, y_grid, Z, cmap=hm_cmap, vmin=0, vmax=1,
                shading='nearest', rasterized=True, zorder=2)
ax_c.plot(ce_vals, sub_p.pc50.values, color=COL_MEDIAN, lw=0.8, zorder=4)
ax_c.fill_between(ce_vals, sub_p.pc25.values, sub_p.pc75.values,
                  color=COL_IQR, alpha=0.20, lw=0, zorder=3)
style_ax(ax_c)

# δ¹⁸O on c — with axis label
m_i3 = iso_data.ce >= 1680
ax_c_r = ax_c.twinx()
ax_c_r.plot(iso_data.loc[m_i3, 'ce'].values,
            iso_data.loc[m_i3, 'd18O'].values,
            color=COL_D18O_C, lw=0.4, alpha=0.5, zorder=1)
ax_c_r.set_ylim(*D18O_COMMON)
ax_c_r.invert_yaxis()
ax_c_r.set_ylabel(r'$\delta^{18}$O (‰)', fontsize=5.5, color='#78909C')
ax_c_r.tick_params(axis='y', labelsize=5, colors='#78909C',
                   width=0.3, length=2, direction='in')
ax_c_r.spines['right'].set_linewidth(0.3)

# Baseflow reference + decline trend
ax_c.axhline(BASEFLOW_DRIP, color=COL_BASEFLOW, ls=':', lw=0.4,
             alpha=0.5, zorder=1)
m_dec = sub_p['ce'] >= 1850
if m_dec.sum() > 3:
    sl = linregress(sub_p.loc[m_dec, 'ce'].values,
                    sub_p.loc[m_dec, 'pc50'].values)
    ax_c.plot([1850, 2000],
              [sl.slope * 1850 + sl.intercept,
               sl.slope * 2000 + sl.intercept],
              color=COL_TREND, lw=0.6, ls='--', zorder=5)

# Historical event markers
marker_events = [
    (1788, '^', COL_MARKER),    # i   — 1788 Yangtze flood
    (1815, 'D', '#B71C1C'),     # ii  — Tambora
    (1870, '^', COL_MARKER),    # iii — 1870 Yangtze megaflood (largest in >800 yr; Zhang et al. 2022)
    (1931, '^', COL_MARKER),    # iv  — 1931 Yangtze flood (record Hankou stage; Zhang et al. 2022)
]
numerals = ['i', 'ii', 'iii', 'iv']
for idx_m, (ce_x, mkr, col) in enumerate(marker_events):
    dr_y = float(f_dr_ce(ce_x))
    ax_c.plot(ce_x, dr_y, marker=mkr, color=col, markersize=6,
              markeredgecolor='white', markeredgewidth=0.35, zorder=8)
    offset = 4 if dr_y < 75 else -6
    ax_c.text(ce_x, dr_y + offset, numerals[idx_m], fontsize=4.5,
              ha='center', va='bottom' if offset > 0 else 'top',
              color='#455A64', fontstyle='italic')

ax_c.set_xlim(1680, 2010)
ax_c.set_ylim(0, 90)
ax_c.set_yticks([0, 20, 40, 60, 80]); ax_c.tick_params(axis="y", labelsize=6)
ax_c.set_xlabel('Year CE', fontsize=6)
ax_c.text(0.04, 0.95, 'c', transform=ax_c.transAxes,
          fontsize=9, fontweight='bold', va='top')

# ══════════════════════════════════════════════════════════════════════
# SAVE
# ══════════════════════════════════════════════════════════════════════
for fmt in ('pdf', 'png', 'eps'):
    fig.savefig(f'ED_Fig_event_zooms.{fmt}', dpi=600, bbox_inches='tight')
    print(f"  Saved → ED_Fig_event_zooms.{fmt}")
plt.close(fig)
print("Done.")
