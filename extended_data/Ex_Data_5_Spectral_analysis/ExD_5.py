#!/usr/bin/env python3
"""
ED_Fig_spectral_analysis.py
============================
Extended Data Figure: Spectral analysis of drip rate and δ¹⁸O variability.

Three-panel figure:
  (a) Lomb-Scargle periodogram of detrended drip rate (n ≈ 585)
  (b) Lomb-Scargle periodogram of detrended δ¹⁸O (n ≈ 1,223)
  (c) Cross-spectral coherence + phase (Welch, 20-yr interpolated grid)

Key result: the two proxies share NO significant periodicities at the 5%
false-alarm probability level, confirming spectral independence.

Inputs:
  - drip_rate_summary.csv    (depth-domain BayProX output: depth, pc50)
  - HS4_age_depth.csv        (31 U-Th tie points for depth → age)
  - Drip_rate.xlsx            (sheet '3.Isotopes': depth, δ¹⁸O)

Outputs:
  ED_Fig_spectral_analysis.png (300 dpi)
  ED_Fig_spectral_analysis.pdf (vector)

Requirements:
  pip install astropy scipy matplotlib pandas openpyxl

Usage:
  python ED_Fig_spectral_analysis.py
  python ED_Fig_spectral_analysis.py --summary path/to/drip_rate_summary.csv \
                                      --age_depth path/to/HS4_age_depth.csv \
                                      --isotopes path/to/Drip_rate.xlsx \
                                      --output ED_Fig_spectral_analysis
"""

import argparse
def _panel_label(ax, letter, x=0.02, y=0.96):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight='bold', va='top', ha='left', zorder=20)

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', 'manuscript_figures'))
from ngeo_style import (apply_style, style_ax, panel_label as _ngeo_panel_label, DOUBLE_COL,
                        COL_BG_900, COL_BG_600, COL_BG_400, COL_BG_300,
                        COL_TEAL_800, COL_TEAL_600, COL_NI, COL_CO,
                        COL_BROWN_500, COL_DORANGE, COL_RED_900)
apply_style()
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
from scipy.signal import csd, welch, detrend, find_peaks as fp
from astropy.timeseries import LombScargle
import warnings
warnings.filterwarnings('ignore')

# ── Arguments ────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description='ED spectral analysis figure')
parser.add_argument('--summary', default='drip_rate_summary.csv',
                    help='Path to drip_rate_summary.csv')
parser.add_argument('--age_depth', default='HS4_age_depth.csv',
                    help='Path to HS4_age_depth.csv (31 tie points)')
parser.add_argument('--isotopes', default='Drip_rate.xlsx',
                    help='Path to Drip_rate.xlsx (sheet 3.Isotopes)')
parser.add_argument('--iso_sheet', default='3.Isotopes',
                    help='Sheet name for δ¹⁸O data')
parser.add_argument('--output', default='ED_Fig_spectral_analysis',
                    help='Output filename stem (no extension)')
args = parser.parse_args()

# Nature Communications house style comes from ngeo_style.apply_style() (above).

# ── Load age model ───────────────────────────────────────────────────
ad = pd.read_csv(args.age_depth)
ad.columns = [c.strip() for c in ad.columns]
d_col = [c for c in ad.columns if 'dist' in c.lower() or 'depth' in c.lower()][0]
a_col = [c for c in ad.columns if 'age' in c.lower() and 'error' not in c.lower()][0]
ad_d = pd.to_numeric(ad[d_col], errors='coerce').values
ad_a = pd.to_numeric(ad[a_col], errors='coerce').values
v = np.isfinite(ad_d) & np.isfinite(ad_a)
depth2age = interp1d(ad_d[v], ad_a[v], kind='linear',
                     bounds_error=False, fill_value='extrapolate')

# ── Load drip rate ───────────────────────────────────────────────────
dr = pd.read_csv(args.summary)
dr.columns = [c.lower() for c in dr.columns]
dr['age'] = depth2age(dr['depth'].values)
dr = dr[np.isfinite(dr['age']) & (dr['pc50'] > 0)].sort_values('age')
print(f"DR: {len(dr)} pts, median spacing {dr.age.diff().median():.1f} yr")

# ── Load δ¹⁸O ───────────────────────────────────────────────────────
iso_raw = pd.read_excel(args.isotopes, sheet_name=args.iso_sheet, header=None)
iso_data = iso_raw.iloc[13:, :2].copy()
iso_data.columns = ['depth', 'd18O']
iso_data = iso_data.apply(pd.to_numeric, errors='coerce').dropna()
iso_data['age'] = depth2age(iso_data['depth'].values)
iso_data = iso_data[np.isfinite(iso_data['age'])].sort_values('age')
print(f"δ¹⁸O: {len(iso_data)} pts, median spacing {iso_data.age.diff().median():.1f} yr")

# ── Detrend ──────────────────────────────────────────────────────────
dr_ages = dr['age'].values
dr_vals = detrend(dr['pc50'].values)
iso_ages = iso_data['age'].values
iso_vals = detrend(iso_data['d18O'].values)

# ── Lomb-Scargle periodograms ────────────────────────────────────────
min_freq, max_freq = 1 / 5000, 1 / 40
freqs = np.linspace(min_freq, max_freq, 2000)
periods = 1 / freqs

print("Computing Lomb-Scargle periodograms ...")

ls_dr = LombScargle(dr_ages, dr_vals)
power_dr = ls_dr.power(freqs)
fap_dr = ls_dr.false_alarm_level([0.01, 0.05, 0.10])
print(f"  DR FAP levels (1%, 5%, 10%): {fap_dr}")

ls_iso = LombScargle(iso_ages, iso_vals)
power_iso = ls_iso.power(freqs)
fap_iso = ls_iso.false_alarm_level([0.01, 0.05, 0.10])
print(f"  δ¹⁸O FAP levels (1%, 5%, 10%): {fap_iso}")

# ── Find significant peaks ──────────────────────────────────────────
def find_pks(freqs, power, fap5):
    peaks, _ = fp(power, height=fap5, distance=20)
    if len(peaks) == 0:
        return []
    idx = np.argsort(power[peaks])[::-1]
    return [{'period': 1 / freqs[p], 'power': power[p]} for p in peaks[idx][:10]]

dr_peaks = find_pks(freqs, power_dr, fap_dr[1])
iso_peaks = find_pks(freqs, power_iso, fap_iso[1])

print(f"\nDR significant peaks: {['{:.0f} yr'.format(p['period']) for p in dr_peaks]}")
print(f"δ¹⁸O significant peaks: {['{:.0f} yr'.format(p['period']) for p in iso_peaks]}")

# ── Known climate periodicities ──────────────────────────────────────
KNOWN = {
    'PDO\n(~60 yr)': 60,
    'AMO\n(~65 yr)': 65,
    'Gleissberg\n(~88 yr)': 88,
    'de Vries\n(~210 yr)': 210,
    'Eddy\n(~1000 yr)': 1000,
    'Bond\n(~1500 yr)': 1500,
}

# ── Cross-spectral coherence ────────────────────────────────────────
print("Computing cross-spectral coherence ...")
reg_step = 20  # yr
age_min = max(dr_ages.min(), iso_ages.min())
age_max = min(dr_ages.max(), iso_ages.max())
reg_ages = np.arange(np.ceil(age_min), np.floor(age_max), reg_step)

f_dr_i = interp1d(dr_ages, dr_vals, bounds_error=False, fill_value=np.nan)
f_iso_i = interp1d(iso_ages, iso_vals, bounds_error=False, fill_value=np.nan)
dr_reg = f_dr_i(reg_ages)
iso_reg = f_iso_i(reg_ages)
valid = np.isfinite(dr_reg) & np.isfinite(iso_reg)
dr_reg = dr_reg[valid]
iso_reg = iso_reg[valid]
reg_ages = reg_ages[valid]

nperseg = min(256, len(dr_reg) // 2)
f_w, Pxy = csd(dr_reg, iso_reg, fs=1 / reg_step, nperseg=nperseg,
               noverlap=nperseg // 2)
f_w, Pxx = welch(dr_reg, fs=1 / reg_step, nperseg=nperseg,
                 noverlap=nperseg // 2)
f_w, Pyy = welch(iso_reg, fs=1 / reg_step, nperseg=nperseg,
                 noverlap=nperseg // 2)
coherence = np.abs(Pxy) ** 2 / (Pxx * Pyy + 1e-30)
phase = np.angle(Pxy, deg=True)
per_w = 1 / f_w[1:]  # skip DC
sig_level = 1 / np.sqrt(nperseg)

# ── Figure ───────────────────────────────────────────────────────────
fig = plt.figure(figsize=(DOUBLE_COL, DOUBLE_COL * 0.85))
gs = gridspec.GridSpec(3, 1, figure=fig, hspace=0.32,
                       left=0.09, right=0.83, top=0.97, bottom=0.07)

# ── Panels A & B: Periodograms ───────────────────────────────────────
for panel_idx, (panel_label, power, fap, peaks, color, title) in enumerate([
    ('a', power_dr, fap_dr, dr_peaks, COL_NI, 'Drip rate periodogram (detrended)'),
    ('b', power_iso, fap_iso, iso_peaks, COL_BROWN_500, 'δ¹⁸O periodogram (detrended)'),
]):
    ax = fig.add_subplot(gs[panel_idx])
    ax.plot(periods, power, color=color, lw=0.6, zorder=3)
    ax.axhline(fap[0], color=COL_RED_900, ls='--', lw=0.4, label='1% FAP')
    ax.axhline(fap[1], color=COL_DORANGE, ls='--', lw=0.4, label='5% FAP')
    ax.axhline(fap[2], color=COL_BG_600, ls=':', lw=0.3, label='10% FAP')

    # Known periodicities
    for name, per in KNOWN.items():
        if 40 < per < 5000:
            ax.axvline(per, color=COL_BG_300, ls=':', lw=0.3)
            ax.text(per, ax.get_ylim()[1] * 0.95 if ax.get_ylim()[1] > 0 else 0.05,
                    name, fontsize=4.5, rotation=90, va='top', ha='right', color=COL_BG_600)

    # Annotate significant peaks
    for p in peaks[:5]:
        ax.annotate(f'{p["period"]:.0f} yr', xy=(p['period'], p['power']),
                    xytext=(5, 5), textcoords='offset points', fontsize=5,
                    arrowprops=dict(arrowstyle='->', lw=0.3, color=COL_BG_900),
                    color=COL_BG_900)

    ax.set_xscale('log')
    ax.set_xlabel('Period (yr)')
    ax.set_ylabel('Lomb-Scargle power')
    _ngeo_panel_label(ax, panel_label)
    style_ax(ax)
    ax.legend(loc='center left', bbox_to_anchor=(1.005, 0.5), frameon=False, fontsize=5, handlelength=1.6, borderaxespad=0.0)
    ax.set_xlim(40, 5000)

# ── Panel C: Cross-spectral coherence ────────────────────────────────
ax = fig.add_subplot(gs[2])
ax2 = ax.twinx()

ax.plot(per_w, coherence[1:], color=COL_NI, lw=0.7, zorder=3, label='Coherence')
ax.axhline(sig_level, color=COL_BG_600, ls=':', lw=0.3,
           label=f'~95% significance ({sig_level:.2f})')
ax.set_ylabel('Coherence²', color=COL_NI)
ax.tick_params(axis='y', colors='#2C7BB6')
ax.set_ylim(0, 1)

ax2.scatter(per_w, phase[1:], s=2, color=COL_BROWN_500, alpha=0.5, zorder=2,
            label='Phase (°)')
ax2.set_ylabel('Phase (°)', color=COL_BROWN_500)
ax2.tick_params(axis='y', colors='#D7191C')
ax2.set_ylim(-180, 180)

for name, per in KNOWN.items():
    if 40 < per < 5000:
        ax.axvline(per, color=COL_BG_300, ls=':', lw=0.3)

ax.set_xscale('log')
ax.set_xlabel('Period (yr)')
_ngeo_panel_label(ax, 'c')
style_ax(ax); style_ax(ax2)
h1, l1 = ax.get_legend_handles_labels()
h2, l2 = ax2.get_legend_handles_labels()
ax.legend(h1 + h2, l1 + l2, loc='center left', bbox_to_anchor=(1.06, 0.5), frameon=False, fontsize=5, handlelength=1.6, borderaxespad=0.0)
ax.set_xlim(40, 5000)

# ── Save ─────────────────────────────────────────────────────────────
for fmt in ('pdf', 'png', 'eps'):
    fig.savefig(f'{args.output}.{fmt}', dpi=600, bbox_inches='tight')
plt.close()
print(f'\nSaved {args.output}.png (300 dpi) and .pdf')

# ── Summary ──────────────────────────────────────────────────────────
print("\n=== SUMMARY ===")
dr_per = set([round(p['period'] / 10) * 10 for p in dr_peaks])
iso_per = set([round(p['period'] / 10) * 10 for p in iso_peaks])
shared = dr_per & iso_per
print(f"Shared significant periods: {sorted(shared) if shared else 'None at 5% level'}")
print(f"DR-only: {sorted(dr_per - iso_per)}")
print(f"δ¹⁸O-only: {sorted(iso_per - dr_per)}")