"""
drip_rate_82ka_fine_analysis.py
================================
Fine-grained statistical analysis of the 8.2 ka event in the HS4
drip rate reconstruction.

Uses the depth-domain best-estimate reconstruction (drip_rate_summary.csv)
mapped to age via the HS4 chronology (HS4_age_depth.csv / age_model.csv),
consistent with all other figures in the submission.

Tests:
  1. Tight event windows (matched to Owen et al. 2016 ~160 yr duration)
     against an early-Holocene baseline, using detrended pc50 values.
  2. Progressive 100-yr windows at 50-yr steps across 7,800–9,200 yr BP.
  3. Multi-scale rolling MWU analysis (50, 100, 150, 200 yr windows).
  4. Local flanking-baseline comparison — event core vs immediate
     pre- and post-event neighbours.

All tests are detrended (linear orbital trend removed) and use Mann–Whitney
U with rank-biserial effect size.

Inputs:
  - drip_rate_summary.csv    (depth-domain BayProX output: depth, pc25, pc50, pc75)
  - HS4_age_depth.csv        (31 U-Th tie points) OR age_model.csv (interpolated)

Outputs:
  Console statistics
  Fig_82ka_rolling_analysis.pdf / .png

Usage:
  python drip_rate_82ka_fine_analysis.py
  python drip_rate_82ka_fine_analysis.py --summary drip_rate_summary.csv \
                                          --age_model age_model.csv
"""

import argparse
import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ks_2samp, linregress
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', 'manuscript_figures'))
from ngeo_style import (apply_style, style_ax, panel_label, DOUBLE_COL,
                        COL_BG_900, COL_BG_600, COL_BG_400, COL_BG_300, COL_BG_200,
                        COL_TEAL_800, COL_TEAL_600, COL_TEAL_200, COL_BG_700,
                        COL_NI, COL_CO, COL_BROWN_500, COL_DORANGE, COL_RED_900)
apply_style()
import matplotlib.ticker as ticker
import warnings
warnings.filterwarnings('ignore')

# ── Nature Geoscience style ──────────────────────────────────────────
plt.style.use('default')
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Helvetica', 'Arial', 'DejaVu Sans'],
    'font.size': 7,
    'pdf.fonttype': 42,
    'svg.fonttype': 'none',
    'axes.linewidth': 0.5,
    'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
    'xtick.major.size': 3, 'ytick.major.size': 3,
    'xtick.minor.size': 1.5, 'ytick.minor.size': 1.5,
    'xtick.direction': 'in', 'ytick.direction': 'in',
    'axes.labelsize': 7, 'axes.labelweight': 'bold',
    'xtick.labelsize': 6, 'ytick.labelsize': 6,
    'legend.fontsize': 5.5,
})

# ── Arguments ────────────────────────────────────────────────────────
parser = argparse.ArgumentParser(description='8.2 ka fine-grained event analysis')
parser.add_argument('--summary', default='drip_rate_summary.csv',
                    help='Path to drip_rate_summary.csv (depth-domain)')
parser.add_argument('--age_model', default='age_model.csv',
                    help='Path to age_model.csv (interpolated depth → age)')
parser.add_argument('--tiepoints', default='HS4_age_depth.csv',
                    help='Path to HS4_age_depth.csv (fallback if age_model absent)')
parser.add_argument('--output', default='Fig_82ka_rolling_analysis',
                    help='Output filename stem')
args = parser.parse_args()

# ═══════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════

# Analysis range (yr BP)
RANGE_LO, RANGE_HI = 7500, 9500

# Pre-event baseline for fixed-reference rolling test
BASELINE_LO, BASELINE_HI = 8800, 9300

# Tight event windows (matched to Owen et al. 2016 ~160 yr duration)
TIGHT_WINDOWS = {
    'Far pre-event':   (9000, 9300),
    'Near pre-event':  (8280, 8400),
    'Event onset':     (8200, 8280),
    'Event core':      (8100, 8250),
    'Event peak':      (8160, 8220),
    'Recovery':        (7950, 8100),
    'Post-recovery':   (7800, 7950),
}

# Rolling window sizes (yr)
ROLLING_WINDOWS = [50, 100, 150, 200]

# Sub-windows for progressive resolution test (50-yr steps)
FINE_STEPS = np.arange(7800, 9200, 50)

# ═══════════════════════════════════════════════════════════════════════
# LOAD DATA
# ═══════════════════════════════════════════════════════════════════════

print("=" * 68)
print("8.2 ka FINE-GRAINED EVENT ANALYSIS")
print("  Data: drip_rate_summary.csv on HS4 chronology")
print("=" * 68)

# ── Load drip rate summary ───────────────────────────────────────────
df = pd.read_csv(args.summary)
df.columns = [c.strip().lower() for c in df.columns]

# Find columns
def find_col(df, candidates):
    for c in candidates:
        if c in df.columns:
            return c
    return None

col_depth = find_col(df, ['depth', 'depth_mm', 'depth (mm)'])
col_med = find_col(df, ['pc50', 'p50', 'median', 'drip_rate_median', 'q50'])
col_p25 = find_col(df, ['pc25', 'p25', 'q25'])
col_p75 = find_col(df, ['pc75', 'p75', 'q75'])
col_age = find_col(df, ['age', 'age_bp', 'age_ybp'])

if col_med is None:
    raise KeyError(f"Cannot find median column. Available: {list(df.columns)}")

# ── Resolve ages ─────────────────────────────────────────────────────
if col_age is not None:
    print(f"  Using age column '{col_age}' from summary")
    df['age'] = df[col_age]
else:
    import os
    # Try age_model.csv first, then tiepoints
    am_path = args.age_model if os.path.isfile(args.age_model) else None
    tp_path = args.tiepoints if os.path.isfile(args.tiepoints) else None

    if am_path:
        am = pd.read_csv(am_path)
        am_d = am.iloc[:, 0].values
        am_a = am.iloc[:, 1].values
        print(f"  Loaded age model: {am_path} ({len(am)} pts)")
    elif tp_path:
        tp = pd.read_csv(tp_path)
        tp.columns = [c.strip() for c in tp.columns]
        d_col = [c for c in tp.columns if 'dist' in c.lower() or 'depth' in c.lower()][0]
        a_col = [c for c in tp.columns if 'age' in c.lower() and 'error' not in c.lower()][0]
        am_d = pd.to_numeric(tp[d_col], errors='coerce').values
        am_a = pd.to_numeric(tp[a_col], errors='coerce').values
        v = np.isfinite(am_d) & np.isfinite(am_a)
        am_d, am_a = am_d[v], am_a[v]
        print(f"  Loaded tie points: {tp_path} ({len(am_d)} pts)")
    else:
        raise FileNotFoundError(
            "No age model found. Provide --age_model or --tiepoints.")

    depth2age = interp1d(am_d, am_a, kind='linear',
                         bounds_error=False, fill_value='extrapolate')
    df['age'] = depth2age(df[col_depth].values)

# ── Clean and sort ───────────────────────────────────────────────────
df = df[np.isfinite(df['age']) & np.isfinite(df[col_med]) & (df[col_med] > 0)]
df = df.sort_values('age').reset_index(drop=True)

rec_min, rec_max = df.age.min(), df.age.max()
print(f"\n  Record: {rec_min:.0f} to {rec_max:.0f} yr BP  ({len(df)} pts)")

# Subset to analysis range
df_sub = df[(df.age >= RANGE_LO) & (df.age <= RANGE_HI)].copy()
ages_sorted = np.sort(df_sub.age.values)
med_spacing = np.median(np.diff(ages_sorted))
print(f"  Median sample spacing: {med_spacing:.1f} yr")
print(f"  Points in analysis range: {len(df_sub)}")


# ═══════════════════════════════════════════════════════════════════════
# DETRENDING
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "-" * 68)
print("DETRENDING (orbital removal)")
print("-" * 68)

# Full-record linear trend
sl_full = linregress(df.age.values, df[col_med].values)
print(f"\n  Full-record linear trend:")
print(f"    Slope = {sl_full.slope*1000:.3f} drips min⁻¹ kyr⁻¹")
print(f"    r² = {sl_full.rvalue**2:.3f}")

# Detrend within analysis range
grand_mean = df[col_med].mean()
trend_sub = sl_full.slope * df_sub.age.values + sl_full.intercept
df_sub['med_detrend'] = df_sub[col_med].values - trend_sub + grand_mean

print(f"  Detrended {len(df_sub)} points (re-centred on grand mean {grand_mean:.2f})")


# ═══════════════════════════════════════════════════════════════════════
# HELPERS
# ═══════════════════════════════════════════════════════════════════════

def window_vals(lo, hi, detrended=False):
    """Best-estimate pc50 values in [lo, hi] yr BP."""
    mask = (df_sub.age >= lo) & (df_sub.age <= hi)
    col = 'med_detrend' if detrended else col_med
    return df_sub.loc[mask, col].values

def rank_biserial(u, n1, n2):
    return 1.0 - (2.0 * u) / (n1 * n2)

def two_sample_stats(vals_a, vals_b):
    """Returns dict of MWU p, KS p, rank-biserial, delta-median, medians."""
    if len(vals_a) < 2 or len(vals_b) < 2:
        return None
    u, p_mw = mannwhitneyu(vals_a, vals_b, alternative='two-sided')
    d_ks, p_ks = ks_2samp(vals_a, vals_b)
    rb = rank_biserial(u, len(vals_a), len(vals_b))
    return {
        'p_mw': p_mw, 'p_ks': p_ks, 'rb': rb,
        'med_a': np.median(vals_a), 'med_b': np.median(vals_b),
        'dm': np.median(vals_b) - np.median(vals_a),
        'n_a': len(vals_a), 'n_b': len(vals_b),
    }


# ═══════════════════════════════════════════════════════════════════════
# 1. TIGHT EVENT WINDOWS
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 68)
print("1. TIGHT EVENT WINDOWS (raw and detrended)")
print("=" * 68)

bl_raw = window_vals(BASELINE_LO, BASELINE_HI, detrended=False)
bl_dt  = window_vals(BASELINE_LO, BASELINE_HI, detrended=True)

print(f"\n  Baseline ({BASELINE_LO}–{BASELINE_HI} yr BP):")
print(f"    Raw:       median = {np.median(bl_raw):.2f} drips min⁻¹  (n = {len(bl_raw)})")
print(f"    Detrended: median = {np.median(bl_dt):.2f} drips min⁻¹")

print(f"\n  {'Window':<22} {'Ages':>14}  {'n':>4}  "
      f"{'Med(raw)':>9} {'Med(dt)':>9}  "
      f"{'MWU p(raw)':>10} {'MWU p(dt)':>10}  "
      f"{'r_RB(raw)':>9} {'r_RB(dt)':>9}  "
      f"{'Δmed(raw)':>10} {'Δmed(dt)':>10}")
print("  " + "-" * 140)

for wname, (wlo, whi) in TIGHT_WINDOWS.items():
    vals_raw = window_vals(wlo, whi, detrended=False)
    vals_dt  = window_vals(wlo, whi, detrended=True)
    n_pts    = len(vals_raw)

    s_raw = two_sample_stats(bl_raw, vals_raw)
    s_dt  = two_sample_stats(bl_dt, vals_dt)

    if s_raw and s_dt:
        print(f"  {wname:<22} {wlo:>5}–{whi:<5}  {n_pts:>4}  "
              f"{np.median(vals_raw):>9.2f} {np.median(vals_dt):>9.2f}  "
              f"{s_raw['p_mw']:>10.4f} {s_dt['p_mw']:>10.4f}  "
              f"{s_raw['rb']:>9.3f} {s_dt['rb']:>9.3f}  "
              f"{s_raw['dm']:>10.2f} {s_dt['dm']:>10.2f}")
    else:
        print(f"  {wname:<22} {wlo:>5}–{whi:<5}  {n_pts:>4}  insufficient data")


# ═══════════════════════════════════════════════════════════════════════
# 2. PROGRESSIVE FINE-RESOLUTION TEST (50-yr steps)
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 68)
print("2. PROGRESSIVE 100-yr WINDOWS (50-yr steps)")
print("=" * 68)

print(f"\n  Each 100-yr window tested against baseline ({BASELINE_LO}–{BASELINE_HI})")
print(f"\n  {'Centre':>7} {'Window':>14}  {'n':>4}  "
      f"{'Med(raw)':>9} {'Med(dt)':>9}  "
      f"{'p_MWU(dt)':>10}  {'r_RB(dt)':>9}  {'Δmed(dt)':>10}  {'Sig':>4}")
print("  " + "-" * 95)

fine_results = []
for centre in FINE_STEPS:
    wlo, whi = centre - 50, centre + 50
    vals_dt  = window_vals(wlo, whi, detrended=True)
    vals_raw = window_vals(wlo, whi, detrended=False)
    n_pts = len(vals_raw)

    s_dt  = two_sample_stats(bl_dt, vals_dt)
    s_raw = two_sample_stats(bl_raw, vals_raw)

    if s_dt and n_pts >= 2:
        sig = '*' if s_dt['p_mw'] < 0.05 else ''
        if s_dt['p_mw'] < 0.01: sig = '**'
        if s_dt['p_mw'] < 0.001: sig = '***'

        print(f"  {centre:>7.0f} {wlo:>5.0f}–{whi:<5.0f}  {n_pts:>4}  "
              f"{np.median(vals_raw):>9.2f} {np.median(vals_dt):>9.2f}  "
              f"{s_dt['p_mw']:>10.4f}  {s_dt['rb']:>9.3f}  "
              f"{s_dt['dm']:>10.2f}  {sig:>4}")

        fine_results.append({
            'centre': centre, 'lo': wlo, 'hi': whi, 'n': n_pts,
            'med_raw': np.median(vals_raw), 'med_dt': np.median(vals_dt),
            'p_mw_dt': s_dt['p_mw'], 'rb_dt': s_dt['rb'], 'dm_dt': s_dt['dm'],
            'p_mw_raw': s_raw['p_mw'] if s_raw else np.nan,
        })

fine_df = pd.DataFrame(fine_results)


# ═══════════════════════════════════════════════════════════════════════
# 3. ROLLING WINDOW ANALYSIS (multi-scale)
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 68)
print("3. ROLLING WINDOW ANALYSIS (multi-scale, detrended)")
print("=" * 68)

rolling_results = {}

for win_size in ROLLING_WINDOWS:
    hw = win_size / 2
    step = max(10, int(med_spacing))
    centres = np.arange(RANGE_LO + hw, RANGE_HI - hw, step)

    res = []
    for c in centres:
        wlo, whi = c - hw, c + hw
        vals_dt  = window_vals(wlo, whi, detrended=True)
        vals_raw = window_vals(wlo, whi, detrended=False)
        s = two_sample_stats(bl_dt, vals_dt)
        if s and len(vals_dt) >= 2:
            res.append({
                'centre_ka': c / 1000,
                'med_raw': np.median(vals_raw),
                'med_dt': np.median(vals_dt),
                'p_mw': s['p_mw'],
                'rb': s['rb'],
                'dm': s['dm'],
                'n': len(vals_dt),
            })

    rdf = pd.DataFrame(res)
    rolling_results[win_size] = rdf
    print(f"\n  Window = {win_size} yr ({len(rdf)} steps):")
    near82 = rdf[(rdf.centre_ka >= 7.9) & (rdf.centre_ka <= 8.5)]
    if len(near82) > 0:
        min_row = near82.loc[near82.p_mw.idxmin()]
        print(f"    Most significant near 8.2 ka: centre = {min_row.centre_ka:.2f} ka, "
              f"p = {min_row.p_mw:.4f}, r_RB = {min_row.rb:.3f}, Δmed = {min_row.dm:.2f}")
        sig_82 = near82[near82.p_mw < 0.05]
        print(f"    Steps with p < 0.05 in 7.9–8.5 ka: {len(sig_82)}/{len(near82)}")


# ═══════════════════════════════════════════════════════════════════════
# 4. LOCAL FLANKING-BASELINE COMPARISON
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 68)
print("4. LOCAL FLANKING-BASELINE COMPARISON (detrended)")
print("=" * 68)

# Event core vs combined pre + post flanks (immediate neighbours only)
ev_core = window_vals(8100, 8250, detrended=True)
pre_flank = window_vals(8250, 8400, detrended=True)
post_flank = window_vals(7950, 8100, detrended=True)
flanks = np.concatenate([pre_flank, post_flank])

print(f"\n  Event core (8,100–8,250 yr BP): n = {len(ev_core)}, "
      f"median = {np.median(ev_core):.2f} drips min⁻¹")
print(f"  Pre-flank  (8,250–8,400 yr BP): n = {len(pre_flank)}")
print(f"  Post-flank (7,950–8,100 yr BP): n = {len(post_flank)}")
print(f"  Combined flanks:                 n = {len(flanks)}, "
      f"median = {np.median(flanks):.2f} drips min⁻¹")

if len(ev_core) >= 2 and len(flanks) >= 2:
    delta = np.median(ev_core) - np.median(flanks)
    pct = delta / np.median(flanks) * 100
    u, p = mannwhitneyu(flanks, ev_core, alternative='two-sided')
    rb = rank_biserial(u, len(flanks), len(ev_core))
    print(f"\n  Δ median = {delta:.2f} drips min⁻¹ ({pct:+.1f}%)")
    print(f"  Mann–Whitney U p = {p:.4f}")
    print(f"  Rank-biserial r  = {rb:.3f}")

    if p > 0.05:
        print(f"\n  → No statistically detectable drought signal at 8.2 ka")
    else:
        print(f"\n  → Significant at p < 0.05")
else:
    print("\n  Insufficient data for flanking test")


# ═══════════════════════════════════════════════════════════════════════
# 5. FIGURE
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 68)
print("5. GENERATING FIGURE")
print("=" * 68)

fig, axes = plt.subplots(4, 1, figsize=(DOUBLE_COL, DOUBLE_COL * 1.15), sharex=True,
                          gridspec_kw={'hspace': 0.08,
                                       'height_ratios': [1.2, 0.8, 0.8, 0.6]})

# Event band
for ax in axes:
    ax.axvspan(8.1, 8.3, color=COL_BG_200, alpha=0.5, lw=0, zorder=0)
    ax.axvline(8.2, color=COL_BG_600, lw=0.3, ls=':')

# ── Panel a: Drip rate (raw + detrended) ──
ax = axes[0]
if col_p25 in df_sub.columns and col_p75 in df_sub.columns:
    ax.fill_between(df_sub.age / 1000, df_sub[col_p25], df_sub[col_p75],
                    color=COL_TEAL_200, alpha=0.7, lw=0, zorder=1)
ax.plot(df_sub.age / 1000, df_sub[col_med], color=COL_NI, lw=0.8, zorder=3,
        label='Median drip rate (raw)')
ax.plot(df_sub.age / 1000, df_sub.med_detrend, color=COL_BROWN_500, lw=0.8,
        ls='--', zorder=3, label='Detrended')
# Trend line
x_trend = np.array([RANGE_LO, RANGE_HI]) / 1000
y_trend = sl_full.slope * np.array([RANGE_LO, RANGE_HI]) + sl_full.intercept
ax.plot(x_trend, y_trend, color=COL_BG_600, lw=0.5, ls=':', zorder=2, label='Linear trend')

ax.set_ylabel('Drip rate\n(drips min$^{-1}$)')
ax.legend(loc='center left', bbox_to_anchor=(1.005, 0.5), frameon=False, fontsize=5.5, handlelength=1.6, borderaxespad=0.0)
panel_label(ax, 'a')

# ── Panel b: Rolling median (multi-scale) ──
ax = axes[1]
ls_styles = ['-', '--', '-.', ':']
_win_colors = [COL_NI, COL_BROWN_500, COL_TEAL_600, COL_BG_700]
for i, win_size in enumerate(ROLLING_WINDOWS):
    rdf = rolling_results[win_size]
    ax.plot(rdf.centre_ka, rdf.med_dt, color=_win_colors[i], lw=0.8, ls=ls_styles[i],
            label=f'{win_size}-yr window', zorder=3)

bl_med_val = np.median(window_vals(BASELINE_LO, BASELINE_HI, detrended=True))
ax.axhline(bl_med_val, color=COL_BG_600, lw=0.3, ls=':')
ax.set_ylabel('Detrended\nmedian')
ax.legend(loc='center left', bbox_to_anchor=(1.005, 0.5), frameon=False, fontsize=5.5, handlelength=1.6, borderaxespad=0.0)
panel_label(ax, 'b')

# ── Panel c: Rolling MWU p-value ──
ax = axes[2]
for i, win_size in enumerate(ROLLING_WINDOWS):
    rdf = rolling_results[win_size]
    ax.semilogy(rdf.centre_ka, rdf.p_mw, color=_win_colors[i], lw=0.8, ls=ls_styles[i],
                label=f'{win_size}-yr', zorder=3)

ax.axhline(0.05, color=COL_BG_900, lw=0.4, ls='--', label='p = 0.05')
ax.axhline(0.01, color=COL_BG_600, lw=0.3, ls=':', label='p = 0.01')
ax.set_ylabel('MWU p-value\n(vs baseline)')
ax.set_ylim(1e-4, 1.1)
ax.legend(loc='center left', bbox_to_anchor=(1.005, 0.5), frameon=False, fontsize=5.5, handlelength=1.6, borderaxespad=0.0)
panel_label(ax, 'c')

# ── Panel d: Rolling effect size (rank-biserial) ──
ax = axes[3]
for i, win_size in enumerate(ROLLING_WINDOWS):
    rdf = rolling_results[win_size]
    ax.plot(rdf.centre_ka, rdf.rb, color=_win_colors[i], lw=0.8, ls=ls_styles[i],
            label=f'{win_size}-yr', zorder=3)

ax.axhline(0, color=COL_BG_600, lw=0.3, ls='-')
ax.axhline(0.3, color=COL_DORANGE, lw=0.4, ls=':',
           label='|r| = 0.3 (moderate)')
ax.axhline(-0.3, color=COL_DORANGE, lw=0.4, ls=':')
ax.set_ylabel('Rank-biserial\neffect size')
ax.set_xlabel('Age (ka BP)')
ax.legend(loc='center left', bbox_to_anchor=(1.005, 0.5), frameon=False, fontsize=5.5, handlelength=1.6, borderaxespad=0.0)
panel_label(ax, 'd')

# Shared x
for ax in axes: style_ax(ax)
plt.subplots_adjust(right=0.80)
axes[-1].set_xlim(RANGE_HI / 1000, RANGE_LO / 1000)
axes[-1].xaxis.set_major_locator(ticker.MultipleLocator(0.2))
axes[-1].xaxis.set_minor_locator(ticker.MultipleLocator(0.1))

for fmt in ('pdf', 'png', 'eps'):
    fig.savefig(f'{args.output}.{fmt}', dpi=600, bbox_inches='tight',
                facecolor='white')
    print(f"  Saved → {args.output}.{fmt}")

plt.close(fig)


# ═══════════════════════════════════════════════════════════════════════
# SUMMARY
# ═══════════════════════════════════════════════════════════════════════

print("\n" + "=" * 68)
print("SUMMARY")
print("=" * 68)

print(f"""
  Data: drip_rate_summary.csv ({len(df)} pts) on HS4 chronology
  Age model: consistent with all main-text and ED figures

  The coarse 500-yr window test (9000–9500 vs 8000–8500) in the original
  stationarity script compared apples to oranges — the orbital trend alone
  drives a ~{sl_full.slope * 1500:.1f} drips min⁻¹ difference over that
  interval, swamping any transient event signal.

  The fine-grained, detrended analysis on the best-estimate series shows:
""")

if len(ev_core) >= 2 and len(flanks) >= 2:
    print(f"  - Event core (8,100–8,250): n = {len(ev_core)}, "
          f"median = {np.median(ev_core):.2f} drips min⁻¹")
    print(f"  - Combined flanks:          n = {len(flanks)}, "
          f"median = {np.median(flanks):.2f} drips min⁻¹")
    print(f"  - Δ = {delta:.2f} drips min⁻¹ ({pct:+.1f}%)")
    print(f"  - Mann–Whitney U p = {p:.4f}")
    print(f"  - Rank-biserial r  = {rb:.3f}")
    if p > 0.05:
        print(f"\n  → No statistically detectable drought signal at 8.2 ka")

print("""
  Interpretation: the kinetic proxy does not record a statistically
  significant drip rate reduction at 8.2 ka. The δ¹³C excursion at this
  interval supports enhanced ventilation and prior carbonate precipitation
  rather than a rainfall deficit as the driver of the δ¹⁸O anomaly.
""")

print("Done.\n")
