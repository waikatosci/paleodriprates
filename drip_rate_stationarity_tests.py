"""
drip_rate_stationarity_tests.py
================================
Statistical tests for event significance and stationarity in the HS4
Holocene drip rate reconstruction (Hartland et al., Nature Geoscience).

Two inputs:
  1. CSV  — Monte Carlo realisations of drip rate (drips min-1).
            Age column (any capitalisation) + one column per realisation.
  2. XLSX — Full-posterior summary with columns for age, median, p25, p75
            (column names auto-detected — capitalisation does not matter).

Tests performed:
  1. 5.2 ka aridity pulse  — IQR non-overlap + MWU/KS on realisation pools
  2. 8.2 ka stationarity   — MWU/KS + IQR overlap confirming no event signal
  3. Post-1820 CE surge    — MWU/KS vs late-Holocene baseline
  4. Long-term Holocene trend  — Mann-Kendall + Theil-Sen
  5. Augmented Dickey-Fuller stationarity (requires statsmodels)
"""

import numpy as np
import pandas as pd
from scipy.stats import mannwhitneyu, ks_2samp, norm
import warnings
warnings.filterwarnings('ignore')

try:
    from statsmodels.tsa.stattools import adfuller
    HAS_ADF = True
except ImportError:
    HAS_ADF = False

# ---------------------------------------------------------------------------
# *** SET PATHS HERE ***
# ---------------------------------------------------------------------------
CSV_PATH  = r'C:\Users\hartlana\OneDrive - lincolnagritech.co.nz\Documents\Python\N Geosci - Drip rates\Main\Fig 5 Record\Drip_rate_realisations_HS4.csv'
XLSX_PATH = r'C:\Users\hartlana\OneDrive - lincolnagritech.co.nz\Documents\Python\N Geosci - Drip rates\Main\Fig 5 Record\Drip_rate_data_HS4.xlsx'

# ---------------------------------------------------------------------------
# Load data — normalise ALL column names to lowercase on load
# so the rest of the script never worries about capitalisation
# ---------------------------------------------------------------------------

# Realisations CSV
df_r = pd.read_csv(CSV_PATH)
df_r.columns = [c.lower() for c in df_r.columns]   # 'Age' -> 'age', etc.
real_cols = [c for c in df_r.columns if c != 'age']
n_real = len(real_cols)

# Full-posterior XLSX
df_m = pd.read_excel(XLSX_PATH)
df_m.columns = [c.lower() for c in df_m.columns]   # 'Age'->'age', 'DR_med'->'dr_med'

# Accept flexible naming in XLSX: e.g. 'dr_med', 'median', 'pc25', 'p25', etc.
# Rename to standard names used throughout this script.
col_renames = {}
for col in df_m.columns:
    if col == 'age':
        pass  # already correct
    elif 'med' in col:
        col_renames[col] = 'dr_med'
    elif 'pc25' in col or ('25' in col and col != 'age'):
        col_renames[col] = 'dr_pc25'
    elif 'pc75' in col or ('75' in col and col != 'age'):
        col_renames[col] = 'dr_pc75'
df_m = df_m.rename(columns=col_renames)
df_m = df_m.sort_values('age').reset_index(drop=True)

# Safety check
required = {'age', 'dr_med', 'dr_pc25', 'dr_pc75'}
missing_cols = required - set(df_m.columns)
if missing_cols:
    raise ValueError(
        f"Could not map required columns {missing_cols} from XLSX.\n"
        f"Columns found after lowercasing: {df_m.columns.tolist()}\n"
        f"Edit the col_renames block above to match your file."
    )

print("=" * 68)
print("HS4 DRIP RATE - STATISTICAL TESTS")
print("Hartland et al. (Nature Geoscience)")
print("=" * 68)
print(f"  Realisations : {n_real} realisations, {len(df_r)} age pts, "
      f"{df_r.age.min():.0f} to {df_r.age.max():.0f} yr BP")
print(f"  Full posterior: {len(df_m)} age pts, "
      f"{df_m.age.min():.0f} to {df_m.age.max():.0f} yr BP")
print()

# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------

def window_r(lo, hi):
    """All realisation values in age window [lo, hi] BP."""
    m = (df_r['age'] >= lo) & (df_r['age'] <= hi)
    return df_r.loc[m, real_cols].values.flatten(), m.sum()

def window_m(lo, hi):
    """Rows of full-posterior summary in age window [lo, hi] BP."""
    return df_m[(df_m['age'] >= lo) & (df_m['age'] <= hi)].copy()

def rank_biserial(u, n1, n2):
    """Effect size from Mann-Whitney U."""
    return 1.0 - (2.0 * u) / (n1 * n2)

def bootstrap_median_diff(a, b, n_boot=10_000, seed=42, max_n=50_000):
    """Percentile bootstrap 95% CI on median(b) - median(a).

    If either array exceeds max_n values it is randomly subsampled before
    bootstrapping to avoid prohibitive memory use when pooled realisations
    run into the hundreds of thousands.
    """
    rng = np.random.default_rng(seed)
    if len(a) > max_n:
        a = rng.choice(a, max_n, replace=False)
    if len(b) > max_n:
        b = rng.choice(b, max_n, replace=False)
    d = [np.median(rng.choice(b, len(b), replace=True)) -
         np.median(rng.choice(a, len(a), replace=True))
         for _ in range(n_boot)]
    return np.percentile(d, [2.5, 97.5])

def per_real_mw(lo_a, hi_a, lo_b, hi_b):
    """Mann-Whitney p-value for each realisation individually.

    Prints progress every 100 realisations so the user can confirm
    the loop is running rather than hung.
    """
    ma = (df_r['age'] >= lo_a) & (df_r['age'] <= hi_a)
    mb = (df_r['age'] >= lo_b) & (df_r['age'] <= hi_b)
    ps = []
    n = len(real_cols)
    for i, col in enumerate(real_cols):
        if i % 100 == 0:
            print(f"    ... per-realisation MWU: {i}/{n}", flush=True)
        a = df_r.loc[ma, col].values
        b = df_r.loc[mb, col].values
        if len(a) >= 2 and len(b) >= 2:
            _, p = mannwhitneyu(a, b, alternative='two-sided')
            ps.append(p)
    print(f"    ... per-realisation MWU: {n}/{n} done", flush=True)
    return np.array(ps)

def iqr_sep(ev_lo, ev_hi, bg_lo, bg_hi):
    """
    Compare event vs background IQRs using the full posterior.
    Returns: overlap(bool), gap(float, + means event below bg), ev_iqr, bg_iqr.
    """
    ev = window_m(ev_lo, ev_hi)
    bg = window_m(bg_lo, bg_hi)
    e_p25 = ev['dr_pc25'].min()
    e_p75 = ev['dr_pc75'].max()
    b_p25 = bg['dr_pc25'].median()
    b_p75 = bg['dr_pc75'].median()
    overlaps = (e_p75 >= b_p25) and (b_p75 >= e_p25)
    gap = b_p25 - e_p75
    return overlaps, gap, (e_p25, e_p75), (b_p25, b_p75)

def mannkendall(x):
    """Mann-Kendall tau, Wald Z, two-tailed p."""
    n = len(x)
    s = sum(np.sign(x[j] - x[i])
            for i in range(n - 1) for j in range(i + 1, n))
    var_s = n * (n - 1) * (2 * n + 5) / 18.0
    z = (s - np.sign(s)) / np.sqrt(var_s)
    p = 2.0 * (1.0 - norm.cdf(abs(z)))
    tau = s / (0.5 * n * (n - 1))
    return tau, z, p

def theilsen(x, y):
    """Theil-Sen median slope."""
    slopes = [(y[j] - y[i]) / (x[j] - x[i])
              for i in range(len(x)) for j in range(i + 1, len(x))
              if x[j] != x[i]]
    return np.median(slopes)

def two_sample(label, vals_a, vals_b, la='A', lb='B'):
    """Run MWU + KS + bootstrap and print a standard block. Returns key stats."""
    u, p_mw = mannwhitneyu(vals_a, vals_b, alternative='two-sided')
    d, p_ks = ks_2samp(vals_a, vals_b)
    rb = rank_biserial(u, len(vals_a), len(vals_b))
    dm = np.median(vals_b) - np.median(vals_a)
    ci = bootstrap_median_diff(vals_a, vals_b)
    print(f"\n  {label}")
    print(f"    n({la}) = {len(vals_a)},  n({lb}) = {len(vals_b)}")
    print(f"    Medians: {la} = {np.median(vals_a):.2f},  "
          f"{lb} = {np.median(vals_b):.2f}  drips min-1")
    print(f"    Mann-Whitney U = {u:.0f},  p = {p_mw:.4f}")
    print(f"    KS D = {d:.4f},  p = {p_ks:.4f}")
    print(f"    Rank-biserial r = {rb:.3f}  (|r|>0.3 moderate, >0.5 large)")
    print(f"    Delta-median ({lb}-{la}) = {dm:.2f} drips min-1")
    print(f"    Bootstrap 95% CI: [{ci[0]:.2f}, {ci[1]:.2f}]")
    return p_mw, p_ks, rb, dm, ci

results = {}   # collects results for summary table


# ===========================================================================
# TEST 1: 5.2 ka ARIDITY PULSE
# ===========================================================================
print("-" * 68)
print("TEST 1: 5.2 ka ARIDITY PULSE")
print("-" * 68)

BG52_LO, BG52_HI   = 5500, 6500   # pre-event background
EV52_LO, EV52_HI   = 4600, 4700   # event minimum (IQR assessment)
EV52B_LO, EV52B_HI = 4500, 5200   # broader window for realisation-level tests

bg52 = window_m(BG52_LO, BG52_HI)
ev52 = window_m(EV52_LO, EV52_HI)

print(f"\n  Background ({BG52_LO}-{BG52_HI} BP, n={len(bg52)} pts):")
print(f"    Median = {bg52.dr_med.median():.2f}  "
      f"IQR [{bg52.dr_pc25.median():.2f}, {bg52.dr_pc75.median():.2f}] drips min-1")

print(f"\n  Event minimum ({EV52_LO}-{EV52_HI} BP, n={len(ev52)} pts):")
for _, row in ev52.iterrows():
    print(f"    Age {row.age:.0f} BP:  "
          f"dr_med = {row.dr_med:.2f}  [{row.dr_pc25:.2f}, {row.dr_pc75:.2f}]")

ovlp52, gap52, ei52, bi52 = iqr_sep(EV52_LO, EV52_HI, BG52_LO, BG52_HI)
print(f"\n  Full-posterior IQR overlap:")
print(f"    Event IQR     : [{ei52[0]:.2f}, {ei52[1]:.2f}]")
print(f"    Background IQR: [{bi52[0]:.2f}, {bi52[1]:.2f}]")
if not ovlp52:
    print(f"    Overlap: NO — event p75 is {gap52:.2f} drips min-1 below background p25")
    print(f"    => Event minimum ENTIRELY BELOW background IQR (strong separation)")
else:
    print(f"    Overlap: YES")

vals_bg52, _ = window_r(BG52_LO,  BG52_HI)
vals_ev52, _ = window_r(EV52B_LO, EV52B_HI)
p52_mw, p52_ks, rb52, dm52, ci52 = two_sample(
    f"Background ({BG52_LO}-{BG52_HI} BP) vs event ({EV52B_LO}-{EV52B_HI} BP)",
    vals_bg52, vals_ev52, la='bg', lb='event')

ps52 = per_real_mw(BG52_LO, BG52_HI, EV52B_LO, EV52B_HI)
print(f"\n  Per-realisation MWU ({n_real} realisations):")
print(f"    Median p = {np.median(ps52):.4f},  "
      f"range [{ps52.min():.4f}, {ps52.max():.4f}]")
print(f"    Fraction p < 0.05: {(ps52 < 0.05).mean():.2f}  "
      f"({(ps52 < 0.05).sum()}/{len(ps52)})")

results['5.2ka'] = (p52_mw, p52_ks, rb52, dm52, not ovlp52)


# ===========================================================================
# TEST 2: 8.2 ka STATIONARITY
# ===========================================================================
print()
print("-" * 68)
print("TEST 2: 8.2 ka STATIONARITY")
print("-" * 68)

BG82_LO,   BG82_HI   = 9000, 9500
EV82_LO,   EV82_HI   = 8000, 8500
POST82_LO, POST82_HI = 7500, 8000

rec_max = df_m.age.max()

if rec_max < EV82_LO:
    print(f"\n  Record ends at {rec_max:.0f} BP — 8.2 ka window is outside the data.")
    print("  Generate realisations covering 7,000-9,500 BP to enable this test.")
    results['8.2ka'] = None
else:
    bg82 = window_m(BG82_LO,  BG82_HI)
    ev82 = window_m(EV82_LO,  EV82_HI)
    po82 = window_m(POST82_LO, POST82_HI)

    print(f"\n  Pre-event reference ({BG82_LO}-{BG82_HI} BP, n={len(bg82)} pts):")
    print(f"    Median = {bg82.dr_med.median():.2f}  "
          f"IQR [{bg82.dr_pc25.median():.2f}, {bg82.dr_pc75.median():.2f}]")
    print(f"\n  Event window ({EV82_LO}-{EV82_HI} BP, n={len(ev82)} pts):")
    print(f"    Median = {ev82.dr_med.median():.2f}  "
          f"IQR [{ev82.dr_pc25.median():.2f}, {ev82.dr_pc75.median():.2f}]")
    print(f"\n  Post-event reference ({POST82_LO}-{POST82_HI} BP, n={len(po82)} pts):")
    print(f"    Median = {po82.dr_med.median():.2f}  "
          f"IQR [{po82.dr_pc25.median():.2f}, {po82.dr_pc75.median():.2f}]")

    ovlp82, gap82, ei82, bi82 = iqr_sep(EV82_LO, EV82_HI, BG82_LO, BG82_HI)
    print(f"\n  Full-posterior IQR overlap (event vs pre-event):")
    print(f"    Event IQR  : [{ei82[0]:.2f}, {ei82[1]:.2f}]")
    print(f"    Pre-event  : [{bi82[0]:.2f}, {bi82[1]:.2f}]")
    if ovlp82:
        print("    Overlap: YES — consistent with stationarity (no event signal)")
    else:
        print(f"    Overlap: NO — separation of {abs(gap82):.2f} drips min-1")

    vals_bg82, _ = window_r(BG82_LO,  BG82_HI)
    vals_ev82, _ = window_r(EV82_LO,  EV82_HI)
    vals_po82, _ = window_r(POST82_LO, POST82_HI)

    p82_mw, p82_ks, rb82, dm82, ci82 = two_sample(
        f"Pre-event ({BG82_LO}-{BG82_HI} BP) vs event ({EV82_LO}-{EV82_HI} BP)",
        vals_bg82, vals_ev82, la='pre', lb='event')
    print("    (Stationarity expected: p > 0.05, small |r|, CI bracketing zero)")

    if len(vals_po82) >= 2:
        two_sample(
            f"Post-event ({POST82_LO}-{POST82_HI} BP) vs event ({EV82_LO}-{EV82_HI} BP)",
            vals_po82, vals_ev82, la='post', lb='event')

    ps82 = per_real_mw(BG82_LO, BG82_HI, EV82_LO, EV82_HI)
    if len(ps82) > 0:
        print(f"\n  Per-realisation MWU ({n_real} realisations):")
        print(f"    Median p = {np.median(ps82):.4f},  "
              f"range [{ps82.min():.4f}, {ps82.max():.4f}]")
        print(f"    Fraction p < 0.05: {(ps82 < 0.05).mean():.2f}  "
              f"({(ps82 < 0.05).sum()}/{len(ps82)})")

    results['8.2ka'] = (p82_mw, p82_ks, rb82, dm82, ovlp82)


# ===========================================================================
# TEST 3: POST-1820 CE SURGE
# ===========================================================================
print()
print("-" * 68)
print("TEST 3: POST-1820 CE PRECIPITATION SURGE")
print("-" * 68)

LH_LO,  LH_HI  = 500, 2000
MOD_LO, MOD_HI = -400, 50     # catches all CE points back to ~1550 CE

lh_sub  = window_m(LH_LO,  LH_HI)
mod_sub = window_m(MOD_LO, MOD_HI)

print(f"\n  Late-Holocene baseline ({LH_LO}-{LH_HI} BP, n={len(lh_sub)} pts):")
print(f"    Median = {lh_sub.dr_med.median():.2f}  "
      f"IQR [{lh_sub.dr_pc25.median():.2f}, {lh_sub.dr_pc75.median():.2f}]")
print(f"\n  Post-1820 CE window ({MOD_LO} to {MOD_HI} BP, n={len(mod_sub)} pts):")
if len(mod_sub) > 0:
    for _, row in mod_sub.iterrows():
        print(f"    Age {row.age:.0f} BP:  "
              f"dr_med = {row.dr_med:.2f}  [{row.dr_pc25:.2f}, {row.dr_pc75:.2f}]")
else:
    print("    No data — check that ages extend to negative BP values")

vals_lh,  _ = window_r(LH_LO,  LH_HI)
vals_mod, _ = window_r(MOD_LO, MOD_HI)

if len(vals_mod) >= 4:
    ovlp_mod, gap_mod, ei_mod, bi_mod = iqr_sep(MOD_LO, MOD_HI, LH_LO, LH_HI)
    print(f"\n  Full-posterior IQR overlap (post-1820 vs baseline):")
    print(f"    Post-1820 IQR : [{ei_mod[0]:.2f}, {ei_mod[1]:.2f}]")
    print(f"    Baseline IQR  : [{bi_mod[0]:.2f}, {bi_mod[1]:.2f}]")
    direction = 'ABOVE' if gap_mod < 0 else 'BELOW'
    if not ovlp_mod:
        print(f"    Overlap: NO — post-1820 is {abs(gap_mod):.2f} drips min-1 "
              f"{direction} baseline IQR")
    else:
        print(f"    Overlap: YES")

    p_mod_mw, p_mod_ks, rb_mod, dm_mod, ci_mod = two_sample(
        "Late-Holocene baseline vs post-1820 CE",
        vals_lh, vals_mod, la='LH', lb='post-1820')

    ps_mod = per_real_mw(LH_LO, LH_HI, MOD_LO, MOD_HI)
    if len(ps_mod) > 0:
        print(f"\n  Per-realisation MWU:")
        print(f"    Median p = {np.median(ps_mod):.4f},  "
              f"range [{ps_mod.min():.4f}, {ps_mod.max():.4f}]")
        print(f"    Fraction p < 0.05: {(ps_mod < 0.05).mean():.2f}  "
              f"({(ps_mod < 0.05).sum()}/{len(ps_mod)})")

    results['surge'] = (p_mod_mw, p_mod_ks, rb_mod, dm_mod)
else:
    print(f"\n  Only {len(vals_mod)} realisation values in post-1820 window.")
    print("  Need at least 4 for a two-sample test.")
    results['surge'] = None


# ===========================================================================
# TEST 4: HOLOCENE TREND — MANN-KENDALL + THEIL-SEN
# ===========================================================================
print()
print("-" * 68)
print("TEST 4: LONG-TERM HOLOCENE TREND (MANN-KENDALL + THEIL-SEN)")
print("-" * 68)

ages_mk = df_m['age'].values
med_mk  = df_m['dr_med'].values

tau, z_mk, p_mk = mannkendall(med_mk)
slope = theilsen(ages_mk, med_mk)

print(f"\n  Full series: n = {len(med_mk)},  "
      f"ages {ages_mk.min():.0f} to {ages_mk.max():.0f} yr BP")
print(f"  Mann-Kendall tau = {tau:.4f},  Z = {z_mk:.3f},  p = {p_mk:.4f}")
print(f"  Theil-Sen slope  = {slope:.5f} drips min-1 yr-1")
print(f"                   = {slope * 1000:.3f} drips min-1 kyr-1")
if slope > 0:
    print("  Direction: positive => older = higher drip rates (Holocene drying)")
else:
    print("  Direction: negative => older = lower drip rates")
print(f"  Implied change over {ages_mk.max() - ages_mk.min():.0f} yr: "
      f"{slope * (ages_mk.max() - ages_mk.min()):.2f} drips min-1")

taus_r, ps_r = [], []
for col in real_cols:
    t, _, p = mannkendall(df_r[col].values)
    taus_r.append(t); ps_r.append(p)

print(f"\n  Per-realisation MK ({n_real} realisations):")
print(f"    Median tau = {np.median(taus_r):.4f},  "
      f"range [{min(taus_r):.4f}, {max(taus_r):.4f}]")
print(f"    Median p   = {np.median(ps_r):.4f}")
print(f"    Fraction p < 0.05: {(np.array(ps_r) < 0.05).mean():.2f}  "
      f"({(np.array(ps_r) < 0.05).sum()}/{n_real})")

early_sub = window_m(5500, 6500)
late_sub  = window_m(0, 2000)
vals_e, _ = window_r(5500, 6500)
vals_l, _ = window_r(0, 2000)

if len(vals_e) >= 2 and len(vals_l) >= 2:
    dm_hl = late_sub.dr_med.median() - early_sub.dr_med.median()
    ci_hl = bootstrap_median_diff(vals_e, vals_l)
    u_hl, p_hl = mannwhitneyu(vals_e, vals_l, alternative='two-sided')
    _, p_hl_ks = ks_2samp(vals_e, vals_l)
    rb_hl = rank_biserial(u_hl, len(vals_e), len(vals_l))

    print(f"\n  Magnitude — early (5,500-6,500 BP) vs late (0-2,000 BP):")
    print(f"    Early  = {early_sub.dr_med.median():.2f}  "
          f"[{early_sub.dr_pc25.median():.2f}, {early_sub.dr_pc75.median():.2f}]")
    print(f"    Late   = {late_sub.dr_med.median():.2f}  "
          f"[{late_sub.dr_pc25.median():.2f}, {late_sub.dr_pc75.median():.2f}]")
    print(f"    Delta-median (late - early) = {dm_hl:.2f} drips min-1")
    print(f"    Bootstrap 95% CI: [{ci_hl[0]:.2f}, {ci_hl[1]:.2f}]")
    print(f"    Mann-Whitney p = {p_hl:.4f},  KS p = {p_hl_ks:.4f}")
    results['trend'] = (p_hl, p_hl_ks, rb_hl, dm_hl, tau, p_mk, slope)


# ===========================================================================
# TEST 5: AUGMENTED DICKEY-FULLER
# ===========================================================================
print()
print("-" * 68)
print("TEST 5: AUGMENTED DICKEY-FULLER STATIONARITY")
print("-" * 68)

if HAS_ADF:
    adf_stat, adf_p, lags, nobs, crit, _ = adfuller(med_mk, autolag='AIC')
    print(f"\n  ADF statistic = {adf_stat:.4f},  p = {adf_p:.4f}")
    print(f"  Lags used = {lags},  n obs = {nobs}")
    print(f"  Critical values:  1% = {crit['1%']:.3f},  "
          f"5% = {crit['5%']:.3f},  10% = {crit['10%']:.3f}")
    if adf_p < 0.05:
        print("  => Reject unit root (p < 0.05): series is stationary around its trend")
    else:
        print("  => Cannot reject unit root (p >= 0.05): consistent with non-stationarity")
    print()
    print("  NOTE: ADF tests for a random-walk unit root, not flatness.")
    print("  The orbitally-forced Holocene decline makes this series trend-")
    print("  stationary — physically expected for a forced climate record.")
else:
    print("\n  statsmodels not installed — skipping ADF test.")
    print("  Install with:  conda install statsmodels  or  pip install statsmodels")


# ===========================================================================
# SUMMARY TABLE
# ===========================================================================
print()
print("=" * 68)
print("SUMMARY")
print("=" * 68)

def sig(p):
    if p is None: return '  n/a'
    if p < 0.001: return '***'
    if p < 0.01:  return ' **'
    if p < 0.05:  return '  *'
    return '  ns'

print(f"\n  {'Test':<40} {'MW p':>8}  {'KS p':>8}  {'r_RB':>6}  {'Dm':>7}  {''}")
print(f"  {'-'*40} {'-'*8}  {'-'*8}  {'-'*6}  {'-'*7}  ---")

if results.get('5.2ka'):
    p1, k1, r1, d1, sep1 = results['5.2ka']
    print(f"  {'5.2ka: background vs event':<40} "
          f"{p1:>8.4f}  {k1:>8.4f}  {r1:>6.3f}  {d1:>7.2f}  {sig(p1)}")
    print(f"    IQR non-overlap (full posterior): "
          f"{'CONFIRMED' if sep1 else 'not confirmed'}")

if results.get('8.2ka'):
    p2, k2, r2, d2, o2 = results['8.2ka']
    print(f"  {'8.2ka: pre-event vs event':<40} "
          f"{p2:>8.4f}  {k2:>8.4f}  {r2:>6.3f}  {d2:>7.2f}  {sig(p2)}")
    print(f"    IQR overlap (stationarity): "
          f"{'YES — no event signal' if o2 else 'NO — anomaly present'}")
elif results.get('8.2ka') is None:
    print(f"  {'8.2ka':<40} {'outside record range':>40}")

if results.get('surge'):
    p3, k3, r3, d3 = results['surge']
    print(f"  {'Post-1820 CE surge':<40} "
          f"{p3:>8.4f}  {k3:>8.4f}  {r3:>6.3f}  {d3:>7.2f}  {sig(p3)}")

if results.get('trend'):
    p4, k4, r4, d4, tau4, pmk4, sl4 = results['trend']
    print(f"  {'Holocene: early vs late':<40} "
          f"{p4:>8.4f}  {k4:>8.4f}  {r4:>6.3f}  {d4:>7.2f}  {sig(p4)}")
    print(f"  {'Mann-Kendall tau / p':<40} "
          f"{tau4:>8.4f}  {'p='+f'{pmk4:.4f}':>8}  "
          f"  slope = {sl4*1000:.3f} drips min-1 kyr-1")

print()
print("  *** p<0.001  ** p<0.01  * p<0.05  ns not significant")
print("  Dm  = Delta-median (group B - group A) in drips min-1")
print("  r_RB = rank-biserial effect size (|r|>0.3 moderate, >0.5 large)")
print()
