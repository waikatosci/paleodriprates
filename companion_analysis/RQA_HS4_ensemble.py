"""
Recurrence Quantification Analysis — HS4 Speleothem Study
==========================================================
Author  : Adam Hartland (ensemble RQA framework)
Figure  : Fig. 5 + Supplementary RP figure

INPUT FILES  (all produced by Drip_rate_parallel_fr.py)
--------------------------------------------------------
drip_rate_realisations.csv   — age + r0...r999  (drip rate, drips min-1)
d18O_realisations.csv        — age + r0...r999  (d18O, per mil VPDB)
Drip_rate.xlsx               — summary statistics
    sheet '5.OutDripRate'    : age, pc05, pc10, pc25, med, pc75, pc90, pc95
    sheet '6.OutIsotope'     : age, median

OUTPUT FILES
------------
rqa_ensemble_results.csv     — per-window ensemble statistics + global tau/m
Fig5_RQA_HS4.pdf/.png/.eps   — main figure
FigS1_RP_HS4.pdf/.png        — supplementary recurrence plot figure

RQA PARAMETER RATIONALE
------------------------
Time delay (tau):
    Estimated ONCE from the full median time series of each proxy via the
    first minimum of Average Mutual Information (Fraser & Swinney 1986,
    Phys Rev A 33:1134). This global tau is then applied uniformly to all
    windows and all realisations. Per-window estimation on 100-point
    segments of palaeoclimate proxy data is insufficiently stable;
    a globally-estimated tau ensures a consistent phase-space
    representation throughout the record.

Embedding dimension (m):
    Estimated ONCE from the full median time series via Cao's method
    (Cao 1997, Physica D 110:43-50), bounded to m in [2, 5].
    Same rationale as tau — global estimation is more stable and
    produces a consistent attractor reconstruction.

Recurrence threshold (epsilon):
    Fixed Recurrence Rate (RR = 0.05) following Kraemer et al. (2018,
    Chaos 28:085720) and Marwan et al. (2007, Phys Rep 438:237-329).
    Fixed RR makes DET and TRANS directly comparable across windows and
    realisations even when signal amplitude varies.

Theiler window (w):
    w = tau x (m-1), excluding autocorrelated neighbours along the
    line of identity (Theiler 1986; Marwan 2010, arXiv:1007.2215).

Minimum diagonal line length (l_min = 2):
    Standard (Zbilut & Webber 1992, Phys Lett A 171:199).

ENSEMBLE STRATEGY
-----------------
RQA is run independently on each of the 1000 Monte-Carlo realisations
for both proxies. Ensemble median and 5th/95th percentile envelopes of
DET and TRANS are reported at each window centre, propagating the full
model uncertainty into the dynamical measures.
"""

# ---------------------------------------------------------------------------
# IMPORTS
# ---------------------------------------------------------------------------
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import interp1d
from scipy.spatial.distance import pdist, squareform
from scipy.ndimage import gaussian_filter1d
from itertools import groupby

warnings.filterwarnings('ignore')

# ---------------------------------------------------------------------------
# MATPLOTLIB STYLE  (Nature Geoscience)
# ---------------------------------------------------------------------------
plt.style.use('default')
plt.rcParams.update({
    'font.family'       : 'sans-serif',
    'font.sans-serif'   : ['Helvetica', 'Arial'],
    'font.size'         : 7,
    'pdf.fonttype'      : 42,
    'svg.fonttype'      : 'none',
    'axes.linewidth'    : 0.5,
    'xtick.major.width' : 0.5,
    'ytick.major.width' : 0.5,
    'xtick.major.size'  : 2.5,
    'ytick.major.size'  : 2.5,
    'legend.fontsize'   : 6,
    'axes.labelsize'    : 8,
    'axes.labelweight'  : 'bold',
    'xtick.labelsize'   : 7,
    'ytick.labelsize'   : 7,
})

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------
DRIP_REAL_FILE = 'drip_rate_realisations.csv'
D18O_REAL_FILE = 'd18O_realisations.csv'
DRIP_XLSX      = 'Drip_rate.xlsx'
DRIP_SHEET     = '5.OutDripRate'
D18O_SHEET     = '6.OutIsotope'

OUTPUT_CSV     = 'rqa_ensemble_results.csv'
OUTPUT_PDF     = 'Fig5_RQA_HS4.pdf'
OUTPUT_PNG     = 'Fig5_RQA_HS4.png'
OUTPUT_EPS     = 'Fig5_RQA_HS4.eps'
SUPP_PDF       = 'FigS1_RP_HS4.pdf'
SUPP_PNG       = 'FigS1_RP_HS4.png'

# RQA parameters
WINDOW_SIZE   = 100     # points
WINDOW_STEP   = 5       # points
TARGET_RR     = 0.05    # fixed recurrence rate (Kraemer et al. 2018)
MIN_LINE      = 2       # minimum diagonal line length
MAX_M         = 5       # upper bound on embedding dimension (Cao's method)
SMOOTH_SIGMA  = 2.0     # Gaussian smoothing sigma (points) for display
DET_THRESHOLD = 0.4     # reference line on RQA panels
MIN_SUSTAIN   = 80      # min consecutive windows below threshold for marker
N_YBINS       = 300     # PDF heatmap y-resolution

# Age of interest for RP supplementary annotation (yr BP)
AGE_MARKER_KA = 5200

# Set to integer (e.g. 200) for quick test; None = all 1000
N_REALISATIONS = None

# ---------------------------------------------------------------------------
# PARAMETER OVERRIDES
# ---------------------------------------------------------------------------
# By default (None) tau and m are estimated from the full median series via
# AMI and Cao's method.  Set integer values here to override for BOTH proxies,
# e.g. to compare data-driven vs. literature-informed choices before committing
# to a full run.
#
# Suggested comparisons:
#   Data-driven (estimated):  FORCE_TAU = None,  FORCE_M = None
#   Conservative literature:  FORCE_TAU = 1,     FORCE_M = 3
#   Current estimated result: FORCE_TAU = 1,     FORCE_M = 5
#
# Theiler window is always recomputed as max(1, tau*(m-1)) regardless.
FORCE_TAU = None   # e.g. 1  -- set to integer to override AMI estimation
FORCE_M   = None   # e.g. 3  -- set to integer to override Cao estimation

# ---------------------------------------------------------------------------
# PARAMETER ESTIMATION  (global, from median series)
# ---------------------------------------------------------------------------

def ami_time_delay(x, max_tau=30):
    """
    First minimum of Average Mutual Information (Fraser & Swinney 1986).
    Applied to the FULL median series for stable global tau estimation.
    Falls back to tau=1 if AMI is monotone.
    """
    n        = len(x)
    prev_ami = np.inf
    best_tau = 1
    for tau in range(1, min(max_tau + 1, n // 4)):
        x1   = x[:-tau]
        x2   = x[tau:]
        bins = max(10, int(np.sqrt(len(x1))))
        h2d, _, _ = np.histogram2d(x1, x2, bins=bins, density=True)
        hx   = h2d.sum(axis=1)
        hy   = h2d.sum(axis=0)
        mask = (h2d > 0) & (np.outer(hx, hy) > 0)
        ami  = np.sum(h2d[mask] * np.log(
               h2d[mask] / np.outer(hx, hy)[mask]))
        if ami < prev_ami:
            prev_ami = ami
        else:
            best_tau = max(1, tau - 1)
            break
    return best_tau


def cao_embedding_dim(x, tau=1, max_m=5):
    """
    Cao's method (Cao 1997). Applied to FULL median series.
    Returns m in [2, max_m].
    """
    n   = len(x)
    E1s = []
    for m in range(1, max_m + 1):
        N = n - m * tau
        if N < 10:
            break
        vecs = np.array([x[i:i + m * tau:tau] for i in range(N)])
        if len(vecs) < 4:
            break
        dists = squareform(pdist(vecs, metric='chebyshev'))
        np.fill_diagonal(dists, np.inf)
        nn_idx = np.argmin(dists, axis=1)
        N_ext  = n - (m + 1) * tau
        if N_ext < 4:
            break
        ratios = []
        for i in range(min(N, N_ext)):
            j   = nn_idx[i]
            d_m = dists[i, j]
            if j < N_ext and np.isfinite(d_m) and d_m > 1e-12:
                if i + m * tau < n and j + m * tau < n:
                    ratios.append(
                        abs(x[i + m * tau] - x[j + m * tau]) / d_m)
        if ratios:
            E1s.append(np.mean(ratios))
    if len(E1s) < 2:
        return 2
    for i in range(1, len(E1s)):
        if E1s[i - 1] > 0 and abs(E1s[i] / E1s[i - 1] - 1.0) < 0.1:
            return min(i + 1, max_m)
    return max_m


def estimate_global_params(median_series, label='proxy'):
    """
    Estimate tau and m once from the full normalised median series.
    Respects FORCE_TAU / FORCE_M overrides if set.
    Prints and returns (tau, m, theiler).
    """
    x = np.asarray(median_series, dtype=float)
    x = x[np.isfinite(x)]
    x = (x - np.mean(x)) / np.std(x)

    if FORCE_TAU is not None:
        tau = int(FORCE_TAU)
        print(f"  {label}: tau={tau} (FORCED)")
    else:
        tau = ami_time_delay(x, max_tau=30)
        print(f"  {label}: tau={tau} (AMI estimated)")

    if FORCE_M is not None:
        m = int(FORCE_M)
        print(f"  {label}: m={m} (FORCED)")
    else:
        m = cao_embedding_dim(x, tau=tau, max_m=MAX_M)
        print(f"  {label}: m={m} (Cao estimated)")

    theiler = max(1, tau * (m - 1))
    print(f"  {label}: Theiler window={theiler}")
    return tau, m, theiler


# ---------------------------------------------------------------------------
# RQA COMPUTATION  (fixed tau and m passed in)
# ---------------------------------------------------------------------------

def build_recurrence_matrix(x, tau, m, target_rr=TARGET_RR, theiler=None):
    """
    Build the recurrence matrix for a normalised series x with fixed tau, m.
    Returns RM (int8 array, LOI neighbourhood zeroed) or None on failure.
    """
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 20 or np.std(x) < 1e-10:
        return None

    x  = (x - np.mean(x)) / np.std(x)
    N  = len(x) - (m - 1) * tau
    if N < 10:
        return None

    PS = np.array([x[i:i + m * tau:tau] for i in range(N)])
    D  = squareform(pdist(PS, metric='euclidean'))
    n  = len(PS)

    upper_tri = D[np.triu_indices(n, k=1)]
    if len(upper_tri) == 0:
        return None
    eps = np.percentile(upper_tri, target_rr * 100)
    RM  = (D <= eps).astype(np.int8)

    # Theiler window
    if theiler is None:
        theiler = max(1, tau * (m - 1))
    for k in range(-theiler, theiler + 1):
        idx  = np.arange(max(0, k), min(n, n + k))
        jdx  = idx - k
        mask = (idx >= 0) & (idx < n) & (jdx >= 0) & (jdx < n)
        RM[idx[mask], jdx[mask]] = 0

    return RM


def rqa_from_matrix(RM):
    """
    Compute DET and TRANS from a pre-built recurrence matrix.
    Returns (det, trans) or (nan, nan).
    """
    if RM is None:
        return np.nan, np.nan

    total_rec = int(RM.sum())
    if total_rec == 0:
        return 0.0, 0.0

    n = RM.shape[0]

    # Determinism
    diag_pts = 0
    for offset in range(1, n):
        diag = np.diagonal(RM, offset=offset)
        runs = [len(list(g)) for k, g in groupby(diag) if k == 1]
        diag_pts += 2 * sum(r for r in runs if r >= MIN_LINE)
    det = float(np.clip(diag_pts / total_rec, 0.0, 1.0))

    # Transitivity (Donner et al. 2010)
    RM_f  = RM.astype(float)
    k_vec = RM_f.sum(axis=1)
    tri   = np.trace(RM_f @ RM_f @ RM_f)
    denom = np.sum(k_vec * (k_vec - 1.0))
    trans = float(np.clip(tri / denom, 0.0, 1.0)) if denom > 0 else 0.0

    return det, trans


def rqa_window(data, tau, m, theiler):
    """Convenience wrapper: build RM then compute DET/TRANS."""
    RM = build_recurrence_matrix(data, tau, m, theiler=theiler)
    return rqa_from_matrix(RM)


# ---------------------------------------------------------------------------
# WINDOWED ENSEMBLE RQA
# ---------------------------------------------------------------------------

def run_ensemble_rqa(real_df, tau, m, theiler, label='proxy'):
    """
    Windowed RQA across all realisations using fixed global tau, m, theiler.
    """
    ages      = real_df['age'].values
    real_cols = [c for c in real_df.columns if c.startswith('r')]
    if N_REALISATIONS is not None:
        real_cols = real_cols[:N_REALISATIONS]
    n_real    = len(real_cols)
    n_pts     = len(ages)
    n_windows = max(0, (n_pts - WINDOW_SIZE) // WINDOW_STEP + 1)

    det_ens   = np.full((n_windows, n_real), np.nan)
    trans_ens = np.full((n_windows, n_real), np.nan)
    win_ages  = np.full(n_windows, np.nan)

    print(f"\n  Ensemble RQA: {label}  "
          f"({n_real} realisations x {n_windows} windows, "
          f"tau={tau}, m={m})")

    for w in range(n_windows):
        i0 = w * WINDOW_STEP
        i1 = i0 + WINDOW_SIZE
        win_ages[w] = np.mean(ages[i0:i1])
        for r, col in enumerate(real_cols):
            det_ens[w, r], trans_ens[w, r] = \
                rqa_window(real_df[col].values[i0:i1], tau, m, theiler)
        if (w + 1) % 25 == 0 or w == n_windows - 1:
            sys.stdout.write(
                f"\r    {w+1}/{n_windows}  "
                f"({100*(w+1)/n_windows:.0f}%)")
            sys.stdout.flush()
    print()
    return win_ages, det_ens, trans_ens


def ensemble_stats(det_ens, trans_ens):
    with warnings.catch_warnings():
        warnings.simplefilter('ignore')
        dm   = np.nanmedian(det_ens,   axis=1)
        dp5  = np.nanpercentile(det_ens,   5, axis=1)
        dp95 = np.nanpercentile(det_ens,  95, axis=1)
        tm   = np.nanmedian(trans_ens, axis=1)
        tp5  = np.nanpercentile(trans_ens,  5, axis=1)
        tp95 = np.nanpercentile(trans_ens, 95, axis=1)
    return dm, dp5, dp95, tm, tp5, tp95


def interp_smooth(win_ages, vals, target_ages):
    ok = np.isfinite(vals)
    if ok.sum() < 2:
        return np.full_like(target_ages, np.nan, dtype=float)
    f = interp1d(win_ages[ok], vals[ok], kind='linear',
                 bounds_error=False,
                 fill_value=(vals[ok][0], vals[ok][-1]))
    v = f(target_ages)
    return gaussian_filter1d(np.where(np.isfinite(v), v, 0.0),
                             sigma=SMOOTH_SIGMA)


# ---------------------------------------------------------------------------
# DATA LOADING
# ---------------------------------------------------------------------------

def _find_data_start(raw):
    for i, row in raw.iterrows():
        try:
            v = float(row.iloc[0])
            if v == v:  # not NaN
                return i
        except (TypeError, ValueError):
            continue
    raise ValueError("Could not find numeric data block")


def load_drip_summary(xlsx, sheet):
    raw = pd.read_excel(xlsx, sheet_name=sheet, header=None)
    ds  = _find_data_start(raw)
    df  = raw.iloc[ds:, :8].copy().reset_index(drop=True)
    df.columns = ['age','pc05','pc10','pc25','med','pc75','pc90','pc95']
    df = df.apply(pd.to_numeric, errors='coerce')
    return df.dropna(subset=['age']).sort_values('age').reset_index(drop=True)


def load_d18o_summary(xlsx, sheet):
    raw = pd.read_excel(xlsx, sheet_name=sheet, header=None)
    ds  = _find_data_start(raw)
    df  = raw.iloc[ds:, :2].copy().reset_index(drop=True)
    df.columns = ['age', 'median']
    df = df.apply(pd.to_numeric, errors='coerce')
    return df.dropna(subset=['age']).sort_values('age').reset_index(drop=True)


def load_realisations(filepath):
    df = pd.read_csv(filepath)
    return df.sort_values('age').reset_index(drop=True)


# ---------------------------------------------------------------------------
# PDF HEATMAP BUILDER
# ---------------------------------------------------------------------------

def build_pdf_matrix(real_df, age_grid, y_edges, n_real=None):
    real_cols = [c for c in real_df.columns if c.startswith('r')]
    if n_real is not None:
        real_cols = real_cols[:n_real]

    mat = np.full((len(age_grid), len(real_cols)), np.nan)
    for ri, col in enumerate(real_cols):
        f = interp1d(real_df['age'].values, real_df[col].values,
                     kind='linear', bounds_error=False,
                     fill_value=(real_df[col].values[0],
                                 real_df[col].values[-1]))
        mat[:, ri] = f(age_grid)

    pdf = np.zeros((len(age_grid), len(y_edges) - 1))
    for ti in range(len(age_grid)):
        vals = mat[ti, :]
        vals = vals[np.isfinite(vals)]
        if len(vals) > 0:
            h, _ = np.histogram(vals, bins=y_edges)
            mx = h.max()
            pdf[ti, :] = h / mx if mx > 0 else h
    return pdf


# ---------------------------------------------------------------------------
# SUPPLEMENTARY RECURRENCE PLOT FIGURE
# ---------------------------------------------------------------------------

def make_rp_figure(drip_real, d18o_real,
                   drip_tau, drip_m, drip_theiler,
                   d18o_tau, d18o_m, d18o_theiler,
                   age_marker=AGE_MARKER_KA):
    """
    Full-record recurrence plots for median drip rate and median d18O,
    computed on the ensemble median series. Saves SUPP_PDF and SUPP_PNG.
    """
    print("\n  Building supplementary recurrence plot figure ...")

    def median_series(real_df):
        real_cols = [c for c in real_df.columns if c.startswith('r')]
        ages = real_df['age'].values
        med  = np.nanmedian(real_df[real_cols].values, axis=1)
        return ages, med

    drip_ages, drip_med = median_series(drip_real)
    d18o_ages, d18o_med = median_series(d18o_real)

    # Build full-record RMs from median series
    drip_RM = build_recurrence_matrix(drip_med, drip_tau, drip_m,
                                      theiler=drip_theiler)
    d18o_RM = build_recurrence_matrix(d18o_med, d18o_tau, d18o_m,
                                      theiler=d18o_theiler)

    fig, axes = plt.subplots(1, 2, figsize=(8.5, 4.0),
                             gridspec_kw={'wspace': 0.35,
                                          'left': 0.10, 'right': 0.96,
                                          'top': 0.90, 'bottom': 0.13})

    titles = [r'$\delta^{18}$O', 'Drip rate']
    RMs    = [d18o_RM, drip_RM]
    ages_l = [d18o_ages, drip_ages]

    for ax, RM, ages, title in zip(axes, RMs, ages_l, titles):
        if RM is None:
            ax.text(0.5, 0.5, 'Insufficient data',
                    ha='center', va='center', transform=ax.transAxes)
            continue

        n = RM.shape[0]
        # Map matrix indices to ages
        # The series is sorted ascending (youngest→oldest index),
        # so we flip for display (oldest at left/bottom)
        age_min = ages.min()
        age_max = ages.max()

        ax.imshow(RM, cmap='binary', origin='lower', aspect='auto',
                  extent=[age_min, age_max, age_min, age_max],
                  interpolation='none')

        # 5.2 ka annotation lines
        ax.axvline(age_marker, color='red', lw=0.8, ls='--', alpha=0.8)
        ax.axhline(age_marker, color='red', lw=0.8, ls='--', alpha=0.8)
        ax.text(age_marker + 100, age_max * 0.97,
                f'{age_marker//1000:.1f} ka',
                color='red', fontsize=6, va='top')

        ax.set_xlabel('Age (yr BP)', fontsize=8, fontweight='bold')
        ax.set_ylabel('Age (yr BP)', fontsize=8, fontweight='bold')
        ax.set_title(title, fontsize=9, fontweight='bold', pad=6)
        ax.tick_params(axis='both', which='major',
                       width=0.5, length=2.5, direction='in')

        # RQA stats for annotation
        det, trans = rqa_from_matrix(RM)
        ax.text(0.02, 0.98,
                f'DET={det:.3f}   TRANS={trans:.3f}',
                transform=ax.transAxes, fontsize=6,
                va='top', ha='left',
                bbox=dict(boxstyle='round,pad=0.3', fc='white',
                          ec='grey', alpha=0.8))

    fig.suptitle('Recurrence plots — HS4 stalagmite, Heshang Cave\n'
                 r'RR = 5%, fixed $\tau$ and $m$ from full median series',
                 fontsize=8, y=0.98)

    for path in (SUPP_PDF, SUPP_PNG):
        fig.savefig(path, dpi=300, bbox_inches='tight')
        print(f"  Saved -> {path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# PLOTTING HELPERS (shared with plot script)
# ---------------------------------------------------------------------------

def draw_sustained_markers(plot_age, det_arr, ax_rqa, ax_pdf,
                            threshold=DET_THRESHOLD,
                            min_sustain=MIN_SUSTAIN):
    below    = (det_arr < threshold).astype(int)
    in_run   = False
    run_start = None
    for k in range(len(below)):
        if below[k] == 1 and not in_run:
            run_start = k
            in_run    = True
        elif below[k] == 0 and in_run:
            if (k - run_start) >= min_sustain:
                age_c = plot_age[run_start]
                ax_rqa.axvline(age_c, color='black',
                               lw=0.8, ls='--', alpha=0.7)
                ax_pdf.axvline(age_c, color='white',
                               lw=0.8, ls='--', alpha=0.7)
            in_run = False


# ---------------------------------------------------------------------------
# MAIN FIGURE
# ---------------------------------------------------------------------------

def make_main_figure(plot_age, age_hi, age_lo,
                     d18o_det_p, d18o_det_lo, d18o_det_hi,
                     d18o_tran_p, d18o_tran_lo, d18o_tran_hi,
                     drip_det_p, drip_det_lo, drip_det_hi,
                     drip_tran_p, drip_tran_lo, drip_tran_hi,
                     d18o_pdf, d18o_med_grid, y_lo, y_hi,
                     drip_pdf, drip_med_grid,
                     drip_pc25_grid, drip_pc75_grid):

    def R(a): return a[::-1]

    fig = plt.figure(figsize=(8.5, 7.0))
    gs  = gridspec.GridSpec(4, 1, figure=fig, hspace=0.22,
                            left=0.12, right=0.88, top=0.96, bottom=0.09)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax3 = fig.add_subplot(gs[2], sharex=ax1)
    ax4 = fig.add_subplot(gs[3], sharex=ax1)

    for ax in (ax1, ax2, ax3, ax4):
        ax.spines['right'].set_visible(True)
        ax.spines['right'].set_linewidth(0.5)
        ax.tick_params(axis='both', which='major',
                       width=0.5, length=2.5, direction='in', right=True)

    # Panel A: RQA d18O
    ax1.fill_between(plot_age, R(d18o_det_lo), R(d18o_det_hi),
                     color='dimgrey', alpha=0.22, lw=0, label='DET 5-95%ile')
    ax1.fill_between(plot_age, R(d18o_tran_lo), R(d18o_tran_hi),
                     color='darkgreen', alpha=0.18, lw=0,
                     label='TRANS 5-95%ile')
    ax1.plot(plot_age, R(d18o_det_p),
             color='dimgrey',   lw=0.9, label='Determinism (DET)')
    ax1.plot(plot_age, R(d18o_tran_p),
             color='darkgreen', lw=0.9, label='Transitivity (TRANS)')
    ax1.axhline(DET_THRESHOLD, color='black', lw=0.5, ls=':', alpha=0.6)
    ax1.set_ylabel(r'RQA ($\delta^{18}$O)', labelpad=8)
    ax1.set_ylim(0, 1)
    ax1.legend(loc='lower right', frameon=True, framealpha=0.9,
               ncol=2, fontsize=6)
    ax1.text(0.01, 0.97, 'A', transform=ax1.transAxes,
             fontsize=10, fontweight='bold', va='top')

    # Panel B: d18O PDF
    img_ext_d18o = [plot_age[0], plot_age[-1], y_lo, y_hi]
    im2 = ax2.imshow(R(d18o_pdf).T,
                     extent=img_ext_d18o, aspect='auto',
                     cmap='cividis', origin='lower', vmin=0, vmax=1)
    ax2.plot(plot_age, R(d18o_med_grid), color='red', lw=0.8, label='Median')
    ax2.set_ylabel(r'$\delta^{18}$O (‰ VPDB)', labelpad=8)
    ax2.set_ylim(y_hi, y_lo)   # inverted: more negative at top
    ax2.legend(loc='upper right', frameon=True, framealpha=0.9, fontsize=6)
    ax2.text(0.01, 0.97, 'B', transform=ax2.transAxes,
             fontsize=10, fontweight='bold', va='top')

    # Panel C: RQA drip rate
    ax3.fill_between(plot_age, R(drip_det_lo), R(drip_det_hi),
                     color='dimgrey', alpha=0.22, lw=0, label='DET 5-95%ile')
    ax3.fill_between(plot_age, R(drip_tran_lo), R(drip_tran_hi),
                     color='darkgreen', alpha=0.18, lw=0,
                     label='TRANS 5-95%ile')
    ax3.plot(plot_age, R(drip_det_p),
             color='dimgrey',   lw=0.9, label='Determinism (DET)')
    ax3.plot(plot_age, R(drip_tran_p),
             color='darkgreen', lw=0.9, label='Transitivity (TRANS)')
    ax3.axhline(DET_THRESHOLD, color='black', lw=0.5, ls=':', alpha=0.6)
    ax3.set_ylabel('RQA (drip rate)', labelpad=8)
    ax3.set_ylim(0, 1)
    ax3.legend(loc='lower right', frameon=True, framealpha=0.9,
               ncol=2, fontsize=6)
    ax3.text(0.01, 0.97, 'C', transform=ax3.transAxes,
             fontsize=10, fontweight='bold', va='top')

    # Panel D: Drip rate PDF
    img_ext_drip = [plot_age[0], plot_age[-1], 0, 40]
    im4 = ax4.imshow(R(drip_pdf).T,
                     extent=img_ext_drip, aspect='auto',
                     cmap='cividis', origin='lower', vmin=0, vmax=1)
    ax4.plot(plot_age, R(drip_med_grid),
             color='red', lw=0.8, label='Median')
    ax4.plot(plot_age, R(drip_pc25_grid),
             color='red', ls='--', lw=0.4)
    ax4.plot(plot_age, R(drip_pc75_grid),
             color='red', ls='--', lw=0.4)
    ax4.set_ylabel('Drip rate (min$^{-1}$)', labelpad=8)
    ax4.set_ylim(0, 40)
    ax4.set_xlabel('Age (yr BP)')
    ax4.legend(loc='upper right', frameon=True, framealpha=0.9, fontsize=6)
    ax4.text(0.01, 0.97, 'D', transform=ax4.transAxes,
             fontsize=10, fontweight='bold', va='top')

    # Sustained DET transition markers
    draw_sustained_markers(plot_age, R(d18o_det_p), ax1, ax2)
    draw_sustained_markers(plot_age, R(drip_det_p), ax3, ax4)

    # Shared x-axis
    for ax in (ax1, ax2, ax3, ax4):
        ax.set_xlim(age_hi, age_lo)
    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)

    # Shared colorbar
    cbar_ax = fig.add_axes([0.905, 0.09, 0.013, 0.87])
    cbar    = fig.colorbar(im4, cax=cbar_ax, ticks=[0, 0.5, 1])
    cbar.set_label('Normalised PDF', fontsize=6, labelpad=4)
    cbar.ax.tick_params(labelsize=6)

    for path in (OUTPUT_PDF, OUTPUT_PNG, OUTPUT_EPS):
        fig.savefig(path, dpi=300, bbox_inches='tight')
        print(f"  Saved -> {path}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# MAIN
# ---------------------------------------------------------------------------

def main():
    print("=" * 65)
    print(" HS4 Ensemble RQA  —  Fig. 5")
    print("=" * 65)

    # Load data
    print(f"\nLoading {DRIP_REAL_FILE}")
    drip_real = load_realisations(DRIP_REAL_FILE)
    print(f"  {drip_real.shape[0]} time steps, "
          f"{drip_real.shape[1]-1} realisations | "
          f"age {drip_real['age'].min():.0f}-{drip_real['age'].max():.0f} BP")

    print(f"Loading {D18O_REAL_FILE}")
    d18o_real = load_realisations(D18O_REAL_FILE)
    print(f"  {d18o_real.shape[0]} time steps, "
          f"{d18o_real.shape[1]-1} realisations | "
          f"age {d18o_real['age'].min():.0f}-{d18o_real['age'].max():.0f} BP")

    print(f"\nLoading summary from {DRIP_XLSX}")
    drip_sum = load_drip_summary(DRIP_XLSX, DRIP_SHEET)
    d18o_sum = load_d18o_summary(DRIP_XLSX, D18O_SHEET)

    # Common age range
    age_lo = max(drip_real['age'].min(), d18o_real['age'].min(),
                 drip_sum['age'].min(), d18o_sum['age'].min())
    age_hi = min(drip_real['age'].max(), d18o_real['age'].max(),
                 drip_sum['age'].max(), d18o_sum['age'].max())

    def trim(df):
        return df[(df['age'] >= age_lo) & (df['age'] <= age_hi)]\
               .reset_index(drop=True)

    drip_real = trim(drip_real)
    d18o_real = trim(d18o_real)
    drip_sum  = trim(drip_sum)
    d18o_sum  = trim(d18o_sum)
    print(f"  Common age range: {age_lo:.0f}-{age_hi:.0f} yr BP")

    # ── Global tau and m from median series ───────────────────────────────
    print("\nEstimating global embedding parameters from median series ...")
    drip_cols = [c for c in drip_real.columns if c.startswith('r')]
    d18o_cols = [c for c in d18o_real.columns if c.startswith('r')]

    drip_median_full = np.nanmedian(drip_real[drip_cols].values, axis=1)
    d18o_median_full = np.nanmedian(d18o_real[d18o_cols].values, axis=1)

    drip_tau, drip_m, drip_theiler = estimate_global_params(
        drip_median_full, label='drip rate')
    d18o_tau, d18o_m, d18o_theiler = estimate_global_params(
        d18o_median_full, label='d18O')

    # ── Ensemble RQA ──────────────────────────────────────────────────────
    drip_wages, drip_det_e, drip_trans_e = run_ensemble_rqa(
        drip_real, drip_tau, drip_m, drip_theiler, label='drip rate')
    d18o_wages, d18o_det_e, d18o_trans_e = run_ensemble_rqa(
        d18o_real, d18o_tau, d18o_m, d18o_theiler, label='d18O')

    # ── Ensemble statistics ───────────────────────────────────────────────
    drip_dm, drip_dp5, drip_dp95, drip_tm, drip_tp5, drip_tp95 = \
        ensemble_stats(drip_det_e, drip_trans_e)
    d18o_dm, d18o_dp5, d18o_dp95, d18o_tm, d18o_tp5, d18o_tp95 = \
        ensemble_stats(d18o_det_e, d18o_trans_e)

    print(f"\n  Drip DET  range: "
          f"{np.nanmin(drip_dm):.3f}-{np.nanmax(drip_dm):.3f}")
    print(f"  d18O DET  range: "
          f"{np.nanmin(d18o_dm):.3f}-{np.nanmax(d18o_dm):.3f}")

    # ── Save CSV ──────────────────────────────────────────────────────────
    rqa_df = pd.DataFrame({
        'age_window'      : drip_wages,
        'drip_DET_med'    : drip_dm,   'drip_DET_p5'   : drip_dp5,
        'drip_DET_p95'    : drip_dp95, 'drip_TRANS_med': drip_tm,
        'drip_TRANS_p5'   : drip_tp5,  'drip_TRANS_p95': drip_tp95,
        'd18O_DET_med'    : d18o_dm,   'd18O_DET_p5'   : d18o_dp5,
        'd18O_DET_p95'    : d18o_dp95, 'd18O_TRANS_med': d18o_tm,
        'd18O_TRANS_p5'   : d18o_tp5,  'd18O_TRANS_p95': d18o_tp95,
    })
    rqa_df.to_csv(OUTPUT_CSV, index=False)

    # Save parameters as a small metadata CSV alongside
    params_df = pd.DataFrame({
        'proxy'  : ['drip_rate', 'd18O'],
        'tau'    : [drip_tau,    d18o_tau],
        'm'      : [drip_m,      d18o_m],
        'theiler': [drip_theiler, d18o_theiler],
        'RR'     : [TARGET_RR,   TARGET_RR],
        'l_min'  : [MIN_LINE,    MIN_LINE],
        'window' : [WINDOW_SIZE, WINDOW_SIZE],
        'step'   : [WINDOW_STEP, WINDOW_STEP],
    })
    params_df.to_csv('rqa_parameters.csv', index=False)
    print(f"  RQA statistics -> {OUTPUT_CSV}")
    print(f"  RQA parameters -> rqa_parameters.csv")

    # ── Interpolate + smooth onto summary age grid ────────────────────────
    age_grid = drip_sum['age'].values

    def IS(wages, vals):
        return interp_smooth(wages, vals, age_grid)

    drip_det_p   = IS(drip_wages, drip_dm)
    drip_det_lo  = IS(drip_wages, drip_dp5)
    drip_det_hi  = IS(drip_wages, drip_dp95)
    drip_tran_p  = IS(drip_wages, drip_tm)
    drip_tran_lo = IS(drip_wages, drip_tp5)
    drip_tran_hi = IS(drip_wages, drip_tp95)

    d18o_det_p   = IS(d18o_wages, d18o_dm)
    d18o_det_lo  = IS(d18o_wages, d18o_dp5)
    d18o_det_hi  = IS(d18o_wages, d18o_dp95)
    d18o_tran_p  = IS(d18o_wages, d18o_tm)
    d18o_tran_lo = IS(d18o_wages, d18o_tp5)
    d18o_tran_hi = IS(d18o_wages, d18o_tp95)

    # ── PDF heatmaps ──────────────────────────────────────────────────────
    print("\nBuilding PDF heatmaps ...")

    drip_y_edges = np.linspace(0, 40, N_YBINS + 1)
    drip_pdf     = build_pdf_matrix(drip_real, age_grid, drip_y_edges,
                                    n_real=N_REALISATIONS)

    f_med  = interp1d(drip_sum['age'].values, drip_sum['med'].values,
                      kind='linear', bounds_error=False,
                      fill_value=(drip_sum['med'].values[0],
                                  drip_sum['med'].values[-1]))
    f_pc25 = interp1d(drip_sum['age'].values, drip_sum['pc25'].values,
                      kind='linear', bounds_error=False,
                      fill_value=(drip_sum['pc25'].values[0],
                                  drip_sum['pc25'].values[-1]))
    f_pc75 = interp1d(drip_sum['age'].values, drip_sum['pc75'].values,
                      kind='linear', bounds_error=False,
                      fill_value=(drip_sum['pc75'].values[0],
                                  drip_sum['pc75'].values[-1]))
    drip_med_grid  = f_med(age_grid)
    drip_pc25_grid = f_pc25(age_grid)
    drip_pc75_grid = f_pc75(age_grid)

    # d18O y range from realisations
    d18o_real_cols = [c for c in d18o_real.columns if c.startswith('r')]
    sample_vals    = d18o_real[d18o_real_cols].values.ravel()
    sample_vals    = sample_vals[np.isfinite(sample_vals)]
    y_lo = np.floor(np.percentile(sample_vals, 0.5) * 2) / 2
    y_hi = np.ceil( np.percentile(sample_vals, 99.5) * 2) / 2

    d18o_y_edges = np.linspace(y_lo, y_hi, N_YBINS + 1)
    d18o_pdf     = build_pdf_matrix(d18o_real, age_grid, d18o_y_edges,
                                    n_real=N_REALISATIONS)

    f_d18o_med = interp1d(d18o_sum['age'].values, d18o_sum['median'].values,
                          kind='linear', bounds_error=False,
                          fill_value=(d18o_sum['median'].values[0],
                                      d18o_sum['median'].values[-1]))
    d18o_med_grid = f_d18o_med(age_grid)

    # ── Plot arrays ───────────────────────────────────────────────────────
    plot_age = age_grid[::-1]   # oldest LEFT

    # ── Main figure ───────────────────────────────────────────────────────
    print("\nRendering main figure ...")
    make_main_figure(
        plot_age, age_hi, age_lo,
        d18o_det_p, d18o_det_lo, d18o_det_hi,
        d18o_tran_p, d18o_tran_lo, d18o_tran_hi,
        drip_det_p, drip_det_lo, drip_det_hi,
        drip_tran_p, drip_tran_lo, drip_tran_hi,
        d18o_pdf, d18o_med_grid, y_lo, y_hi,
        drip_pdf, drip_med_grid, drip_pc25_grid, drip_pc75_grid,
    )

    # ── Supplementary RP figure ───────────────────────────────────────────
    make_rp_figure(
        drip_real, d18o_real,
        drip_tau, drip_m, drip_theiler,
        d18o_tau, d18o_m, d18o_theiler,
    )

    print("\nComplete.\n")


if __name__ == '__main__':
    main()