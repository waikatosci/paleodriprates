"""
FigS_RQA_HS4.py
================
Supplementary / Extended Data figure: Recurrence Quantification Analysis
for the HS4 stalagmite drip rate and δ¹⁸O records.

Four panels (equal height):
  A  RQA metrics (DET, TRANS) for δ¹⁸O
  B  δ¹⁸O PDF heatmap
  C  RQA metrics (DET, TRANS) for drip rate
  D  Drip rate PDF heatmap

Plus a supplementary recurrence-plot figure (2 panels).

Adapted from the original Fig5_RQA_HS4_plot.py — all data loading and
computation logic preserved; only the role changes from main figure to
supplementary.

Required input files:
    rqa_ensemble_results.csv    — output of Fig5_RQA_HS4_ensemble.py
    drip_rate_realisations.csv  — output of Drip_rate_parallel_fr.py
    d18O_realisations.csv       — output of Drip_rate_parallel_fr.py
    Drip_rate.xlsx              — sheets 5.OutDripRate, 6.OutIsotope

Outputs:
    FigS_RQA_HS4.pdf / .png / .eps
    FigS_RP_HS4.pdf  / .png
"""

import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter1d
from scipy.spatial.distance import pdist, squareform
from itertools import groupby

warnings.filterwarnings('ignore')

# ═══════════════════════════════════════════════════════════════════════
# CONFIGURATION
# ═══════════════════════════════════════════════════════════════════════
RQA_CSV        = 'rqa_ensemble_results.csv'
DRIP_REAL_FILE = 'drip_rate_realisations.csv'
D18O_REAL_FILE = 'd18O_realisations.csv'
DRIP_XLSX      = 'Drip_rate.xlsx'
DRIP_SHEET     = '5.OutDripRate'
D18O_SHEET     = '6.OutIsotope'

OUTPUT_PDF     = 'FigS_RQA_HS4.pdf'
OUTPUT_PNG     = 'FigS_RQA_HS4.png'
OUTPUT_EPS     = 'FigS_RQA_HS4.eps'
SUPP_PDF       = 'FigS_RP_HS4.pdf'
SUPP_PNG       = 'FigS_RP_HS4.png'
PARAMS_CSV     = 'rqa_parameters.csv'

SMOOTH_SIGMA     = 2.0
TRANS_THRESHOLD  = 0.5
N_YBINS          = 300
N_REALISATIONS   = None

# Event windows
EVENTS = {
    '8.2 ka':     {'centre': 8200, 'hw': 250, 'color': '#4A90D9', 'alpha': 0.10},
    '5.2 ka':     {'centre': 5200, 'hw': 200, 'color': '#E07B39', 'alpha': 0.10},
    'Post-1800':  {'centre':  65,  'hw':  65, 'color': '#C0392B', 'alpha': 0.10},
}

# ═══════════════════════════════════════════════════════════════════════
# STYLE
# ═══════════════════════════════════════════════════════════════════════
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


# ═══════════════════════════════════════════════════════════════════════
# HELPERS (preserved from original)
# ═══════════════════════════════════════════════════════════════════════

def interp_smooth(src_ages, vals, target_ages):
    ok = np.isfinite(vals)
    if ok.sum() < 2:
        return np.full_like(target_ages, np.nan, dtype=float)
    f = interp1d(src_ages[ok], vals[ok], kind='linear',
                 bounds_error=False,
                 fill_value=(vals[ok][0], vals[ok][-1]))
    v = f(target_ages)
    return gaussian_filter1d(np.where(np.isfinite(v), v, 0.0),
                             sigma=SMOOTH_SIGMA)


def _find_data_start(raw):
    for i, row in raw.iterrows():
        try:
            v = float(row.iloc[0])
            if not (v != v):
                return i
        except (TypeError, ValueError):
            continue
    raise ValueError("Could not find numeric data block")


def load_drip_summary(xlsx, sheet):
    raw = pd.read_excel(xlsx, sheet_name=sheet, header=None)
    ds = _find_data_start(raw)
    df = raw.iloc[ds:, :8].copy().reset_index(drop=True)
    df.columns = ['age', 'pc05', 'pc10', 'pc25', 'med', 'pc75', 'pc90', 'pc95']
    df = df.apply(pd.to_numeric, errors='coerce')
    return df.dropna(subset=['age']).sort_values('age').reset_index(drop=True)


def load_d18o_summary(xlsx, sheet):
    raw = pd.read_excel(xlsx, sheet_name=sheet, header=None)
    ds = _find_data_start(raw)
    df = raw.iloc[ds:, :2].copy().reset_index(drop=True)
    df.columns = ['age', 'median']
    df = df.apply(pd.to_numeric, errors='coerce')
    return df.dropna(subset=['age']).sort_values('age').reset_index(drop=True)


def build_pdf_matrix(real_df, age_grid, y_edges, n_real=None):
    real_cols = [c for c in real_df.columns if c.startswith('r')]
    if n_real is not None:
        real_cols = real_cols[:n_real]
    mat = np.full((len(age_grid), len(real_cols)), np.nan)
    for ri, col in enumerate(real_cols):
        f = interp1d(real_df['age'].values, real_df[col].values,
                     kind='linear', bounds_error=False,
                     fill_value=(real_df[col].values[0], real_df[col].values[-1]))
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


def draw_trans_excursion_markers(plot_age, trans_arr, ax_rqa, ax_pdf,
                                  threshold=TRANS_THRESHOLD):
    above = trans_arr > threshold
    in_run = False
    run_start = None
    for k in range(len(above)):
        if above[k] and not in_run:
            run_start = k
            in_run = True
        elif not above[k] and in_run:
            seg = trans_arr[run_start:k]
            peak = run_start + int(np.argmax(seg))
            age_c = plot_age[peak]
            ax_rqa.axvline(age_c, color='black', lw=0.8, ls='--', alpha=0.7)
            ax_pdf.axvline(age_c, color='white', lw=0.8, ls='--', alpha=0.7)
            in_run = False
    if in_run:
        seg = trans_arr[run_start:]
        peak = run_start + int(np.argmax(seg))
        age_c = plot_age[peak]
        ax_rqa.axvline(age_c, color='black', lw=0.8, ls='--', alpha=0.7)
        ax_pdf.axvline(age_c, color='white', lw=0.8, ls='--', alpha=0.7)


def draw_event_bands(axes, events):
    for name, ev in events.items():
        lo = ev['centre'] - ev['hw']
        hi = ev['centre'] + ev['hw']
        for ax in axes:
            ax.axvspan(lo, hi, color=ev['color'], alpha=ev['alpha'],
                       zorder=0, lw=0)


# ═══════════════════════════════════════════════════════════════════════
# RECURRENCE PLOT HELPERS
# ═══════════════════════════════════════════════════════════════════════

def build_recurrence_matrix(x, tau, m, target_rr=0.05, theiler=None):
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    if len(x) < 20 or np.std(x) < 1e-10:
        return None
    x = (x - np.mean(x)) / np.std(x)
    N = len(x) - (m - 1) * tau
    if N < 10:
        return None
    PS = np.array([x[i:i + m * tau:tau] for i in range(N)])
    D = squareform(pdist(PS, metric='euclidean'))
    n = len(PS)
    upper_tri = D[np.triu_indices(n, k=1)]
    if len(upper_tri) == 0:
        return None
    eps = np.percentile(upper_tri, target_rr * 100)
    RM = (D <= eps).astype(np.int8)
    if theiler is None:
        theiler = max(1, tau * (m - 1))
    for k in range(-theiler, theiler + 1):
        idx = np.arange(max(0, k), min(n, n + k))
        jdx = idx - k
        mask = (idx >= 0) & (idx < n) & (jdx >= 0) & (jdx < n)
        RM[idx[mask], jdx[mask]] = 0
    return RM


def rqa_from_matrix(RM, min_line=2):
    if RM is None:
        return np.nan, np.nan
    total_rec = int(RM.sum())
    if total_rec == 0:
        return 0.0, 0.0
    n = RM.shape[0]
    diag_pts = 0
    for offset in range(1, n):
        diag = np.diagonal(RM, offset=offset)
        runs = [len(list(g)) for k, g in groupby(diag) if k == 1]
        diag_pts += 2 * sum(r for r in runs if r >= min_line)
    det = float(np.clip(diag_pts / total_rec, 0.0, 1.0))
    RM_f = RM.astype(float)
    k_vec = RM_f.sum(axis=1)
    tri = np.trace(RM_f @ RM_f @ RM_f)
    denom = np.sum(k_vec * (k_vec - 1.0))
    trans = float(np.clip(tri / denom, 0.0, 1.0)) if denom > 0 else 0.0
    return det, trans


# ═══════════════════════════════════════════════════════════════════════
# MAIN: 4-PANEL RQA FIGURE
# ═══════════════════════════════════════════════════════════════════════

def main():
    print("=" * 60)
    print("  Supplementary RQA Figure — HS4")
    print("=" * 60)

    # ── Load data ─────────────────────────────────────────────────────
    print("\nLoading data ...")
    rqa = pd.read_csv(RQA_CSV)

    try:
        params = pd.read_csv(PARAMS_CSV)
        drip_row = params[params['proxy'] == 'drip_rate'].iloc[0]
        d18o_row = params[params['proxy'] == 'd18O'].iloc[0]
        print(f"  Drip: tau={int(drip_row['tau'])}, m={int(drip_row['m'])}")
        print(f"  d18O: tau={int(d18o_row['tau'])}, m={int(d18o_row['m'])}")
    except Exception as e:
        print(f"  Warning: could not load {PARAMS_CSV} ({e})")
        params = None

    win_ages = rqa['age_window'].values

    drip_sum = load_drip_summary(DRIP_XLSX, DRIP_SHEET)
    d18o_sum = load_d18o_summary(DRIP_XLSX, D18O_SHEET)
    drip_real = pd.read_csv(DRIP_REAL_FILE).sort_values('age').reset_index(drop=True)
    d18o_real = pd.read_csv(D18O_REAL_FILE).sort_values('age').reset_index(drop=True)

    # Common age range
    age_lo = max(win_ages.min(), drip_real['age'].min(), d18o_real['age'].min(),
                 drip_sum['age'].min(), d18o_sum['age'].min())
    age_hi = min(win_ages.max(), drip_real['age'].max(), d18o_real['age'].max(),
                 drip_sum['age'].max(), d18o_sum['age'].max())

    def trim(df):
        return df[(df['age'] >= age_lo) & (df['age'] <= age_hi)].reset_index(drop=True)

    drip_real = trim(drip_real)
    d18o_real = trim(d18o_real)
    drip_sum = trim(drip_sum)
    d18o_sum = trim(d18o_sum)
    age_grid = drip_sum['age'].values

    # ── RQA curves ────────────────────────────────────────────────────
    def IS(vals):
        return interp_smooth(win_ages, vals, age_grid)

    drip_det_p   = IS(rqa['drip_DET_med'].values)
    drip_det_lo  = IS(rqa['drip_DET_p5'].values)
    drip_det_hi  = IS(rqa['drip_DET_p95'].values)
    drip_tran_p  = IS(rqa['drip_TRANS_med'].values)
    drip_tran_lo = IS(rqa['drip_TRANS_p5'].values)
    drip_tran_hi = IS(rqa['drip_TRANS_p95'].values)

    d18o_det_p   = IS(rqa['d18O_DET_med'].values)
    d18o_det_lo  = IS(rqa['d18O_DET_p5'].values)
    d18o_det_hi  = IS(rqa['d18O_DET_p95'].values)
    d18o_tran_p  = IS(rqa['d18O_TRANS_med'].values)
    d18o_tran_lo = IS(rqa['d18O_TRANS_p5'].values)
    d18o_tran_hi = IS(rqa['d18O_TRANS_p95'].values)

    # ── PDF heatmaps ──────────────────────────────────────────────────
    print("Building PDF heatmaps ...")

    drip_y_edges = np.linspace(0, 40, N_YBINS + 1)
    drip_pdf = build_pdf_matrix(drip_real, age_grid, drip_y_edges, N_REALISATIONS)

    f_med  = interp1d(drip_sum['age'], drip_sum['med'], bounds_error=False,
                      fill_value=(drip_sum['med'].values[0], drip_sum['med'].values[-1]))
    f_pc25 = interp1d(drip_sum['age'], drip_sum['pc25'], bounds_error=False,
                      fill_value=(drip_sum['pc25'].values[0], drip_sum['pc25'].values[-1]))
    f_pc75 = interp1d(drip_sum['age'], drip_sum['pc75'], bounds_error=False,
                      fill_value=(drip_sum['pc75'].values[0], drip_sum['pc75'].values[-1]))
    drip_med_grid  = f_med(age_grid)
    drip_pc25_grid = f_pc25(age_grid)
    drip_pc75_grid = f_pc75(age_grid)

    d18o_real_cols = [c for c in d18o_real.columns if c.startswith('r')]
    if N_REALISATIONS is not None:
        d18o_real_cols = d18o_real_cols[:N_REALISATIONS]
    sample_vals = d18o_real[d18o_real_cols].values.ravel()
    sample_vals = sample_vals[np.isfinite(sample_vals)]
    y_lo = np.floor(np.percentile(sample_vals, 0.5) * 2) / 2
    y_hi = np.ceil(np.percentile(sample_vals, 99.5) * 2) / 2

    d18o_y_edges = np.linspace(y_lo, y_hi, N_YBINS + 1)
    d18o_pdf = build_pdf_matrix(d18o_real, age_grid, d18o_y_edges, N_REALISATIONS)

    f_d18o_med = interp1d(d18o_sum['age'], d18o_sum['median'], bounds_error=False,
                          fill_value=(d18o_sum['median'].values[0],
                                      d18o_sum['median'].values[-1]))
    d18o_med_grid = f_d18o_med(age_grid)

    # ── Reverse for plotting ──────────────────────────────────────────
    def R(a): return a[::-1]
    plot_age = age_grid[::-1]

    img_ext_drip = [plot_age[0], plot_age[-1], 0, 40]
    img_ext_d18o = [plot_age[0], plot_age[-1], y_lo, y_hi]

    # ── FIGURE ────────────────────────────────────────────────────────
    print("Rendering supplementary RQA figure ...")
    fig = plt.figure(figsize=(8.5, 7.5))
    gs = gridspec.GridSpec(4, 1, figure=fig, hspace=0.22,
                           left=0.12, right=0.88, top=0.96, bottom=0.08)
    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1], sharex=ax1)
    ax3 = fig.add_subplot(gs[2], sharex=ax1)
    ax4 = fig.add_subplot(gs[3], sharex=ax1)
    all_axes = [ax1, ax2, ax3, ax4]

    for ax in all_axes:
        ax.spines['right'].set_visible(True)
        ax.spines['right'].set_linewidth(0.5)
        ax.tick_params(axis='both', which='major',
                       width=0.5, length=2.5, direction='in', right=True)

    # Event bands
    draw_event_bands(all_axes, EVENTS)

    # Panel A: RQA δ¹⁸O
    ax1.fill_between(plot_age, R(d18o_det_lo), R(d18o_det_hi),
                     color='dimgrey', alpha=0.22, lw=0, label='DET 5–95%ile')
    ax1.fill_between(plot_age, R(d18o_tran_lo), R(d18o_tran_hi),
                     color='darkgreen', alpha=0.18, lw=0, label='TRANS 5–95%ile')
    ax1.plot(plot_age, R(d18o_det_p),
             color='dimgrey', lw=0.9, label='Determinism (DET)')
    ax1.plot(plot_age, R(d18o_tran_p),
             color='darkgreen', lw=0.9, label='Transitivity (TRANS)')
    ax1.axhline(TRANS_THRESHOLD, color='darkgreen', lw=0.5, ls=':', alpha=0.6)
    ax1.set_ylabel(r'RQA ($\delta^{18}$O)', labelpad=8)
    ax1.set_ylim(0, 1)
    ax1.legend(loc='lower right', frameon=True, framealpha=0.9, ncol=2, fontsize=6)
    ax1.text(0.01, 0.97, 'a', transform=ax1.transAxes,
             fontsize=10, fontweight='bold', va='top')

    # Panel B: δ¹⁸O PDF
    im2 = ax2.imshow(R(d18o_pdf).T, extent=img_ext_d18o, aspect='auto',
                     cmap='Blues', origin='lower', vmin=0, vmax=1)
    ax2.plot(plot_age, R(d18o_med_grid), color='red', lw=0.8, label='Median')
    ax2.set_ylabel(r'$\delta^{18}$O (‰ VPDB)', labelpad=8)
    ax2.set_ylim(y_hi, y_lo)
    ax2.legend(loc='upper right', frameon=True, framealpha=0.9, fontsize=6)
    ax2.text(0.01, 0.97, 'b', transform=ax2.transAxes,
             fontsize=10, fontweight='bold', va='top')

    # Panel C: RQA drip rate
    ax3.fill_between(plot_age, R(drip_det_lo), R(drip_det_hi),
                     color='dimgrey', alpha=0.22, lw=0, label='DET 5–95%ile')
    ax3.fill_between(plot_age, R(drip_tran_lo), R(drip_tran_hi),
                     color='darkgreen', alpha=0.18, lw=0, label='TRANS 5–95%ile')
    ax3.plot(plot_age, R(drip_det_p),
             color='dimgrey', lw=0.9, label='Determinism (DET)')
    ax3.plot(plot_age, R(drip_tran_p),
             color='darkgreen', lw=0.9, label='Transitivity (TRANS)')
    ax3.axhline(TRANS_THRESHOLD, color='darkgreen', lw=0.5, ls=':', alpha=0.6)
    ax3.set_ylabel('RQA (drip rate)', labelpad=8)
    ax3.set_ylim(0, 1)
    ax3.legend(loc='lower right', frameon=True, framealpha=0.9, ncol=2, fontsize=6)
    ax3.text(0.01, 0.97, 'c', transform=ax3.transAxes,
             fontsize=10, fontweight='bold', va='top')

    # Panel D: Drip rate PDF
    im4 = ax4.imshow(R(drip_pdf).T, extent=img_ext_drip, aspect='auto',
                     cmap='YlGn', origin='lower', vmin=0, vmax=1)
    ax4.plot(plot_age, R(drip_med_grid), color='#1A1A1A', lw=0.9, label='Median')
    ax4.plot(plot_age, R(drip_pc25_grid), color='#2C3E50', ls='--', lw=0.4)
    ax4.plot(plot_age, R(drip_pc75_grid), color='#2C3E50', ls='--', lw=0.4)
    ax4.set_ylabel('Drip rate (min$^{-1}$)', labelpad=8)
    ax4.set_ylim(0, 40)
    ax4.set_xlabel('Age (yr BP)')
    ax4.legend(loc='upper right', frameon=True, framealpha=0.9, fontsize=6)
    ax4.text(0.01, 0.97, 'd', transform=ax4.transAxes,
             fontsize=10, fontweight='bold', va='top')

    # TRANS excursion markers — drip rate panels only
    draw_trans_excursion_markers(plot_age, R(drip_tran_p), ax3, ax4)

    # Shared x-axis
    for ax in all_axes:
        ax.set_xlim(plot_age[0], plot_age[-1])
        ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=12))

    plt.setp(ax1.get_xticklabels(), visible=False)
    plt.setp(ax2.get_xticklabels(), visible=False)
    plt.setp(ax3.get_xticklabels(), visible=False)

    # Colorbar
    cbar_ax = fig.add_axes([0.905, 0.08, 0.013, 0.88])
    cbar = fig.colorbar(im4, cax=cbar_ax, ticks=[0, 0.5, 1])
    cbar.set_label('Normalised PDF', fontsize=6, labelpad=4)
    cbar.ax.tick_params(labelsize=6)

    for path in (OUTPUT_PDF, OUTPUT_PNG, OUTPUT_EPS):
        fig.savefig(path, dpi=300, bbox_inches='tight')
        print(f"  Saved → {path}")
    plt.close(fig)

    # ── Supplementary recurrence plot figure ──────────────────────────
    if params is not None:
        make_rp_figure(drip_real, d18o_real, params)

    print("\nDone.\n")


# ═══════════════════════════════════════════════════════════════════════
# SUPPLEMENTARY RECURRENCE PLOT FIGURE
# ═══════════════════════════════════════════════════════════════════════

def make_rp_figure(drip_real, d18o_real, params, age_marker=5200):
    print("\nBuilding supplementary recurrence plot figure ...")

    def median_series(real_df):
        real_cols = [c for c in real_df.columns if c.startswith('r')]
        return real_df['age'].values, np.nanmedian(real_df[real_cols].values, axis=1)

    drip_ages, drip_med = median_series(drip_real)
    d18o_ages, d18o_med = median_series(d18o_real)

    drip_row = params[params['proxy'] == 'drip_rate'].iloc[0]
    d18o_row = params[params['proxy'] == 'd18O'].iloc[0]

    drip_RM = build_recurrence_matrix(
        drip_med, int(drip_row['tau']), int(drip_row['m']),
        target_rr=float(drip_row['RR']), theiler=int(drip_row['theiler']))
    d18o_RM = build_recurrence_matrix(
        d18o_med, int(d18o_row['tau']), int(d18o_row['m']),
        target_rr=float(d18o_row['RR']), theiler=int(d18o_row['theiler']))

    fig, axes = plt.subplots(1, 2, figsize=(8.5, 4.0),
                             gridspec_kw={'wspace': 0.35,
                                          'left': 0.10, 'right': 0.96,
                                          'top': 0.88, 'bottom': 0.13})

    for ax, RM, ages, title in zip(
            axes, [d18o_RM, drip_RM], [d18o_ages, drip_ages],
            [r'$\delta^{18}$O', 'Drip rate']):

        if RM is None:
            ax.text(0.5, 0.5, 'Insufficient data',
                    ha='center', va='center', transform=ax.transAxes)
            continue

        age_min, age_max = ages.min(), ages.max()
        ax.imshow(RM, cmap='binary', origin='lower', aspect='auto',
                  extent=[age_min, age_max, age_min, age_max],
                  interpolation='none')
        ax.axvline(age_marker, color='red', lw=0.8, ls='--', alpha=0.8)
        ax.axhline(age_marker, color='red', lw=0.8, ls='--', alpha=0.8)
        ax.text(age_marker + 80, age_max * 0.97,
                f'{age_marker/1000:.1f} ka',
                color='red', fontsize=6, va='top')
        ax.set_xlabel('Age (yr BP)', fontsize=8, fontweight='bold')
        ax.set_ylabel('Age (yr BP)', fontsize=8, fontweight='bold')
        ax.text(0.02, 0.97, title.split()[0] if title and title[0].islower() and len(title.split()[0])<=2 else '', transform=ax.transAxes, fontsize=9, fontweight='bold', va='top')
        ax.tick_params(width=0.5, length=2.5, direction='in')

        det, trans = rqa_from_matrix(RM)
        ax.text(0.02, 0.98, f'DET={det:.3f}   TRANS={trans:.3f}',
                transform=ax.transAxes, fontsize=6, va='top', ha='left',
                bbox=dict(boxstyle='round,pad=0.3', fc='white',
                          ec='grey', alpha=0.8))


    for path in (SUPP_PDF, SUPP_PNG):
        fig.savefig(path, dpi=300, bbox_inches='tight')
        print(f"  Saved → {path}")
    plt.close(fig)


if __name__ == '__main__':
    main()
