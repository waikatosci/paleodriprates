#!/usr/bin/env python3
"""
dripwater_source_stability.py — distributional analysis of the HS4 dripwater
source composition underpinning the static-median treatment of X_a and Y_a.

NCOMMS-26-041445-T companion analysis (Reviewer 1, source-stationarity query).

The kinetic inversion fixes the aqueous source composition at the observed
medians of the 2007-2016 HS4 dripwater monitoring series:
Ca 67,437 ppb; Ni 4.370 ppb; Co 0.460 ppb (n = 92) — the exact values carried
by the Dr Paleo canonical parameter set. This script quantifies how variable
that source composition actually was over the eight-year window, and documents
the distributional facts the defence rests on:

  1. Ca is pinned at a saturation ceiling: GSD 1.065 (+/-6 %), total range
     1.44x, LEFT-skewed — fixing Ca at its median is near-exact by observation.
  2. The trace metals occupy bounded multiplicative bands (Co 3.4x, Ni 5.8x);
     lognormality is REJECTED for all series (Shapiro-Wilk on logs, p < 0.05) —
     the distributions are bounded/compressed, not heavy-tailed.
  3. Me/Ca ratio stability derives from a near-constant denominator, NOT from
     covariance cancellation: r(ln Ni, ln Ca) = -0.11, r(ln Co, ln Ca) = +0.06,
     Var[ln(Ni/Ca)] ~ 1.1 x Var[ln Ni]. The metals share a common carrier
     instead: r(ln Ni, ln Co) = +0.73.

Inputs   ../dripwater/HS4_dripwater_canonical.csv   (92 samples, 2007-2016)
Outputs  ../manuscript_figures/output/TableS_dripwater_distributions.csv
         ../manuscript_figures/output/FigS_dripwater_distributions.{png,pdf}

The companion propagation of this dispersion through the kinetic inversion is
source_variation_propagation.py.
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "manuscript_figures"))
from ngeo_style import panel_label, apply_style, COL_TEAL_800, COL_TEAL_600, COL_BROWN_500, COL_BG_900  # noqa: E402

import matplotlib.pyplot as plt  # noqa: E402

DW_CSV = os.path.join(ROOT, "dripwater", "HS4_dripwater_canonical.csv")
OUTDIR = os.path.join(ROOT, "manuscript_figures", "output")

SERIES = [
    ("Ca (ppb)",            "Ca_ppb"),
    ("Co (ppb)",            "Co_ppb"),
    ("Ni (ppb)",            "Ni_ppb"),
    ("Ni/Ca (umol/mol)",    "NiCa_umol_per_mol"),
    ("Co/Ca (umol/mol)",    "CoCa_umol_per_mol"),
]


def dist_row(label, x):
    x = pd.Series(x).dropna().astype(float)
    ln = np.log(x)
    sw = stats.shapiro(ln)
    return dict(
        variable=label,
        n=len(x),
        median=x.median(),
        q25=x.quantile(0.25),
        q75=x.quantile(0.75),
        GSD=float(np.exp(ln.std(ddof=1))),
        CV_pct=100.0 * x.std(ddof=1) / x.mean(),
        fold_range=float(x.max() / x.min()),
        skew_raw=float(stats.skew(x)),
        shapiro_p_log=float(sw.pvalue),
    )


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    dw = pd.read_csv(DW_CSV, parse_dates=["date"])

    # ── Table ────────────────────────────────────────────────────────────
    tab = pd.DataFrame([dist_row(lab, dw[col]) for lab, col in SERIES])
    tab_path = os.path.join(OUTDIR, "TableS_dripwater_distributions.csv")
    tab.round(4).to_csv(tab_path, index=False)

    lnNi, lnCo, lnCa = np.log(dw.Ni_ppb), np.log(dw.Co_ppb), np.log(dw.Ca_ppb)
    r_ni_ca, r_co_ca = lnNi.corr(lnCa), lnCo.corr(lnCa)
    r_ni_co = lnNi.corr(lnCo)
    var_infl = np.log(dw.NiCa_umol_per_mol).var(ddof=1) / lnNi.var(ddof=1)

    print("HS4 dripwater source-composition summary "
          f"(n = {len(dw)}, {dw.date.min().date()} to {dw.date.max().date()})")
    print(tab.round(3).to_string(index=False))
    print(f"\nr(ln Ni, ln Ca) = {r_ni_ca:+.2f}   r(ln Co, ln Ca) = {r_co_ca:+.2f}   "
          f"r(ln Ni, ln Co) = {r_ni_co:+.2f}")
    print(f"Var[ln(Ni/Ca)] / Var[ln Ni] = {var_infl:.2f} "
          "(Ca-normalisation cancels essentially nothing; the denominator is simply constant)")
    print("Shapiro-Wilk on logs rejects lognormality for all five series (p < 0.05): "
          "bounded/compressed distributions, not heavy tails.")
    med = dict(Ca=dw.Ca_ppb.median(), Ni=dw.Ni_ppb.median(), Co=dw.Co_ppb.median())
    print(f"Medians carried into the Dr Paleo canonical parameter set: "
          f"Ca {med['Ca']:.1f} ppb, Ni {med['Ni']:.6f} ppb, Co {med['Co']:.6f} ppb")

    # ── Figure ───────────────────────────────────────────────────────────
    apply_style()
    fig, axs = plt.subplots(2, 2, figsize=(7.09, 4.6))

    # (a) aligned marginals: ln(x/median), one strip per series
    ax = axs[0, 0]
    rng = np.random.default_rng(42)
    for i, (lab, col) in enumerate(SERIES):
        x = dw[col].dropna().astype(float)
        z = np.log(x / x.median())
        jitter = rng.uniform(-0.16, 0.16, len(z))
        c = COL_BG_900 if col == "Ca_ppb" else (COL_TEAL_800 if "Ni" in col else COL_BROWN_500)
        ax.scatter(i + jitter, z, s=2.5, alpha=0.45, lw=0, color=c, rasterized=True)
        q = np.log(np.quantile(x, [0.25, 0.75]) / x.median())
        ax.plot([i - 0.24, i + 0.24], [q[0]] * 2, color=c, lw=0.6)
        ax.plot([i - 0.24, i + 0.24], [q[1]] * 2, color=c, lw=0.6)
        ax.plot([i - 0.3, i + 0.3], [0, 0], color=c, lw=1.0)
    ax.axhline(0, color="0.85", lw=0.4, zorder=0)
    ax.set_xticks(range(len(SERIES)))
    ax.set_xticklabels(["Ca", "Co", "Ni", "Ni/Ca", "Co/Ca"])
    ax.set_xlim(-0.9, len(SERIES) - 0.5)
    ax.set_ylabel("ln(value / median)")
    panel_label(ax, "a")

    # (b) skew vs dispersion — the anti-"skew constrains range" panel
    ax = axs[0, 1]
    for lab, col in SERIES:
        x = dw[col].dropna().astype(float)
        gsd = np.exp(np.log(x).std(ddof=1))
        sk = stats.skew(x)
        c = COL_BG_900 if col == "Ca_ppb" else (COL_TEAL_800 if "Ni" in col else COL_BROWN_500)
        ax.scatter(sk, gsd, s=14, color=c, zorder=3)
        # labels to the left of the marker so nothing runs past the axes frame
        ax.annotate(lab.split(" ")[0], (sk, gsd), textcoords="offset points",
                    xytext=(-4, -1), fontsize=5.5, color=c, ha="right", va="center")
    ax.axvline(0, color="0.85", lw=0.4, zorder=0)
    ax.margins(x=0.18, y=0.15)
    ax.set_xlabel("Skewness (raw values)")
    ax.set_ylabel("Geometric SD")
    panel_label(ax, "b")

    # (c) Me-Ca decoupling
    ax = axs[1, 0]
    ax.scatter(lnCa, lnNi, s=4, color=COL_TEAL_800, alpha=0.6, lw=0,
               label=f"Ni  (r = {r_ni_ca:+.2f})", rasterized=True)
    ax.scatter(lnCa, lnCo + np.log(dw.Ni_ppb.median() / dw.Co_ppb.median()),
               s=4, color=COL_BROWN_500, alpha=0.6, lw=0,
               label=f"Co (offset; r = {r_co_ca:+.2f})", rasterized=True)
    ax.set_xlabel("ln Ca (ppb)")
    ax.set_ylabel("ln Me (ppb, Co offset to Ni)")
    ax.legend(frameon=False, loc="lower left")
    panel_label(ax, "c")

    # (d) Ni-Co coupling — the shared carrier
    ax = axs[1, 1]
    ax.scatter(lnCo, lnNi, s=5, color=COL_TEAL_600, alpha=0.7, lw=0, rasterized=True)
    b, a = np.polyfit(lnCo, lnNi, 1)
    xs = np.linspace(lnCo.min(), lnCo.max(), 10)
    ax.plot(xs, a + b * xs, color=COL_BG_900, lw=0.7, ls="--")
    ax.text(0.12, 0.92, f"r = {r_ni_co:+.2f}", transform=ax.transAxes, fontsize=6, color=COL_BG_900)
    ax.set_xlabel("ln Co (ppb)")
    ax.set_ylabel("ln Ni (ppb)")
    panel_label(ax, "d")

    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUTDIR, f"FigS_dripwater_distributions.{ext}"))
    print(f"\nwrote {tab_path}")
    print(f"wrote {os.path.join(OUTDIR, 'FigS_dripwater_distributions.png')} (+pdf)")


if __name__ == "__main__":
    main()
