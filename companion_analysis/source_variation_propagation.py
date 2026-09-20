#!/usr/bin/env python3
"""
source_variation_propagation.py — how much drip-rate uncertainty does the
static-median source composition actually inject?

NCOMMS-26-041445-T companion analysis (Reviewer 1, source-stationarity query;
Reviewer 3, uncertainty budget). Companion to dripwater_source_stability.py.

The kinetic inversion holds the aqueous source composition (X_a, Y_a) fixed at
the 2007-2016 dripwater medians. If the true source ratio at deposition time
differed from the median by d = Delta ln(Me/Ca)_source, the inversion
misattributes that offset to drip rate. Because the forward model is

    X_s = K0 * Phi(V),   K0 prop. (Me/Ca)_source,   Phi(V) = 1 - nS*E1(V),

the induced drip-rate error at fixed measured X_s is

    Delta ln V = - d / S(V),      S(V) = d ln Phi / d ln V,

so 1/|S(V)| is the DAMPING FACTOR of the kinetic transfer: where the proxy is
responsive (kinetic window), source noise is damped; where the curve flattens,
the proxy saturates and no reconstruction is attempted anyway.

Two source-dispersion estimates are propagated:
  (a) CONSERVATIVE — the entire observed dripwater ratio dispersion treated as
      source variation (upper bound: part of that spread is real signal);
  (b) PARTITIONED — the residual dispersion after removing the component of
      dripwater Me/Ca that covaries with observed drip rate (86 paired samples).

Propagation is a joint empirical bootstrap (N = 10,000): paired (Ni, Co)
deviations resampled from the observed 92-sample record — preserving the
bounded, non-lognormal shapes and the Ni-Co carrier covariance (r = +0.73) —
mapped through the inversion at reference drip rates, and combined as the
precision-weighted two-proxy estimate (weights prop. S_i^2, matching the joint
posterior product used in the reconstruction; the equal-weight variant is also
reported and is nearly identical).

Inputs   ../dripwater/HS4_dripwater_canonical.csv
Outputs  ../manuscript_figures/output/TableS_source_propagation.csv
         ../manuscript_figures/output/FigS_source_propagation.{png,pdf}
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "manuscript_figures"))
from ngeo_style import (panel_label, label_bubble, apply_style, COL_NI, COL_CO, COL_BG_900, COL_BG_300, COL_BG_600,  # noqa: E402
                        COL_TEAL_200, COL_DORANGE, DOUBLE_COL)

import matplotlib.pyplot as plt  # noqa: E402

DW_CSV = os.path.join(ROOT, "dripwater", "HS4_dripwater_canonical.csv")
OUTDIR = os.path.join(ROOT, "manuscript_figures", "output")

# ── Canonical kinetic parameters (Dr Paleo production set, 2026-07-20) ────
PARAMS = {
    "Ni": dict(kd_mn=-3.8372, kd_sd=1.282549830161864, F=0.01),
    "Co": dict(kd_mn=-5.4171, kd_sd=1.282549830161864, F=0.01),
}
RATIO_COL = {"Ni": "NiCa_umol_per_mol", "Co": "CoCa_umol_per_mol"}

# Record anchors (canonical native run 20260720_202730)
V_MODERN = 15.62        # top-of-record pc50, drips/min (16.66-target refit)
V_FLOOR = 2.81          # lowest resolved point (157.01 cm; 568-point record, 2026-09-21)
LN_SIGNAL = np.log(V_MODERN / V_FLOOR)   # ln-range of the native record

V_REFS = [V_FLOOR, 2.0, 5.0, 10.0, V_MODERN, 20.0]
N_BOOT = 10_000
K_HALFWIDTH_SD = 12.0
K_NPTS = 4001


def phi_of_V(V_drips_min, kd_mn, kd_sd, F):
    """Forward kinetic transfer Phi(V) = 1 - nS*E1(V); V in drips/min."""
    V = np.atleast_1d(np.asarray(V_drips_min, dtype=float))
    nS = 1.0 - F
    k = np.linspace(kd_mn - K_HALFWIDTH_SD * kd_sd,
                    kd_mn + K_HALFWIDTH_SD * kd_sd, K_NPTS)
    gk = np.exp(-(k - kd_mn) ** 2 / (2 * kd_sd ** 2)) / (kd_sd * np.sqrt(2 * np.pi))
    w = np.gradient(k)
    tau = 60.0 / V                                   # seconds between drips
    E1 = np.einsum("k,vk->v", gk * w, np.exp(-np.exp(k)[None, :] * tau[:, None]))
    return 1.0 - nS * E1


def slope_S(V, kd_mn, kd_sd, F, eps=1e-4):
    """S(V) = d ln Phi / d ln V (central difference in ln V)."""
    lo = phi_of_V(V * np.exp(-eps), kd_mn, kd_sd, F)
    hi = phi_of_V(V * np.exp(+eps), kd_mn, kd_sd, F)
    return (np.log(hi) - np.log(lo)) / (2 * eps)


def main():
    os.makedirs(OUTDIR, exist_ok=True)
    dw = pd.read_csv(DW_CSV, parse_dates=["date"])
    rng = np.random.default_rng(42)

    # ── Source deviations, both estimates ────────────────────────────────
    dev_a = {}   # conservative: ln(x) - ln(median), full n = 92
    dev_b = {}   # partitioned: residuals of ln(x) ~ ln(DR), paired n = 86
    paired = dw.dropna(subset=["DR_drops_per_min"])
    lnDR = np.log(paired.DR_drops_per_min.values)
    for m in ("Ni", "Co"):
        x = dw[RATIO_COL[m]].astype(float)
        dev_a[m] = np.log(x.values) - np.log(x.median())
        y = np.log(paired[RATIO_COL[m]].astype(float).values)
        b1, b0 = np.polyfit(lnDR, y, 1)
        dev_b[m] = y - (b0 + b1 * lnDR)
        print(f"{m}: sigma_ln(source) conservative = {dev_a[m].std(ddof=1):.3f} "
              f"(GSD {np.exp(dev_a[m].std(ddof=1)):.3f}); "
              f"partitioned residual = {dev_b[m].std(ddof=1):.3f} "
              f"(slope vs ln DR = {b1:+.2f})")
    r_pair = np.corrcoef(dev_a["Ni"], dev_a["Co"])[0, 1]
    print(f"carrier covariance preserved in bootstrap: r(dNi, dCo) = {r_pair:+.2f}\n")

    # ── Transfer slopes ──────────────────────────────────────────────────
    Vgrid = np.geomspace(0.3, 60, 400)
    S = {m: slope_S(Vgrid, **PARAMS[m]) for m in ("Ni", "Co")}
    S_at = {m: {v: float(slope_S(np.array([v]), **PARAMS[m])[0]) for v in V_REFS}
            for m in ("Ni", "Co")}

    # ── Bootstrap propagation ────────────────────────────────────────────
    rows = []
    for tag, dev in (("conservative", dev_a), ("partitioned", dev_b)):
        n = len(dev["Ni"])
        idx = rng.integers(0, n, N_BOOT)
        dN, dC = dev["Ni"][idx], dev["Co"][idx]           # paired draws
        for v in V_REFS:
            sN, sC = S_at["Ni"][v], S_at["Co"][v]
            d_joint_pw = -(sN * dN + sC * dC) / (sN ** 2 + sC ** 2)
            d_joint_eq = -0.5 * (dN / sN + dC / sC)
            for scheme, d in (("precision", d_joint_pw), ("equal", d_joint_eq)):
                q = np.quantile(d, [0.025, 0.16, 0.84, 0.975])
                rows.append(dict(
                    estimate=tag, weighting=scheme, V_ref=v,
                    S_Ni=sN, S_Co=sC,
                    sigma_lnV=d.std(ddof=1),
                    q16_84_pct=100 * (np.exp(q[2]) - np.exp(q[1])) / 2,
                    env68_lo=v * np.exp(q[1]), env68_hi=v * np.exp(q[2]),
                    env95_lo=v * np.exp(q[0]), env95_hi=v * np.exp(q[3]),
                    frac_of_signal=d.std(ddof=1) / LN_SIGNAL,
                ))
        # single-proxy disagreement bound (source noise alone)
        for v in (V_FLOOR, V_MODERN):
            sN, sC = S_at["Ni"][v], S_at["Co"][v]
            dis = dN / sN - dC / sC
            print(f"[{tag}] predicted single-proxy disagreement from source noise "
                  f"at V = {v:.2f}: RMS = {dis.std(ddof=1):.3f} ln-units "
                  f"({100 * (np.exp(dis.std(ddof=1)) - 1):.0f} %)")

    tab = pd.DataFrame(rows)
    tab_path = os.path.join(OUTDIR, "TableS_source_propagation.csv")
    tab.round(4).to_csv(tab_path, index=False)

    hd = tab[(tab.estimate == "conservative") & (tab.weighting == "precision")]
    hd = hd.set_index("V_ref")
    print("\nHEADLINE (conservative, precision-weighted joint):")
    for v in V_REFS:
        r = hd.loc[v]
        print(f"  V = {v:6.2f} drips/min: sigma_lnV = {r.sigma_lnV:.3f} "
              f"(+/-{100 * (np.exp(r.sigma_lnV) - 1):5.1f} %), 68 % envelope "
              f"{r.env68_lo:5.2f}-{r.env68_hi:5.2f}, = {100 * r.frac_of_signal:.1f} % "
              f"of the record ln-range")
    worst = hd.loc[V_REFS].sigma_lnV.max()
    print(f"\n  Worst-case across anchors: sigma_lnV = {worst:.3f} "
          f"= {100 * worst / LN_SIGNAL:.0f} % of the native record ln-range "
          f"(ln {V_MODERN}/{V_FLOOR} = {LN_SIGNAL:.2f}); the 5.2 ka megadrought "
          f"collapse is ~{LN_SIGNAL / worst:.0f}x the worst-case source-induced spread.")

    # ── Figure ───────────────────────────────────────────────────────────
    apply_style()
    fig, axs = plt.subplots(1, 2, figsize=(DOUBLE_COL, 2.5))

    # (a) transfer curves + mapped source band at V = 14.14
    ax = axs[0]
    for m, c in (("Ni", COL_NI), ("Co", COL_CO)):
        ax.plot(Vgrid, phi_of_V(Vgrid, **PARAMS[m]), color=c, lw=1.0, label=m)
    m, c = "Ni", COL_NI
    sig = dev_a[m].std(ddof=1)
    phi_ref = phi_of_V(np.array([V_MODERN]), **PARAMS[m])[0]
    band = (phi_ref * np.exp(-sig), phi_ref * np.exp(+sig))
    # invert band edges numerically on the grid
    phiN = phi_of_V(Vgrid, **PARAMS["Ni"])
    v_lo = np.interp(-band[1], -phiN, Vgrid)   # phi decreasing in V
    v_hi = np.interp(-band[0], -phiN, Vgrid)
    ax.axhspan(band[0], band[1], color=COL_TEAL_200, alpha=0.6, lw=0, zorder=0)
    ax.axvspan(v_lo, v_hi, color=COL_TEAL_200, alpha=0.9, lw=0, zorder=0)
    ax.plot([V_MODERN], [phi_ref], "o", ms=3, color=COL_BG_900, zorder=5)
    label_bubble(ax, f"source band (GSD {np.exp(sig):.2f})",
                 xy=(1.3, 0.1544), xytext=(0.5935, 0.3191), va="center",
                 color=COL_NI, bubble_ec=COL_NI, arrow_color=COL_NI)
    label_bubble(ax, f"mapped drip-rate band\n{v_lo:.1f}\u2013{v_hi:.1f} drips min$^{{-1}}$",
                 xy=(12.35, 0.5162), xytext=(3.407, 0.7823), va="center",
                 color=COL_NI, bubble_ec=COL_NI, arrow_color=COL_NI)
    ax.set_ylim(0, 1.0)
    ax.set_xscale("log")
    ax.set_xlabel("Drip rate V (drips min$^{-1}$)")
    ax.set_ylabel(r"Kinetic transfer  $\Phi(V)$")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2, fontsize=5.5)
    panel_label(ax, "a")

    # (b) induced sigma_lnV(V), both estimates, joint
    ax = axs[1]
    sig_a = {m: dev_a[m].std(ddof=1) for m in ("Ni", "Co")}
    sig_b = {m: dev_b[m].std(ddof=1) for m in ("Ni", "Co")}
    cov_a = np.cov(dev_a["Ni"], dev_a["Co"])[0, 1]
    cov_b = np.cov(dev_b["Ni"], dev_b["Co"])[0, 1]
    for (sig, cov, ls, lab) in ((sig_a, cov_a, "-", "conservative"),
                                (sig_b, cov_b, "--", "partitioned")):
        num = (S["Ni"] * sig["Ni"]) ** 2 + (S["Co"] * sig["Co"]) ** 2 \
            + 2 * S["Ni"] * S["Co"] * cov
        s_joint = np.sqrt(num) / (S["Ni"] ** 2 + S["Co"] ** 2)
        ax.plot(Vgrid, 100 * (np.exp(s_joint) - 1), color=COL_BG_900, ls=ls,
                lw=1.0, label=f"joint ({lab})")
    for m, c in (("Ni", COL_NI), ("Co", COL_CO)):
        ax.plot(Vgrid, 100 * (np.exp(sig_a[m] / np.abs(S[m])) - 1), color=c,
                lw=0.7, alpha=0.8, label=f"{m} only (conservative)")
    for v in (V_FLOOR, V_MODERN):
        ax.axvline(v, color=COL_BG_300, lw=0.5, ls=":")
    label_bubble(ax, "drought floor", xy=(1.236, 255.1), xytext=(2.669, 255), va="center", color=COL_BG_600, bubble_ec=COL_BG_600, arrow_color=COL_BG_600)
    label_bubble(ax, "modern",        xy=(16.07, 224.3), xytext=(28.27, 387.4), va="center", color=COL_BG_600, bubble_ec=COL_BG_600, arrow_color=COL_BG_600)
    ax.axhline(100 * (np.exp(LN_SIGNAL) - 1), color=COL_DORANGE, lw=0.6)
    label_bubble(ax, "record signal range",
                 xy=(3, 100 * (np.exp(LN_SIGNAL) - 1)),
                 xytext=(6.907, 568.7), va="center",
                 color=COL_DORANGE, bubble_ec=COL_DORANGE, arrow_color=COL_DORANGE)
    ax.set_ylim(10, 3000)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Drip rate V (drips min$^{-1}$)")
    ax.set_ylabel("Source-induced DR uncertainty (%)")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.20), ncol=2, fontsize=5.5)
    panel_label(ax, "b")

    fig.tight_layout()
    fig.subplots_adjust(bottom=0.24)
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUTDIR, f"FigS_source_propagation.{ext}"))
    print(f"\nwrote {tab_path}")
    print(f"wrote {os.path.join(OUTDIR, 'FigS_source_propagation.png')} (+pdf)")


if __name__ == "__main__":
    main()
