#!/usr/bin/env python3
"""
make_supp_figs_1_5.py -- regenerate Supplementary Figures 1-5 of Hartland et al.
(NCOMMS-26-041445-T) in the manuscript house style (ngeo_style: bare bold
panel letters, no titles), from the released data and the locked calibration.

  SF1  Ni-Co covariance by epoch, exponential fit (outliers of SM3 excluded)
  SF2  Sensitivity of the kinetic transfer Phi(V) to the fast fraction lambda_F,
       with the induced shift of the Holocene drip-rate record (inset)
  SF3  Recurrence plots of the median d18O and drip-rate series (SM6)
  SF4  Modelled kinetic isotope fractionation of the labile pool (SM8, eq. 18)
  SF5  Native-resolution drip-rate record with the Ni and Co series

Inputs  ../manuscript_figures/external/HS4_TE_multielement.csv
        ../manuscript_figures/external/drip_rate_summary_hr_censored.csv
        ../manuscript_figures/external/drip_rate_summary_ap.csv
        ../extended_data/Ex_Data_2_RQA/Drip_rate.xlsx (sheet 6.OutIsotope)
        ../extended_data/Ex_Data_2_RQA/rqa_parameters.csv
Outputs ./output/SuppFig{1..5}.{png,pdf}
"""
import os
import sys

import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(ROOT, "manuscript_figures"))
sys.path.insert(0, os.path.join(ROOT, "companion_analysis"))
sys.path.insert(0, os.path.join(ROOT, "extended_data", "Ex_Data_2_RQA"))
from ngeo_style import (apply_style, style_ax, panel_label, SINGLE_COL, DOUBLE_COL,  # noqa: E402
                        COL_NI, COL_CO, COL_BG_900, COL_BG_400, COL_BG_300, COL_BG_200,
                        COL_TEAL_600, COL_TEAL_200, COL_BROWN_500, COL_DORANGE,
                        COL_GREEN_800, COL_RED_900)
from source_variation_propagation import phi_of_V, PARAMS                       # noqa: E402
import matplotlib.pyplot as plt                                                 # noqa: E402

EXT = os.path.join(ROOT, "manuscript_figures", "external")
OUT = os.path.join(HERE, "output")
os.makedirs(OUT, exist_ok=True)
OUTLIER_DEPTHS = (7.1, 111.1)          # the two SM3 outliers (cm; C-992, B-422)
KD_NI, KD_CO = PARAMS["Ni"]["kd_mn"], PARAMS["Co"]["kd_mn"]
SIG = PARAMS["Ni"]["kd_sd"]


def save(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(OUT, f"{name}.{ext}"))
    plt.close(fig)
    print("wrote", name)


def load_te():
    te = pd.read_csv(os.path.join(EXT, "HS4_TE_multielement.csv"))
    te = te.dropna(subset=["age_yBP"]).sort_values("depth_cm").reset_index(drop=True)
    return te


# ── SF1: Ni-Co by epoch ─────────────────────────────────────────────────
def sf1(te):
    d = te[~np.isclose(te.depth_cm.values[:, None], OUTLIER_DEPTHS, atol=0.1).any(1)].dropna(subset=["Ni", "Co"])
    epochs = [("Early Holocene (9.5-6 ka)", 6000, 20000, COL_BG_400),
              ("5.2 ka event (6-5 ka)", 5000, 6000, COL_DORANGE),
              ("Late Holocene (5-0.2 ka)", 200, 5000, COL_TEAL_600),
              ("Last 200 years", -100, 200, COL_BG_900)]
    f = lambda x, a, b: a * np.exp(b * x)
    (a, b), _ = curve_fit(f, d.Ni, d.Co, p0=(0.1, 0.3))
    apply_style()
    fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL * 0.9))
    for lab, lo, hi, c in epochs:
        s = d[(d.age_yBP > lo) & (d.age_yBP <= hi)]
        ax.scatter(s.Ni, s.Co, s=7, color=c, alpha=0.75, lw=0, label=lab)
    xs = np.linspace(0.5, 11.5, 200)
    ax.plot(xs, f(xs, a, b), color=COL_BG_900, lw=0.8, ls="--",
            label=f"Co = {a:.2f} exp({b:.2f} Ni)")
    ax.set_xlabel("Ni (ppm)")
    ax.set_ylabel("Co (ppm)")
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 4)
    ax.legend(frameon=False, loc="upper left")
    style_ax(ax)
    fig.tight_layout()
    save(fig, "SuppFig1")
    print(f"  SF1 fit: Co = {a:.3f} exp({b:.3f} Ni), n = {len(d)}")


# ── SF2: lambda_F sensitivity ────────────────────────────────────────────
def sf2():
    ap = pd.read_csv(os.path.join(EXT, "drip_rate_summary_ap.csv"))
    ap = ap[ap.age >= 0]
    V = np.geomspace(0.1, 60, 400)
    Fs = [(1e-6, "$\\lambda_F$ = 10$^{-6}$", COL_TEAL_600),
          (0.01, "$\\lambda_F$ = 0.01 (adopted)", COL_BG_900),
          (0.05, "$\\lambda_F$ = 0.05", COL_DORANGE)]
    apply_style()
    fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL * 0.8))
    for F, lab, c in Fs:
        ax.plot(V, phi_of_V(V, KD_NI, SIG, F), color=c, lw=0.9, label=lab)
    ax.set_xscale("log")
    ax.set_xlabel("Drip rate (drips min$^{-1}$)")
    ax.set_ylabel("Kinetic transfer $\\Phi$ (Ni)")
    ax.set_ylim(0, 1)
    ax.legend(frameon=False, loc="lower left")
    style_ax(ax)
    # inset: record with +/-50 % lambda_F band. At fixed measured Phi, the
    # drip rate that lambda_F' implies is Phi^{-1}(Phi(V0; 0.01); lambda_F').
    ins = ax.inset_axes([0.53, 0.55, 0.44, 0.40])
    grid = np.geomspace(0.05, 200, 2000)
    phi0 = phi_of_V(ap.pc50.values, KD_NI, SIG, 0.01)
    band = []
    for F in (0.005, 0.015):
        phiF = phi_of_V(grid, KD_NI, SIG, F)
        band.append(np.interp(-phi0, -phiF, grid))      # Phi decreasing in V
    lo, hi = np.minimum(*band), np.maximum(*band)
    ka = ap.age.values / 1000
    ins.fill_between(ka, lo, hi, color=COL_BG_200, lw=0, label="$\\lambda_F$ 0.005-0.015")
    ins.plot(ka, ap.pc50, color=COL_BG_900, lw=0.5, label="median (0.01)")
    ins.set_xlim(9.7, 0)
    ins.set_ylim(0, 35)
    ins.set_xlabel("Age (ka BP)", fontsize=5)
    ins.set_ylabel("Drip rate", fontsize=5)
    ins.tick_params(labelsize=4.5)
    ins.legend(frameon=False, fontsize=4.2, loc="upper right")
    style_ax(ins)
    fig.tight_layout()
    save(fig, "SuppFig2")
    m = (ap.age >= 5100) & (ap.age <= 5200)
    rel = 100 * (hi - lo) / (2 * ap.pc50.values)
    print(f"  SF2 +/-50 % lambda_F band: median half-width {np.median(rel):.0f} %, "
          f"5.1-5.2 ka {np.median(rel[m]):.0f} %")


# ── SF3: recurrence plots ────────────────────────────────────────────────
def sf3():
    from RQA_HS4_ensemble import (build_recurrence_matrix, rqa_from_matrix,
                                  load_d18o_summary)
    ap = pd.read_csv(os.path.join(EXT, "drip_rate_summary_ap.csv"))
    xl = os.path.join(ROOT, "extended_data", "Ex_Data_2_RQA", "Drip_rate.xlsx")
    d18 = load_d18o_summary(xl, "6.OutIsotope")
    par = pd.read_csv(os.path.join(ROOT, "extended_data", "Ex_Data_2_RQA",
                                   "rqa_parameters.csv")).set_index("proxy")
    lo = max(ap.age.min(), d18.age.min()); hi = min(ap.age.max(), d18.age.max())
    ap = ap[(ap.age >= lo) & (ap.age <= hi)]; d18 = d18[(d18.age >= lo) & (d18.age <= hi)]
    series = [("$\\delta^{18}$O", d18.age.values, d18["median"].values, par.loc["d18O"]),
              ("Drip rate", ap.age.values, ap.pc50.values, par.loc["drip_rate"])]
    apply_style()
    fig, axs = plt.subplots(1, 2, figsize=(DOUBLE_COL, DOUBLE_COL * 0.5))
    for ax, (lab, ages, x, p), L in zip(axs, series, "ab"):
        RM = build_recurrence_matrix(x, int(p.tau), int(p.m), theiler=int(p.theiler))
        det, trans = rqa_from_matrix(RM)
        # display at <= 1500 px per side (block maximum), the matrix itself is unchanged
        n = RM.shape[0]; k = int(np.ceil(n / 1500)); npad = k * int(np.ceil(n / k))
        D = np.zeros((npad, npad), bool); D[:n, :n] = RM.astype(bool)
        D = D.reshape(npad // k, k, npad // k, k).any(axis=(1, 3))
        del RM
        ax.imshow(D, cmap="binary", origin="lower", aspect="equal",
                  extent=[ages.min(), ages.max(), ages.min(), ages.max()],
                  interpolation="none")
        for fn in (ax.axvline, ax.axhline):
            fn(5200, color=COL_DORANGE, lw=0.6, ls="--")
        ax.text(0.98, 0.03, f"{lab}\nDET = {det:.3f}   TRANS = {trans:.3f}",
                transform=ax.transAxes, fontsize=5.5, ha="right", va="bottom",
                bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=COL_BG_300, lw=0.3))
        ax.set_xlabel("Age (yr BP)"); ax.set_ylabel("Age (yr BP)")
        panel_label(ax, L)
        style_ax(ax)
        print(f"  SF3 {lab}: tau {int(p.tau)} m {int(p.m)} DET {det:.3f} TRANS {trans:.3f}")
    fig.tight_layout()
    save(fig, "SuppFig3")


# ── SF4: kinetic isotope fractionation (eq. 18) ─────────────────────────
def sf4():
    alpha, inert = 0.9995, 0.2
    metals = [("$\\delta^{60}$Ni", 0.0, KD_NI, COL_NI), ("$\\delta^{65}$Cu", -0.3, KD_NI, COL_BROWN_500)]

    def delta(V, d0, kd, fi):
        f = (1 - fi) * phi_of_V(V, kd, SIG, 0.01)
        return d0 + 1000 * (1 - alpha) * np.log(1 - f)

    apply_style()
    fig, ax = plt.subplots(figsize=(SINGLE_COL, SINGLE_COL * 0.8))
    V = np.linspace(0.02, 30, 600)
    for lab, d0, kd, c in metals:
        ax.plot(V, delta(V, d0, kd, 0.0), color=c, lw=0.9, label=f"{lab}, fully labile")
        ax.plot(V, delta(V, d0, kd, inert), color=c, lw=0.9, ls="--",
                label=f"{lab}, inert fraction {inert}")
    ax.set_xlim(0, 30); ax.set_ylim(-1.4, 0.1)
    ax.set_xlabel("Drip rate (drips min$^{-1}$)")
    ax.set_ylabel("$\\delta_{labile}$ (‰)")
    top = ax.secondary_xaxis("top", functions=(lambda v: v, lambda v: v))
    ticks = [2, 5, 10, 20, 30]
    top.set_xticks(ticks); top.set_xticklabels([f"{60 / t:.0f}" for t in ticks])
    top.set_xlabel("Residence time $\\tau$ (s)")
    ax.legend(frameon=False, loc="lower left", bbox_to_anchor=(0.06, 0.02), fontsize=4.8)
    style_ax(ax); style_ax(top)
    ins = ax.inset_axes([0.56, 0.14, 0.41, 0.40])
    Vl = np.geomspace(0.01, 30, 400)
    for lab, d0, kd, c in metals:
        ins.plot(Vl, delta(Vl, d0, kd, 0.0), color=c, lw=0.7)
        ins.plot(Vl, delta(Vl, d0, kd, inert), color=c, lw=0.7, ls="--")
    ins.set_xscale("log"); ins.set_ylim(-1.4, 0.1)
    ins.set_xlabel("Drip rate (log)", fontsize=5); ins.tick_params(labelsize=4.5)
    style_ax(ins)
    fig.tight_layout()
    save(fig, "SuppFig4")


# ── SF5: native record + Ni, Co ──────────────────────────────────────────
def sf5(te):
    hr = pd.read_csv(os.path.join(EXT, "drip_rate_summary_hr_censored.csv"))
    m = pd.merge_asof(te.sort_values("depth_cm"), hr.sort_values("depth"),
                      left_on="depth_cm", right_on="depth", tolerance=0.02, direction="nearest")
    m = m.dropna(subset=["pc50"]).sort_values("age_yBP")
    ka = m.age_yBP / 1000
    apply_style()
    fig, axs = plt.subplots(3, 1, figsize=(DOUBLE_COL, DOUBLE_COL * 0.62), sharex=True)
    ax = axs[0]
    ok = m.censored == 0
    ax.fill_between(ka[ok], m.pc25[ok], m.pc75[ok], color=COL_TEAL_200, lw=0, label="25th-75th percentile")
    ax.plot(ka[ok], m.pc50[ok], color=COL_NI, lw=0.6, label="median")
    if (~ok).any():
        ax.plot(ka[~ok], np.ones((~ok).sum()), "v", ms=3, color=COL_DORANGE, lw=0,
                label="censored (≤ 1 drip min$^{-1}$)")
    ax.set_ylabel("Drip rate (drips min$^{-1}$)")
    ax.legend(frameon=False, loc="upper right", ncol=3)
    axs[1].plot(ka, m.Ni, color=COL_NI, lw=0.6)
    axs[1].set_ylabel("Ni (ppm)")
    axs[2].plot(ka, m.Co, color=COL_CO, lw=0.6)
    axs[2].set_ylabel("Co (ppm)")
    axs[2].set_xlabel("Age (ka BP)")
    axs[2].set_xlim(9.7, 0)
    for ax, L in zip(axs, "abc"):
        panel_label(ax, L); style_ax(ax)
    fig.tight_layout(h_pad=0.6)
    save(fig, "SuppFig5")
    print(f"  SF5: {len(m)} native points, {(~ok).sum()} censored")


if __name__ == "__main__":
    te = load_te()
    sf1(te); sf2(); sf3(); sf4(); sf5(te)
