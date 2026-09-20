"""
make_fig_sigma_limits.py -- Supplementary Figure 10: the width of the dissociation-
rate distribution (sigma) is bounded below by the kinetics, cannot be resolved from
a drip-rate series near that bound, and only deepens the reconstructed droughts
when it is increased. Numerical demonstration replacing the analytical derivation
of the earlier S14.2 (Gumbel-window algebra), which is retained in the git history.

Panels
  a  Fraction dissociated within the drip interval for a single rate constant
     (sigma -> 0) and for the adopted population width (sigma = pi/sqrt(6)),
     against drip rate. Even one rate constant responds over a factor ~22 in
     drip rate between 10 % and 90 % dissociated.
  b  Width of the drip-rate response (SD of the sensitivity weighting in ln V,
     computed numerically from the forward model) against the population width
     sigma. The response never narrows below the single-rate floor and follows
     the quadrature sum; near the floor the response is insensitive to sigma,
     so sigma cannot be recovered from a drip-rate calibration.
  c  Drip rate inferred for the drought minimum (157.01 cm, A-880) from Ni,
     from Co and jointly, against sigma, with mu re-anchored to 16.66 drips/min
     for every sigma. All three fall monotonically as sigma increases, so the
     reported drought magnitudes are lower bounds.
  (the former panel d, single-proxy posteriors for a censored point, was
     removed 2026-09-21 with the exclusion of the second-laboratory samples;
     no point in the 568-point record falls below the joint resolution limit).

Outputs ../manuscript_figures/output/FigS_sigma_limits.{png,pdf}
"""
import os
import sys
import importlib.util

import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "..", "manuscript_figures")
sys.path.insert(0, FIGDIR)
from ngeo_style import (apply_style, style_ax, panel_label, label_bubble,  # noqa: E402
                        COL_NI, COL_CO, COL_MARKER, COL_BG_900, COL_BG_300, COL_BG_600,
                        COL_TEAL_600, DOUBLE_COL, SINGLE_COL)
spec = importlib.util.spec_from_file_location("ck", os.path.join(HERE, "calibrate_kd.py"))
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
TRZ = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

SIG_FLOOR = np.pi / np.sqrt(6)          # 1.2825
SIG_CAL = 1.39                          # width returned by the modern calibration
TARGET = 16.66                          # 2004-2023 annual-baseflow mean, drips/min
CAL_MED = {"Ni": 3.6382, "Co": 0.3121}  # HS4 calcite medians over the calibration window (ppm)
PT_MIN = {"Ni": 7.4381, "Co": 1.4920}    # 157.01 cm (A-880), drought minimum of the 568-point record (2026-09-21)
# PT_CEN (former censored point at 156.23 cm) removed 2026-09-21: that sample came from the excluded second-laboratory run
REL = 0.08                               # illustrative relative concentration uncertainty
V = np.logspace(np.log10(2e-3), np.log10(60), 4000); lnV = np.log(V)


def E1(V_sec, mu, sd, nF=0.01):
    """Vectorised form of calibrate_kd._E1 (same grid and weights)."""
    V_sec = np.atleast_1d(np.asarray(V_sec, float))
    k = np.linspace(mu - mod.KLIM * sd, mu + mod.KLIM * sd, mod.KRES)
    krsw = 0.5 * np.r_[k[1] - k[0], k[2:] - k[:-2], k[-1] - k[-2]]
    g = np.exp(-(k - mu) ** 2 / (2 * sd ** 2)) / (sd * np.sqrt(2 * np.pi)) * krsw
    out = np.empty(len(V_sec))
    for i0 in range(0, len(V_sec), 500):
        blk = V_sec[i0:i0 + 500]
        out[i0:i0 + 500] = (np.exp(-np.exp(k)[None, :] / blk[:, None]) * g[None, :]).sum(1)
    return out


def h(V, metal, mu, sigma):
    m = mod.METALS[metal]
    Xa = (1.0 - m["inertF"]) * m["aq_ppb"] / (1e6 * m["mw"])
    Ya = mod.CA_PPB / (1e6 * 40.078)
    K0 = (m["Kp"] * mod.Y_S) * (Xa / Ya); nS = 1.0 - 0.01
    return K0 * (1.0 - nS * E1(np.atleast_1d(V) / 60.0, mu, sigma, 0.01)) * 1e3 * m["mw"]


def anchor(metal, sigma):
    return brentq(lambda mu: h(TARGET, metal, mu, sigma)[0] - CAL_MED[metal], -15, 5)


def posterior(metal, sigma, obs):
    mu = anchor(metal, sigma)
    hv = h(V, metal, mu, sigma)
    p = np.exp(-((hv - obs[metal]) ** 2) / (2 * (REL * obs[metal]) ** 2))
    a = TRZ(p, lnV)
    return p / a if a > 0 else p


def response_width(sigma, metal="Ni"):
    """SD (in ln V) of the normalised sensitivity weighting -dPhi/dlnV."""
    mu = anchor(metal, sigma)
    Vw = np.logspace(-7, 8, 6000); lw = np.log(Vw)     # wide grid: the weighting has long tails
    phi = h(Vw, metal, mu, sigma) / (mod.METALS[metal]["Kp"] * mod.Y_S)  # shape only matters
    w = -np.gradient(phi, lw); w = np.clip(w, 0, None); w /= TRZ(w, lw)
    m1 = TRZ(w * lw, lw)
    return np.sqrt(TRZ(w * (lw - m1) ** 2, lw))


def mode(p):
    return V[np.argmax(p)]


if __name__ == "__main__":
    apply_style()
    fig, axs = plt.subplots(1, 3, figsize=(DOUBLE_COL, DOUBLE_COL * 0.34))

    # (a) single rate constant vs population
    ax = axs[0]
    tau = 60.0 / V
    kd = np.exp(anchor("Ni", 1e-3))          # a single rate constant anchored the same way
    R1 = 1 - np.exp(-kd * tau)
    muP = anchor("Ni", SIG_FLOOR)
    RP = 1 - E1(V / 60.0, muP, SIG_FLOOR, 0.01)
    ax.plot(V, R1, color=COL_BG_900, lw=1.0, label="single rate constant")
    ax.plot(V, RP, color=COL_NI, lw=1.0, label="population, $\\sigma$ = $\\pi/\\sqrt{6}$")
    v10, v90 = np.interp([0.9, 0.1], R1[::-1], V[::-1])
    ax.axvspan(v10, v90, color=COL_BG_300, alpha=0.25, lw=0)
    label_bubble(ax, f"10-90 % dissociated:\nfactor {v90 / v10:.0f} in drip rate",
                 xy=(1.953, 0.2884), xytext=(0.05986, 0.5674), va="center",
                 color=COL_BG_900, bubble_ec=COL_BG_900, arrow_color=COL_BG_900)
    ax.set_xscale("log"); ax.set_xlim(2e-3, 60); ax.set_ylim(0, 1.02)
    ax.set_xlabel("Drip rate (drips min$^{-1}$)")
    ax.set_ylabel("Fraction dissociated in drip interval")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, fontsize=5.2, columnspacing=1.0, handlelength=1.4)
    panel_label(ax, "a")
    print(f"(a) single-k 10-90 % span: {v10:.3f}-{v90:.2f} drips/min, factor {v90 / v10:.1f}")

    # (b) response width vs sigma
    ax = axs[1]
    sig = np.concatenate([[0.02, 0.05, 0.1], np.linspace(0.2, 4.0, 20)])
    wid = np.array([response_width(s) for s in sig])
    ax.plot(sig, np.sqrt(sig ** 2 + SIG_FLOOR ** 2), color=COL_BG_300, lw=0.8,
            label="quadrature sum $\\sqrt{\\sigma^2 + (\\pi/\\sqrt{6})^2}$")
    ax.plot(sig, wid, "o", ms=2.6, color=COL_NI, lw=0, label="forward model (numerical)")
    ax.axhline(SIG_FLOOR, color=COL_BG_900, lw=0.5, ls=":")
    ax.axvline(SIG_FLOOR, color=COL_BG_900, lw=0.5, ls=":")
    label_bubble(ax, "single-rate\nfloor $\\pi/\\sqrt{6}$",
                 xy=(2.785, SIG_FLOOR), xytext=(3.361, 2.184), va="center",
                 color=COL_BG_900, bubble_ec=COL_BG_900, arrow_color=COL_BG_900)
    ax.axvline(SIG_CAL, color=COL_MARKER, lw=0.6, ls="--")
    label_bubble(ax, f"modern calibration\n$\\sigma$ = {SIG_CAL}",
                 xy=(SIG_CAL, 2.5), xytext=(2.431, 3.869), va="center",
                 color=COL_MARKER, bubble_ec=COL_MARKER, arrow_color=COL_MARKER)
    ax.set_xlim(0, 4.05); ax.set_ylim(0, 4.4)
    ax.set_xlabel("Population width $\\sigma$ (ln $k_d$ units)")
    ax.set_ylabel("Width of drip-rate response (ln units)")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=2, fontsize=5.2, columnspacing=1.0, handlelength=1.4)
    panel_label(ax, "b")
    i0 = np.argmin(abs(sig - 0.02)); ic = np.argmin(abs(sig - SIG_FLOOR))
    print(f"(b) response width at sigma 0.02: {wid[i0]:.3f}; at sigma 1.28: {wid[ic]:.3f} "
          f"(quadrature {np.sqrt(sig[ic]**2 + SIG_FLOOR**2):.3f}); at sigma 4: {wid[-1]:.3f}")

    # (c) inferred drip rate at the resolved minimum vs sigma
    ax = axs[2]
    sig_c = np.linspace(0.3, 4.0, 16)
    vN, vC, vJ = [], [], []
    for s in sig_c:
        pN = posterior("Ni", s, PT_MIN); pC = posterior("Co", s, PT_MIN)
        vN.append(mode(pN)); vC.append(mode(pC)); vJ.append(mode(np.sqrt(pN * pC)))
    ax.plot(sig_c, vN, color=COL_NI, lw=0.9, label="Ni alone")
    ax.plot(sig_c, vC, color=COL_CO, lw=0.9, label="Co alone")
    ax.plot(sig_c, vJ, color=COL_BG_900, lw=1.2, label="joint")
    ax.axvline(SIG_FLOOR, color=COL_BG_900, lw=0.5, ls=":")
    ax.axhline(TARGET, color=COL_BG_300, lw=0.5, ls="--")
    label_bubble(ax, "modern baseflow", xy=(2.006, 15.22), xytext=(2.961, 5.115), va="center",
                 color=COL_BG_600, bubble_ec=COL_BG_600, arrow_color=COL_BG_600)
    ax.set_yscale("log"); ax.set_xlim(0, 4.05)
    ax.set_xlabel("Population width $\\sigma$ (ln $k_d$ units)")
    ax.set_ylabel("Inferred drip rate, 157.01 cm\n(drips min$^{-1}$)")
    ax.legend(frameon=False, loc="upper center", bbox_to_anchor=(0.5, -0.22), ncol=3, fontsize=5.2, columnspacing=1.0, handlelength=1.4)
    panel_label(ax, "c")
    j = np.argmin(abs(sig_c - SIG_FLOOR))
    print(f"(c) 157.01 cm at sigma {sig_c[j]:.2f}: Ni {vN[j]:.2f} Co {vC[j]:.2f} joint {vJ[j]:.2f}; "
          f"at sigma 4.0: Ni {vN[-1]:.3f} Co {vC[-1]:.3f} joint {vJ[-1]:.3f}")

    for ax in axs:
        style_ax(ax)
    fig.tight_layout(w_pad=0.6)
    fig.subplots_adjust(bottom=0.20, left=0.055, right=0.985, wspace=0.34)
    for fmt in ("png", "pdf"):
        fig.savefig(os.path.join(FIGDIR, "output", f"FigS_sigma_limits.{fmt}"))
    print("wrote FigS_sigma_limits.png/.pdf")
