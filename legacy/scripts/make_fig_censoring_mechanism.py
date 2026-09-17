"""
make_fig_censoring_mechanism.py — SI figure: why the drought-minimum points censor.

Mechanism (illustrative). The joint drip-rate posterior is the geometric mean of
the single-proxy Ni and Co posteriors, so it has appreciable mass only where the
two posteriors overlap. For the 156.23 cm drought point (Ni = 10.7, Co = 3.73 ppm)
the two proxies indicate strongly discordant drip rates and their posteriors are
essentially disjoint at every physically admissible dissociation-rate width: the
geometric-mean joint is < 0.02% of the single-proxy peak height at the operating
floor width and remains so even at CLE-scale width. Broadening sigma does not
reconcile the proxies — it slides both peaks further apart and deeper, toward the
drip-stopped limit. There is therefore no width that yields a usable joint
estimate; at the operating floor the honest result is "below resolution", and the
point is reported as censored (<= 1 drip/min). All annotation is in the caption.

Computed from the forward model at the real 156.23 cm Ni/Co concentrations, mu
re-anchored per sigma (calibration median -> 16.66 drips/min). The 8% relative
concentration uncertainty is illustrative; the (non-)overlap behaviour is robust
to it. Illustrates the coupled-inversion mechanism, not Dr Paleo's exact posteriors.
"""
import sys
import os
import numpy as np
from scipy.optimize import brentq
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

FIGDIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "manuscript_figures")
sys.path.insert(0, FIGDIR)
from ngeo_style import apply_style, style_ax, COL_NI, COL_CO, COL_MARKER, DOUBLE_COL, MM
apply_style()

import importlib.util
spec = importlib.util.spec_from_file_location(
    "ck", os.path.join(os.path.dirname(os.path.abspath(__file__)), "calibrate_kd.py"))
mod = importlib.util.module_from_spec(spec); spec.loader.exec_module(mod)
TRZ = np.trapezoid if hasattr(np, "trapezoid") else np.trapz

def h(V, metal, mu, sigma):
    m = mod.METALS[metal]
    Xa = (1.0 - m["inertF"]) * m["aq_ppb"] / (1e6 * m["mw"])
    Ya = mod.CA_PPB / (1e6 * 40.078)
    K0 = (m["Kp"] * mod.Y_S) * (Xa / Ya); nS = 1.0 - 0.01
    return K0 * (1.0 - nS * mod._E1(V / 60.0, mu, sigma, 0.01)) * 1e3 * m["mw"]

CAL_MED = {"Ni": 3.6382, "Co": 0.3121}; TARGET = 16.66
def anchor(metal, sigma):
    return brentq(lambda mu: h(TARGET, metal, mu, sigma) - CAL_MED[metal], -15, 5)

V = np.logspace(np.log10(2e-3), np.log10(60), 5000); lnV = np.log(V)
OBS = {"Ni": 10.7, "Co": 3.73}
REL = 0.08
FLOOR = 1.0

def posterior(metal, sigma):
    mu = anchor(metal, sigma)
    hv = np.array([h(v, metal, mu, sigma) for v in V])
    p = np.exp(-((hv - OBS[metal]) ** 2) / (2 * (REL * OBS[metal]) ** 2))
    a = TRZ(p, lnV)
    return p / a if a > 0 else p

SIG_A = np.pi / np.sqrt(6)
SIG_B = 4.5
MAG = 50.0   # overlap magnification for visibility (stated in caption)

fig, axes = plt.subplots(1, 2, figsize=(DOUBLE_COL, 60 * MM), sharey=True)

for ax, sigma, lab in zip(axes, [SIG_A, SIG_B], ["a", "b"]):
    style_ax(ax)
    pN = posterior("Ni", sigma); pC = posterior("Co", sigma)
    hmax = max(pN.max(), pC.max()); gm = np.sqrt(pN * pC)
    vN, vC = V[pN.argmax()], V[pC.argmax()]

    ax.plot(V, pN / hmax, color=COL_NI, lw=1.2)
    ax.plot(V, pC / hmax, color=COL_CO, lw=1.2)
    # The true geometric-mean joint is < 0.02% of the peak height and is invisible at honest
    # scale — that vanishing overlap IS the censoring. Instead of a magnified fill that would
    # overstate it, shade the drip-rate no-man's-land between the two proxies' 5-95% ranges,
    # where a joint would have to sit but essentially no mass exists.
    def band(p):
        c = np.cumsum(p * np.gradient(lnV)); c /= c[-1]
        return np.exp(np.interp(0.05, c, lnV)), np.exp(np.interp(0.95, c, lnV))
    nlo, nhi = band(pN); clo, chi = band(pC)
    inner_lo, inner_hi = min(chi, nhi), max(clo, nlo)
    if inner_lo < inner_hi:
        ax.axvspan(inner_lo, inner_hi, color=COL_MARKER, alpha=0.10, lw=0)
        ax.text(np.sqrt(inner_lo * inner_hi), 0.32, "no overlap", rotation=90,
                fontsize=5.6, color=COL_MARKER, ha="center", va="center", alpha=0.95)

    ax.axvline(FLOOR, ls=":", lw=0.7, color="#9aa0a6")
    ax.set_xscale("log"); ax.set_xlim(V.min(), 60); ax.set_ylim(0, 1.1)
    ax.set_xlabel(r"Inferred drip rate (drips min$^{-1}$)")
    if lab == "a":
        ax.set_ylabel("Single-proxy posterior (norm.)")
    ax.text(0.025, 0.96, lab, transform=ax.transAxes, fontsize=10,
            fontweight="bold", va="top", ha="left")
    ax.annotate("Ni", xy=(vN, 1.0), xytext=(vN, 1.045), fontsize=6.2, color=COL_NI,
                ha="center", va="bottom", clip_on=False)
    ax.annotate("Co", xy=(vC, 1.0), xytext=(vC, 1.045), fontsize=6.2, color=COL_CO,
                ha="center", va="bottom", clip_on=False)

fig.tight_layout(w_pad=1.3)
for fmt in ("pdf", "png"):
    fig.savefig(f"{os.path.dirname(os.path.abspath(__file__))}/../manuscript_figures/output/FigS_censoring_mechanism.{fmt}", dpi=600, bbox_inches="tight")
plt.close(fig)
for sigma in [SIG_A, SIG_B]:
    pN = posterior("Ni", sigma); pC = posterior("Co", sigma)
    hmax = max(pN.max(), pC.max()); gm = np.sqrt(pN * pC)
    print(f"sigma={sigma:.3f}: Ni_peak={V[pN.argmax()]:.3f} Co_peak={V[pC.argmax()]:.4f} "
          f"joint_peak={100*gm.max()/hmax:.4f}% of single-proxy peak")
