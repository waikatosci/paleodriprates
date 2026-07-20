#!/usr/bin/env python3
"""
calibrate_kd.py
===============
Calibration of the dissociation-rate distribution for the Co/Ni kinetic
drip-rate proxy (Hartland et al., NCOMMS-26-041445-T).

WHAT IS FITTED
    mu  (= mean of ln k_d), one value per metal, is fitted so that the forward
        model applied to the median HS4 metal concentration reproduces the
        observed modern drip-rate level (cave-monitoring mean).

WHAT IS NOT FITTED
    sigma (= width of the ln k_d distribution) is FIXED by the dissociation-rate
        sampling window, sigma = pi/sqrt(6) ~= 1.283 (Supplementary Methods 1).
        It is a kinetic constant, identical for every metal, and is NOT tuned to
        any calibration target. (The legacy runs used a rounded value, 1.385.)

    lambda_F (= fast/instantly-labile fraction, nF) is FIXED at 0.01, matching the
        published inversion (Drip_rate.xlsx) and MS Fig. 4 / Supplementary Table S1.

ERROR PROPAGATION (feeds Reviewer 3, Major 1 — uncertainty at each stage)
    mu uncertainty is obtained by Monte Carlo: the aqueous medians (Ni, Co, Ca),
    the partition coefficients (K_p), the target drip level, and the HS4 metal
    median (bootstrap) are each resampled within their uncertainties; mu is
    refitted on every draw; the spread of the resulting mu ensemble is reported.

MODES
    faithfulness : sigma=1.385, target=16.0 drips/min (2005-2015 monitoring)
                   -> regression test; should reproduce published mu
                      (mu_Ni ~= -3.91, mu_Co ~= -5.57).
    refit        : sigma=1.283, target=15.09 drips/min (2004-2023 monitoring)
                   -> revised calibration for the R2.1 response.

OUTPUT
    kd_calibration_parameters.csv
"""
import argparse, numpy as np, pandas as pd, warnings
from scipy.optimize import brentq
warnings.filterwarnings("ignore")

# ------------------------------------------------------------------ constants
VMIN, VMAX, VRES = 0.02, 100.0, 5000     # drip-rate grid (drips/min), per params.py
KRES, KLIM       = 5000, 20.0            # ln k_d integration grid
Y_S              = 400.0 / 40.078        # mol Ca per kg calcite (40 wt% Ca)
SIGMA_KINETIC    = np.pi / np.sqrt(6.0)  # = 1.2825 : Gumbel-window SD (fixed sigma)

# per-metal recipe (aqueous ppb medians, inert fraction lambda_I, K_p, mol wt)
METALS = {
    "Ni": dict(Kp=1.1, aq_ppb=4.370, aq_cv=0.20, inertF=0.10, mw=58.693),
    "Co": dict(Kp=4.4, aq_ppb=0.4605, aq_cv=0.25, inertF=0.40, mw=58.9332),
}
CA_PPB, CA_CV = 67437.0, 0.10            # aqueous Ca median (ppb) and its CV
KP_CV         = 0.15                     # partition-coefficient uncertainty
CALIB_WIN     = (1953, 2012)            # growth window over which metal median is taken


# ------------------------------------------------------------ forward model
def forward_V(Xs_mol, Kp, Xa_mol, Ya_mol, mu, sd, nF=0.01):
    """Vectorised inverse of the kinetic model: drip rate (drips/min) for an
    array of calcite metal concentrations Xs_mol (mol/L).  Matches
    model.dr_timeseries to 3 dp but avoids per-point root finding."""
    nS = 1.0 - nF
    K0 = (Kp * Y_S) * (Xa_mol / Ya_mol)
    Vsec = np.linspace(VMIN, VMAX, VRES) / 60.0
    k = np.linspace(mu - KLIM * sd, mu + KLIM * sd, KRES)
    krsw = 0.5 * np.r_[k[1] - k[0], k[2:] - k[:-2], k[-1] - k[-2]]
    integ = ((1.0 / (sd * np.sqrt(2 * np.pi)))
             * np.exp(-np.exp(k[None, :]) / Vsec[:, None])
             * np.exp(-(k[None, :] - mu) ** 2 / (2 * sd ** 2)))
    Xs_model = K0 * (1.0 - nS * (integ * krsw[None, :]).sum(1))
    o = np.argsort(Xs_model)
    return np.interp(np.atleast_1d(Xs_mol), Xs_model[o], (Vsec * 60.0)[o])


def _E1(V_sec, mu, sd, nF=0.01):
    """Expected surviving (bound) fraction E1 = int N(k;mu,sd) exp(-e^k/V) dk.
    1-D quadrature over ln k_d only — the fast kernel behind the fit."""
    k = np.linspace(mu - KLIM * sd, mu + KLIM * sd, KRES)
    krsw = 0.5 * np.r_[k[1] - k[0], k[2:] - k[:-2], k[-1] - k[-2]]
    integ = ((1.0 / (sd * np.sqrt(2 * np.pi)))
             * np.exp(-np.exp(k) / V_sec)
             * np.exp(-(k - mu) ** 2 / (2 * sd ** 2)))
    return float((integ * krsw).sum())


def fit_mu(metal_median_ppm, Kp, aq_ppb, inertF, mw, target_drip, sigma, nF=0.01):
    """Return mu such that the forward model at the median metal reproduces the
    target drip level.  Solved directly at V = target via the E1 relation
    Xs = K0 (1 - n_S E1), so only a 1-D integral is evaluated per iteration."""
    Xa = (1.0 - inertF) * aq_ppb / (1e6 * mw)
    Ya = CA_PPB_LOCAL[0] / (1e6 * 40.078)
    K0 = (Kp * Y_S) * (Xa / Ya)
    Xs = metal_median_ppm / (1e3 * mw)
    nS = 1.0 - nF
    target_E1 = (1.0 - Xs / K0) / nS          # required surviving fraction at V=target
    if not (0.0 < target_E1 < 1.0):
        raise ValueError("target outside model range")
    Vsec = target_drip / 60.0
    g = lambda mu: _E1(Vsec, mu, sigma, nF) - target_E1   # E1 decreasing in mu
    return brentq(g, -12.0, 4.0, xtol=1e-5)


# a tiny mutable holder so fit_mu can see the current Ca draw during MC
CA_PPB_LOCAL = [CA_PPB]


# ------------------------------------------------------------------- loader
def load_metal_medians(te_path):
    raw = pd.read_excel(te_path, sheet_name="HS4-Wuhan", header=None)
    df = raw.iloc[4:].copy()
    df.columns = ["age", "x1", "Co", "CoCorr", "Ni", "Cu", "x6", "x7"]
    for c in df.columns:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["age"])
    df["year"] = 2000 - df["age"]               # 'yr before 2000' -> calendar year
    w = df[(df.year >= CALIB_WIN[0]) & (df.year <= CALIB_WIN[1])]
    return {"Ni": w["Ni"].dropna().values, "Co": w["Co"].dropna().values}


# ---------------------------------------------------------------- MC engine
def calibrate(metal_samples, target_drip, target_sd, sigma, n_mc=400, seed=0):
    rng = np.random.default_rng(seed)
    out = {}
    for el, cfg in METALS.items():
        s = metal_samples[el]
        med0 = np.median(s)
        # point estimate at nominal inputs
        CA_PPB_LOCAL[0] = CA_PPB
        mu0 = fit_mu(med0, cfg["Kp"], cfg["aq_ppb"], cfg["inertF"], cfg["mw"],
                     target_drip, sigma)
        # Monte-Carlo propagation
        draws = np.empty(n_mc)
        for i in range(n_mc):
            med  = np.median(rng.choice(s, size=len(s), replace=True))   # metal bootstrap
            aq   = cfg["aq_ppb"] * np.exp(rng.normal(0, cfg["aq_cv"]))    # aqueous metal
            kp   = cfg["Kp"]     * np.exp(rng.normal(0, KP_CV))           # partition coeff
            CA_PPB_LOCAL[0] = CA_PPB * np.exp(rng.normal(0, CA_CV))       # aqueous Ca
            tgt  = target_drip + rng.normal(0, target_sd)                 # drip level
            try:
                draws[i] = fit_mu(med, kp, aq, cfg["inertF"], cfg["mw"], tgt, sigma)
            except ValueError:
                draws[i] = np.nan
        draws = draws[np.isfinite(draws)]
        out[el] = dict(mu=mu0, mu_sd=np.std(draws, ddof=1),
                       mu_lo=np.percentile(draws, 2.5),
                       mu_hi=np.percentile(draws, 97.5),
                       metal_median_ppm=med0, n_samples=len(s))
    CA_PPB_LOCAL[0] = CA_PPB
    return out


# ----------------------------------------------------------------- reporting
def run(mode, te_path):
    if mode == "faithfulness":
        sigma, target, target_sd = 1.385, 16.0, 0.62      # 2005-2015 monitoring
        published = {"Ni": -3.91, "Co": -5.57}
    else:  # refit
        sigma, target, target_sd = SIGMA_KINETIC, 15.09, 0.60   # 2004-2023 monitoring
        published = {"Ni": None, "Co": None}

    samples = load_metal_medians(te_path)
    res = calibrate(samples, target, target_sd, sigma)

    rows = []
    print(f"\n=== mode={mode}  sigma={sigma:.4f}  target={target} drips/min ===")
    for el in ("Ni", "Co"):
        r = res[el]
        pub = published[el]
        tag = f"(published {pub:+.2f}, delta {r['mu']-pub:+.3f})" if pub else ""
        print(f"  mu_{el} = {r['mu']:+.3f} +/- {r['mu_sd']:.3f}  "
              f"[95% {r['mu_lo']:+.3f}, {r['mu_hi']:+.3f}]  {tag}")
        rows.append(dict(mode=mode, metal=el, parameter="mu_ln_kd",
                         value=round(r["mu"], 4), sd=round(r["mu_sd"], 4),
                         ci95_lo=round(r["mu_lo"], 4), ci95_hi=round(r["mu_hi"], 4),
                         sigma_ln_kd=round(sigma, 4),
                         sigma_source=("pi/sqrt(6) kinetic window (fixed)"
                                       if abs(sigma - SIGMA_KINETIC) < 1e-3
                                       else "legacy fixed 1.385"),
                         target_drip_min=target,
                         metal_median_ppm=round(r["metal_median_ppm"], 4),
                         n_metal_samples=r["n_samples"],
                         published_mu=("" if pub is None else pub)))
    return rows


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--te", default="../Calibration_data_updated/TE_dat_HS4_HS6_age.xls")
    ap.add_argument("--out", default="kd_calibration_parameters.csv")
    args = ap.parse_args()
    all_rows = run("faithfulness", args.te) + run("refit", args.te)
    pd.DataFrame(all_rows).to_csv(args.out, index=False)
    print(f"\nwrote {args.out}")
