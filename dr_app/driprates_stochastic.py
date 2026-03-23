"""
Stochastic driprates wrapper
=============================

Drop-in replacement for drip_rate_util.driprates() that adds support
for lognormal concentration priors. Import as:

    from driprates_stochastic import driprates

When TE dict contains 'ca_prior' and/or 'aq_prior' keys (instances of
ConcentrationPrior), the function samples N concentration pairs and
marginalises. Otherwise it delegates to the original scalar path.

The function signature is backward-compatible: all existing call sites
work unchanged.
"""

import numpy as np

import model
from concentration_prior import ConcentrationPrior, N_CONC_SAMPLES
from model_stochastic import dr_pdfseries_stochastic


def driprates(Kd_mn, Kd_sd, K_e, ConcAq=1.0, TE=None,
              calib=False, pbar_on=False, n_conc=None):
    """
    Returns drip rate PDFs for a given choice of mean Kd.

    Extended from drip_rate_util.driprates to support stochastic
    aqueous concentrations via ConcentrationPrior objects in the
    TE dict ('ca_prior', 'aq_prior' keys).

    Parameters
    ----------
    Kd_mn, Kd_sd : float
        Mean and std dev of ln(Kd) — OMC dissociation rate distribution.
    K_e : float
        OMC dissociation efficiency multiplier.
    ConcAq : float
        Legacy aqueous concentration parameter (unused in normal path).
    TE : dict
        Trace element parameter dictionary. Must contain 'PDist',
        'mol_wt', 'Kp', 'F', 'InertF', 'aq_conc', 'ca_conc', etc.
        Optionally: 'ca_prior' and 'aq_prior' (ConcentrationPrior).
    calib : bool
        If True, restrict to calibration period (1953–2012).
    pbar_on : bool
        Show progress bar.
    n_conc : int or None
        Number of concentration samples for stochastic mode.
        Defaults to N_CONC_SAMPLES (50).

    Returns
    -------
    V_pdf : ndarray, shape (VRES, T_steps)
        Drip rate probability density at each timestep.
    age : ndarray, shape (T_steps,)
        Age axis in years BP.
    V_span : ndarray, shape (VRES,)
        Drip rate axis in drips per minute.
    """
    n_conc = n_conc or N_CONC_SAMPLES

    # ── Detect stochastic mode ───────────────────────────────────────
    ca_prior = TE.get('ca_prior') if TE else None
    aq_prior = TE.get('aq_prior') if TE else None

    if ca_prior is None and aq_prior is None:
        # Scalar path — delegate to original implementation
        return _driprates_scalar(Kd_mn, Kd_sd, K_e, ConcAq, TE,
                                 calib, pbar_on)

    # ── Stochastic path ──────────────────────────────────────────────

    # 1. Get proxy PDist (call te_pdfseries once, with proxyspan protection)
    _proxyspan_backup = TE['PDist'].proxyspan.copy()
    age, Xs_pdf = model.te_pdfseries(TE)
    TE['PDist'].proxyspan = _proxyspan_backup      # restore after mutation
    Xs_pdf = np.array(Xs_pdf, dtype="object")

    # Calibration period restriction (same logic as original)
    if calib:
        Yi, Yf = 1953, 2012
        age = (1950. - age).astype("int")
        age = age[::-1]
        Xs_pdf = Xs_pdf[::-1]
        i = (age >= Yi) * (age <= Yf)
        age, Xs_pdf = age[i], Xs_pdf[i]

    # 2. Sample concentration pairs
    rng = np.random.default_rng()       # thread-safe fresh generator

    if ca_prior is not None:
        ca_ppb = ca_prior.sample(n_conc, rng)
    else:
        ca_ppb = np.full(n_conc, float(TE['ca_conc']))

    if aq_prior is not None:
        aq_ppb = aq_prior.sample(n_conc, rng)
    else:
        aq_ppb = np.full(n_conc, float(TE['aq_conc']))

    # Convert ppb → mol/L (same conversion as original driprates)
    mol_wt = TE['mol_wt']
    aq_mol = aq_ppb / (1E6 * mol_wt)     # TE in mol/L
    ca_mol = ca_ppb / (1E6 * 40.078)     # Ca in mol/L

    # 3. Model parameters
    Kp     = TE['Kp']
    if Kp < 0:
        T = TE['Temp_C']
        Kp = model.kp_theory(T, TE['elem'])

    inertF = TE['InertF']
    etaF   = TE['F']

    # 4. Run stochastic forward model (precomputes E1/E2 once)
    V_pdf, V_span = dr_pdfseries_stochastic(
        Xs_pdf, Kp, Kd_mn, Kd_sd, K_e, etaF,
        ca_mol, aq_mol, inertF,
    )

    return V_pdf, age, V_span


def _driprates_scalar(Kd_mn, Kd_sd, K_e, ConcAq, TE, calib, pbar_on):
    """
    Original scalar driprates path — reproduces drip_rate_util.driprates
    exactly, but with proxyspan protection.
    """
    # Protect proxyspan from in-place mutation in te_pdfseries
    _proxyspan_backup = TE['PDist'].proxyspan.copy()
    age, Xs_pdf = model.te_pdfseries(TE)
    TE['PDist'].proxyspan = _proxyspan_backup
    Xs_pdf = np.array(Xs_pdf, dtype="object")

    if calib:
        Yi, Yf = 1953, 2012
        age = (1950. - age).astype("int")
        age = age[::-1]
        Xs_pdf = Xs_pdf[::-1]
        i = (age >= Yi) * (age <= Yf)
        age, Xs_pdf = age[i], Xs_pdf[i]

    inertF = TE['InertF']
    Xa = (1.0 - inertF) * TE['aq_conc']
    Ya = TE['ca_conc']
    Xa /= (1E6 * TE['mol_wt'])
    Ya /= (1E6 * 40.078)

    etaF = TE['F']

    Kp = TE['Kp']
    if Kp < 0:
        T = TE['Temp_C']
        Kp = model.kp_theory(T, TE['elem'])

    V_pdf, V_span = model.dr_pdfseries(
        Xs_pdf, Xa, Ya, Kp, Kd_mn, Kd_sd, K_e, etaF, pbar=pbar_on)

    return V_pdf, age, V_span
