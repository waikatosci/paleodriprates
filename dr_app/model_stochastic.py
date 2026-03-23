"""
Stochastic forward model for drip rate PDF with concentration uncertainty
========================================================================

Companion to model.py — does NOT modify the original.

The key optimization: E_1 and E_2 (the expensive k-integration) depend
only on the Kd distribution and the V grid, NOT on the aqueous
concentrations. Since h(V) = K_0 × (1 − nS × E_1) is linear in K_0,
we precompute E_1 and E_2 once, then for each concentration sample we
just rescale by that sample's K_0 and re-evaluate the interp1d PDFs.

This reduces cost from O(N × KRES × VRES) to
    O(KRES × VRES) + O(N × VRES × T_steps)
— a ~100× speedup over naively calling dr_pdfseries N times.

Usage (from driprates_stochastic.py):
    from model_stochastic import dr_pdfseries_stochastic
    V_pdf, V_span_min = dr_pdfseries_stochastic(
        Xs_pdf, Kp, Kd_mn, Kd_sd, K_e, n_F,
        ca_samples_mol, aq_samples_mol, inertF)
"""

import numpy as np
np.seterr(over="ignore")    # match model.py

from params import VMIN, VMAX, VRES, KRES, KLIM


def dr_pdfseries_stochastic(
    Xs_pdf,                    # list of interp1d functions (T_steps,)
    Kp: float,                 # partition coefficient
    Kd_mn: float,              # mean of ln(Kd)
    Kd_sd: float,              # std dev of ln(Kd)
    K_e: float,                # OMC dissociation efficiency multiplier
    n_F: float,                # fast fraction (XF)
    ca_samples_mol: np.ndarray,  # shape (N,) — Ca in mol/L per sample
    aq_samples_mol: np.ndarray,  # shape (N,) — TE in mol/L per sample
    inertF: float = 0.0,        # inert fraction (for Xa = (1-InertF)*aq)
) -> tuple:
    """
    Vectorised forward model over N sampled concentration pairs.

    Returns
    -------
    V_pdf : ndarray, shape (VRES, T_steps)
        Marginalised (geometric mean over N samples) drip rate PDF.
    V_span_min : ndarray, shape (VRES,)
        Drip rate axis in drips per minute.
    """
    N_samples = len(ca_samples_mol)
    T_steps   = len(Xs_pdf)

    # ── Kd integration grid (computed once) ──────────────────────────
    k_mu, k_sd = Kd_mn, Kd_sd
    k_span = np.linspace(k_mu - KLIM * k_sd, k_mu + KLIM * k_sd, KRES)
    k_rsw  = 0.5 * np.r_[
        k_span[1] - k_span[0],
        k_span[2:] - k_span[:-2],
        k_span[-1] - k_span[-2],
    ]

    # ── V grid (drips per second internally, converted at output) ────
    V_span = np.linspace(VMIN, VMAX, VRES) / 60.0
    V_rsw  = 0.5 * np.r_[
        V_span[1] - V_span[0],
        V_span[2:] - V_span[:-2],
        V_span[-1] - V_span[-2],
    ]

    # ── Precompute E_1(v) and E_2(v) — independent of concentrations ─
    nk, nv = KRES, VRES
    k_arr = k_span.repeat(nv).reshape(nk, nv).T       # (nv, nk)
    v_arr = V_span.repeat(nk).reshape(nv, nk)          # (nv, nk)

    # Gaussian envelope
    gauss = np.exp(-(k_arr - k_mu) ** 2 / (2. * k_sd ** 2))
    gauss /= (k_sd * np.sqrt(2. * np.pi))

    # exp(−exp(k)/v)
    exp_kv = np.exp(-np.exp(k_arr) / v_arr)

    # E_1(v) = ∫ gauss(k) × exp(−exp(k)/v) dk
    E_1 = (exp_kv * gauss * k_rsw).sum(axis=1)         # (nv,)

    # E_2(v) = ∫ gauss(k) × (exp(k)/v²) × exp(−exp(k)/v) dk
    E_2 = (exp_kv * gauss * np.exp(k_arr) / v_arr ** 2 * k_rsw).sum(axis=1)

    # ── Constants ────────────────────────────────────────────────────
    Y_s = 400.0 / 40.078                               # mol Ca per kg calcite
    n_S = 1.0 - n_F

    # ── Accumulator for log-PDF (geometric mean marginalisation) ─────
    log_pdf_sum = np.zeros((VRES, T_steps))
    n_valid     = np.zeros(T_steps, dtype=int)

    for s in range(N_samples):
        # Per-sample K_0
        Xa_s = (1.0 - inertF) * aq_samples_mol[s]
        Ya_s = ca_samples_mol[s]
        if Ya_s <= 0 or Xa_s <= 0:
            continue
        K_0 = (Kp * Y_s) * (Xa_s / Ya_s) * K_e

        # h(v) and h'(v) for this K_0 — just scalar rescaling
        h_v  = K_0 * (1.0 - n_S * E_1)                # (nv,)
        hp_v = -K_0 * n_S * E_2                        # (nv,)

        # Per-timestep PDF evaluation
        V_pdf_s = np.zeros((VRES, T_steps))
        for t in range(T_steps):
            pdf = -Xs_pdf[t](h_v) * hp_v
            pdf_sum = np.nansum(pdf * V_rsw)
            if pdf_sum > 0 and np.isfinite(pdf_sum):
                pdf /= pdf_sum
                V_pdf_s[:, t] = pdf

        # Accumulate log for geometric mean (skip zeros/NaN)
        mask = V_pdf_s > 0
        with np.errstate(divide='ignore', invalid='ignore'):
            log_pdf = np.where(mask, np.log(V_pdf_s), 0.0)
        log_pdf_sum += log_pdf
        # Track which timesteps got valid contributions
        valid_cols = np.any(mask, axis=0)
        n_valid += valid_cols.astype(int)

    # ── Geometric mean marginalisation ───────────────────────────────
    # Where n_valid > 0, divide accumulated log by count; else leave as 0
    with np.errstate(divide='ignore', invalid='ignore'):
        n_safe = np.where(n_valid > 0, n_valid, 1)
        V_pdf = np.exp(log_pdf_sum / n_safe)

    # Renormalise each timestep
    C = (V_pdf.T * V_rsw).sum(axis=1)
    C[(C == 0) | ~np.isfinite(C)] = 1.0
    V_pdf = V_pdf / C

    return V_pdf, V_span * 60.0
