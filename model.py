"""
CORE FUNCTIONS THAT MODEL DRIPRATES BASED ON TRACE ELEMENT CONCENTRATIONS
=========================================================================

Let,
    V   := driprate,
    X_s := concentration of trace metal in calcite,
    Y_s := concentration of calcium in calcite,
    X_a := concentration of trace element in aqueous solution,
    Y_a := concentration of calcium in aqueous solution,
    n_S := slow (kinetically fractionated) fraction of trace element,
    K_p := partition coefficient of trace element into calcite, and
    K_d := decay rate of trace elements from humic acid complexes in solution.

We have two main relations,

    V   = g(X_s), and
    X_s = h(V),

such that h == g^-1, as g is strictly decreasing in X_s.

The second relation, X_s = h(V), is expressed as :

    X_s = h(V)  = K_0 (1 - n_S E_1),

where,

    K_0 = (K_p Y_s X_a) / Y_a, and
    E_1 = E_fk[exp(- e^k/v)]
        = \int (1/(k_sd\sqrt(2pi))) exp(- (k - k_mu) / (2k_sd^2))
                exp(- e^k/v) dk;

    k = np.log(K_d),
    k_sd = np.log(1 + (K_d_sd / K_d_mu) ^2), and
    k_mu = np.log(K_d_mu / \sqrt(1 + (K_d_sd / K_d_mu)^2); and

    K_d_mu, and K_d_sd are the mean and standard deviations of decay rates of
    the population of decay rates of humic acid complexes in the solution.

Now, if we know the marginal probabilities p(x_s|t) of the trace element
concentrations, we can derive the marginal probability p(v|t),

    p(v|t) = - p(h(v)|t) h'(v),

where,

    h'(v) = - K_0 n_S E_2

    E_2 = E_fk[(e^k/v^2) exp(- e^k / v)]
        = \int  (1/(k_sd\sqrt(2pi))) exp(- (k - k_mu) / (2k_sd^2))
                (e^k/v^2) exp(- e^k/v) dk

"""
# Created: Mon Jul 15, 2019  01:58pm
# Last modified: Mon Nov 02, 2020  03:48pm
# Copyright: Bedartha Goswami <bedartha.goswami@uni-tuebingen.de>


import numpy as np
np.seterr(over="ignore")    # ignore overflow error for np.exp; integrand_1
import pandas as pd

from scipy.integrate import quad
from scipy.optimize import root_scalar
from scipy.optimize import minimize

import utils


from params import DATPATH, FIGPATH
from params import OUTLIER_WINDOW_SIZE as OWS
from params import BAYPROX_SAMPLING_RES as SRES
from params import VMIN, VMAX, VRES
from params import KRES, KLIM


def te_timeseries(TE, series):
    """
    Returns the trace element time series estimated using BayProX.
    """
    FIN = DATPATH + "output/002_bayprox_%s_ows%d_sres%d_proxyrecord.csv" \
                    % (TE, OWS, SRES)
    DAT = pd.read_csv(FIN, delimiter=",")
    age = DAT["Age (yBP)"].to_numpy()
    Xs = DAT["%s (ppb) %s" % (TE, series)].to_numpy()
    molwt = {
                "Ni": 58.693,
                "Co": 58.933,
            }
    Xs /= (1E3 * molwt[TE])        # convert to mol / L

    return age, Xs


def aqueous_concentrations(TE, mode="median"):
    """
    Returns aqueous concentrations from monitoring data.
    """
    DAT = np.genfromtxt(DATPATH + "processed/Monitoring_NiCoCa.csv",
                        delimiter=",", names=True, usecols=(2, 3, 4),
                        missing_values=("NA"), filling_values=np.nan,
                        autostrip=True)
    # set variables
    Xa = DAT["%s_ppb" % TE]
    Ya = DAT["Ca_ppb"]
    # get the quantities in mol / L and drip rate in drips / min
    # x ppb => Z g solute in 10^9 g solution
    #       => Z / 10^6 g solute in 1000 g solution ~ 1L solution
    #       => Z / (10^6 x Atomic wt.) mol in 1L solution
    if TE=='Ni':
        Xa /= (1E6 * 58.693)        # := atomic wt of nickel
    elif TE=='Co':
        Xa /= (1E6 * 58.9332)        # := atomic wt of nickel
    else:
        Xa /= (1E6 * 58.693)        # := atomic wt of nickel
    Ya /= (1E6 * 40.078)


    out = {
           "mean": (np.mean(Xa), np.mean(Ya)),
           "median": (np.median(Xa), np.median(Ya)),
           "min": (np.min(Xa), np.min(Ya)),
           "max": (np.max(Xa), np.max(Ya)),
           "none": (Xa, Ya),
          }

    return out[mode]



def cave_temperature(mode="median"):
    """
    Returns the cave temperature for Heshang from monitoring data.
    """
    FN = DATPATH + "processed/.older/hs4_monitoring_climatic_observables.csv"
    T = np.genfromtxt(FN,
                      delimiter=",", skip_header=2, usecols=(1), dtype=float)
    T += 273.15

    out = {
           "mean": np.mean(T),
           "median": np.median(T),
           "min": np.min(T),
           "max": np.max(T),
          }

    return out[mode]


def kp_theory(T, TE):
    """
    Returns theoretical partition coefficient as per Wang and Xu (2001).
    """
    R = 8.314472            # kg m^2 / s^2 / K / mol := J / K / mol
    R *= 0.0002390057       # convert R to units of kcal / K / mol
    a_MHX = 0.968           # dimensionless
    b_MHX = 75.168          # kcal / mol / A
    # # the following are rounded off values that Wang and Xu use in their
    # # estimates
    # a_MHX = 0.968           # dimensionless
    # b_MHX = 75.170          # kcal / mol / A

    # Ca
    r_Ca = 1                # A / r_Ca
    delG_f_Ca = -132.12     # kcal / mol
    delG_n_Ca = -10.83      # kcal / mol

    params = {
                "Cu": (
                        0.73,           # r_        A / r_Ca
                        15.55,          # delG_f_   kcal / mol
                        160.39          # delG_n_   kcal / mol
                    ),
                "Ni": (
                        0.70,           # r_        A / r_Ca
                        -10.90,         # delG_f_   kcal / mol
                        136.86          # delG_n_   kcal / mol
                    ),
                "Co": (
                        0.735,          # r_        A / r_Ca
                        -13.00,         # delG_f_   kcal / mol
                        131.36          # delG_n_   kcal / mol
                    ),
            }
    r_, delG_f_, delG_n_ = params[TE]

    # calculate the three main terms
    term1 = a_MHX * (delG_n_ - delG_n_Ca)
    term2 = b_MHX * (r_ - r_Ca)
    term3 = delG_f_ - delG_f_Ca

    # calculate log(Kp)
    log_Kp = term1 + term2 - term3
    log_Kp /= (-2.303  * R * T)
    val = np.power(10., log_Kp)
    # val /= 3.0
    if TE=="Co":
        val = 4.4
    elif TE=="Ni":
        val = 1.1

    return val


def dr_timeseries(X_s, X_a, Y_a, K_p, K_d_mu, K_d_sd, n_F=0.001, pbar=False):
    """
    Returns the driprate for given trace element concentration in calcite.
    """
    def __integrand_1(k, v):
        """
        The integrand inside expectation E1 above.

        Variables which are used from parent function namespace:
            k_mu, k_sd
        """
        term0 = - np.exp(k) / v
        # term0 = - np.exp(k) / (v/60.0)
        term1 = np.exp(term0)
        term2 = np.exp(- (k - k_mu) ** 2 / (2. * k_sd ** 2))
        const = 1. / (k_sd * np.sqrt(2. * np.pi))
        return const * term1 * term2


    def __expectation_1(v):
        """
        Expectation E_1 = E_fk[exp(-k/v)] above.

        """
        k_span = np.linspace(k_mu - KLIM * k_sd, k_mu + KLIM * k_sd, KRES)
        k_rsw = 0.5 * np.r_[
                            k_span[1] - k_span[0],
                            k_span[2:] - k_span[:-2],
                            k_span[-1] - k_span[-2]
                           ]
        f_k = __integrand_1(k_span, v)
        return (f_k * k_rsw).sum()


    def __residual(v):
        """
        Function for which we need to find roots.

        Variables which are used from parent function namespace:
            x_s_CURR, K_0, n_S
        """
        return x_s_CURR - K_0 * (1. - n_S * __expectation_1(v))

    N = X_s.shape[0]
    if np.isscalar(X_a):
        X_a = X_a.repeat(N)
    if np.isscalar(Y_a):
        Y_a = Y_a.repeat(N)
    assert X_a.shape == Y_a.shape, "Aqueous concentrations are different sizes"

    # for the solid Ca in calcite, it is known that Ca is 40% by weight in
    # calcite, so ---
    #       40 g of Ca in 100 g of calcite
    #       400 g of Ca in 1000 g of calcite equiv to 1 L of calcite
    #       400 / 40.078 mol of Ca in 1 kg of calcite
    Y_s = 400. / 40.078
    n_S = 1. - n_F
    # k_sd = np.log(1. + (K_d_sd / K_d_mu) ** 2)
    # k_mu = np.log(K_d_mu / np.sqrt(1. + (K_d_sd / K_d_mu) ** 2))
    k_sd = K_d_sd
    k_mu = K_d_mu

    V = np.zeros(N)
    pb = utils._progressbar_start(max_value=N, pbar_on=pbar)
    for i in range(N):
        x_s_CURR = X_s[i]
        K_0 = (K_p * Y_s) * (X_a[i] / Y_a[i])
        try:
            sol = root_scalar(f=__residual,
                              method="brentq", bracket=[0.001, 100.],
                              x0=13., xtol=1E-3, rtol=1E-3,
                              )
        except ValueError:
            V[i] = np.nan
        else:
            V[i] = sol.root * 60.
        utils._progressbar_update(pb, i)
    utils._progressbar_finish(pb)

    return V



def te_pdfseries(TE):
    """
    Returns the trace element time series estimated using BayProX.
    """
    # FIN = DATPATH + "output/002_bayprox_%s_ows%d_sres%d_proxyrecord.npz" \
    #                 % (TE, OWS, SRES)
    # DAT = np.load(FIN, allow_pickle=True)

    # PDist = DAT["PDist"].item()
    PDist = TE['PDist']
    age = PDist.calbp
    X_s = PDist.proxyspan
    pdfmat = PDist.pdfmat
    molwt = {
                "Ni": 58.693,
                "Co": 58.933,
            }
    X_s /= (1E3 * TE['mol_wt'])        # convert to mol / L

    # use scipy.interpolate.interp1d to estimate PDF series
    from scipy.interpolate import interp1d
    pdfseries = []
    for i in range(pdfmat.shape[1]):
        pdfseries.append(interp1d(X_s, pdfmat[:, i],
                                  bounds_error=False, fill_value=0.
                                  )
                        )

    return age, pdfseries


def dr_pdfseries(X_s_pdf, X_a, Y_a, K_p, K_d_mu, K_d_sd, K_e,
                 n_F=0.001, pbar=False):
    """
    Returns the pdf series of driprate. 
    Evaluates rho(v|t) = rho(h(v)|t)*h'(v) = rho(Xs|t)*h'(v)
    X_s_pdf is the PDF of the trace element vs time t
    X_a is the avarage trace element concentration in drip water
    Y_a is the average calcite concentration in drip water
    K_p is the partitioning coefficient
    K_d_mu is the average of LN(kinetic rate)
    K_d_sd is the standard deviation of LN(kinetic rate)
    """
    def h(v, k, k_rsw):
        """
        Returns trace element as a function of driprate.
        """
        f_k = __integrand_1(k, v)
        E_1 = (f_k * k_rsw).sum(axis=1)
        return K_0 * (1. - n_S * E_1)


    def __integrand_1(k, v):
        """
        The integrand inside expectation E1 above.
        """
        term0 = - np.exp(k) / v
        # term0 = - np.exp(k) / (v/60)
        term1 = np.exp(term0)
        term2 = np.exp(- (k - k_mu) ** 2 / (2. * k_sd ** 2))
        const = 1. / (k_sd * np.sqrt(2. * np.pi))
        return const * term1 * term2


    def hprime(v, k, k_rsw):
        """
        Returns trace element as a function of driprate.
        """
        f_k = __integrand_2(k, v)
        E_2 = (f_k * k_rsw).sum(axis=1)
        return - K_0 * n_S * E_2


    def __integrand_2(k, v):
        """
        The integrand inside expectation E2 above.
        """
        term0 = - np.exp(k) / v
        # term0 = - np.exp(k) / (v/60.0)
        term1 = np.exp(term0)
        term2 = np.exp(- (k - k_mu) ** 2 / (2. * k_sd ** 2))
        term3 = np.exp(k) / v ** 2
        # term3 = np.exp(k) / (v/60) ** 2
        const = 1. / (k_sd * np.sqrt(2. * np.pi))
        return const * term1 * term2 * term3

    # X_a *= 1.5
    # K_p *= 3.0
    # K_e /= 3.0
    N = len(X_s_pdf)
    if np.isscalar(X_a):
        X_a = X_a.repeat(N)
    if np.isscalar(Y_a):
        Y_a = Y_a.repeat(N)
    assert X_a.shape == Y_a.shape, "Aqueous concentrations are different sizes"

    # for the solid Ca in calcite, it is known that Ca is 40% by weight in
    # calcite, so ---
    #       40 g of Ca in 100 g of calcite
    #       400 g of Ca in 1000 g of calcite equiv to 1 L of calcite
    #       400 / 40.078 mol of Ca in 1 kg of calcite
    Y_s = 400. / 40.078
    n_S = 1. - n_F
    # k_sd = np.log(1. + (K_d_sd / K_d_mu) ** 2)
    # k_mu = np.log(K_d_mu / np.sqrt(1. + (K_d_sd / K_d_mu) ** 2))
    k_sd = K_d_sd
    k_mu = K_d_mu

    # create the arrays for the integrations with respect to k
    k_span = np.linspace(k_mu - KLIM * k_sd, k_mu + KLIM * k_sd, KRES)
    k_rsw = 0.5 * np.r_[
                        k_span[1] - k_span[0],
                        k_span[2:] - k_span[:-2],
                        k_span[-1] - k_span[-2]
                       ]
    V_span = np.linspace(VMIN, VMAX, VRES) / 60.
    nk, nv = k_span.shape[0], V_span.shape[0]
    k_arr = k_span.repeat(nv).reshape(nk, nv).T
    v_arr = V_span.repeat(nk).reshape(nv, nk)

    # estimate K_0 for each time instant and evaluate h(v) and h_prime(v)
    # beforehand using the array formulations
    K_0 = np.zeros(N)
    for i in range(N):
        K_0 = (K_p * Y_s) * (X_a[i] / Y_a[i]) * K_e
    h_v_arr = h(v_arr, k_arr, k_rsw)
    h_prime_v_arr = hprime(v_arr, k_arr, k_rsw)

    # final loop to get the pdf matrix entries
    V_rsw = 0.5 * np.r_[
                        V_span[1] - V_span[0],
                        V_span[2:] - V_span[:-2],
                        V_span[-1] - V_span[-2]
                       ]
    V_pdf = np.zeros((VRES, N))
    pb = utils._progressbar_start(max_value=N, pbar_on=pbar)
    # N is the number of t evaluated
    for i in range(N):
        pdf = - X_s_pdf[i](h_v_arr) * h_prime_v_arr
        # pdf /= (pdf * V_rsw).sum()
        pdf /= np.nansum(pdf * V_rsw)
        V_pdf[:, i] = pdf
        utils._progressbar_update(pb, i)
    utils._progressbar_finish(pb)

    return V_pdf, V_span * 60.
