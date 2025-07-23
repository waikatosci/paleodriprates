# Various functions to calculate drip rate using trace element concentrations
#
import sys, os, openpyxl
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from scipy.interpolate import interp1d
from scipy.stats import skew, pareto, lognorm

from params import DATPATH, FIGPATH
from params import AXLABFS, TIKLABFS
from params import MAJTIKSZ, MINTIKSZ


from params import OUTLIER_WINDOW_SIZE as WS
from params import BAYPROX_SAMPLING_RES as sampling_res
# from params import CALAGE_MIN
# from params import CALAGE_MAX
# from params import CALAGE_ARR as calage
from params import BAYPROX_MAX_DEGREE as max_degree
from params import BAYPROX_SOLVER_METHOD as method

import model

from utils import _progressbar_start, _progressbar_update, _progressbar_finish

import bayprox
dat = bayprox.data
ad = bayprox.agedepth
prx = bayprox.proxyrecord



def detect_outliers(name, x, y, winsize=11, zero_tol=1E-3, pdf_tol=1E1):
    """
    Removes outliers from specified trace element data series.
    """
    absdev_scaled = __like_hampel(y, winsize=winsize, zero_tol=zero_tol)

    # fit a lognormal distribution to obtained sample of absolute deviations
    # from local medians
    absdev_scaled += 1          # add 1 to sample values to allow logarithm
    hc, bc, be = __loghist(absdev_scaled)
    span = np.logspace(np.log10(bc.min()), np.log10(bc.max()), 100)

    fit = pareto.fit(absdev_scaled)
    pdf = pareto.pdf(span, fit[0], fit[1], fit[2])
    pdf_bc = pareto.pdf(bc, fit[0], fit[1], fit[2])


    # find out where estimated histogram densities are greater than the
    # predicted PDF values from the lognormal fit by a factor of PDF_TOL
    i = np.where(hc / pdf_bc >= pdf_tol)[0]

    # if there are bin centers where estimated histogram densities are bigger
    # than predicted PDF values then get the corresponding bin edges as upper
    # and lower thresholds and use them to identify the anomalies
    j = []
    if len(i) > 0:
        for k in range(len(i)):
            thr_lo, thr_hi = be[i[k]], be[i[k] + 1]
            m = (absdev_scaled >= thr_lo) * (absdev_scaled <= thr_hi)
            j.extend(np.where(m)[0])

    return j



def __like_hampel(y, winsize=11, zero_tol=1E-3):
    """
    Returns the absolute deviations from local median in units of MAD
    """
    winsize = int(winsize)

    n = len(y)
    absdev_scaled = np.zeros(n)
    k = int(winsize / 2)
    for i in range(n):
        if i < k:
            y_i = y[:winsize]
        elif (i >= k) * (i < (n - k)):
            y_i = y[i-k:i+k+1]
        elif i >= (n - k):
            y_i = y[winsize:]
        md = np.median(y_i)
        absdev = np.abs(y_i - md)
        mad = np.median(absdev)
        if (mad > zero_tol):
            absdev_scaled[i] = absdev[k] / mad
        else:
            absdev_scaled[i] = 0.
    return absdev_scaled



def __doane(arr):
    """
    Returns the number of bins according to Doane's formula.

    More info:
        https://en.wikipedia.org/wiki/Histogram#Number_of_bins_and_width
    """
    n = float(len(arr))
    g1 = skew(arr)
    sig_g1 = np.sqrt((6. * (n - 2)) / ((n + 1) * (n + 3)))
    nbins = int(np.ceil(1. + np.log2(n) + np.log2(1 + np.abs(g1) / sig_g1)))

    return nbins



def __loghist(arr):
    """
    Returns the histogram counts on a logarithmic binning.
    """
    nbins = __doane(np.log( arr ))
    bins = np.logspace(np.log10(arr.min()),
                       np.log10(arr.max() + 0.01),
                       nbins + 1)
    hc, be = np.histogram(arr, bins=bins, density=True)
    bc = 0.5 * (be[1:] + be[:-1])
    return hc, bc, be



def preremvar_optimize_residual(pre_remvar, *args):
    """
    Returns the residual for a given choice of prior remainder variance.
    """
    DT, PD, calage, max_degree, method, verb = args
    pre_remvar *= DT.ageerror.std()
    DWF = ad.DWF(DT)
    DWF.set_data(PD.depth,
                 max_degree, pre_remvar, verb,
                 cal_age_min=None, cal_age_max=None, calBP_step=None,
                 nsteps=None
                 )
    regressionfunc = DWF.set_rmagemodel([max_degree, pre_remvar, verb])
    rmagemod = DWF.get_rmagemodel(regressionfunc, PD.depth, method)
    ae, rae = DT.ageerror, rmagemod[1]
    u = ae.mean() - rae.mean()
    v = ae.std() - rae.std()
    residual = np.sqrt(u ** 2 + v ** 2)
    return residual



def get_cdfmat(pdfmat, var_span, verbose=False, pbar=False):
    """
    Returns Cumulative Distribution Functions from given PDFs.

    pdfmat.shape = (nv, nt)
    """
    pdfmat = pdfmat.T
    nt = pdfmat.shape[0]
    bj = 0.5 * np.r_[
                     var_span[1] - var_span[0],
                     var_span[2:] - var_span[:-2],
                     var_span[-1] - var_span[-2]
                    ]                               # Riemann sum width
    cdfmat = np.zeros(pdfmat.shape)
    for i in range(nt):
        pdf = pdfmat[i]
        cdfmat[i] = np.cumsum(pdf * bj)

    return cdfmat.T



def get_proxy_percentile(pc, cdfmat, var_span, verbose=False):
    """
    Returns the Inter-Quartile Range for the paleo dataset.

    cdfmat.shape = (nv, nt)
    """
    cdfmat = cdfmat.T
    nt = cdfmat.shape[0]
    q = np.zeros(nt)
    for i in range(nt):
        q[i] = np.interp(pc / 100., cdfmat[i], var_span)
    return q



def driprates(Kd_mn, Kd_sd, K_e, ConcAq=1.0, TE=None,
                calib=False, pbar_on=False):
    """
    Returns the driprates for a given choice of mean Kd
    """
    # Kd_mn=Kd_mn1;Kd_sd=Kd_sd1;K_e1;TE=TE1;calib=False;pbar_on=False;ConcAq=1.0
    # concentration of trace element in calcite
    age, Xs_pdf = model.te_pdfseries(TE)
    Xs_pdf = np.array(Xs_pdf, dtype="object")

    # if called by the parameter estimation function
    if calib:
        age = (1950. - age).astype("int")
        age = age[::-1]
        Xs_pdf = Xs_pdf[::-1]
        i = (age >= Yi) * (age <= Yf)
        age, Xs_pdf = age[i], Xs_pdf[i]

    # aqueous concentrations of Ca and trace metal
    # Xa, Ya = model.aqueous_concentrations(TE, "median")
    # Xa, Ya = model.aqueous_concentrations(TE, "median")
    Xa, Ya = TE['aq_conc'],TE['ca_conc']
    Xa /= (1E6 * TE['mol_wt'])
    Ya /= (1E6 * 40.078)
    # Xa = ConcAq

    # slow and fast fractions
    # etaF = 0.001
    etaF = 0.001

    # get partition coefficient
    Kp = TE['Kp']
    if Kp<0:
        T = TE['Temp_C']
        Kp = model.kp_theory(T, TE['elem'])

    # driprates
    # Using the DWF (depth-spanning weight function), we can evaluate the posterior
    # PDF of a proxy variable (e.g., trace element concentrations) at a given time 
    # using the proxy depth data. The PDF represents the relationship between trace 
    # element concentrations and age. For example, 1000 BP may correspond to a 
    # collection of trace element concentrations (e.g., 4.6 ppm with 20% weight, 
    # 4.81 ppm with 15% weight, 5.05 ppm with 35% weight, and 6.1 ppm with 30% weight).

    # V_pdf is a matrix of the PDF of driprate (drips/minute) vs. time 
    # V_span is the array of the drip rates estimated (0.02 to 100 per minute)
 
    V_pdf, V_span = model.dr_pdfseries(Xs_pdf, Xa, Ya, Kp, Kd_mn, Kd_sd, K_e,
                                       etaF, pbar=pbar_on)

    return V_pdf, age, V_span


