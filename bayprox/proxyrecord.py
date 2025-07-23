"""Contains classes and methods related to posterior proxy estimation.

    This module contains the following classes for age-depth estimation:
        1.      ProxyDistributions
        2.      ProxyEstimates
"""

import scipy as sp
from scipy.interpolate import interp1d
from progressbar import ProgressBar

from . import agedepth
from . import data


class ProxyDistributions(agedepth.DWF, data.ProxyDepth):
    """Contains methods to estimate the posterior PDFs and CDFs of proxy."""
    def __init__(self, DWF):
        """Initializes a dummy class instance."""
        self.verbose = DWF.verbose

    def get_pdf(self, DWF, PD, res=1000, sd_mult=3.0, limits=[0., 0.]):
        """Estimates posterior PDF of proxy given DWF & proxy-depth data."""
        if sum(limits) == 0.:
            limits = ProxyDistributions.get_limits(self, DWF, PD, sd_mult)
        self.res = res
        
        # create the distance array --> POI := Point Of Interest
        span = sp.linspace(limits[0], limits[1], res)
        self.proxyspan = span
        pobs = PD.proxy
        diff = sp.array([span[i] - pobs for i in range(res)])
        sign = sp.diff(sp.sign(diff), axis=1)
        indices = [sp.where(row != 0)[0] for row in sign]

        # main time loop
        rsw = DWF.get_riemannsum_intervalwidth() # Riemann sum width
        dwf_f = DWF.final
        rsw_pxy = 0.5 * sp.r_[
                                span[1] - span[0],
                                span[2:] - span[:-2],
                                span[-1] - span[-2]
                            ]
        
        pdfmat = sp.zeros((res, DWF.n_timepts))
        if self.verbose == 1:
            print("proxy pdf ...")
            pbar = ProgressBar(maxval=DWF.n_timepts).start()
        for t in range(DWF.n_timepts):
            wts = dwf_f[t] * rsw
            weighted_dwf_t = interp1d(PD.depth, wts)
            # proxy density in depth is Delta function
            if PD.proxyerror == 0:
                pp = sp.zeros(res)
                k = 0
                for i in range(res):
                    idx = indices[i]
                    d_ = 0.5 * (PD.depth[idx] + PD.depth[idx + 1])
                    if len(idx) == 0:
                        k += 1
                        pp[i] = 0.
                    else:
                        pp[i] = weighted_dwf_t(d_).sum() / wts.sum()
            else:   # proxy density is a Gaussian
                pass
            pdfmat[:, t] = pp / (pp * rsw_pxy).sum()
            if self.verbose == 1: pbar.update(t)
        if self.verbose == 1: pbar.finish()
        self.pdfmat = pdfmat
        cdfmat = self.get_cdf(pdfmat)
        self.cdfmat = cdfmat
        self.calbp = DWF.calbp
        return pdfmat

    def get_cdf(self, pdfmat):
        """Estimates posterior CDF of proxy given posterior PDF."""
        span = self.proxyspan
        rsw_pxy = 0.5 * sp.r_[
                        span[1] - span[0],
                        span[2:] - span[:-2],
                        span[-1] - span[-2]
                ]
        X = (pdfmat.T * rsw_pxy).T
        cdfmat = sp.cumsum(X, axis=0)
        return cdfmat

    def get_limits(self, DWF, PD, sd_mult):
        """Estimate limits of proxy axis within which to estimate PDF & CDF"""
        proxymean = ProxyEstimates.get_mean(DWF, PD)
        proxyvar = ProxyEstimates.get_variance(DWF, PD)
        upper_lim = max(proxymean + sd_mult*sp.sqrt(proxyvar))
        lower_lim = min(proxymean - sd_mult*sp.sqrt(proxyvar)) # errors at times
        if PD.proxyrange[0] != -sp.inf:
            lower_lim = PD.proxyrange[0]
        elif PD.proxyrange[1] != sp.inf:
            upper_lim = PD.proxyrange[1]
        limits = [lower_lim, upper_lim]
        self.proxyrange = limits
        return limits


class ProxyEstimates(object):
    """Contains methods to get mean, median, variance & inter-quantiles."""
    def __init__(self):
        """Initializes a dummy class instance."""
    pass

    @classmethod
    def get_mean(self, DWF, PD):
        """Gets posterior mean from DWF object & ProxyDepth object (PD)."""
        riemannsum_width = DWF.get_riemannsum_intervalwidth()
        dwf = DWF.final
        weighted_dwf = sp.tile(riemannsum_width, (dwf.shape[0], 1)) * dwf
        obs_matrix = sp.tile(PD.proxy, (dwf.shape[0], 1))
        numer = sp.sum(weighted_dwf * obs_matrix, 1)
        denom = sp.sum(weighted_dwf, 1)
        exp_val = numer/denom
        self.mean_value = exp_val
        return exp_val

    @classmethod
    def get_variance(self, DWF, PD, mean=None):
        """Gets variance from DWF object & ProxyDepth object (PD)."""
        if mean is None:
            mean = self.get_mean(DWF, PD)
        mean = sp.reshape(mean, (len(mean), 1))
        riemannsum_width = DWF.get_riemannsum_intervalwidth()
        dwf = DWF.final
        mean_matrix = sp.tile(mean, (1, dwf.shape[1]))
        obs_matrix = sp.tile(PD.proxy, (dwf.shape[0], 1))
        var_innate = sp.square(obs_matrix - mean_matrix)
        err_matrix = sp.tile(PD.proxyerror, (dwf.shape[0], 1))
        var_instru = sp.square(err_matrix)
        weighted_dwf = sp.tile(riemannsum_width, (dwf.shape[0], 1)) * dwf
        numer = sp.sum(weighted_dwf*(var_innate + var_instru), 1)
        denom = sp.sum(weighted_dwf, 1)
        post_var = numer/denom
        self.variance = post_var
        return post_var

    @classmethod
    def get_median(self, cdf, limits, res):
        """Estimates the median based on given CDF."""
        proxygrid = sp.linspace(limits[0], limits[1], res)
        nT = len(cdf)
        median_val = sp.empty((1,nT)).squeeze()
        for t in range(nT):
            cdf_estimate = cdf[t](proxygrid)
            temp = sp.where(cdf_estimate <= 0.5)[0][-1]
            median_val[t] = proxygrid[temp]
        self.median_value = median_val
        return median_val

    @classmethod
    def get_quantiles(self, cdf, limits, res, quantiles):
        """Estimates the specified quantiles based on given CDF."""
        self.quantiles = []
        proxygrid = sp.linspace(limits[0], limits[1], res)
        nT = len(cdf)
        [qLo, qHi] = [sp.empty((1,nT)).squeeze(), sp.empty((1,nT)).squeeze()]
        for t in range(nT):
            cdf_estimate = cdf[t](proxygrid)
            qLo[t] = proxygrid[sp.where(cdf_estimate <= quantiles[0])[0][-1]]
            qHi[t] = proxygrid[sp.where(cdf_estimate <= quantiles[1])[0][-1]]
        self.quantiles.append([qLo, qHi])
        return qLo, qHi
