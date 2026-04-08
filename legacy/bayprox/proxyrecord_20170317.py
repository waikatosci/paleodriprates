"""Contains classes and methods related to posterior proxy estimation.

    This module contains the following classes for age-depth estimation:
        1.      ProxyDistributions
        2.      ProxyEstimates
"""

import scipy as sp
import numpy as np
import agedepth
import data
from scipy.interpolate import interp1d
from progressbar import ProgressBar

class ProxyDistributions(agedepth.DWF, data.ProxyDepth):
    """Contains methods to estimate the posterior PDFs and CDFs of proxy."""
    def __init__(self, DWF):
        """Initializes a dummy class instance."""
        self.verbose = DWF.verbose

    def get_pdf(self, DWF, PD, res=1000):
        """Estimates posterior PDF of proxy given DWF & proxy-depth data."""
        limits = ProxyDistributions.get_limits(self, DWF, PD)
        self.res = res
        # create the distance array --> POI := Point Of Interest
        proxydensity_grid = np.linspace(limits[0], limits[1], res)
        self.proxyspan = proxydensity_grid
        proxydensity_grid = np.reshape(proxydensity_grid, (res, 1))
        proxy_poigrid = np.tile(proxydensity_grid, (1, DWF.n_depthpts))
        proxy_obsgrid = np.tile(PD.proxy, (res, 1))
        proxy_distgrid = proxy_obsgrid - proxy_poigrid
        depth_grid = np.tile(PD.depth, (res, 1))
        # main time loop
        riemannsum_width = DWF.get_riemannsum_intervalwidth()
        dwf_f = DWF.final
        weighted_dwf = np.tile(riemannsum_width, (dwf_f.shape[0], 1)) * dwf_f
        pdfmat = np.zeros((len(self.proxyspan), DWF.n_timepts))
        if self.verbose == 1:
            print("proxy pdf ...")
            pbar = ProgressBar(maxval=DWF.n_timepts).start()
        for t in range(DWF.n_timepts):
            dwf = interp1d(PD.depth, weighted_dwf[t,:])
            if PD.proxyerror == 0: # proxy density in depth is Delta function
                row, col = np.where(np.diff(np.sign(proxy_distgrid)) != 0)
                d1, d2 = depth_grid[row, col], depth_grid[row, col+1]
                p1, p2 = proxy_obsgrid[row, col], proxy_obsgrid[row, col+1]
                q = proxy_poigrid[row, col]
                d0 = d1 + ((q-p1)/(p2-p1)) * (d2-d1)
                pp = np.zeros(proxy_obsgrid.shape)
                prob_list = dwf(d0)
                pp[row, col] = prob_list
            else:   # proxy density is a Gaussian
                proxy_vargrid = np.tile(PD.proxyerror, (res, 1))
                gt = np.exp(-0.5 * (np.square(proxy_distgrid/proxy_vargrid)))
                co = 1.0/np.sqrt(2.0 * np.pi * proxy_vargrid)
                pp = co * gt    # co:= constant; gt := Gaussian term
            XXX = np.sum(weighted_dwf[t,:], axis=0)
            YYY = np.sum(pp*weighted_dwf[t,:], axis=1)
            # import sys
            # sys.exit()
            # pp = np.sum(pp, axis=1)/np.sum(weighted_dwf[t,:], axis=0)
            pp = YYY / XXX #np.sum(pp, axis=1)/np.sum(weighted_dwf[t,:], axis=0)
            # pp = pp/sum(pp)
            #pdf.append(interp1d(proxydensity_grid[:,0], pp))
            pdfmat[:, t] = pp
            if self.verbose == 1: pbar.update(t)
        if self.verbose == 1: pbar.finish()
        self.pdfmat = pdfmat
        return pdfmat

    def get_cdf(self, pdfmat):
        """Estimates posterior CDF of proxy given posterior PDF."""
        return np.cumsum(pdfmat, axis=1)
        # cdf = []
        # nT = len(pdf)
        # proxydensity_grid = self.proxyspan
        # if self.verbose == 1:
        #     print("proxy cdf ...")
        #     pbar = ProgressBar(maxval=nT).start()
        # for t in range(nT):
        #     prob_dens_func = pdf[t]
        #     prob_dens_est = prob_dens_func(proxydensity_grid)
        #     cdf_est = np.cumsum(prob_dens_est)
        #     cdf.append(interp1d(proxydensity_grid, cdf_est))
        #     if self.verbose == 1: pbar.update(t)
        # if self.verbose == 1: pbar.finish()
        # self.cdf = cdf
        # return cdf

    def get_limits(self, DWF, PD):
        """Estimate limits of proxy axis within which to estimate PDF & CDF"""
        proxymean = ProxyEstimates.get_mean(DWF, PD)
        proxyvar = ProxyEstimates.get_variance(DWF, PD)
        upper_lim = max(proxymean + 2.5*np.sqrt(proxyvar))
        lower_lim = min(proxymean - 2.5*np.sqrt(proxyvar)) # errors at times
        if PD.proxyrange[0] != -np.inf:
            lower_lim = PD.proxyrange[0]
        elif PD.proxyrange[1] != np.inf:
            upper_lim = PD.proxyrange[1]
        limits = np.array([lower_lim, upper_lim])
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
        weighted_dwf = np.tile(riemannsum_width, (dwf.shape[0], 1)) * dwf
        obs_matrix = np.tile(PD.proxy, (dwf.shape[0], 1))
        numer = np.sum(weighted_dwf * obs_matrix, 1)
        denom = np.sum(weighted_dwf, 1)
        exp_val = numer/denom
        self.mean_value = exp_val
        return exp_val

    @classmethod
    def get_variance(self, DWF, PD, mean=None):
        """Gets variance from DWF object & ProxyDepth object (PD)."""
        if mean is None:
            mean = self.get_mean(DWF, PD)
        mean = np.reshape(mean, (len(mean), 1))
        riemannsum_width = DWF.get_riemannsum_intervalwidth()
        dwf = DWF.final
        mean_matrix = np.tile(mean, (1, dwf.shape[1]))
        obs_matrix = np.tile(PD.proxy, (dwf.shape[0], 1))
        var_innate = np.square(obs_matrix - mean_matrix)
        err_matrix = np.tile(PD.proxyerror, (dwf.shape[0], 1))
        var_instru = np.square(err_matrix)
        weighted_dwf = np.tile(riemannsum_width, (dwf.shape[0], 1)) * dwf
        numer = np.sum(weighted_dwf*(var_innate + var_instru), 1)
        denom = np.sum(weighted_dwf, 1)
        post_var = numer/denom
        self.variance = post_var
        return post_var

    @classmethod
    def get_median(self, cdf, limits, res):
        """Estimates the median based on given CDF."""
        proxygrid = np.linspace(limits[0], limits[1], res)
        nT = len(cdf)
        median_val = np.empty((1,nT)).squeeze()
        for t in range(nT):
            cdf_estimate = cdf[t](proxygrid)
            temp = np.where(cdf_estimate <= 0.5)[0][-1]
            median_val[t] = proxygrid[temp]
        self.median_value = median_val
        return median_val

    @classmethod
    def get_quantiles(self, cdf, limits, res, quantiles):
        """Estimates the specified quantiles based on given CDF."""
        self.quantiles = []
        proxygrid = np.linspace(limits[0], limits[1], res)
        nT = len(cdf)
        [qLo, qHi] = [np.empty((1,nT)).squeeze(), np.empty((1,nT)).squeeze()]
        for t in range(nT):
            cdf_estimate = cdf[t](proxygrid)
            qLo[t] = proxygrid[np.where(cdf_estimate <= quantiles[0])[0][-1]]
            qHi[t] = proxygrid[np.where(cdf_estimate <= quantiles[1])[0][-1]]
        self.quantiles.append([qLo, qHi])
        return qLo, qHi
