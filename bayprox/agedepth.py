"""Contains classes and methods related to calibration & DWF.

    This module contains the following classes for age-depth estimation:
        1.      Calibration
        2.      DWF
"""

import scipy as sp
from scipy.interpolate import interp1d


from . import data
from .other.motabar import function as motabar
from progressbar import ProgressBar


class Calibration(object):
    """Loads specified IntCal & post-bomb calibration curves.

        Methods:
        get_calcurve    loads specified calibration curve (default = IntCal13)
        get_postbomb    loads post-bomb calibration curve (default = NH1)
        update          updates information from data.DatingTable instance

        Example usage:
        >> import agedepth
        >> X = agedepth.Caibration()
        >> intcal09 = X.get_intcal()
        >> postbomb_NH2 = X.get_postbomb(zone='NH2')
        >> X = X.update(D)  # where D is an data.DatingTable
        """
    def __init__(self):
        pass

    @classmethod
    def get_calcurve(self, ver='13', cal_type='intcal'):
        """Loads specified version of the IntCal calibration curve.

        Usage: If X is an agedepth.Calibration instance, then
        >> intcal09 = X.get_intcal(ver='09')
        NOTE:   Default version is IntCal09. And right now, this is the
        only version available. The data is as given by Reimer et al. in
        P.J.REIMER ET AL.,
        INTCAL09 AND MARINE09 RADIOCARBON AGE CALIBRATION CURVES,
        0-50,000 YEARS CAL BP,
        RADIOCARBON, Vol 51, Nr 4, 2009, pp 1111-1150
        """
        dirpath = data.__file__[:-9]
        # if 'intcal' in dir(self):
        #     ver = self.intcal
        fpath = '/other/calibration/'
        fname = dirpath + fpath + cal_type + ver +'.14c'
        calcurve = sp.flipud(sp.genfromtxt(fname, delimiter=','))
        return calcurve

    @classmethod
    def get_postbomb(self, zone=None):
        """Loads specified version of the post-bomb calibration curve.

        Usage: If X is an agedepth.Calibration instance, then to
               get the calibration info for Southern Hemisphere (SH)-
        >> postbomb_SH = X.get_postbomb(zone='SH')
        NOTE:   This calibration data is available only till 2000.
                The options for the various zones are:
                Not needed for current data : None (DEFAULT)
                Northern Hemisphere Zone 1 : 'NH1'
                Northern Hemisphere Zone 2 : 'NH2'
                Northern Hemisphere Zone 3 : 'NH3'
                Southern Hemisphere        : 'SH'
        The data for these curves are as provided by Hua and Barbetti in
        Q. HUA AND M. BARBETTI,
        REVIEW OF TROPOSPHERIC BOMB 14C DATA FOR CARBON CYCLE MODELING
        AND AGE CALIBRATION PURPOSES,
        RADIOCARBON, Vol 46, Nr 3, 2004, pp 1273-1298
        """
        dirpath = data.__file__[:-9]
        if 'postbomb' in dir(self):
            zone = self.postbomb
        if zone:
            fpath = '/other/calibration/postbomb/'
            fname = dirpath+fpath+zone+'.csv'
            postbomb = sp.genfromtxt(fname, delimiter=',')
        else:
            postbomb = None
        return postbomb

    @classmethod
    def update(self, datingtable):
        """Update calibration info according to data.DatingTable instance.

        Usage: If X is an agedepth.Calibration instance and D is an
               data.DatingTable instance, then X can update/load the
               calibration info stored in D as:
        >> X = X.update(D)
        Thereafter, X.intcal is the same as D.intcal_version; and
        X.postbomb is the same as D.postbomb_zone.
        NOTE: After X.update(), calls to get_intcal() and get_postbomb()
        will automatically use X.intcal and X.postbomb as the respective
        options. To load a calibration curve other than the ones
        specified in D (or X, after update), one has to use either a
        new agedepth.Calibration instance or change the data in the
        attributes 'intcal' and 'postbomb' respectively in X.
        """
        if datingtable.info.cal_reqd == 'Yes':
            datingtable.calibration(pb_zone=datingtable.info.postbomb_zone)
            self.cal_version = datingtable.info.cal_version
            self.cal_type = datingtable.info.cal_type
            self.postbomb = datingtable.info.postbomb_zone


class DWF(data.DatingTable, Calibration):
    def __init__(self, datingtable):
        data.DatingTable.__init__(self,
                                   datingtable.depth,
                                   datingtable.age,
                                   datingtable.ageerror,
                                   datingtable.info.datingmethod,
                                   datingtable.info)
        Calibration.update(datingtable)

    def set_rmagemodel(self, regression_params):
        """Sets the RM Age Model using MoTaBaR regression."""
        if self.verbose == 1: print("RM age model...")
        [max_degree, pre_remvar, verbose] = regression_params
        regression = motabar.Function(p=max_degree, verbosity=verbose)
        regression.set_priors(sr=pre_remvar)
        regression.set_data(self.depth, self.age, sy=self.ageerror)
        return regression

    def get_rmagemodel(self, regressionfunc, eval_depth=None, method="LUpinv"):
        """Gets output of RM Age Model at specified depths."""
        if eval_depth is not None:
            eval_depth = self.proxydepth
        estimate, precision = regressionfunc(eval_depth,
                                             method=method,
                                             verbosity=self.verbose)
        rm_mean, rm_dev = estimate[:,0], sp.sqrt(precision[:,0,0]**(-1))
        self.rmagemod = [eval_depth, rm_mean, rm_dev]
        return rm_mean, rm_dev

    def get_riemannsum_intervalwidth(self):
        """Estimates size of interval width for Riemann sum used in DWF."""
        depth = self.proxydepth
        width = depth[2:] - depth[:-2]
        width = sp.hstack((depth[1]-depth[0], width, depth[-1]-depth[-2]))
        width = abs(0.5*width)
        return width

    def get_cal_axis(self, rm_age_top, rm_age_bott, calBP_step):
        """Gets calendar age axis for DWF with given RM ages & stepsize."""
        [intcal_curve, postbomb_curve] = self.calib_curves
        if postbomb_curve is not None:
            turn_pt = sp.where(postbomb_curve[:,1] ==
                               min(postbomb_curve[:,1]))[0]
            postbomb_curve = postbomb_curve[0:turn_pt]
            calib_dat = sp.vstack((sp.flipud(postbomb_curve), intcal_curve))
        else:
            calib_dat = intcal_curve
        cal_rufly = interp1d(calib_dat[:,1], calib_dat[:,0])
        [cal_min, cal_max] = [cal_rufly(rm_age_top), cal_rufly(rm_age_bott)]
        cal_axis = sp.arange(cal_min, cal_max, calBP_step)
        return cal_axis

    def get_cal_dat(self, cal_ax):
        """Gets the calibration function for specified calendar age axis."""
        if self.info.cal_reqd == 'Yes':
            [intcal_curve, postbomb_curve] = self.calib_curves
            if postbomb_curve is not None:
                calib_dat = sp.vstack((sp.flipud(postbomb_curve), intcal_curve))
            else:
                calib_dat = intcal_curve
            cal_func = interp1d(calib_dat[:,0], calib_dat[:,1])
            cal_errfunc = interp1d(calib_dat[:,0], calib_dat[:,2])
        elif self.info.cal_reqd == 'No':
            cal_func = lambda x: x
            cal_errfunc = lambda x: sp.zeros(x.shape)
        rm_age, rm_age_err = cal_func(cal_ax), cal_errfunc(cal_ax)
        self.calcurve = [cal_ax, rm_age, rm_age_err]
        return rm_age, rm_age_err

    def get_dwf_pre(self,
                    rm_calage, rm_calage_err,
                    rm_pdepth, rm_pdepth_err,
                    riemannsum_width):
        """Estimates the initial non-monotonic DWF."""
        self.n_timepts, self.n_depthpts = len(rm_calage), len(rm_pdepth)
        dwf_pre = sp.zeros((self.n_timepts, self.n_depthpts))
        if self.verbose == 1:
            print("DWF pre...")
            pbar = ProgressBar(maxval=self.n_timepts).start()
        for t in range(self.n_timepts):
            numer = sp.square(rm_pdepth - rm_calage[t])
            denom = sp.square(rm_calage_err[t]) + sp.square(rm_pdepth_err)
            const = 1. / sp.sqrt(denom)
            dwf_pre[t,:] = const * sp.exp((-.5) * numer/denom)
            denom = sum(riemannsum_width*dwf_pre[t,:])
            dwf_pre[t,:] = dwf_pre[t,:]/denom
            if self.verbose == 1: pbar.update(t)
        if self.verbose == 1: pbar.finish()
        return dwf_pre

    def get_cdwf_pre(self, dwf_pre, riemannsum_width):
        """Estimates the initial non-monotonic CDWF from given DWF_pre."""
        dwf_pre = dwf_pre * riemannsum_width.reshape((1,-1))
        dwf_pre /= dwf_pre.sum(axis=1).reshape((-1,1))       # normalization
        [nT, nZ] = [self.n_timepts, self.n_depthpts]
        cdwf_pre = sp.zeros((nT, nZ))
        if self.verbose == 1:
            print("CDWF pre...")
            pbar = ProgressBar(maxval=nZ).start()
        for depth in range(nZ):
            cdwf_pre[:,depth] = dwf_pre[:,:depth+1].sum(axis=1)
            if self.verbose == 1: pbar.update(depth)
        if self.verbose == 1: pbar.finish()
        cdwf_pre /= cdwf_pre[:,-1].reshape((-1,1))           # normalization
        return cdwf_pre

    def initialize_cdwf_post(self, cdwf_pre):
        """Initialize monotonic CDWF with pragmatic min/max choice."""
        cdwf_mon_ini = sp.zeros(cdwf_pre.shape)
        if self.verbose == 1:
            print("Initialize CDWF post...")
            pbar = ProgressBar(maxval=cdwf_pre.shape[0]).start()
        for t in range(cdwf_pre.shape[0]):
            cdwf_mon_ini[t,:] = ( cdwf_pre[:t+1,:].min(axis=0)
                                + cdwf_pre[t:,:].max(axis=0) )/2.
            if self.verbose == 1: pbar.update(t)
        if self.verbose == 1: pbar.finish()
        return cdwf_mon_ini

    def workspace_mat(self, cdwf_mon_ini):
        """Get workspace matrix with buffer rows/columns for shifting."""
        [nT, nZ] = [self.n_timepts, self.n_depthpts]
        vert_buff = sp.inf*sp.ones((nT, 1))
        hori_buff = sp.inf*sp.ones((1, nZ+2))
        work = sp.hstack((-vert_buff, cdwf_mon_ini, vert_buff))
        work = sp.vstack((hori_buff, work, -hori_buff))
        return work

    def get_cdwf_post(self, cdwf_pre, cdwf_mon_ini, work_mat):
        """Estimates final monotonic CDWF from initial monotonic CDWF."""
        spmin = sp.minimum
        spmax = sp.maximum
        [nT, nZ] = [self.n_timepts, self.n_depthpts]
        nsteps = self.relax_dyn_params
        dtau = 10./nsteps 					# 10/steps!
        cdwf_post = cdwf_mon_ini
        if self.verbose == 1:
            print("CDWF post...")
            pbar = ProgressBar(maxval=nsteps).start()
        for it in range(nsteps):
            dPhi = dtau*(cdwf_pre - cdwf_post)
            cand = cdwf_post + dPhi
            cdwf_post = work_mat[1:nT+1, 1:nZ+1] \
                      = sp.where(dPhi > 0,
                                 spmin(spmin(cand, work_mat[1:nT+1, 2:]),
                                       work_mat[:nT, 1:nZ+1]
                                       ),
                                 spmax(spmax(cand, work_mat[1:nT+1, :nZ]),
                                       work_mat[2:, 1:nZ+1]
                                       )
                                 )
            if self.verbose == 1: pbar.update(it)
        if self.verbose == 1: pbar.finish()
        return cdwf_post

    def smooth_cdwf(self, cdwf_post):
        """Final smoothing to repair small errors in CDWF_post estimate."""
        [nT, nZ] = [self.n_timepts, self.n_depthpts]
        cdwf_smooth = sp.zeros((nT,nZ))
        if self.verbose == 1:
            print("Smooth CDWF...")
            pbar = ProgressBar(maxval=nT).start()
        for i in range(nT):
            cdwf_smooth[i,:] = ( cdwf_post[:i+1,:].min(axis=0)
                               + cdwf_post[i:,:].max(axis=0) )/2.
            if self.verbose == 1: pbar.update(i)
        if self.verbose == 1: pbar.finish()
        return cdwf_smooth

    def get_dwf_post(self, cdwf_post_smooth):
        """Get DWF_post from final smoothened CDWF_post."""
        [nT, nZ] = [self.n_timepts, self.n_depthpts]
        dwf_post = sp.zeros((nT,nZ))
        dwf_post[:,0] = cdwf_post_smooth[:,0]
        dwf_post[:,1:] = cdwf_post_smooth[:,1:] - cdwf_post_smooth[:,:-1]
        return dwf_post

    def set_data(self,
                 proxydepth,
                 max_degree, pre_remvar, verbose,
                 cal_age_min, cal_age_max, calBP_step,
                 nsteps):
        """Sets various data/paremeter values for estimation of DWF."""
        self.proxydepth = proxydepth
        self.regression_params = [max_degree, pre_remvar, verbose]
        self.relax_dyn_params = nsteps
        self.cal_age_lims = [cal_age_min, cal_age_max, calBP_step]
        self.verbose = verbose
        if self.info.cal_reqd == 'Yes':
            cal_crv = Calibration.get_calcurve(ver=self.cal_version,
                                               cal_type=self.cal_type)
            cal_crv = cal_crv[:, 0:3]
            pb_crv = Calibration.get_postbomb(zone=self.postbomb)
            self.calib_curves = [cal_crv, pb_crv]
        else:
            self.calib_curves = None, None

    def __call__(self,
                 proxydepth,
                 max_degree=2, pre_remvar=None, method="LUpinv", verbose=0,
                 cal_age_min=None, cal_age_max=None, calBP_step=10,
                 calBP_axis=None,
                 nsteps=1000):
        # proxydepth=PD.depth; calBP_axis=calage; verbose=0; nsteps=1000; calBP_step=10
        self.verbose = verbose
        DWF.set_data(self,
                     proxydepth,
                     max_degree, pre_remvar, verbose,
                     cal_age_min, cal_age_max, calBP_step,
                     nsteps)
        rm_age_fun = DWF.set_rmagemodel(self, self.regression_params)
        rm_age_mod = DWF.get_rmagemodel(self,
                                            rm_age_fun,
                                            self.proxydepth,
                                            method=method)
        rm_pdepth, rm_pdepth_err = rm_age_mod
        ## TODO: put in a check here to see if rm_pdeptherr has imaginary values
        ##       this happens if the remainder variance is None and the age
        ##       measurement errors are very low. Then exit the run and prompt
        ##       the user to put in a suitable pre_remvar e.g. 1, 5, 10.
        if calBP_axis is None:
            if self.info.cal_reqd == 'Yes':
                if (cal_age_min is not None) and (cal_age_max is not None):
                    calBP_axis = sp.arange(cal_age_min, cal_age_max, calBP_step)
                else:
                    # TODO: Fix the error with interpolation out of range
                    #print rm_pdepth[0], rm_pdepth[-1],
                    calBP_axis = DWF.get_cal_axis(self,
                                                  rm_pdepth[0],
                                                  rm_pdepth[-1],
                                                  calBP_step)
            elif self.info.cal_reqd == 'No':
                if self.cal_age_lims[0] is None and self.cal_age_lims[1] is None:
                    calBP_axis = sp.arange(min(self.age - self.ageerror),
                                        max(self.age + self.ageerror),
                                        calBP_step)
                else:
                    calBP_axis = sp.arange(self.cal_age_lims[0], #calagemin
                                        self.cal_age_lims[1], #calagemax
                                        self.cal_age_lims[2], #calbpstep
                                        )
        rm_calage, rm_calage_err = DWF.get_cal_dat(self, calBP_axis)
        self.calbp = calBP_axis
        riemannsum_width = DWF.get_riemannsum_intervalwidth(self)
        dwf_pre = DWF.get_dwf_pre(self,
                                  rm_calage, rm_calage_err,
                                  rm_pdepth, rm_pdepth_err,
                                  riemannsum_width)
        cdwf_pre = DWF.get_cdwf_pre(self, dwf_pre, riemannsum_width)
        cdwf_mon_ini = DWF.initialize_cdwf_post(self, cdwf_pre)
        work = DWF.workspace_mat(self, cdwf_mon_ini)
        cdwf_post = DWF.get_cdwf_post(self, cdwf_pre, cdwf_mon_ini, work)
        cdwf_post = DWF.smooth_cdwf(self, cdwf_post)
        dwf_post = DWF.get_dwf_post(self, cdwf_post)
        self.final = dwf_post
        self.output = {
                        "CDWF_final": cdwf_post,
                        "DWF_final": self.final,
                        "calbp": self.calbp,
                        "RM_agemod": rm_age_mod,
                        "Regression_params": self.regression_params,
                        "Relax_Dyn_params": self.relax_dyn_params,
                }
#        self.output = [self.final, self.calbp, rm_age_mod,
#                       self.regression_params, self.relax_dyn_params]
        return self.output
