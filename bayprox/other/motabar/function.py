from numpy import array, ceil, cumsum, diag, diagonal, double, dot, \
    eye, fromfunction, inf, \
    ndarray, pi, product, repeat, sort, sqrt, sum, vstack, where, zeros
from numpy.linalg import inv#, pinv
from scipy.linalg import lstsq, pinv
from scipy.linalg.matfuncs import sqrtm
from scipy.special import factorial
import scipy.sparse
from scipy.sparse import csc_matrix, spdiags, lil_matrix
from scipy.sparse.linalg import factorized
from scipy.optimize import brent #fmin
from .utils import *
import progressbar

class Function(object):
    """
    A MoTaBar regression/interpolation function for an unknown real-valued 
    function  f  of one or more real-valued arguments.
    """

    # names:
    # errx, erry: gamma, epsilon in article
    # tP, tmu: tilde P, tilde mu in article
    
    
    # PROPERTIES:
    
    @property
    def verbosity(self): 
        """(int >= 0)  
        Verbosity level.
        """
        return self._verbosity
    @verbosity.setter
    def verbosity(self, verbosity): self._verbosity = verbosity

    @property
    def d(self): 
        """(read-only, int > 0)  
        Dimension (= number of arguments of  f).
        Set during instantiation.
        """
        return self._d

    @property
    def p(self): 
        """(read-only, int > 0)  
        Degree (derivatives of smaller order will be estimated).
        Set during instantiation.
        """
        return self._p

    @property
    def m(self): 
        """(read-only, int > 0)  
        Number of individual derivatives that will be estimated.
        """
        return self._m

    @property
    def nr(self): 
        """(read-only, int > 0)  
        Number of remainders in Taylor polynomial.
        """
        return self._nr

    @property
    def alphas_phi(self): 
        """(read-only, array(m,d))  
        Multiindices of derivatives to estimate.
        """ 
        return self._alphas_phi
    
    @property
    def alphas_psi(self): 
        """(read-only, array(nr,d))  
        Multiindices of remainders.
        """
        return self._alphas_psi

#    @property
#    def abs_alphas(self): 
#        """(read-only, array(m))  
#        Absolute values of multiindices of estimate derivatives.
#        """
#        return self._abs_alphas_phi
#
#    @property
#    def abs_alphas_r(self): 
#        """(read-only, array(nr))  
#        Absolute values of multiindices of remainders.
#        """
#        return self._abs_alphas_psi

    @property
    def standard_corr_phi(self): 
        """(read-only, array(m,m))  
        Standard prior correlation matrix for estimated derivatives.
        """
        return self._standard_corr_phi

    @property
    def mu_phi(self): 
        """(read-only, array(m))
        Prior mean for estimated derivatives. Set using set_priors().
        """
        return self._mu_phi

    @property
    def mu_psi_i(self): 
        """(read-only, array(nr))  
        Mean of p-th order derivatives. Set using set_priors().
        """
        return self._mu_psi_i

    @property
    def Sigma_phi(self): 
        """(read-only, array(m,m))
        Prior covariance matrix for estimated derivatives.
        Set using set_priors().
        """
        return self._Sigma_phi

    @property
    def Sigma_psi_i(self): 
        """(read-only, csc_matrix(nr,nr))  
        Covariance matrix for p-th order derivatives.
        Set using set_priors().
        """
        return self._Sigma_psi_i

    @property
    def Sigma_gamma_i(self): 
        """(read-only, csc_matrix(d,d))  
        Current default covariance matrix of each observation's argument error.
        Set using set_default_errors().
        """
        return self._Sigma_gamma_i

    @property
    def sigma_eps_i(self): 
        """(read-only, float >= 0)  
        Current default value error standard deviation for future added data.
        Set using set_default_errors().
        """
        return self._sigma_eps_i

    @property
    def corr_gamma_ij(self): 
        """(read-only, csc_matrix(d,d))  
        Current default correlation matrix of argument errors of different 
        observations. Set using set_default_errors().
        """
        return self._corr_gamma_ij

    @property
    def corr_eps_ij(self): 
        """(read-only, float)  
        Current default correlation of value errors of different observations.
        Set using set_default_errors().
        """
        return self._corr_eps_ij

    @property
    def preprocessed(self): 
        """(read-only, bool)  
        Whether data has already been preprocessed for faster calls. 
        If False, preprocessing is done automatically with the next call.   
        """
        return self._preprocessed

    
    # CONSTRUCTOR:
    
    def __init__(self, d=1, p=2, verbosity=0):
        """
        Initialize a MoTaBar regression/interpolation function.
        
        @type d: int > 0
        @arg d:  Dimension (= number of arguments of  f). Default: 1.
        
        @type p: int > 0
        @arg p:  Degree (derivatives of smaller order will be estimated).
                 Default: 2.
                 
        @type verbosity: int > 0
        @arg verbosity:  Verbosity level.
                         Default: 0.
        """
        assert isinstance(d,int) and d > 0
        self._d = d       
        assert isinstance(p,int) and p > 0
        self._p = p
        assert isinstance(verbosity,int) and p >= 0
        self._verbosity = verbosity
        
        # prepare list of multi-indices alpha ordered by |alpha|:

        alphas = last_alphas = [zeros(d,int)]
        last_max_k = [0]
        for abs_a in range(1,p+1):
            new_alphas = []
            new_max_k = []
            for a in range(len(last_alphas)):
                alpha = last_alphas[a]
                for k in range(last_max_k[a],d):
                    new_alpha = 1*alpha
                    new_alpha[k] += 1
                    new_alphas.append(new_alpha)
                    new_max_k.append(k)
            alphas += new_alphas
            last_alphas = new_alphas
            last_max_k = new_max_k
                
        nr = self._nr = len(new_alphas)
        m = self._m = len(alphas) - nr
        self._alphas = array(alphas)
        self._alphas_phi = self._alphas[:m,:]
        self._alphas_psi = array(new_alphas)
        self._abs_alphas_phi = self._alphas_phi.sum(axis=1)
        self._abs_alphas_psi = self._alphas_psi.sum(axis=1)
        self._alpha_fac_invs = 1. / product(factorial(self._alphas),axis=1)\
                                                .reshape((1,1,m+nr))

        # initialize rest:

        self._I_phi = scipy.sparse.eye(m,m)
        self._mu_phi = None
        self._mu_psi_i = None
        self._Sigma_phi = None
        self._P_phi = None
        self._P_phi_mu_phi = None
        self._Sigma_psi_i = None
        self._sigma2_psi_i = None
        self._correlated_psi = None
        
        self._Sigma_gamma_i = zeros((d,d))
        self._sigma_eps_i = 0
        self._corr_gamma_ij = zeros((d,d))
        self._corr_eps_ij = 0
        
        self.clear_data(verbosity=verbosity)
        self.set_default_errors(verbosity=verbosity)
        self.set_priors(verbosity=verbosity)

    # METHODS:

    def set_priors(self, mf=0, mr=0, 
                   sf=inf, sr=None, s0=None, omega=None, l=None, 
                   corr1=0, corr2=0,  
                   verbosity=None):
        """
        Set prior distributions for phi and psi

        If covariances are not given, a standard correlation structure is used:
        The correlation between the alpha-th and beta-th derivatives is
        corr1 * corr2**abs(|alpha|-|beta|-1) if |alpha|-|beta| is odd
        and     corr2**abs(|alpha|-|beta|)   if it is even.
        
        s0 and omega or l can be used to generate standard deviations of the 
        form  s0 * omega**alpha  or  s0 * (2*pi/l)**alpha 
        
        @type mf: 1d array or float
        @arg mf:  Prior mean of phi. Default: 0. 

        @type mr: 1d array or float
        @arg mr:  Common prior mean of all psi[i]. Default: 0. 

        @type sf: 1d or 2d array, or inf
        @arg sf:  If 1d, prior standard deviations of phi. 
                  If 2d, prior covariance matrix of phi.
                  If inf, a noninformative prior for phi is used.
                  Default: inf

        @type sr: 1d or 2d array, float, or None
        @arg sr:  If float, common prior standard devs. of all psi[i,k].
                  If 1d, common prior standard devs. of all psi[i].
                  If 2d, common prior covariance matrix of all psi[i].
                  If None, will be estimated from the data (see article). 
                  Default: None 

        @type s0: float > 0
        @arg s0:  Prior standard deviation of f(x). 
        
        @type omega: 1d array or float > 0
        @arg omega:  Optional "typical angular frequencies" in each argument.
                     Default: 2*pi / l
        
        @type l: 1d array or float > 0
        @arg l:  Optional "typical length scale" in each argument.
                      
        @type corr1,corr2: float in [-1,1]
        @arg corr1,corr2:  Correlation structure. 
                           For f(x)=exp(ax), corr1=corr2>0 is suitable.
                           For f(x)=sin(ax), corr1=0 and corr2<0 is suitable.
                           Default: corr1=corr2=0
                           USE WITH CARE! 

        @type verbosity: int > 0
        @arg verbosity:  Optional verbosity level.
                         Default: as set by class property "verbosity".
        """
        # TODO: allow functions of x for parameters!
        # rename parameters to match article's notation:
        mu_phi = mf
        mu_psi_i = mr
        Sigma_phi = sf
        Sigma_psi_i = sr
        
        d = self._d
        m = self._m
        nr = self._nr
        p = self._p
        mu_phi = double_vector(mu_phi,m)
        mu_psi_i = double_vector(mu_psi_i,nr)
        
        if l is not None: 
            assert all(l > 0), "Typical lengths l must be positive."
            omega = 2*pi / l
        if s0 is not None and omega is not None:
            self._verbose(
             "Using s0 and omega to set/overwrite sf and sr\n"
             +"  to values suitable for a harmonic sin(omega*x)*s0*sqrt(2).")
            s0 = double(s0)
            assert isinstance(s0,double) and s0 > 0, \
                "s0 must be a positive real number."
            omega = double_vector(omega,d)
            assert all(omega > 0), \
                "Typical angular frequencies omega must be positive."
            o = omega.reshape((1,d))
            Sigma_phi = s0 * product(o**self._alphas_phi,axis=1)
            Sigma_psi_i = s0 * product(o**self._alphas_psi,axis=1)
            
        if isinstance(Sigma_phi,float) and Sigma_phi == inf:
            self._verbose(
                "Using a noninformative prior for estimated derivatives.")
            Sigma_phi = diag(repeat(inf,p))
            P_phi = zeros((p,p))
        else:
            Sigma_phi = array(double(Sigma_phi))
            assert isinstance(Sigma_phi,ndarray)
            if Sigma_phi.shape == (m,): # not elif!
                assert all(Sigma_phi > 0), \
                    "Standard deviations in sf must be positive."
                self._verbose(
                 "Assuming standard correlations for estimated derivatives.")
                corr_phi = zeros((m,m))
                for a in range(m):
                    for b in range(m):
                        dabs = abs(self._abs_alphas_phi[a] - self._abs_alphas_phi[b])
                        if dabs % 2 == 0:
                            corr_phi[a,b] = corr2**dabs
                        else:
                            corr_phi[a,b] = corr1 * corr2**(dabs-1)
                p_phi = 1. / Sigma_phi
                Sigma_phi = Sigma_phi.reshape((m,1)) \
                                * Sigma_phi.reshape((1,m)) \
                                * corr_phi
                P_phi = p_phi.reshape((m,1)) * p_phi.reshape((1,m)) \
                                * pinv(corr_phi)
            else:
                Sigma_phi = Sigma_phi.reshape((m,m))
                P_phi = inv(Sigma_phi)
            #print Sigma_phi,P_phi
            assert is_covmatrix(Sigma_phi,m), \
             "Covariance matrix sf is not a non-negative definite square matrix."  
        
        if Sigma_psi_i is None:
            self._Sigma_psi_i = None
            self._sigma2_psi_i = None
            self._correlated_psi = False           
        else:
            Sigma_psi_i = array(double(Sigma_psi_i))
            assert isinstance(Sigma_psi_i,ndarray)
            if Sigma_psi_i.size == 1: 
                assert Sigma_psi_i >= 0, \
                    "Standard deviations in sr must be non-negative."
                self._verbose(
                 "Assuming equal std.dev. for all p-th order derivatives.")
                Sigma_psi_i = repeat(Sigma_psi_i,nr)
            if Sigma_psi_i.shape == (nr,): # not elif!
                self._verbose(
                 "Assuming uncorrelated p-th order derivatives.")
                sigma2_psi_i = Sigma_psi_i**2
                Sigma_psi_i = diag(sigma2_psi_i)
                correlated_psi = False
            else:
                sigma2_psi_i = diag(Sigma_psi_i)
                correlated_psi = True
            assert is_covmatrix(Sigma_psi_i,nr), \
             "Covariance matrix sr is not a non-negative definite square matrix." 
            self._Sigma_psi_i = csc_matrix(Sigma_psi_i)
            self._sigma2_psi_i = sigma2_psi_i
            self._correlated_psi = correlated_psi 

        self._mu_phi = mu_phi.reshape((-1,1))
        self._mu_psi_i = mu_psi_i.reshape((-1,1))
        self._Sigma_phi = Sigma_phi
        self._P_phi = P_phi
        self._P_phi_mu_phi = dot(P_phi,mu_phi).reshape((-1,1))
        
        self._verbose("Sigma_phi="+str(Sigma_phi),level=2)
        

    def set_default_errors(self, sx=0, sy=0, cx=0, cy=0, verbosity=None):
        """
        Set default error distributions for all subsequently added data.
        
        Standard deviations, correlations, or covariances not specified 
        are taken as 0.
        
        @type sx: 1d or 2d array or float #or function(x_i)
        @arg sx:  If 1d, standard deviations of each errx[i].
                  If 2d, covariance matrix of each errx[i].
                  Default: 0

        @type sy: float #or function(x_i)
        @arg sy:  Standard deviation of each erry[i].
                  Default: 0

        @type cx: 2d array or float #or function(x_i,x_j)
        @arg cx:  Correlation matrix between each pair errx[i],errx[j] (i!=j).
                  Default: 0 

        @type cy: float #or function(x_i,x_j)
        @arg cy:  Correlation between each pair erry[i],erry[j] (i!=j).
                  Default: 0 

        @type verbosity: int > 0
        @arg verbosity:  Optional verbosity level.
                         Default: as set by class property "verbosity".
        """
        # rename parameters to match article's notation:
        Sigma_gamma_i = sx
        sigma_eps_i = sy
        corr_gamma_ij = cx
        corr_eps_ij = cy
        
        d = self._d
        
        Sigma_gamma_i = array(double(Sigma_gamma_i))
        assert isinstance(Sigma_gamma_i,ndarray)
        if Sigma_gamma_i.size == 1: 
            assert Sigma_gamma_i >= 0, \
                "Standard deviations in sx must be non-negative."
            self._verbose(
             "Assuming equal error std.dev. for all argument components.")
            Sigma_gamma_i = repeat(Sigma_gamma_i,d)
        if Sigma_gamma_i.shape == (d,): # not elif!
            self._verbose(
             "Assuming uncorrelated components of argument error.")
            Sigma_gamma_i = diag(Sigma_gamma_i**2)
        assert is_covmatrix(Sigma_gamma_i,d) 
        
        sigma_eps_i = double(sigma_eps_i)
        assert isinstance(sigma_eps_i,double) and sigma_eps_i >= 0, \
                "Standard deviation sy must be a non-negative real number."
        
        corr_gamma_ij = array(double(corr_gamma_ij))
        assert isinstance(corr_gamma_ij,ndarray) \
            and all(-1 <= corr_gamma_ij <= 1), \
            "Correlations in cx must be between -1 and 1."
        if corr_gamma_ij.size == 1: 
            self._verbose(
             "Assuming equal inter-observation correlation of argument " \
             + "errors of the same component.")
            corr_gamma_ij = repeat(corr_gamma_ij,d)
        if corr_gamma_ij.shape == (d,): # not elif!
            self._verbose(
             "Assuming no inter-observation correlation of argument " \
             + "errors of different components.")
            corr_gamma_ij = diag(corr_gamma_ij)
        assert corr_gamma_ij.shape == (d,d) 
        
        # verify that Sigma_gamma_i and corr_gamma_ij are compatible:
        S = eye(2*d)
        S[:d,:d] = S[d:,d:] = Sigma_gamma_i
        S[:d,d:] = S[d:,:d] = corr_gamma_ij \
                                * diagonal(Sigma_gamma_i).reshape((1,d))
        assert is_covmatrix(S), \
            "sx and cx are incompatible."
        
        corr_eps_ij = double(corr_eps_ij)
        assert isinstance(corr_eps_ij,double) and -1 <= corr_eps_ij <= 1, \
            "Correlation cy must be between -1 and 1."
        
        self._Sigma_gamma_i = csc_matrix(Sigma_gamma_i)      
        self._sigma_eps_i = sigma_eps_i
        self._corr_gamma_ij = csc_matrix(corr_gamma_ij)
        self._corr_eps_ij = corr_eps_ij

        # TODO: init P_gamma_i etc.?

    
    def set_data(self, x, y, **kwargs):
        """
        Set some data. Same as clear_data() followed by add_data(...).
        """
        self.clear_data()
        self.add_data(x, y, **kwargs)
        
            
    def add_data(self, x, y, sx=None, sy=None, verbosity=None):
        """
        Add some data.
        
        @type x: float or 1d or 2d array
        @arg x:  Argument measurement(s).
                 If d=1 and float, x is one measurement.
                 If d>1 and 1d, x[:] is one measurement.
                 If d=1 and 1d, each x[i] is one measurement.
                 If d>1 and 2d, each x[i,:] is one measurement.

        @type y: float or 1d array
        @arg y:  Corresponding value measurement(s).
        
        @type sx: 1d-4d array or float
        @arg sx:  If float, common standard deviation of all errx[i,k].
                  If 1d, common list of standard devs. of all errx[i].
                  If 2d, common covariance matrix of all errx[i].
                  If 3d, Sigma_gamma[i,:,:] = cov. matrix of errx[i].
                  If 4d, Sigma_gamma[i,j,:,:] = cov.m. of errx[i],errx[j].
                  Default: As specified in set_default_errors().
                        
        @type sy: 1d or 2d array or float
        @arg sy:  If float, common standard deviation of all erry[i].
                  If 1d, individual standard devs. for each erry[i],
                  using correlation from set_default_errors().
                  If 2d, covariance matrix of the pairs erry[i],erry[j].
                  Default: As specified in set_default_errors().

        @type verbosity: int > 0
        @arg verbosity:  Optional verbosity level.
                         Default: as set by class property "verbosity".
        """
        # rename parameters to match article's notation:
        Sigma_gamma = sx
        Sigma_eps = sy 
        
        d = self._d
        
        x = array(double(x))
        assert isinstance(x,ndarray) and x.size % d == 0
        n = int(x.size / d)
        y = array(double(y))
        assert isinstance(y,ndarray) and y.size == n, \
            "x and y have incompatible shapes." 

        if Sigma_eps is None: 
            self._verbose(
             "Using default value error distribution.")
            Sigma_eps = self._sigma_eps_i
        Sigma_eps = array(double(Sigma_eps))
        assert isinstance(Sigma_eps,ndarray)
        if Sigma_eps.size == 1: # not elif! 
            assert Sigma_eps >= 0, \
                "Standard deviations in sy must be non-negative."
            self._verbose(
             "Assuming equal std.dev. for all value errors.")
            Sigma_eps = repeat(Sigma_eps,n)
        if Sigma_eps.shape == (n,): # not elif!
            self._verbose(
             "Assuming default correlation of value errors.")
            sigma_eps = Sigma_eps.copy()
            Sigma_eps = sigma_eps.reshape((1,n)) \
                            * sigma_eps.reshape((n,1)) \
                            * (eye(n) + (1-eye(n)) * self._corr_eps_ij)
        else:
            sigma_eps = sqrt(diagonal(Sigma_eps))
        assert is_covmatrix(Sigma_eps,n), \
         "Covariance matrix sy is not a non-negative definite square matrix."
        
        # TODO: treat Sigma_gamma!
        
        # store data and sigmas:
        self._xs.append(x.reshape((-1,d)))
        self._ys.append(y.reshape((-1,1)))
        self._Sigma_epss.append(Sigma_eps)
        self._sigma_epss.append(sigma_eps)
        self._preprocessed = False


    def preprocess(self, verbosity=None):
        """
        Preprocess data for faster calls. (Called automatically by calls)
        
        @type verbosity: int > 0
        @arg verbosity:  Optional verbosity level.
                         Default: as set by class property "verbosity".                        
        """
        if self._preprocessed: return

        # compute y, n:
        self._x = vstack(self._xs)
        y = self._y = vstack(self._ys).reshape((-1,1))
        n = self._n = y.size
        
        # compose Sigma_eps: 
        nd = len(self._Sigma_epss)
        offs = zeros(nd+1)
        off = 0
        if self._corr_eps_ij is not None and self._corr_eps_ij > 0:
            sigma_eps = zeros(n)
            for i in range(nd):
                offs[i] = off
                thiss = self._sigma_epss[i]
                thisn = thiss.size
                next_off = off + thisn
                sigma_eps[off:next_off] = thiss
                off = next_off
            offs[nd] = off
            Sigma_eps = zeros((n,n)) + self._corr_eps_ij
            for i in range(nd):
                Sigma_eps[offs[i]:offs[i+1],offs[i]:offs[i+1]] \
                    = self._Sigma_epss[i]
            assert is_covmatrix(Sigma_eps), \
                "Specified value error correlations are incompatible."
        else:
            Sigma_eps = lil_matrix((n,n))
            for i in range(nd):
                offs[i] = off
                thisS = self._Sigma_epss[i]
                thisn = thisS.shape[0]
                next_off = off + thisn
                Sigma_eps[off:next_off,off:next_off] = thisS
                off = next_off
        self._Sigma_eps = Sigma_eps
        self._sigma2_eps = Sigma_eps.diagonal()
        
        self._preprocessed = True


    def __call__(self, x, return_w=False, verbosity=None, obs=inf, 
                 method="LUpinv"):
        """
        Return posterior mean and precision matrix of f and its 
        derivatives of order <p for one or more arguments of interest.
        
        NOTE: does not consider argument error yet!
        
        @type x: float or 1d or 2d array
        @arg x:  Argument(s) of interest.
                 If d=1 and float, x is one argument of interest.
                 If d>1 and 1d, x[:] is one argument of interest.
                 If d=1 and 1d, each x[i] is one argument of interest.
                 If d>1 and 2d, each x[i,:] is one argument of interest.
                 
        @type return_w: bool
        @arg return_w:  If True, return MoTaBaR coefficient vectors as well.
                        Default: False
        
        @type obs: int>0 or float <1
        @arg obs:  If int>0, this many observations will be used in the
                   estimation for each argument of interest.
                   If float<1, as many observations will be used as are 
                   necessary to ensure that the remaining observations have a
                   total weight <obs.
                   Default: inf (all observations will be used)
                   USE WITH CARE!

        @type method: "LU", "LUpinv", "pinvLU", "pinv", or "lstsq"
        @arg method:  Which numerical method to use for solving the involved
                      linear systems of equations. See code for details.
                      Default: "LUpinv"
        
        @type verbosity: int > 0
        @arg verbosity:  Optional verbosity level.
                         Default: as set by class property "verbosity".
                         
        @rtype:  array, array [, array]
        @return: Posterior mean and precision matrix of phi.
                 If x contains more than one argument of interest, the first 
                 index of each return array is its index.
                 If return_w=True, (each row of) the third array contains a
                 the vector of MoTaBaR coefficients used to derive the 
                 posterior mean from the data values.  
        """
        assert self._mu_phi is not None, "Please call set_priors() first."
        self.preprocess(verbosity=verbosity)

        assert isinstance(return_w,bool)
        
        assert obs == inf or (isinstance(obs,int) and obs > 0) or \
                (isinstance(obs,float) and 0 < obs < 1) 
                
        p = self._p
        d = self._d
        m = self._m
        nr = self._nr
        n = self._n

        
        xis = array(double(x)).reshape((-1,d))
        assert isinstance(xis,ndarray) and xis.size % d == 0
        nx = xis.shape[0]
        pbar = (nx > 0) * (verbosity > 0)

        x_n1d = self._x.reshape((n,1,d))
        xis_nx11d = xis.reshape((nx,1,1,d))
        ex = self._alphas.reshape((1,m+nr,d))
        fac = self._alpha_fac_invs.reshape((1,m+nr))
        
        tP_phis = zeros((nx,m,m))
        tmu_phis = zeros((nx,m))
        if return_w:
            ws = zeros((nx,m,n))
            
        # if nx>10:
        if pbar:
            pb = progressbar.ProgressBar().start()
            pbinc = ceil(nx / 200)
            
        for ind in range(nx):
            xi = xis[ind,:]
            self._verbose("x="+str(xi),verbosity=verbosity,level=2)
            # TODO: maybe use scipy's Vandermonde generator instead:
            Xp = product((x_n1d - xis_nx11d[ind,:,:,:])**ex, axis=2) * fac
            Xr = Xp[:,m:]
            Xr2= Xr**2
            sumXr2 = sum(Xr2, axis=1)
                            
            def get_tmu_phi(Sigma_y,return_w=False,return_loss=False):
                w = None

                if obs < n:
                    wgts = 1.0/Sigma_y.diagonal()
                    sortedwgts = sort(wgts)
                    if obs >= 1: 
                        ndropped = n - obs
                    else: 
                        ndropped = min(where(cumsum(sortedwgts) 
                                             >= sum(wgts)*obs)[0])
                    filter = where(wgts >= sortedwgts[ndropped])[0]
                    this_Sigma_y = Sigma_y[filter,:][:,filter]
                    this_Xr = Xr[filter,:]
                    self._verbose("  using ",filter.size,"observations",
                                  verbosity=verbosity,level=2)
                    X = Xp[filter,:][:,:m]
                    y = self._y[filter]
                else:
                    this_Sigma_y = Sigma_y
                    this_Xr = Xr
                    X = Xp[:,:m]
                    y = self._y
                    
                XT = X.T # TODO: remove need for this by replacing dot by other
                y_minus_mu_r_T = y.reshape((1,-1)) - \
                    sum(this_Xr*self._mu_psi_i.reshape((1,nr)), axis=1).reshape((1,-1))
    
                if method == "LU":
                    # use LU decompositions for Sigma_y and tP_phi:
                    solver1 = factorized(this_Sigma_y)
                    WX = msolve(this_Sigma_y, X, solver=solver1)
                    tP_phi = self._P_phi + dot(XT, WX) # dot since dense*dense!
                    solver2 = factorized(tP_phi)              
                    tmu_phi = msolve(tP_phi, 
                        self._P_phi_mu_phi + dot(y_minus_mu_r_T, WX).T, 
                        solver=solver2)
                    if return_w:
                        w = msolve(tP_phi, WX.T, solver=solver2)
                elif method == "LUpinv":
                    # use LU decomposition for Sigma_y 
                    # but a pseudo-inverse for tP_phi:
                    solver1 = factorized(this_Sigma_y)
                    WX = msolve(this_Sigma_y, X, solver=solver1)
                    tP_phi = self._P_phi + dot(XT, WX) # dot since dense*dense!
                    tP_phi_pinv = pinv(tP_phi)
                    tmu_phi = dot(tP_phi_pinv, 
                        self._P_phi_mu_phi + dot(y_minus_mu_r_T, WX).T)
                    if return_w:
                        w = dot(tP_phi_pinv, WX.T)
                elif method == "pinvLU":
                    # use a pseudo-inverse for Sigma_y
                    # but a LU decomposition for tP_phi: 
                    W = pinv(this_Sigma_y.todense())
                    WX = dot(W, X)
                    tP_phi = self._P_phi + dot(XT, WX) # dot since dense*dense!
                    solver2 = factorized(tP_phi)              
                    tmu_phi = msolve(tP_phi, 
                        self._P_phi_mu_phi + dot(y_minus_mu_r_T, WX).T, 
                        solver=solver2)
                    if return_w:
                        w = msolve(tP_phi, WX.T, solver=solver2)
                elif method == "pinv":
                    # use a pseudo-inverse for Sigma_y and tP_phi: 
                    W = pinv(this_Sigma_y.todense())
                    WX = dot(W, X)
                    tP_phi = self._P_phi + dot(XT, WX) # dot since dense*dense!
                    tP_phi_pinv = pinv(tP_phi)
                    tmu_phi = dot(tP_phi_pinv, 
                        self._P_phi_mu_phi + dot(y_minus_mu_r_T, WX).T)
                    if return_w:
                        w = dot(tP_phi_pinv, WX.T)
                elif method == "lstsq":
                    # use least-squares (only applicable if P_phi = 0):
                    assert all(self._P_phi == 0)
                    solver1 = factorized(this_Sigma_y)
                    WX = msolve(this_Sigma_y, X, solver=solver1)
                    tP_phi = self._P_phi + dot(XT, WX) # dot since dense*dense!
                    tP_phi_pinv = pinv(tP_phi)
                    tmu_phi0 = dot(tP_phi_pinv, dot(y_minus_mu_r_T, WX).T)
                    # FIXME: the following doesn't work if xi=x for some i
                    # since then Sigma_y is singular...    
                    sqrt_W = inv(sqrtm(this_Sigma_y.todense()))
                    tmu_phi,_res,_rank,_sigma = lstsq(dot(sqrt_W,X),dot(sqrt_W,y))
                    if return_w:
                        w = dot(tP_phi_pinv, WX.T)
                        
                if return_loss:
                    if not method in ["pinvLU","pinv"]:
                        W = pinv(this_Sigma_y.todense())
                    dev = self._y - dot(Xp[:,:m],tmu_phi)
                    return (n + m - dot(dot(dev.T,W),dev)[0,0])**2
                else:
                    return tmu_phi,tP_phi,w

            if self._correlated_psi:
                raise MotabarError("Sorry, correlated sr not yet implemented!")
            
            elif self._sigma2_psi_i is None:
                if self._d>1:
                    raise MotabarError(
                        "Sorry, estimation of sr not yet implemented for d>1!")
                
                def loss(sigma2_psi_i):
                    Sigma_y = csc_matrix(self._Sigma_eps 
                         + spdiags([sumXr2 * sigma2_psi_i],[0],n,n,format='lil'))
                    return get_tmu_phi(Sigma_y=Sigma_y,return_loss=True) 
                
                sigma2_psi_i = brent(loss,brack=(1e-4,1e-3)) #fmin(loss,1) # TODO: start value!
                        #print sigma2_psi_i
                Sigma_y = csc_matrix(self._Sigma_eps 
                     + spdiags([sumXr2 * sigma2_psi_i],[0],n,n,format='lil'))
                tmu_phi,tP_phi,w = get_tmu_phi(Sigma_y=Sigma_y,return_w=return_w)
                    
            else:
                sigma2_r = sum(Xr2 * self._sigma2_psi_i.reshape((1,nr)), 
                               axis=1)
                Sigma_y = csc_matrix(self._Sigma_eps 
                                     + spdiags([sigma2_r],[0],n,n,format='lil'))
                tmu_phi,tP_phi,w = get_tmu_phi(Sigma_y=Sigma_y,return_w=return_w)
                
            tP_phis[ind,:,:] = tP_phi
            tmu_phis[ind,:] = tmu_phi.flatten()
            if return_w:
                ws[ind,:,:] = w
                
            # if nx>10:
            if pbar:
                if (ind % pbinc) == 0:
                    pb.update(int(100.0 * ind / nx))
        
        # if nx>10: pb.finish()   
        if pbar:
            pb.finish()
         
        # TODO: use errx!
        
        if nx == 1:
            if return_w:
                return tmu_phis.reshape((m,)), tP_phis.reshape((m,m)), \
                        ws.reshape((m,n))
            else:
                return tmu_phis.reshape((m,)), tP_phis.reshape((m,m))
        else:
            if return_w:
                return tmu_phis, tP_phis, ws
            else:
                return tmu_phis, tP_phis


    def clear_data(self, verbosity=None):
        """
        Remove all data added earlier.

        @type verbosity: int > 0
        @arg verbosity:  Optional verbosity level.
                         Default: as set by class property "verbosity".
        """
        self._xs = []
        self._ys = []
        self._Sigma_gammas = []
        self._Sigma_epss = []
        self._sigma_epss = []        
        self._preprocessed = False
        

    def _verbose(self, text, verbosity=None, level=1):
        """
        Print text if verbosity >= level. 
        """
        if verbosity is None: verbosity = self._verbosity
        if verbosity >= level: print(text)
