from numpy import all, array, double, eye, ndarray, repeat, zeros
from numpy.linalg import eigvalsh
from scipy.sparse import lil_matrix
from scipy.sparse.linalg import factorized, minres

class MotabarError (Exception):
    """
    Used for all exceptions raised by the package.
    """
    def __init__(self, value):
        self.value = value
    def __str__(self):
        return repr(self.value)

def double_vector(a, l):
    """
    Assert a is of length l. 
    If a is scalar, transform to vector by repetition.
    """             
    a = array(double(a))
    assert isinstance(a,ndarray)
    if a.size == 1: a = repeat(a,l)
    assert a.shape == (l,)
    return a

def is_covmatrix(S,dim=None):
    """
    Test whether S is symmetric and non-negative definite 
    (and optionally of dimension dim).
    """
    if dim is None: dim = S.shape[0]
    return S.shape == (dim,dim) and all(S == S.T) and min(eigvalsh(S)) > -1e-10  

def msolve(M, R, solver=None):
    """
    Solve M*X = R with sparse M.
    """
    if solver is None: solver = factorized(M)
    n = R.shape[1]
    X = zeros((M.shape[1],n)) 
    for i in range(n):
        X[:,i] = solver(array(R[:,i]).flatten())
    return X

def singsolve(M, R):
    """
    EXPERIMENTAL: Solve M*X = R with sparse near-singular M.
    """
    n = R.shape[1]
    X = zeros((M.shape[1],n)) 
    for i in range(n):
        X[:,i],res = minres(M, array(R[:,i]).flatten())
        if res != 0: print("singsolve did not converge, rc =", res)
    return X

def sparseinv(M):
    """
    Invert sparse M.
    """
    # TODO: remove need to transpose!
    solve_M = factorized(M)
    n = M.shape[0]
    I = eye(n)
    XT = lil_matrix((n,n)) 
    for i in range(n):
        XT[i,:] = solve_M(I[:,i])
    return XT.T
