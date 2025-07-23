import function
import numpy
from numpy import *
from pylab import *

def small_example():
    """run a very small test example."""

    x = linspace(0,20,1000)
    xs = linspace(0,20,100)
    sy = 1.0
    ys = sin(xs) + numpy.random.normal(size=xs.size)*sy
    beta = 2
    
    def plotconf(x,y,P,color=None):
        """plot confidence bands."""
        conf = 2.0*(0.5/P)**0.5
        upper = y+conf
        lower = y-conf
        if color is None: color = gca().lines[-1]._color
        x = x.tolist()
        xr = x[:]
        xr.reverse()
        upper = upper.tolist()
        lower = lower.tolist()
        plot(x,upper,":",color=color,zorder=-5)
        plot(x,lower,":",color=color,zorder=-5)
        lower.reverse()
        fill(x+xr,upper+lower,alpha=0.1,ec=color,fc=color,zorder=-10)  
    
    p=8
    f = function.Function(d=1,p=p,verbosity=1)
    f.add_data(xs,ys,sx=0,sy=sy)
    
    figure()
    plot(x,sin(x),"k--",label="true function")
    plot(xs,ys,"k+",label="sampled points")
    for c in [0.0,0.9,0.999,0.99999]:
        corr = zeros((p,p))
        for a in range(p):
            for b in range(p):
                if a == b: corr[a,b] = 1.
                elif (a-b) % 4 == 0: corr[a,b] = c
                elif (a-b) % 2 == 0: corr[a,b] = -c
        
        f.set_priors(sf=beta*corr,sr=beta)
        phi,Pphi = f(x)
        plot(x,phi[:,0],label="c="+str(c))
        plotconf(x,phi[:,0],Pphi[:,0,0])    
    
    title("noisy sine, p="+str(p)+", priors with varying correlation"); legend(loc=2); show()
    
def doit():
    """run a large suite of examples."""
    
    # simple harmonics:
    
    x = linspace(0,20,1000)
    xs = linspace(0,20,100)
    sy = 1.0 #2.0
    ys = sin(xs) + numpy.random.normal(size=xs.size)*sy
    beta = 2
    
    def plotconf(x,y,P,color=None):
        conf = 2.0*(0.5/P)**0.5
        upper = y+conf
        lower = y-conf
        if color is None: color = gca().lines[-1]._color
        x = x.tolist()
        xr = x[:]
        xr.reverse()
        upper = upper.tolist()
        lower = lower.tolist()
        plot(x,upper,":",color=color,zorder=-5)
        plot(x,lower,":",color=color,zorder=-5)
        lower.reverse()
        fill(x+xr,upper+lower,alpha=0.1,ec=color,fc=color,zorder=-10)
    
    
    for p in [8,40]:
        f = function.Function(d=1,p=p,verbosity=1)
        f.add_data(xs,ys,sx=0,sy=sy)
        
        figure()
        plot(x,sin(x),"k--",label="true function")
        plot(xs,ys,"k+",label="sampled points")
        for c in [0.0,0.9,0.999,0.99999]:
            corr = zeros((p,p))
            for a in range(p):
                for b in range(p):
                    if a == b: corr[a,b] = 1.
                    elif (a-b) % 4 == 0: corr[a,b] = c
                    elif (a-b) % 2 == 0: corr[a,b] = -c
            
            f.set_priors(sf=beta*corr,sr=beta)
            phi,Pphi = f(x)
            plot(x,phi[:,0],label="c="+str(c))
            plotconf(x,phi[:,0],Pphi[:,0,0])
        
        title("noisy sine, p="+str(p)+", priors with varying correlation"); legend(loc=2); show()
    
    
    # two harmonics:
    
    x = linspace(0,40,1000)
    xs = linspace(0,40,100)
    sy = 2.0
    fac = 3
    y = zeros(1000)
    y[:500] = sin(x[:500])
    y[500:] = fac*sin(x[500:])
    ys = zeros(100)
    ys[:50] = sin(xs[:50]) + numpy.random.normal(size=50)*sy
    ys[50:] = fac*sin(xs[50:]) + numpy.random.normal(size=50)*sy
    beta = 2*fac**2
    
    for p in [8,40]:
        f = function.Function(d=1,p=p,verbosity=1)
        f.add_data(xs,ys,sx=0,sy=sy)
        
        figure()
        plot(x,y,"k--",label="true function")
        plot(xs,ys,"k+",label="sampled points")
        for c in [0.0,0.9,0.999,0.99999]:
            corr = zeros((p,p))
            for a in range(p):
                for b in range(p):
                    if a == b: corr[a,b] = 1.
                    elif (a-b) % 4 == 0: corr[a,b] = c
                    elif (a-b) % 2 == 0: corr[a,b] = -c
            
            f.set_priors(sf=beta*corr,sr=beta)
            phi,Pphi = f(x)
            plot(x,phi[:,0],label="c="+str(c))
            plotconf(x,phi[:,0],Pphi[:,0,0])
        
        title("noisy sine with break-point, p="+str(p)+", priors with varying correlation, 95% conf. bounds"); legend(loc=2); show()
    
    
    # simple harmonics with wrong frequency:
    
    x = linspace(0,20,1000)
    xs = linspace(0,20,100)
    sy = 2.0
    ys = sin(2*xs) + numpy.random.normal(size=xs.size)*sy
    
    for p in [8,40]:
        f = function.Function(d=1,p=p,verbosity=1)
        f.add_data(xs,ys,sx=0,sy=sy)
        
        figure()
        plot(x,2*sin(x),"k--",label="true function")
        plot(xs,ys,"k+",label="sampled points")
        for c in [0.0,0.9,0.999,0.99999]:
            corr = zeros((p,p))
            for a in range(p):
                for b in range(p):
                    if a == b: corr[a,b] = 1.
                    elif (a-b) % 4 == 0: corr[a,b] = c
                    elif (a-b) % 2 == 0: corr[a,b] = -c
            
            beta = 2
            f.set_priors(sf=beta*corr,sr=beta)
            phi,Pphi = f(x)
            plot(x,phi[:,0],label="c="+str(c))
        
        title("noisy sine with wrong frequency, p="+str(p)+", priors with varying correlation"); legend(loc=2); show()
    
    
    # response to monomials:
    
    x = linspace(-1,2,100)
    
    # constant signal:
    figure()
    for p in range(1,5):
        xs = linspace(0,1,2*p)
        ys = xs**0
        plot(xs,ys,".")
        f = function.Function(d=1,p=p,verbosity=1)
        f.set_priors(sf=diag(repeat(inf,p)),sr=1e4)
        f.add_data(xs,ys,sx=0,sy=1e-4)
        plot(x,f(x)[0],label="p="+str(p))
    
    title("constant signal"); legend(loc=2); gca().set_ybound(-1,2); show()
    
    # linear signal:
    figure()
    for p in range(2,6):
        xs = linspace(0,1,2*p)
        ys = xs**1
        plot(xs,ys,".")
        f = function.Function(d=1,p=p,verbosity=1)
        f.set_priors(sf=diag(repeat(inf,p)),sr=1e4)
        f.add_data(xs,ys,sx=0,sy=1e-4)
        plot(x,f(x)[0],label="p="+str(p))
    
    title("linear signal"); legend(loc=2); gca().set_ybound(-1,2); show()
    
    # quadratic signal:
    figure()
    for p in range(3,7):
        xs = arange(2*p) #linspace(0,1,2*p)
        ys = xs**2
        x = arange(10*p) / 5.
        plot(xs,ys,".")
        f = function.Function(d=1,p=p,verbosity=1)
        f.set_priors(sf=diag(repeat(inf,p)),sr=1e4)
        f.add_data(xs,ys,sx=0,sy=1e-4)
        plot(x,f(x)[0],label="p="+str(p))
    
    title("quadratic signal"); legend(loc=2); gca().set_ybound(-1,2); show()
    
    # quadratic signal:
    p = 10
    xs = linspace(-1,1,11)
    ys = xs**2
    f = function.Function(d=1,p=p,verbosity=1)
    f.add_data(xs,ys,sx=0,sy=1e-4)
    es = linspace(0,10,100)
    phis = zeros((100,p))
    for i in range(100):
        e = es[i]
        f.set_priors(sf=diag(repeat(10**e,p)),sr=1e4)
        phis[i,:] = f(0)[0]
    
    plot(es,phis)
    title("quadratic signal, proper priors"); legend(loc=2); show()
    
    
    
    # p=1:
    x = linspace(-3,5,200)
    
    f = function.Function(d=1,p=1,verbosity=1)
    figure()
    
    # negligible value error:
    
    f.set_default_errors(sy=1e-2) # smaller leads to numerical problems!
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2])
    
    f.set_priors(sf=[inf],sr=1e4) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi,label="negligible value error, improper prior for phi, large prior for psi -> IDW interpolation")
    
    f.set_priors(sf=[1e4],sr=1e4) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi,label="negligible value error, almost improper priors -> faster decay than IDW")
    
    f.set_priors(sf=[1],sr=1)
    phi,P=f(x); plot(x,phi,label="negligible value error, medium priors -> faster decay than IDW")
    
    f.set_priors(sf=[1e-4],sr=1e-4)
    phi,P=f(x); plot(x,phi,label="negligible value error, sharp priors -> almost zero")
    
    f.set_priors(sf=[1e4],sr=1e-4)
    phi,P=f(x); plot(x,phi,label="negligible value error, almost improper prior for value, sharp prior for slope -> mean")
    
    # medium value error:
    
    f.set_default_errors(sy=1) # smaller leads to numerical problems!
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2])
    
    f.set_priors(sf=[inf],sr=1e4) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi,label="medium value error, improper prior for phi, large prior for psi -> IDW interpolation")
    
    f.set_priors(sf=[1e4],sr=1e4) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi,label="medium value error, almost improper priors -> faster decay than IDW")
    
    f.set_priors(sf=[1],sr=1)
    phi,P=f(x); plot(x,phi,label="medium value error, medium priors -> IDW smoothing")
    
    f.set_priors(sf=[1e-4],sr=1e-4)
    phi,P=f(x); plot(x,phi,label="medium value error, sharp priors -> almost zero")
    
    f.set_priors(sf=[1e4],sr=1e-4)
    phi,P=f(x); plot(x,phi,label="medium value error, almost improper prior for value, sharp prior for slope -> mean")
    
    # large value error:
    
    f.set_default_errors(sy=1e6)
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2])
    
    f.set_priors(sf=[inf],sr=1e8) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi,label="large value error, improper prior for phi, large prior for psi -> IDW interpolation")
    
    f.set_priors(sf=[1e8],sr=1e8) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi,label="large value error, almost improper priors -> faster decay than IDW")
    
    f.set_priors(sf=[1],sr=1)
    phi,P=f(x); plot(x,phi,label="large value error, medium priors -> almost zero")
    
    f.set_priors(sf=[1e-4],sr=1e-4)
    phi,P=f(x); plot(x,phi,label="large value error, sharp priors -> almost zero")
    
    f.set_priors(sf=[1e4],sr=1e-4)
    phi,P=f(x); plot(x,phi,label="large value error, almost improper prior for value, sharp prior for slope -> mean")
    
    title("p=1"); legend(loc=2); show()
    
    
    # p=2:
    
    f = function.Function(d=1,p=2,verbosity=1)
    
    
    # omega=1:
    
    figure()
    
    # negligible value error:
    
    f.set_default_errors(sy=1e-2) # smaller leads to numerical problems!
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2])
    
    f.set_priors(sf=[1e4,1e4],sr=1e4) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi[:,0],label="negligible value error, almost improper priors -> slight overshooting")
    
    f.set_priors(sf=[1,1],sr=1)
    phi,P=f(x); plot(x,phi[:,0],label="negligible value error, medium priors -> slight overshooting")
    
    f.set_priors(sf=[1e-4,1e-4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],label="negligible value error, sharp priors -> almost zero")
    
    f.set_priors(sf=[1e4,1e4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],label="negligible value error, almost improper prior for phi, sharp prior for psi -> linear regression")
    
    # medium value error:
    
    f.set_default_errors(sy=1) # smaller leads to numerical problems!
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2])
    
    f.set_priors(sf=[1e4,1e4],sr=1e4) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi[:,0],label="medium value error, almost improper priors -> slight overshooting")
    
    f.set_priors(sf=[1,1],sr=1)
    phi,P=f(x); plot(x,phi[:,0],label="medium value error, medium priors -> smoothed")
    
    f.set_priors(sf=[1e-4,1e-4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],label="medium value error, sharp priors -> almost zero")
    
    f.set_priors(sf=[1e4,1e4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],label="medium value error, almost improper prior for phi, sharp prior for psi -> linear regression")
    
    # large value error:
    
    f.set_default_errors(sy=1e6)
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2])
    
    f.set_priors(sf=[1e8,1e8],sr=1e8) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi[:,0],label="large value error, almost improper priors -> slight overshooting")
    
    f.set_priors(sf=[1,1],sr=1)
    phi,P=f(x); plot(x,phi[:,0],label="large value error, medium priors -> almost zero")
    
    f.set_priors(sf=[1e-4,1e-4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],label="large value error, sharp priors -> almost zero")
    
    f.set_priors(sf=[1e4,1e4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],label="large value error, almost improper prior for phi, sharp prior for psi -> linear regression")
    
    title("p=2, omega=1"); legend(loc=2); show()
    
    
    # omega=2:
    
    figure()
    
    # negligible value error:
    
    f.set_default_errors(sy=1e-2) # smaller leads to numerical problems!
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2])
    
    f.set_priors(s0=1e4,omega=2) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi[:,0],label="negligible value error, almost improper priors -> less overshooting")
    
    f.set_priors(s0=1,omega=2)
    phi,P=f(x); plot(x,phi[:,0],label="negligible value error, medium priors -> less overshooting")
    
    f.set_priors(s0=1e-4,omega=2)
    phi,P=f(x); plot(x,phi[:,0],label="negligible value error, sharp priors -> almost zero")
    
    f.set_priors(sf=[1e4,4e4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],label="negligible value error, almost improper prior for phi, sharp prior for psi -> linear regression")
    
    # medium value error:
    
    f.set_default_errors(sy=1) # smaller leads to numerical problems!
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2])
    
    f.set_priors(s0=1e4,omega=2) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi[:,0],label="medium value error, almost improper priors -> less overshooting")
    
    f.set_priors(s0=1,omega=2)
    phi,P=f(x); plot(x,phi[:,0],label="medium value error, medium priors -> smoothed")
    
    f.set_priors(s0=1e-4,omega=2)
    phi,P=f(x); plot(x,phi[:,0],label="medium value error, sharp priors -> almost zero")
    
    f.set_priors(sf=[1e4,4e4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],label="medium value error, almost improper prior for phi, sharp prior for psi -> linear regression")
    
    # large value error:
    
    f.set_default_errors(sy=1e6)
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2])
    
    f.set_priors(s0=1e8,omega=2) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi[:,0],label="large value error, almost improper priors -> less overshooting")
    
    f.set_priors(s0=1,omega=2)
    phi,P=f(x); plot(x,phi[:,0],label="large value error, medium priors -> almost zero")
    
    f.set_priors(s0=1e-4,omega=2)
    phi,P=f(x); plot(x,phi[:,0],label="large value error, sharp priors -> almost zero")
    
    f.set_priors(sf=[1e4,4e4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],label="large value error, almost improper prior for phi, sharp prior for psi -> linear regression")
    
    title("p=2, omega=2 -> more localized"); legend(loc=2); show()
    
    
    # omega=1/4:
    
    figure()
    
    # negligible value error:
    
    f.set_default_errors(sy=1e-2) # smaller leads to numerical problems!
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2])
    
    f.set_priors(s0=1e4,omega=.25) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi[:,0],label="negligible value error, almost improper priors -> more overshooting")
    
    f.set_priors(s0=1,omega=.25)
    phi,P=f(x); plot(x,phi[:,0],label="negligible value error, medium priors -> more overshooting")
    
    f.set_priors(s0=1e-4,omega=.25)
    phi,P=f(x); plot(x,phi[:,0],label="negligible value error, sharp priors -> almost zero")
    
    f.set_priors(sf=[1e4,.0625e4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],label="negligible value error, almost improper prior for phi, sharp prior for psi -> linear regression")
    
    # medium value error:
    
    f.set_default_errors(sy=1) # smaller leads to numerical problems!
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2])
    
    f.set_priors(s0=1e4,omega=.25) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi[:,0],label="medium value error, almost improper priors -> more overshooting")
    
    f.set_priors(s0=1,omega=.25)
    phi,P=f(x); plot(x,phi[:,0],label="medium value error, medium priors -> very smooth")
    
    f.set_priors(s0=1e-4,omega=.25)
    phi,P=f(x); plot(x,phi[:,0],label="medium value error, sharp priors -> almost zero")
    
    f.set_priors(sf=[1e4,.0625e4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],label="medium value error, almost improper prior for phi, sharp prior for psi -> linear regression")
    
    # large value error:
    
    f.set_default_errors(sy=1e6)
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2])
    
    f.set_priors(s0=1e8,omega=.25) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi[:,0],label="large value error, almost improper priors -> more overshooting")
    
    f.set_priors(s0=1,omega=.25)
    phi,P=f(x); plot(x,phi[:,0],label="large value error, medium priors -> almost zero")
    
    f.set_priors(s0=1e-4,omega=.25)
    phi,P=f(x); plot(x,phi[:,0],label="large value error, sharp priors -> almost zero")
    
    f.set_priors(sf=[1e4,.0625e4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],label="large value error, almost improper prior for phi, sharp prior for psi -> linear regression")
    
    title("p=2, omega=1/4 -> more globalized"); legend(loc=2); show()
    
    
    # omega=.01:
    
    figure()
    
    # negligible value error:
    
    f.set_default_errors(sy=1e-2) # smaller leads to numerical problems!
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2])
    
    f.set_priors(s0=1e4,omega=.01) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi[:,0],label="negligible value error, almost improper priors -> nearing linear regression")
    
    f.set_priors(s0=1,omega=.01)
    phi,P=f(x); plot(x,phi[:,0],label="negligible value error, medium priors -> nearing linear regression")
    
    f.set_priors(s0=1e-4,omega=.01)
    phi,P=f(x); plot(x,phi[:,0],label="negligible value error, sharp priors -> almost zero")
    
    # medium value error:
    
    f.set_default_errors(sy=1) # smaller leads to numerical problems!
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2])
    
    f.set_priors(s0=1e4,omega=.01) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi[:,0],label="medium value error, almost improper priors -> nearing linear regression")
    
    f.set_priors(s0=1,omega=.01)
    phi,P=f(x); plot(x,phi[:,0],label="medium value error, medium priors -> very smooth")
    
    f.set_priors(s0=1e-4,omega=.01)
    phi,P=f(x); plot(x,phi[:,0],label="medium value error, sharp priors -> almost zero")
    
    # large value error:
    
    f.set_default_errors(sy=1e6)
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2])
    
    f.set_priors(s0=1e8,omega=.01) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi[:,0],label="large value error, almost improper priors -> almost linear regression")
    
    f.set_priors(s0=1,omega=.01)
    phi,P=f(x); plot(x,phi[:,0],label="large value error, medium priors -> almost zero")
    
    f.set_priors(s0=1e-4,omega=.01)
    phi,P=f(x); plot(x,phi[:,0],label="large value error, sharp priors -> almost zero")
    
    title("p=2, omega=.01 -> strongly globalized"); legend(loc=2); show()
    
    
    # omega=1, value error correlation:
    
    figure()
    
    # uncorrelated medium value error:
    
    f.set_default_errors(sy=1) # smaller leads to numerical problems!
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2])
    
    f.set_priors(sf=[1e4,1e4],sr=1e4) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi[:,0],label="uncorrelated value error, almost improper priors -> slight overshooting")
    
    f.set_priors(sf=[1,1],sr=1)
    phi,P=f(x); plot(x,phi[:,0],label="uncorrelated value error, medium priors -> smoothed")
    
    f.set_priors(sf=[1e-4,1e-4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],label="uncorrelated value error, sharp priors -> almost zero")
    
    f.set_priors(sf=[1e4,1e4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],label="uncorrelated value error, almost improper prior for phi, sharp prior for psi -> linear regression")
    
    # correlated medium value error:
    
    f.set_default_errors() # smaller leads to numerical problems!
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2],sy=[[1,.99,.99],[.99,1,.99],[.99,.99,1]])
    
    f.set_priors(sf=[1e4,1e4],sr=1e4) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi[:,0],label="correlated value error, almost improper priors -> slight overshooting")
    
    f.set_priors(sf=[1,1],sr=1)
    phi,P=f(x); plot(x,phi[:,0],label="correlated value error, medium priors -> smoothed")
    
    f.set_priors(sf=[1e-4,1e-4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],label="correlated value error, sharp priors -> almost zero")
    
    f.set_priors(sf=[1e4,1e4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],label="correlated value error, almost improper prior for phi, sharp prior for psi -> linear regression")
    
    # anticorrelated medium value error:
    
    f.set_default_errors() # smaller leads to numerical problems!
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2],sy=[[1,.99,-.99],[.99,1,-.99],[-.99,-.99,1]])
    
    f.set_priors(sf=[1e4,1e4],sr=1e4) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi[:,0],label="anticorrelated value error, almost improper priors -> slight overshooting")
    
    f.set_priors(sf=[1,1],sr=1)
    phi,P=f(x); plot(x,phi[:,0],label="anticorrelated value error, medium priors -> smoothed")
    
    f.set_priors(sf=[1e-4,1e-4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],label="anticorrelated value error, sharp priors -> almost zero")
    
    f.set_priors(sf=[1e4,1e4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],label="anticorrelated value error, almost improper prior for phi, sharp prior for psi -> linear regression")
    
    title("p=1, omega=1, different value error correlation"); legend(loc=2); show()
    
    
    # p=3:
    
    f = function.Function(d=1,p=3,verbosity=1)
    
    
    # omega=1:
    
    figure()
    
    # negligible value error:
    
    f.set_default_errors(sy=1e-2) # smaller leads to numerical problems!
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2])
    
    f.set_priors(s0=1e4,omega=1) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi[:,0],"k-",lw=0.5,alpha=0.9,label="negligible value error, almost improper priors -> more overshooting")
    
    f.set_priors(s0=1,omega=1)
    phi,P=f(x); plot(x,phi[:,0],"b-",lw=0.5,alpha=0.9,label="negligible value error, medium priors -> more overshooting")
    
    f.set_priors(s0=1e-4,omega=1)
    phi,P=f(x); plot(x,phi[:,0],"r-",lw=0.5,alpha=0.9,label="negligible value error, sharp priors -> strange behaviour")
    
    f.set_priors(sf=[1e4,1e4,1e4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],"y-",lw=0.5,alpha=0.9,label="negligible value error, almost improper prior for phi, sharp prior for psi -> quadratic interpolation")
    
    # medium value error:
    
    f.set_default_errors(sy=1) # smaller leads to numerical problems!
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2])
    
    f.set_priors(s0=1e4,omega=1) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi[:,0],"k--",lw=1,alpha=0.9,label="medium value error, almost improper priors -> more overshooting")
    
    f.set_priors(s0=1,omega=1)
    phi,P=f(x); plot(x,phi[:,0],"b--",lw=1,alpha=0.9,label="medium value error, medium priors -> very smooth")
    
    f.set_priors(s0=1e-4,omega=1)
    phi,P=f(x); plot(x,phi[:,0],"r--",lw=1,alpha=0.9,label="medium value error, sharp priors -> almost zero")
    
    f.set_priors(sf=[1e4,1e4,1e4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],"y--",lw=1,alpha=0.9,label="medium value error, almost improper prior for phi, sharp prior for psi -> quadratic interpolation")
    
    # large value error:
    
    f.set_default_errors(sy=1e6)
    f.clear_data()
    f.add_data([-1,0,2],[0,-1,2])
    
    f.set_priors(s0=1e8,omega=1) # larger psi leads to numerical problems!
    phi,P=f(x); plot(x,phi[:,0],"k:",lw=3,alpha=0.9,label="large value error, almost improper priors -> more overshooting")
    
    f.set_priors(s0=1,omega=1)
    phi,P=f(x); plot(x,phi[:,0],"b:",lw=3,alpha=0.9,label="large value error, medium priors -> strange behaviour")
    
    f.set_priors(s0=1e-4,omega=1)
    phi,P=f(x); plot(x,phi[:,0],"r:",lw=3,alpha=0.9,label="large value error, sharp priors -> strange behaviour")
    
    f.set_priors(sf=[1e4,1e4,1e4],sr=1e-4)
    phi,P=f(x); plot(x,phi[:,0],"y:",lw=3,alpha=0.9,label="large value error, almost improper prior for phi, sharp prior for psi -> strange behaviour")
    
    title("p=3, omega=1"); legend(loc=2); show()
    
    
    
    # nearest neighbour:
    
    
    figure()
    x = linspace(-3,5,200)
    for p in range(1,11):
        f = function.Function(d=1,p=p,verbosity=1)
        f.set_priors(sf=diag([inf]+repeat(1e-10,p-1).tolist()),sr=1e10) 
        f.add_data([-1,0,2],[0,-1,2],sx=0,sy=1e-4)
        phi,P=f(x); plot(x,phi[:,0],label="p="+str(p))
    
    title("vanishing slope, increasing p -> converges to nearest neighbours"); legend(loc=2); show()
    
    figure()
    x = linspace(-3,3,200)
    for p in range(2,11):
        f = function.Function(d=1,p=p,verbosity=1)
        f.set_priors(sf=diag([inf,inf]+repeat(1e-10,p-2).tolist()),sr=1e10) 
        f.add_data([-1,0,1],[0,2,1],sx=0,sy=1e-4)
        phi,P=f(x); plot(x,phi[:,0],label="p="+str(p))
    
    title("vanishing curvature, increasing p -> converges to a form of local linear"); legend(loc=2); show()
    
    # sine curves:
    
    xs = linspace(0,4*pi,30)
    ys0 = sin(xs)
    x = linspace(-pi,5*pi,1000)
    
    # p=2, no value error:
    figure()
    f = function.Function(d=1,p=2)
    f.set_default_errors(sy=1e-6)
    f.clear_data()
    f.add_data(xs,ys0)
    plot(xs,ys0,"b.")
    f.set_priors(s0=1,omega=1)
    phi,P=f(x); plot(x,phi,"y-",label="s0=1,omega=1")
    f.set_priors(s0=1,omega=.5)
    phi,P=f(x); plot(x,phi,"r-",label="s0=1,omega=.5")
    f.set_priors(s0=1,omega=.1)
    phi,P=f(x); plot(x,phi,"k-",label="s0=1,omega=.1")
    title("p=2, no value error"); legend(loc=2); show()
    
    # p=3, no value error:
    figure()
    f = function.Function(d=1,p=3)
    f.set_default_errors(sy=1e-6)
    f.clear_data()
    f.add_data(xs,ys0)
    plot(xs,ys0,"b.")
    f.set_priors(s0=1,omega=1)
    phi,P=f(x); plot(x,phi,"y-",label="s0=1,omega=1")
    f.set_priors(s0=1,omega=.5)
    phi,P=f(x); plot(x,phi,"r-",label="s0=1,omega=.5")
    f.set_priors(s0=1,omega=.1)
    phi,P=f(x); plot(x,phi,"k-",label="s0=1,omega=.1")
    title("p=3, no value error"); legend(loc=2); show()
    
    # p=4, no value error:
    figure()
    f = function.Function(d=1,p=4)
    f.set_default_errors(sy=1e-6)
    f.clear_data()
    f.add_data(xs,ys0)
    plot(xs,ys0,"b.")
    f.set_priors(s0=1,omega=1)
    phi,P=f(x); plot(x,phi,"y-",label="s0=1,omega=1")
    f.set_priors(s0=1,omega=.5)
    phi,P=f(x); plot(x,phi,"r-",label="s0=1,omega=.5")
    f.set_priors(s0=1,omega=.1)
    phi,P=f(x); plot(x,phi,"k-",label="s0=1,omega=.1")
    title("p=4, no value error"); legend(loc=2); show()
    
    # p=5, no value error:
    figure()
    f = function.Function(d=1,p=5)
    f.set_default_errors(sy=1e-6)
    f.clear_data()
    f.add_data(xs,ys0)
    plot(xs,ys0,"b.")
    f.set_priors(s0=1,omega=.5)
    phi,P=f(x); plot(x,phi,"r-",label="s0=1,omega=.5")
    f.set_priors(s0=1,omega=.1)
    phi,P=f(x); plot(x,phi,"k-",label="s0=1,omega=.1")
    title("p=5, no value error"); legend(loc=2); show()
    
    # p=5, small value error:
    ys = ys0 + numpy.random.normal(size=ys0.shape)*0.1
    figure()
    f = function.Function(d=1,p=5)
    f.set_default_errors(sy=0.1)
    f.clear_data()
    f.add_data(xs,ys)
    plot(xs,ys0,"y*")
    plot(xs,ys,"b.")
    f.set_priors(s0=1,omega=1)
    phi,P=f(x); plot(x,phi[:,:2],"y-",label="s0=1,omega=1")
    f.set_priors(s0=1,omega=.5)
    phi,P=f(x); plot(x,phi[:,:2],"r-",label="s0=1,omega=.5")
    f.set_priors(s0=1,omega=.1)
    phi,P=f(x); plot(x,phi[:,:2],"k-",label="s0=1,omega=.1")
    title("p=5, small value error of std.dev. 0.1"); legend(loc=2); show()
    
    # p=6, small value error:
    ys = ys0 + numpy.random.normal(size=ys0.shape)*0.1
    figure()
    f = function.Function(d=1,p=6)
    f.set_default_errors(sy=0.1)
    f.clear_data()
    f.add_data(xs,ys)
    plot(xs,ys0,"y*")
    plot(xs,ys,"b.")
    f.set_priors(s0=1,omega=1)
    phi,P=f(x); plot(x,phi[:,:2],"y-",label="s0=1,omega=1")
    f.set_priors(s0=1,omega=.5)
    phi,P=f(x); plot(x,phi[:,:2],"r-",label="s0=1,omega=.5")
    f.set_priors(s0=1,omega=.1)
    phi,P=f(x); plot(x,phi[:,:2],"k-",label="s0=1,omega=.1")
    title("p=6, small value error of std.dev. 0.1"); legend(loc=2); show()
    
    # p=6, medium value error:
    ys = ys0 + numpy.random.normal(size=ys0.shape)*0.5
    figure()
    f = function.Function(d=1,p=6)
    f.set_default_errors(sy=0.5)
    f.clear_data()
    f.add_data(xs,ys)
    plot(xs,ys0,"y*")
    plot(xs,ys,"b.")
    f.set_priors(s0=1,omega=1)
    phi,P=f(x); plot(x,phi[:,:2],"y-",label="s0=1,omega=1")
    f.set_priors(s0=1,omega=.5)
    phi,P=f(x); plot(x,phi[:,:2],"r-",label="s0=1,omega=.5")
    f.set_priors(s0=1,omega=.1)
    phi,P=f(x); plot(x,phi[:,:2],"k-",label="s0=1,omega=.1")
    title("p=6, medium value error of std.dev. 0.5"); legend(loc=2); show()
    
    # medium value error, different p:
    ys = ys0 + numpy.random.normal(size=ys0.shape)*0.5
    figure()
    plot(xs,ys0,"y*")
    plot(xs,ys,"b.")
    f = function.Function(d=1,p=4); f.set_default_errors(sy=0.5)
    f.add_data(xs,ys); f.set_priors(s0=1,omega=1)
    phi,P=f(x); plot(x,phi[:,:2],"y-",label="p=4")
    f = function.Function(d=1,p=6); f.set_default_errors(sy=0.5)
    f.add_data(xs,ys); f.set_priors(s0=1,omega=1)
    phi,P=f(x); plot(x,phi[:,:2],"r-",label="p=6")
    f = function.Function(d=1,p=10); f.set_default_errors(sy=0.5)
    f.add_data(xs,ys); f.set_priors(s0=1,omega=1)
    phi,P=f(x); plot(x,phi[:,:2],"k-",label="p=10")
    f = function.Function(d=1,p=20); f.set_default_errors(sy=0.5)
    title("medium value error of std.dev. 0.5, correct omega"); legend(loc=2); show()
    
    # sum of two sines and an exponential trend:
    x = linspace(0,1,1000)
    xs = numpy.random.uniform(size=100)
    ys = exp(xs) + 0.2*sin(19*xs) + 0.05*sin(97*xs) + 0.01*numpy.random.normal(size=100)
    plot(xs,ys,"k+",ms=10,alpha=1,label="measured points",zorder=10)
    f = function.Function(d=1,p=4)
    f.set_default_errors(sy=0.01)
    f.add_data(xs,ys)
    f.set_priors(s0=2,omega=10)
    phi,P=f(x); plot(x,phi[:,0],"y-",lw=5,alpha=0.4,label="p=4, omega=10")
    f.set_priors(s0=2,omega=100,mf=repeat(1.8,4),mr=1.5)
    phi,P=f(x); plot(x,phi[:,0],"r-",lw=3,alpha=0.5,label="p=4, omega=100, with estimated trend")
    f.set_priors(sf=repeat(inf,4),sr=2*100**4,mr=1.8)
    phi,P=f(x); plot(x,phi[:,0],"b-",lw=3,alpha=0.5,label="p=4, improper prior")
    title("MoTaBaR interpolation of measured data",size=40)
    legend(loc=2); show()
