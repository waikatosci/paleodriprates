import os
import scipy as sp
os.chdir('../')
import motabar.function as motabar
os.chdir('./calibration')
fpath = './hua_barbetti/'
fname = 'hua_barbetti_'
fout = './postbomb/'
pb_zones = ['NH1', 'NH2', 'NH3', 'SH']
timeBPend = -63

# regression paramters
max_degree = 2
pre_remvar = None
verbose = 0

for zone in pb_zones:
    print('Processing post-bomb zone '+zone+'...')
    file = fpath+fname+zone+'.csv'
    postbomb = sp.genfromtxt(file, delimiter=',')
    postbomb = postbomb[3::,:]
    timeAD = postbomb[:,0]
    timeBP = 1950 - timeAD
    f14C = postbomb[:,3]
    f14C_err = postbomb[:,4]
    regression = motabar.Function(p=max_degree, verbosity=verbose)
    regression.set_priors(sr=pre_remvar)
    regression.set_data(timeBP, f14C, sy=f14C_err)
    estimate, precision = regression(timeBPend)
    est_mean = estimate[0]
    est_dev = sp.sqrt(precision[0,0] ** (-1))
    f14C = sp.hstack((f14C, est_mean))
    f14C_err = sp.hstack((f14C_err, est_dev))
    timeBP = sp.hstack((timeBP, timeBPend))
    C14 = -8033*sp.log(f14C)
    C14_err = -8033*(f14C_err/f14C)
    X = sp.vstack((timeBP,C14,C14_err)).T
    sp.savetxt(fout+zone+'.csv',X,delimiter=',')
    print('Zone '+zone+' done! CSV file saved to '+fout+zone+'.csv')
