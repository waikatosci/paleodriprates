fpath_='/Users/bedartha/Dropbox/bayprox/BayProX/'
fpath = '/home/goswami/Dropbox/bayprox/BayProX'
from matplotlib import use
# use('cocoaagg')
import os
import agedepth
import input
import scipy as sp
from scipy.interpolate import interp1d
from pprint import pprint
import proxyrecord
import simulate
import visualize

import matplotlib.pyplot as plt
# AgeDepth, ProxyDepth = simulate_core()
# info = input.SampleInfo('Test123', 'new', 'lake sediment')
# D = input.DatingTable(AgeDepth[:,0], AgeDepth[:,1], AgeDepth[:,2],
#                              'C14', info)
# P = input.ProxyDepth(ProxyDepth[:,0], ProxyDepth[:,1], 0, info)
# D.calibration()
# print dir(D)
# print D.info.name
# print D.info.datingmethod
# print info.datingmethod
# print AgeDepth.shape, ProxyDepth.shape
# # # pprint(dir(AgeDepth))
# # # pprint(dir(ProxyDepth))
# D.calibration()
# dwf = agedepth.DWF(D)
# dwf_final = dwf(P.depth, calBP_step=10, verbose=1)
# print dwf_final.shape
# Proxy_est = proxyrecord.ProxyEstimates()
# print dir(dwf)
# P_dist = proxyrecord.ProxyDistributions()
# pdf = P_dist.get_pdf(dwf, P)
# cdf = P_dist.get_cdf(pdf)
# lims = P_dist.proxyrange
# Med = proxyrecord.ProxyEstimates.get_median(cdf, lims, P_dist.res)
# print Med.shape, max(Med), min(Med)
# quants = proxyrecord.ProxyEstimates.get_quantiles(cdf, lims, P_dist.res, [0.25,0.75])
# print quants[0].shape, quants[1].shape
# print max(quants[0]), min(quants[0])
# print max(quants[1]), min(quants[1])

## get artificial archive & measuremetns
simulateArchive = simulate.Archive()
agemod = simulateArchive.growth(100, 20, 7)
proxymods = simulate.Archive.proxy()
measurements = simulate.Measurements()
age_measurements = measurements.radiometric(simulateArchive,
											sp.array([200, 10, 10, 20, 500]),
											num_samples=5)
proxy_measurements = measurements.proxy(simulateArchive, 0, 20)
deep = sp.linspace(0, simulateArchive.maxdepth, 10)
age = simulateArchive.agemodel(deep)
## input the measurements
info = input.SampleInfo('Test123', 'new', 'lake sediment')
D = input.DatingTable(age_measurements[0], age_measurements[1], age_measurements[2], 'C14', info)
P = input.ProxyDepth(proxy_measurements[0], proxy_measurements[1], 0, info)
D.calibration()
# pprint(age_measurements)
## get age depth DWFs
# print help(agedepth.DWF)
dwf = agedepth.DWF(D)
dwf_final = dwf(P.depth, calBP_step=10, verbose=0)
## get proxy estimates
P_dist = proxyrecord.ProxyDistributions()
Proxy_est = proxyrecord.ProxyEstimates()
pdf = P_dist.get_pdf(dwf, P, res=10)
cdf = P_dist.get_cdf(pdf)
lims = P_dist.proxyrange
Med = proxyrecord.ProxyEstimates.get_median(cdf, lims, P_dist.res)
quants = proxyrecord.ProxyEstimates.get_quantiles(cdf, lims, P_dist.res, [0.25,0.75])

ptype = visualize.PlotType(name='trunc', size="2col", num_rows=2)
ptype.trunc(dwf, D, P, Proxy_est)
# plt.show()
# plt.savefig('./sandbox/test.pdf')
pprint(dir(ptype))

# DWF.caldat -> = [cal_ax, rm_age, rm_age_err]
# DWF.rmagemod -> [eval_depth, rm_mean, rm_dev]
# ProxyDepth.depth/proxy/proxyerror
# DatingTable.depth/age/ageerror
# ProxyEstimates.expected_value/variance/median_value/quantiles