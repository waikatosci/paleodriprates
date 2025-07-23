""" Contains classes and methods that allow you to generate synthetic datasets.

    This module first allows you to generate the growth of an artifical proxy
    archive with its two crucial aspects - the age-depth relation, and the
    proxy-depth data. After creating an articial proxy archive thus, you can
    simulate radiometric measurements on the 'true' age-depth relation as well
    as proxy measurements (sampled with a much higher density than the
    rediometric ages) on the 'true' proxy-depth curve. The knowledge of the
    'true' curves allow you to test the efficacy of the methods used here as
    well as explore the impacts of uncertainties and variability on final
    intepretations.
"""

import scipy as sp
from scipy.interpolate import interp1d
import numpy.random as rand
norm = rand.normal

class Archive():
	"""Contains methods to simulate age-depth & proxy-depth relations."""
	def __init__(self, name='Test123'):
		self.name = name

	@classmethod
	def growth(self, num_layers=100, startRate=20, variability=7):
		"""Returns age-depth relation as a SciPy.Interpolate object"""
		age = sp.zeros(num_layers)
		depth_layers = sp.arange(num_layers)
		growRate = sp.zeros(num_layers)
		for layer in depth_layers[1:]:
			growRate[layer] = growRate[layer-1] + variability*norm()
			while growRate[layer] < 0:
				growRate[layer] = growRate[layer-1] + variability*norm()
			age[layer] = sum(growRate[1:layer] *
							(depth_layers[1:layer] - depth_layers[0:layer-1]))
		agemodel = interp1d(depth_layers, age)
		agemodel_inverse = interp1d(age, depth_layers)
		self.agemodel, self.agemodel_inverse = agemodel, agemodel_inverse
		self.maxage, self.maxdepth = age[-1], depth_layers[-1]
		return agemodel

	@classmethod
	def proxy(self, type='sine', params=[0.001, 0.005]):
		"""Returns proxy-time & proxy-depth as Scipy.Interpolate object."""
		maxAge = self.maxage
		age = sp.arange((maxAge)).squeeze()
		if type == 'sine':
			self.type = type
			w1, w2 = params
			proxy_curve = sp.sin(age*w1*2*sp.pi) + sp.sin(age*w2*2*sp.pi)
		age_to_depth = self.agemodel_inverse
		depth = age_to_depth(age) # PRINTS WARNINGS ABOUT DIVISION: depth[0]=nan
		depth[0] = 0. ## CHECK OUT THIS BUG!! THIS LINE SHOULDN'T BE THERE!
		proxy_vs_depth = interp1d(depth, proxy_curve)
		proxy_vs_time = interp1d(age, proxy_curve)
		self.proxydepth, self.proxy_time = proxy_vs_depth, proxy_vs_time
		return proxy_vs_time, proxy_vs_depth



class Measurements():
	"""Contains methods to simulate age-depth & proxy-depth measurements."""
	def __init__(self, method='C14'):
		self.method = method

	def radiometric(self, Archive, error, num_samples=15, errtype='uniform'):
		"""Returns a set of error-based radiometric measurements."""
		measured_depths = sp.linspace(0, Archive.maxdepth, num_samples)
		if type(error) is int:
			error = error * sp.ones(num_samples)
		if errtype == 'minmax':
			error = sp.linspace(error[0], error[1], num_samples)
		measured_ages = Archive.agemodel(measured_depths) + \
						error * norm(size=(error.shape))
		self.age_measurements = measured_depths, measured_ages, error
		return self.age_measurements

	def proxy(self, Archive, error, num_samples=300):
		"""Returns a set of error-based proxy measurements."""
		measured_depths = sp.linspace(0, Archive.maxdepth-1, num_samples)
		if type(error) is int or float:
			error = error * sp.ones(num_samples)
		proxy = Archive.proxydepth(measured_depths) + \
				error * norm(size=(error.shape))
		self.proxy_measurements = measured_depths, proxy, error
		return self.proxy_measurements