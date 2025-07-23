"""	Contains classes and methods that visualize the data and results.

    This module allows you to visualize the relevant datasets along with the
	results of the analysis in concise, publication-ready format. The standard
	option is to plot the information and results in the typical 2x2 format as
	used in the original BayProX manuscript. You also have options to plot the 
	proxy record, radiometric (RM) age model, calibration curve and proxy-depth
	measurement datasets separately.
"""

import scipy as sp
from copy import deepcopy
from matplotlib.ticker import FixedLocator, FixedFormatter, MaxNLocator
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pprint import pprint
from matplotlib.collections import LineCollection

class Drawers(object):
	"""Contains methods for drawing various plot objects."""
	def __init__(self):
		pass

	def plot_line(self, ax, x, y, params):
		lwidth, lcolor = params['linewidth'], params['linecolor']
		ax.plot(x, y, lw=lwidth, color=lcolor, clip_on=False)

	def plot_fill(self, ax, x, yHi, yLo, params):
	    """plot narrower confidence bounds"""
	    fcolor = params['fillcolor']
	    ax.fill_between(x, yHi, yLo, color=fcolor, zorder=0, clip_on=False)

	def plot_o(self, ax, x, y, err, params):
		"""plot errorbars"""
		mcolor, msize = params['markercolor'], params['markersize']
		mcapsize = params['markercapsize']
		e = ax.errorbar(x, y, yerr=err, fmt='.', c=mcolor, mec=mcolor, mew=1.0,
				 	ms=msize, capsize=mcapsize, zorder=ax.get_zorder()+4)
		for item in e:
			if isinstance(item, list) or isinstance(item, tuple):
				for subitem in item:
					subitem.set_clip_on(False)
			else:
				item.set_clip_on(False)

	def axis_ticks(self, ax, axis, where, tiklocs, params=None):
		if not params:
			size, width, col = 6., 1., 'k'
		else:
			size = params['ticklength']
			width = params['tickwidth']
			color = params['tickcolor']
		vmin, vmax = self.get_conjugate_axis_view_interval(ax, axis)
		v, x = self.parse_location(where, vmin, vmax)
		v *= 0.95
		vrange = abs(vmax - vmin)
		tiklen = self.convert_points_to_datacoords(ax, axis, size, vrange)
		ntiks = len(tiklocs)
		tik_start = self.get_tick_endpoints(axis, tiklocs, [v] * ntiks)
		tik_stop = self.get_tick_endpoints(axis, tiklocs, [v+x*tiklen] * ntiks)
		tiks = zip(tik_start, tik_stop)
		lc = LineCollection(tiks,
							linewidths=width, colors=color, linestyles='solid',
							clip_on=False)
		ax.add_collection(lc)

	def convert_points_to_datacoords(self, ax, axis, size_in_pts, datarange):
		if axis == 'x':
			figheight = plt.rcParams['figure.figsize'][1]
			axheight = ax.get_position().height
			axsize_inches = axheight * figheight
		elif axis == 'y':
			figwidth = plt.rcParams['figure.figsize'][0]
			axwidth = ax.get_position().width
			axsize_inches = axwidth * figwidth			
		points_to_dataunits = (1./72.) * (datarange / axsize_inches)
		size_dataunits = points_to_dataunits * size_in_pts
		return size_dataunits

	def get_tick_endpoints(self, axis, xx, yy):
		tick_endpoints = zip(xx, yy)
		if axis == 'y':
			tick_endpoints = zip(yy, xx)
		return tick_endpoints

	def axis_spines(self, ax, axis, where, tiklocs, params=None):
		if not params:
			width, color = 1., 'k'
		else:
			width, color = params['spinewidth'], params['spinecolor']
		vmin, vmax = self.get_conjugate_axis_view_interval(ax, axis)
		v, x = self.parse_location(where, vmin, vmax)
		v *= 0.95
		xx, yy = (min(tiklocs), max(tiklocs)), (v, v)
		if axis == 'y': xx, yy = yy, xx
		spine = Line2D(xx, yy, lw=width, color=color, clip_on=False)
		ax.add_artist(spine)

	def get_conjugate_axis_view_interval(self, ax, axis):
		axes_names = set(['x', 'y'])
		conj_axis = list(axes_names - set([axis]))[0]
		vmin, vmax = eval('ax.get_' + conj_axis + 'axis().get_view_interval()')
		return vmin, vmax

	def parse_location(self, where, vmin, vmax):
		if where in ['bottom', 'left']:
			v, x = vmin, 1
		elif where in ['top', 'right']:
			v, x = vmax, -1
		return v, x

	def tick_labels(self, ax, axis, where, tiklabs, tiklocs,
					  v_off=None, h_off=None, params=None):
		if not params:
			font, size, color =  'Helvetica', 10., 'k'
		else:
			font = params['ticklabelfont']
			size = params['ticklabelsize']
			color = params['ticklabelcolor']
		ntiks = len(tiklocs)
		vmin, vmax = self.get_conjugate_axis_view_interval(ax, axis)
		v, x = self.parse_location(where, vmin, vmax)
		vrange = abs(vmax - vmin)
		offset = 0.05*self.convert_points_to_datacoords(ax, axis, size, vrange)
		v *= 1.375
		v_off, h_off = offset, 0.
		xx, yy = tiklocs, [v+offset]*ntiks
		h_align, v_align = 'center', 'top'
		if axis == 'y':
			xx, yy = yy, xx
			h_align, v_align = v_align, h_align
			v *= 0.55
		for i in range(ntiks):
			if tiklabs[i] == '-0': tiklabs[i] = '0.0'
			ax.text(xx[i], yy[i], tiklabs[i],
					ha='center', va='center',
					fontname=font, size=size, style='normal', color=color)

class Formatters(object):
	"""Contains methods for formatting the plots."""
	def __init__(self):
		pass

	def axrem(self, ax):
		"""remove frame and axis ticks"""
		# ax.set_frame_on(False)
		ax.tick_params(axis='both',
					   top='off', labeltop='off',
					   bottom='off', labelbottom='off',
					   left='off', labelleft='off',
					   right='off', labelright='off')

	def axrev(self, ax, axis):
		"""reverse axis direction"""
		if axis == 'x':
			ax.set_xlim(ax.get_xlim()[::-1])
		elif axis == 'y':
			ax.set_ylim(ax.get_ylim()[::-1])

	def axview(self, ax, axis='both', Lo=None, Hi=None):
		"""set view intervals"""
		if not (Lo or Hi):
			if axis == 'both':
				for axname in ['x', 'y']:
					Lo, Hi = self.get_axlim(ax, axname)
					if axis == 'x': Lo, Hi = Hi, Lo
					eval('ax.set_' + axname + 'lim(Lo, Hi)')
			else:
				Lo, Hi = self.get_axlim(ax, axis)
				eval('ax.set_' + axis + 'lim(Lo, Hi)')
		else:
			eval('ax.set_' + axis + 'lim(Lo, Hi)')

	def get_axlim(self, ax, axis):
		Lo = eval('ax.dataLim.' + axis + 'min')
		Hi = eval('ax.dataLim.' + axis + 'max')
		frac, datarange = 0.0, abs(Hi - Lo)
		# if ax.is_first_col():# and ax.is_first_row()
		# 	frac = 0.0
		Lo, Hi = Lo - frac*datarange, Hi + frac*datarange
		return Lo, Hi

	def axtiks(self, ax, axis, where, params,
			   num_ticks=None, tikspan=None, spines=True, optimize=True):
		"""set axis ticks and tick labels and draw axis spines if specified"""
		if num_ticks and tikspan: optimize = False
		if tikspan:
			tikspan = [float(i) for i in tikspan]
			Lo, Hi = tikspan
		else:
			data = eval('ax.get_lines()[0].get_' + axis + 'data()')
			Lo, Hi = min(data), max(data)
		if not num_ticks:
			num_ticks = 7
		if optimize:
			Lo, Hi, num_ticks = self.optimize_tiklocs(Lo, Hi)
		tiklocs = self.get_tiklocs(Lo, Hi, num_ticks)
		draw = Drawers()
		draw.axis_ticks(ax, axis, where, tiklocs, params=params)
		tiklocs_formatted = self.get_tiklocs(Lo, Hi, num_ticks, decimals=1)
		tiklabs = [str(loc) for loc in tiklocs_formatted]
		draw.tick_labels(ax, axis, where, tiklabs, tiklocs, params=params)
		if spines:
			draw.axis_spines(ax, axis, where, tiklocs, params=params)

	def optimize_tiklocs(self, Lo, Hi, min_num_tiks=5, max_num_tiks=6):
		"""The scoring is based on awarding positive penalties for bad ticks 
			and rewarding negative penalties for nice ticks. the scheme is:
				* odd endings 	=> +1
				* even endings 	=> 0
				* zero endings	=> -2
			In this scheme, a minimal score will give the best scheme of ticks.
		"""
		datarange = Hi - Lo
		thres = 0.05 * datarange
		endpoint = [self.set_endpoint(Lo, thres), self.set_endpoint(Hi, thres)]
		score, current_score = 0, 0
		numtik_out = min_num_tiks
		# minimize the score wrt num_ticks
		for num_tik in range(min_num_tiks, max_num_tiks+1): 
			tikpos = self.get_tiklocs(endpoint[0], endpoint[1],
									  num_tik, decimals=1)
			last_digit = [int(str(i)[-1]) for i in tikpos]
			for i in [0, -1]:
				if last_digit[i] == 0:
					current_score -= 2
					del last_digit[i]
			current_score += sum([i%2 for i in last_digit])
			if current_score < score:
				score = current_score
				numtik_out = num_tik
		return endpoint[0], endpoint[1], numtik_out

	def set_endpoint(self, endpt, thres):
		maxidx = 1.
		while round(endpt/maxidx) != 0.:
			maxidx = maxidx * 10.
		maxidx = sp.log10(maxidx / 10.)
		for i in sp.arange(1, maxidx+1):
			nu_endpt = round(endpt / 10**i) * 10**i
			dist = abs(nu_endpt - endpt)
			if dist < thres:
				endpt = nu_endpt
		return endpt

	def get_tiklocs(self, Lo, Hi, num_ticks, decimals=None):
		step = abs(Hi - Lo)/(num_ticks - 1)
		tiklocs = sp.arange(Lo, Hi + 0.001, step)
		if decimals:
			if any(tiklocs  > 1000.):
				tiklocs = tiklocs/1000.
			tiklocs = sp.around(tiklocs, decimals=decimals)
		return tiklocs	

class PlotType(Drawers, Formatters):
	"""Contains methods for specifying the plot types."""
	def __init__(self, name, size='2col', num_rows=2, top=0.5):
		self.name = name
		self.size = size
		if name != "trunc":
			num_rows = 1
			axis_heights = self.set_axis_heights(num_rows, top=0.95)
		else:
			axis_heights = self.set_axis_heights(num_rows, top=top)
		def_params = {'subplotlabel': {'font': 'Helvetica',
									   'size': 14,
									   'color': 'black'
									   },
					  'calcurve': {'linewidth': 1,
					  			   'linecolor': 'red',
					  			   'fillcolor': 'lightblue',
					  			   'sigma': 2.,
					  			   'ticklength': 8.,
					  			   'tickwidth': 1.,
					  			   'tickcolor': 'k',
					  			   'ticklabelfont': 'Helvetica',
					  			   'ticklabelsize': 10.,
					  			   'ticklabelcolor': 'k',
					  			   'axislabelfont': 'Helvetica',
					  			   'axislabelsize': 12,
					  			   'axislabelcolor': 'k',
					  			   'spinecolor': 'k',
					  			   'spinewidth':1.
					  			   },
					  'rmagemod': {'linewidth': 1,
					  			   'linecolor': 'red',
					  			   'fillcolor': 'lightblue',
					  			   'sigma': 2.,
					  			   'markercolor': 'black',
					  			   'markersize': 8.,
					  			   'markercapsize': 5.,
					  			   'ticklength': 8.,
					  			   'tickwidth': 1.,
					  			   'tickcolor': 'k',
					  			   'ticklabelfont': 'Helvetica',
					  			   'ticklabelsize': 10.,
					  			   'ticklabelcolor': 'k',
					  			   'axislabelfont': 'Helvetica',
					  			   'axislabelsize': 12,
					  			   'axislabelcolor': 'k',
					  			   'spinecolor': 'k',
					  			   'spinewidth':1.
					  			   },
					  'proxydepth': {'linewidth': [1],
					  			     'linecolor': ['black'],
					  			     'fillcolor': [None],
					  			     'sigma': [None],
					  			     'ticklength': [8.],
					  			   	 'tickwidth': [1.],
						  			 'tickcolor': ['k'],
						  			 'ticklabelfont': ['Helvetica'],
						  			 'ticklabelsize': [10.],
						  			 'ticklabelcolor': ['k'],
						  			 'axislabelfont': ['Helvetica'],
						  			 'axislabelsize': [12],
						  			 'axislabelcolor': ['k'],
						  			 'spinecolor': ['k'],
						  			 'spinewidth':[1.]
					  			   	 },
					  'proxyrecord': {'linewidth': [1],
					  			      'linecolor': ['red'],
					  			      'fillwide': ['lightblue'],
					  			      'fillthin': ['darkblue'],
					  			      'sigma': [2, 4],
					  			      'ticklength': [8.],
					  			   	  'tickwidth': [1.],
						  			  'tickcolor': ['r'],
						  			  'ticklabelfont': ['Helvetica'],
						  			  'ticklabelsize': [10.],
						  			  'ticklabelcolor': ['r'],
						  			  'axislabelfont': ['Helvetica'],
						  			  'axislabelsize': [12],
						  			  'axislabelcolor': ['r'],
						  			  'spinecolor': ['r'],
						  			  'spinewidth':[1.]
					  			   	  },
					  'figure': {'leftmargin': 0.05,
					  			 'bottmargin': 0.100,
					  			 'centmargin': 0.500,
					  			 'inter_col_space': 0.1,
					  			 'inter_row_space': 0.1,
					  			 'colwidth': [0.4, 0.4],
					  			 'rowheight': axis_heights,
					  			 },
					  }
		self.params = {'default': deepcopy(def_params),
					   'current': deepcopy(def_params),
					   }
		if size == "1col":
			self.set_1col_params()
		# TODO: The following TeX params update leads to plt.show() error
		# texparams = {'text.usetex': True,
		# 			 'text.latex.preamble': [
		# 			 r"\usepackage[adobe-utopia, sfscaled=true]{mathdesign}",
		# 			 r"\usepackage{wasysym}", 
		# 			 r"\renewcommand{\sfdefault}{phv}"],
		# 			 'savefig.dpi': 120,
		# 			}
		# plt.rcParams.update(texparams)

	def set_axis_heights(self, num_rows, top=0.5, bottom=0.25):
		self.num_rows = num_rows
		axhite = [top]
		axhite.extend([bottom] * (self.num_rows - 1))
		return axhite

	def set_fig_params(self, **kwargs):
		for key in kwargs:
			obj, prop = self.parse_params(key)
			val = kwargs[key]
			self.params['current'][obj][prop] = val

	def get_fig_params(self, *args):
		for key in args:
			obj, prop = self.parse_params(key)
			val = self.params['current'][obj][prop]
			return val

	def set_1col_params(self):
		axis_heights = self.set_axis_heights(self.num_rows, bottom=0.3)
		params_1col = {'axislabel': {'size': 10},
					   'subplotlabel': {'size': 12},
					   'ticklabel': {'size': 8},
					   'ticksize': {'calage': 8,
					  			   'rmage': 8,
					  			   'depth': 8,
					  			   'proxy': 8},
					   'rmagemod': {'markersize': 5,
					   				'markercapsize': 4,
					   				},
					   'figure': {'leftmargin': 0.000,
					  			  'bottmargin': 0.075,
					  			  'centmargin': 0.510,
					  			  'inter_col_space': 0.05,
					  			  'inter_row_space': 0.05,
					  			  'colwidth': 0.500,
					  			  'rowheight': axis_heights
					  			 }
					  }
		params_2col = self.params['default']
		for key in params_2col.iterkeys():
			if key in params_1col.keys():
				self.params['default'][key].update(params_1col[key])
		self.params['current'] = deepcopy(self.params['default'])

	def parse_params(self, key):
		separator = key.index('_')
		obj = key[0:separator]
		prop = key[separator+1:]
		return obj, prop

	def reset_fig_params(self):
		self.params['current'] = deepcopy(self.params['default'])

	def trunc(self, DWF=[], DatingTable=[], ProxyDepth=[], ProxyEstimates=[],
			  TrueProxy=None,
			  center="median", spread="quantiles", fillslim=False,
			  figsize=[7.48031, 3.93701]):
		"""Track Uncertainty (TrUnc) visualizes all data in a 2-column grid."""
		drawer = Drawers()
		plotparams = self.params['current']
		plotnames = ['calcurve', 'rmagemod', 'proxyrecord', 'proxydepth']
		ncols, nrows = 2, self.num_rows
		prox_ax_hite = 1.37795
		if self.size == "1col":
			figsize = [3.54331, 2.36220]
			prox_ax_hite = 0.70866
		figsize[1] = figsize[1] + (nrows - 2) * prox_ax_hite
		plt.rcParams['figure.figsize'] = figsize
		fig = plt.figure(facecolor='w', edgecolor='w')
		num_plots = ncols * nrows
		axname = []
		x_off = plotparams['figure']['inter_col_space']
		y_off = plotparams['figure']['inter_row_space']
		bm = 1.0
		axes_bgcols = ['w']*4 #['r', 'b', 'g', 'y']
		for i in range(num_plots):
			axname.append('ax' + str(i))
			if i == 0:
				axX, axY = None, None
			elif i == 1:
				axX, axY = None, eval('ax0')
			elif i%2 == 0 and i != 0:
				axX, axY = eval('ax0'), None
			elif i%2 == 1 and i != 1:
				axX, axY = eval('ax1'), eval('ax' + str(i-1))
			exec(axname[i] + ' = fig.add_subplot(nrows,ncols,' + str(i+1) + ','
				 + 'sharex=axX, sharey=axY, axis_bgcolor=axes_bgcols[i], ' + 
				 'clip_on=False)')
			if i%2 == 0:
				lm = plotparams['figure']['leftmargin']
				xh = plotparams['figure']['rowheight'][i/2]
				bm = bm -xh - y_off 
			else:
				lm = plotparams['figure']['centmargin'] + x_off*0.5
			wd = plotparams['figure']['colwidth'][i%2]
			exec(axname[i] + '.set_position([lm,bm,wd,xh])')
		# plot the top row
		Hi_pt, Lo_pt = 0., 0.
		for i in [0,1]:
			ax = eval('ax' + str(i))
			objname = plotnames[i]
			params = deepcopy(plotparams[objname])
			[x, y, y_err] = eval('DWF.' + objname)
			s = params['sigma']			
			errHi, errLo = y + s*y_err, y - s*y_err
			if max(errHi) > Hi_pt: Hi_pt = max(errHi)
			if min(errLo) < Lo_pt: Lo_pt = min(errLo)
			ax = eval('ax' + str(i))
			drawer.plot_line(ax, x, y, params)
			drawer.plot_fill(ax, x, errHi, errLo, params)
			if i == 1:
				drawer.plot_o(ax,
							   DatingTable.depth,
							   DatingTable.age,
							   s * DatingTable.ageerror,
							   params)
				ageerrHi = DatingTable.age + s * DatingTable.ageerror
				ageerrLo = DatingTable.age - s * DatingTable.ageerror
				if max(ageerrHi) > Hi_pt: Hi_pt = max(ageerrHi)
				if min(ageerrLo) < Lo_pt: Lo_pt = min(ageerrLo)
			del params
		# plot the rest of the rows
		for i in range(2, num_plots):
			ax = eval('ax' + str(i))
			objname = plotnames[i]
			params = deepcopy(plotparams[objname])
			if i%2 == 0:	# proxy record i.e. left column
				if TrueProxy is not None:
					drawer.plot_line(ax, TrueProxy[0], TrueProxy[1],
									 {'linewidth': 1., 'linecolor': '0.3'})
				k = i/2 - 1
				for key in params.iterkeys():
				    val = params[key][k]
				    params[key] = val
				x = DWF.calcurve[0]
				y = eval('ProxyEstimates.' + center + '_value')				
				drawer.plot_line(ax, x, y, params)
				conf_bnds = eval('ProxyEstimates.' + spread) ##BUG@VARIANCE!!
				filltypes = ['fillwide']
				if fillslim: filltypes = ['fillwide', 'fillslim']
				prxHi, prxLo = [0.]*nrows, [0.]*nrows
				for bounds, fill in zip(conf_bnds, filltypes):
					params['fillcolor'] = params[fill]
					drawer.plot_fill(ax, x, bounds[1], bounds[0], params)
					if max(bounds[1]) > prxHi[k]: prxHi[k] = max(bounds[1])
					if min(bounds[0]) < prxLo[k]: prxLo[k] = min(bounds[0])
			else:			# proxy-depth i.e. right column
				k = i/2 - 1
				for key in params.iterkeys():
				    val = params[key][k]
				    params[key] = val
				x = ProxyDepth.depth
				y = ProxyDepth.proxy
				yerr = ProxyDepth.proxyerror
				drawer.plot_line(ax, x, y, params)
			del params
		# customize the look of the plot
		format = Formatters()
		# top row
		for i in range(num_plots):
			ax = eval('ax' + str(i))
			format.axrem(ax)
			objname = plotnames[i]
			params = deepcopy(plotparams[objname])
			if i == 0:
				format.axview(ax, axis='x')
				# ax.set_xlim(min(DWF.calcurve[0]), max(DWF.calcurve[0]))
				format.axrev(ax, 'x')
				# format.axtiks(ax, 'x', 'bottom', params)
				# format.axtiks(ax, 'y', 'right', params)
			elif i == 1:
				# ax.set_xlim(min(ProxyDepth.depth), max(ProxyDepth.depth))
				# ax.set_axis_off()
				format.axview(ax, axis='x')
				# print Lo_pt, Hi_pt
				format.axview(ax, axis='y', Lo=Lo_pt, Hi=Hi_pt)
				# format.axtiks(ax, 'x', 'bottom', params)
			elif i%2 == 0:
				k = i/2 - 1
				for key in params.iterkeys():
					val = params[key][k]
					params[key] = val
				format.axrev(ax, 'y')
				format.axview(ax, axis='y', Lo=prxLo[k], Hi=prxHi[k])
				# format.axtiks(ax, 'y', 'right', params)
			del params
			# print ax.dataLim.max
			# print ax.dataLim.min
			# print ax.get_position()
		# pprint(dir(fig))

	def rmagemod():
		pass

	def calcurve():
		pass

	def seacliff():
		pass

	def proxydepth():
		pass

	def proxyrec():
		pass
