#!/usr/local/python
import os
import numpy as np
import pickle
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
import datetime
import sys
from matplotlib.lines import Line2D
from matplotlib import rc
from matplotlib import rcParams
import cPickle as pickle

rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
rcParams['mathtext.default'] = 'regular'
fig_size = plt.rcParams["figure.figsize"]
# Prints: [8.0, 6.0]

# Set figure width to 9 and height to 9
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=False)
rcParams['mathtext.default'] = 'regular'
fig_size = plt.rcParams["figure.figsize"]
# Prints: [8.0, 6.0]

# Set figure width to 9 and height to 9
fig_size[1] = 9
fig_size[0] = 9
size = 20
plt.rcParams["figure.figsize"] = fig_size
matplotlib.rcParams.update({'font.size': 22})
plt.ioff()
print "Current size:", fig_size


def single_uvcov(u, v, freqs, ax,color,telescope_name):
	c = 299792458.
        print('Plotting uv coverage for %s' % telescope_name)
	for i, freqi in enumerate(freqs):
                print('Frequency plotted: %.5e' % freqi)
		fc = freqi/c/1e6
		if i == 0:
			ax.plot(+u*fc, +v*fc, marker='.', ms=0.01, ls='',
					 color=color, mec=color, label=telescope_name,alpha=0.005, rasterized=True)
			ax.plot(-u*fc, -v*fc, marker='.', ms=0.01, ls='',
					 color=color, mec=color,alpha=0.005,rasterized=True)
		ax.plot(+u*fc, +v*fc, marker='.', ms=0.01, ls='',
				 color=color, mec=color,alpha=0.005,rasterized=True)
		ax.plot(-u*fc, -v*fc, marker='.', ms=0.01, ls='',
				 color=color, mec=color,alpha=0.005,rasterized=True)
	#lgnd = ax.legend(numpoints=1, markerscale=6, frameon=False, ncol=2, prop={'size':8})
	return ax


def make_uvcov(msfiles,plotfile):
	#logger.info('Plotting uv-coverage for all sources'.format())
	custom_lines = []
	colors = ['k','r','b','g']
	telescope_names = []
	fig = plt.figure()
	ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
	for i, msfile in enumerate(msfiles):
		arrays = pickle.load( open( msfile, "rb" ) )
		u = arrays['u']
		v = arrays['v']
		freqs = arrays['freqs']
		telescope_name = msfile.split('_uv.pkl')[0]

		#for f in msinfo['sources']['allsources'].split(','):
		#ax1 = fig.add_axes([0.8, 0.1, 0.02, 0.8])
		single_uvcov(u, v, freqs, ax,color=colors[i],telescope_name=telescope_name)
		custom_lines = custom_lines + [Line2D([0], [0], color=colors[i], lw=4)]
		telescope_names = telescope_names + [telescope_name]

	ax.legend(custom_lines, telescope_names)
		#ax1.yaxis.set_label_position("right")
		#ax1.set_ylabel('Frequency [GHz]')
	ax.set_xlabel('V [Mlambda]')
	ax.set_ylabel('U [Mlambda]')
	ax.set_aspect('equal')
	main_lim = np.max(np.abs(ax.get_ylim() + ax.get_xlim()))
	ax.set_xlim(-main_lim, +main_lim)
	ax.set_ylim(-main_lim, +main_lim)
	#ax.legend()
	fig.savefig(plotfile, dpi=1500, bbox_inches='tight')


try:
	plotfile = sys.argv[sys.argv.index('make_uvcoverage_plots_pickle.py')+1]
	msfiles = sys.argv[sys.argv.index('make_uvcoverage_plots_pickle.py')+2:]
except:
	print('Error: Usage python make_uvcoverage_plots.py plotfile <pkl1> <pkl2> etc.')

make_uvcov(msfiles,plotfile)
