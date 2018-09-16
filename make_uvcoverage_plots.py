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

def get_freqs(msfile, allfreqs=False):
	ms.open(msfile)
	axis_info = ms.getdata(['axis_info'],ifraxis=True)
	ms.close()
	if allfreqs:
		channels = axis_info['axis_info']['freq_axis']['chan_freq']
		freqs=np.sort(axis_info['axis_info']['freq_axis']['chan_freq'].flatten())[::len(channels)/4]
	else:
		freqs = axis_info['axis_info']['freq_axis']['chan_freq'].mean(axis=0)
	return freqs

def single_uvcov(u, v, freqs, ax,color,telescope_name):
	c = 299792458.
	for i, freqi in enumerate(freqs):
		fc = freqi/c/1e6
		if i == 0:
			ax.plot(+u*fc, +v*fc, marker='.', ms=0.01, ls='',
					 color=color, mec=color, label=telescope_name,alpha=0.1)
			ax.plot(-u*fc, -v*fc, marker='.', ms=0.01, ls='',
					 color=color, mec=color,alpha=0.1)
		ax.plot(+u*fc, +v*fc, marker='.', ms=0.01, ls='',
				 color=color, mec=color,alpha=0.1)
		ax.plot(-u*fc, -v*fc, marker='.', ms=0.01, ls='',
				 color=color, mec=color,alpha=0.1)
	#lgnd = ax.legend(numpoints=1, markerscale=6, frameon=False, ncol=2, prop={'size':8})
	return ax


def read_uvw(msfile, field):
	ms.open(msfile)
	staql={'field':field, 'spw':'0'}
	ms.msselect(staql)
	uv = ms.getdata(['u', 'v'])
	ms.close()
	u = uv['u']
	v = uv['v']
	return u, v

def make_uvcov(msfiles,save_uv,plotfile):
	#logger.info('Plotting uv-coverage for all sources'.format())
	custom_lines = []
	colors = ['k','r','b','g']
	telescope_names = []
	fig = plt.figure()
	ax = fig.add_axes([0.1, 0.1, 0.8, 0.8])
	for i, msfile in enumerate(msfiles):
		freqs = get_freqs(msfile, allfreqs=True)
		tb.open(msfile+'/FIELD')
		fields_ms = tb.getcol('NAME')
		tb.close()
		tb.open(msfile+'/OBSERVATION')
		telescope_name = tb.getcol('TELESCOPE_NAME')[0]
		print telescope_name
		tb.close()

		#for f in msinfo['sources']['allsources'].split(','):
		#ax1 = fig.add_axes([0.8, 0.1, 0.02, 0.8])
		if save_uv == 'True':
			u_save = np.array([])
			v_save = np.array([])
		for f in fields_ms:
			u, v = read_uvw(msfile, f)
			single_uvcov(u, v, freqs, ax,color=colors[i],telescope_name=telescope_name)
			if save_uv == 'True':
				u_save = np.hstack([u_save,u])
				v_save = np.hstack([v_save,v])
		else:
			print 'Cannot plot uvcov for {0}. Source not in ms.'.format(f)
		if save_uv == 'True':
			uv_dict = {'u':u_save,'v':v_save,'freqs':freqs}
			pickle.dump(uv_dict, open( "%s_uv.pkl"%telescope_name, "wb" ) )
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
	fig.savefig(plotfile, dpi=150, bbox_inches='tight')


try:
	plotfile = sys.argv[sys.argv.index('make_uvcoverage_plots.py')+1]
	save_uv = sys.argv[sys.argv.index('make_uvcoverage_plots.py')+2]
	msfiles = sys.argv[sys.argv.index('make_uvcoverage_plots.py')+3:]
except:
	print('Error: Usage casa -c make_uvcoverage_plots.py save_uv plotfile <ms1> <ms2> etc.')

make_uvcov(msfiles,save_uv,plotfile)
