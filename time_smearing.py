import numpy as np
from scipy.special import erf
from scipy.interpolate import UnivariateSpline
from matplotlib import pyplot as plt

def RdeltatInverse(r, r_hpbw, tau):

  retVar = 1.e0 - 1.22e-9*tau*tau*(r/r_hpbw)**2.

  return retVar

plt.close('all')
pix_scale = 1

'''
# VLA-FIRST values
theta_max_arcmin = 30.83
theta_max_arcsec = theta_max_arcmin*60.
theta_max = theta_max_arcmin/pix_scale
r_max_arcsec = 10.e0
r_max_arcmin = r_max_arcsec*60.
r_max = r_max_arcmin/pix_scale
r_hpbw_arcsec = 5
r_hpbw_arcmin = r_hpbw_arcsec*60.
r_hpbw = r_hpbw_arcmin/pix_scale
delta_nu = 3.e6
nu0 = 1.4e9
tau = 5

'''

# SKA1-MID values
theta_max_arcmin = 20.73
theta_max_arcsec = theta_max_arcmin*60.
theta_max = theta_max_arcmin/pix_scale
r_max_arcsec = 3.e0
r_max_arcmin = r_max_arcsec*60.
r_max = r_max_arcmin/pix_scale
r_hpbw_arcsec = 0.5
r_hpbw_arcmin = r_hpbw_arcsec*60.
r_hpbw = r_hpbw_arcmin/pix_scale
delta_nu = 1.e6
nu0 = 1.4e9
tau = 1

r = np.linspace(-r_max_arcsec, r_max_arcsec, 2048)
times_arr = np.linspace(0.001, 10, 100)
broadening = np.zeros_like(times_arr)
req = np.ones_like(times_arr)*0.00082
req2 = np.ones_like(times_arr)*0.00035

for i,time in enumerate(times_arr):
  broadening[i] = RdeltatInverse(theta_max_arcsec, r_hpbw_arcsec, time)

q = broadening
e = (1.e0 - q)/(1.e0 + q)

nom_spline = UnivariateSpline(times_arr, e, s=0)
ska_nom = nom_spline(0.1)*np.ones_like(times_arr)

plt.figure(1)
plt.loglog(times_arr, e)
plt.plot(times_arr, req, label='$\mathrm{SKA1 \, WL \, requirement}$')
plt.plot(times_arr, req2, label='$\mathrm{SKA2 \, WL \, requirement}$')
plt.plot(times_arr, ska_nom, label='$\mathrm{Correlator \, capability \, (0.1\,seconds)}$')
plt.xlabel('$\mathrm{Gridding \, integration \, time} \, \Delta t \, [\mathrm{seconds}]$')
plt.ylabel('$\mathrm{PSF \, ellipticity} \, e$')
plt.legend(loc='upper left', fontsize='small', frameon=False)
plt.gcf().set_size_inches(4.5, 3.75)
plt.savefig('time_smearing.png', bbox_inches='tight', dpi=160)