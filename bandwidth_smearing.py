'''
Code to calculate beam ellipticity from bandwidth smearing.

- Finds distorted beam shape due to bandwidth smearing
- Fits a spline to beam, finds roots to get FWHM
- Assumes smearing only in tangential direction => calc ellipticity
- Compares to requirements for smearing to be negligible for SKA WL surveys
'''
import numpy as np
from scipy.special import erf
from scipy.interpolate import UnivariateSpline
from scipy.ndimage.filters import gaussian_filter1d
from matplotlib import pyplot as plt

def Beam_SquareBandpassGaussianTaper(r, r_hpbw):
  '''From Bridle, A. H., Schwab, F. R., 1999,
  in Synthesis Imaging in Radio Astronomy II
  pp. 375
  eqn. (18-22)
  '''
  gamma = 2.e0*np.sqrt(np.log(2))
  retVar = np.exp(-gamma*gamma*r*r/(r_hpbw*r_hpbw))

  return retVar

def DistortedBeam_SquareBandpassGaussianTaper(theta, r, r_hpbw,
                                              delta_nu, nu0):
  '''From Bridle, A. H., Schwab, F. R., 1999,
  in Synthesis Imaging in Radio Astronomy II
  pp. 375
  eqn. (18-23)
  '''
  alpha = r/r_hpbw
  beta = (delta_nu/nu0)*(theta/r_hpbw)
  print(delta_nu/1.e6, beta)
  gamma = 2.e0*np.sqrt(np.log(2))

  retVar = (np.sqrt(np.pi)/(2.e0*gamma*beta))*(erf(gamma*(alpha + beta/2.e0)) - 
                                               erf(gamma*(alpha - beta/2.e0)))
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
theta_max_arcmin = 30
theta_max_arcsec = theta_max_arcmin*60.
theta_max = theta_max_arcmin/pix_scale
r_max_arcsec = 3.e0
r_max_arcmin = r_max_arcsec*60.
r_max = r_max_arcmin/pix_scale
r_hpbw_arcsec = 0.5
r_hpbw_arcmin = r_hpbw_arcsec*60.
r_hpbw = r_hpbw_arcmin/pix_scale
delta_nu = 1.e6
nu0 = 950e6
tau = 1

theta = np.linspace(0, theta_max_arcsec, 2048) # from pointing centre
r = np.linspace(-r_max_arcsec, r_max_arcsec, 2048) # distance across beam

# calculate undistorted beam
y1 = Beam_SquareBandpassGaussianTaper(r, r_hpbw_arcsec)
# calculate distorted beam at different distances from pointing centre
y2 = DistortedBeam_SquareBandpassGaussianTaper(1.e-10,
                                               r, r_hpbw_arcsec,
                                               delta_nu, nu0)
y3 = DistortedBeam_SquareBandpassGaussianTaper(theta_max_arcsec/2.e0,
                                               r, r_hpbw_arcsec,
                                               delta_nu, nu0)
y4 = DistortedBeam_SquareBandpassGaussianTaper(theta_max_arcsec,
                                               r, r_hpbw_arcsec,
                                               delta_nu, nu0)

# splines. finding roots of these gives the fwhm of the smeared beam.
sp1 = UnivariateSpline(r, y1 - y1.max()/2, s=0)
sp2 = UnivariateSpline(r, y2, s=0)
sp3 = UnivariateSpline(r, y3, s=0)
sp4 = UnivariateSpline(r, y4, s=0)

print(sum(abs(sp1.roots()))) # should be 0.5

# plot these beams
plt.figure(1)
plt.clf()
plt.plot(r, y2, label='$\\theta = 0\',\, \\beta=0$')
plt.plot(r, y3, label='$\\theta = 10\',\, \\beta=0.8$')
plt.plot(r, y4, label='$\\theta = 20\',\, \\beta=1.78$')
plt.text(-2.5, 0.9, '$\Delta\\nu = 1\mathrm{MHz}$')
plt.xlabel('$r \, [\mathrm{arcsec}]$')
plt.ylabel('$B_{D}(r)$')
plt.legend(loc='upper right', fontsize='small', frameon=False)
plt.gcf().set_size_inches(4.5, 3.75)
plt.savefig('beam_shapes_demo.png', bbox_inches='tight', dpi=160)
'''
# array of different channel widths
delta_nu_arr = np.logspace(3, 6, 20)
distorted_fwhm = np.zeros_like(delta_nu_arr)

# requirements on beam smearing
req = np.ones_like(delta_nu_arr)*0.00082 # GREAT c for SKA1
req2 = np.ones_like(delta_nu_arr)*0.00035 # ...for SKA2
req3 = np.ones_like(delta_nu_arr)*0.0082 # ...for SuperCLASS
req4 = np.ones_like(delta_nu_arr)*0.004 # ...for VLASS-DEEP

# for each channel width, find new width of beam at edge of FoV
for i,nu in enumerate(delta_nu_arr):
  beam = DistortedBeam_SquareBandpassGaussianTaper(theta_max_arcsec,
                                                   r, r_hpbw_arcsec,
                                                   nu, nu0)
  spline = UnivariateSpline(r, beam - beam.max()/2, s=0)
  distorted_fwhm[i] = sum(abs(spline.roots()))

# calculate ellipticity of new beam
# nb assuming smearing in radial / tangential direction only
q = r_hpbw_arcsec/distorted_fwhm
e = (1.e0 - q)/(1.e0 + q)

# find the capability at full SKA correlator resolution
nom_spline = UnivariateSpline(delta_nu_arr, e, s=0)
ska_nom = nom_spline(300.e6/64.e3)*np.ones_like(delta_nu_arr)

# plot all of this info
plt.figure(2)
plt.loglog(delta_nu_arr/1.e6, e)
#plt.plot(delta_nu_arr/1.e6, req3,
#         label='$\mathrm{SuperCLASS \, WL \, requirement}$')
#plt.plot(delta_nu_arr/1.e6, req4,
#         label='$\mathrm{VLASS-DEEP \, WL \, requirement}$')
plt.plot(delta_nu_arr/1.e6, req,
         label='$\mathrm{SKA1 \, WL \, requirement}$')
plt.plot(delta_nu_arr/1.e6, req2,
         label='$\mathrm{SKA2 \, WL \, requirement}$')
plt.plot(delta_nu_arr/1.e6, ska_nom,
         label='$\mathrm{Correlator \, capability \, (64k \, channels)}$')
plt.xlabel('$\mathrm{Gridding \, frequency \, width} \, \Delta\\nu \, \
           [\mathrm{MHz}]$')
plt.ylabel('$\mathrm{PSF \, ellipticity} \, e$')
plt.legend(loc='upper left', fontsize='small', frameon=False)
plt.gcf().set_size_inches(4.5, 3.75)
plt.savefig('bandwidth_smearing.png', bbox_inches='tight', dpi=160)
'''