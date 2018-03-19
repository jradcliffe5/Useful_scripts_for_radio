from astropy.io import fits
from astropy.wcs import WCS
import os
from astropy.coordinates import SkyCoord, Angle
import numpy as np
import matplotlib
from numpy import linspace, meshgrid
from matplotlib.mlab import griddata
import scipy.ndimage as ndimage
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.interpolate import interp2d
from scipy.ndimage.filters import gaussian_filter
import matplotlib.pyplot as plt
from regions import CircleSkyRegion
from matplotlib import *
from matplotlib.patches import Rectangle
from astropy.visualization.wcsaxes import SphericalCircle
import astropy.units as u
rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
from matplotlib import rcParams
rcParams['mathtext.default'] = 'regular'
fig_size = plt.rcParams["figure.figsize"]
# Prints: [8.0, 6.0]
print "Current size:", fig_size
matplotlib.rcParams.update({'font.size': 22})
# Set figure width to 9 and height to 9
fig_size[1] = 9
fig_size[0] = 9
plt.rcParams["figure.figsize"] = fig_size
plt.ioff()
plt.rcParams['axes.axisbelow'] =True
markers = ['v','.','p','*','+']
size =[80,150,70,100,100]
matplotlib.rcParams.update({'font.size':22})
def convertAIPStoPythonImage(filename,outfilename):
    hdu_list = fits.open(filename)

    head = hdu_list['PRIMARY'].header
    head['CRVAL1'] = 360 - (head['CRVAL1']*-1)
    head['NAXIS'] = 2
    del head['NAXIS3']
    del head['NAXIS4']
    del head['CTYPE3'],  head['CRVAL3'], head['CDELT3'], head['CRPIX3'], head['CROTA3']
    del head['CTYPE4'], head['CRVAL4'], head['CDELT4'], head['CRPIX4'], head['CROTA4']
    #if filename.endswith('NA.fits') == False:
    #head.set('BMIN', float(hdu_list[1].data.field('BMIN')))
    #head.set('BMAJ', float(hdu_list[1].data.field('BMAJ')))
    #head.set('BPA', float(hdu_list[1].data.field('BPA')))
    image_data = hdu_list[0].data
    image_data = image_data[0,0,:,:]
    hdu = fits.PrimaryHDU(image_data, header=head)
    hdulist = fits.HDUList([hdu])
    hdulist.writeto(outfilename,overwrite=True)
    return outfilename

def grid(x, y, z, resX=1000, resY=1000):
    "Convert 3 column data to matplotlib grid"
    xi = linspace(min(x), max(x), resX)
    yi = linspace(min(y), max(y), resY)
    Z = griddata(x, y, z, xi, yi,interp='linear')
    Z2 = ndimage.gaussian_filter(Z, sigma=1.0, order=1)
    #Z = interp2d(x, y, z, kind='cubic')
    #Z2 = Z(xi,yi)
    X, Y = meshgrid(xi, yi)
    return X, Y, Z

## Create array of ra, dec, & intensity
RA_c=[]
DEC_c = []
RA_a=[]
RA_s=[]
RA_l=[]
RA_b= []
DEC_a=[]
DEC_l=[]
DEC_b=[]
DEC_s=[]
rms = []
os.system('rm *Py.fits')
print len(os.listdir('./'))

for file in os.listdir('./'):
    if file.endswith('.fits'):
        hdu = fits.open(file)
        if file.startswith('HDFC'):
            print file
            RA_c.append(360+float(hdu[0].header['CRVAL1']))
            DEC_c.append(float(hdu[0].header['CRVAL2']))
        if file.startswith('HDFA'):
            RA_a.append(360+float(hdu[0].header['CRVAL1']))
            DEC_a.append(float(hdu[0].header['CRVAL2']))
        if file.startswith('HDFS'):
            RA_s.append(360+float(hdu[0].header['CRVAL1']))
            DEC_s.append(float(hdu[0].header['CRVAL2']))
        if file.startswith('HDFL'):
            RA_l.append(360+float(hdu[0].header['CRVAL1']))
            DEC_l.append(float(hdu[0].header['CRVAL2']))
        if file.startswith('HDFB'):
            RA_b.append(360+float(hdu[0].header['CRVAL1']))
            DEC_b.append(float(hdu[0].header['CRVAL2']))
            print file
        #print 360+hdu[0].header['CRVAL1']
        #RA.append(360+float(hdu[0].header['CRVAL1']))
        #DEC.append(float(hdu[0].header['CRVAL2']))
        #image_data = hdu[0].data
        #image_data = image_data[0,0,:,:]
        #rms.append(np.std(image_data[154:518,700:900]))
print len(DEC_c)
#print RA, DEC, rms
#print rms
#print np.nonzero(rms > 10E-6)
#print rms[np.argmax(rms>100E-6)]
#print np.mean(rms), np.std(rms)
c = SkyCoord(RA_c,DEC_c,unit='deg',frame='icrs')
a = SkyCoord(RA_a,DEC_a,unit='deg',frame='icrs')
s = SkyCoord(RA_s,DEC_s,unit='deg',frame='icrs')
l = SkyCoord(RA_l,DEC_l,unit='deg',frame='icrs')
b = SkyCoord(RA_b,DEC_b,unit='deg',frame='icrs')

### Create axes & scatter plot
filename = 'HDFC0155_NA_PBCOR_IM.fits' ## input fits file for correct wcs definition
outfilename = 'HDFC0155_NA_PBCOR_IMPy.fits'
convertAIPStoPythonImage(filename,outfilename)
hdu = fits.open(outfilename)
wcs = WCS(hdu[0].header)
fig = plt.figure(1)
ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)
#ax1 = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)
lon = ax.coords['ra']
lat = ax.coords['dec']
#lon.grid(color='k',alpha=0.5,linestyle='dashed')
#lat.grid(color='k',alpha=0.5,linestyle='dashed')
lon.set_major_formatter('hh:mm')
lat.set_major_formatter('dd:mm')
lon.set_axislabel('Right Ascension (J2000)', minpad=1.5)
lat.set_axislabel('Declination (J2000)', minpad=1)
ax.set_xlim(-1100000-50000, 1100000-50000)
ax.set_ylim(-1100000, 1100000)
#ax1.set_xlim(-500000, 500000)HDFC0155_PBCOR_IM.fits
#ax1.set_ylim(-500000, 500000)X, Y, Z = grid(c.ra.degree, c.dec.degree, rms)'magma',vmin=np.amin(rms), vmax=50E-6, alpha=1)
submm = ax.scatter(s.ra.degree, s.dec.degree, transform=ax.get_transform('icrs'),s=size[0], marker=markers[0],color='y',alpha=1,label=r'Sub-mm')
eMERGE=ax.scatter(c.ra.degree, c.dec.degree, transform=ax.get_transform('icrs'),s=size[1], marker=markers[1],color='r',alpha=1,label=r'eMERGE')
annulus=ax.scatter(a.ra.degree, a.dec.degree, transform=ax.get_transform('icrs'),s=size[2], marker=markers[2],color='b',alpha=1,label=r'Annulus')
legacy=ax.scatter(l.ra.degree, l.dec.degree, transform=ax.get_transform('icrs'),s=size[3], marker=markers[3],color='g',alpha=1,label=r'Legacy')
bright=ax.scatter(b.ra.degree, b.dec.degree, transform=ax.get_transform('icrs'),s=size[4], marker=markers[4],color='k',alpha=1,label=r'Bright')
ax.grid(linestyle='dashed',alpha=0.8)

leg = ax.legend((eMERGE,submm,legacy,annulus,bright),('eMERGE','Sub-mm','Legacy','Annulus','Bright'),loc=3,prop={'size':16})
leg.get_frame().set_facecolor('white')
fig.savefig('positions.pdf',bbox_inches='tight',dpi=fig.dpi,format="pdf")

#ax.legend_.remove()
ax.grid(linestyle='dashed',alpha=0.8)
#print np.amin(rms), np.amax(rms)
fig.savefig('positions_grid.pdf',bbox_inches='tight',dpi=fig.dpi,format="pdf")
#plt.show()
plt.clf()
### Create axes & scatter plot
filename = 'HDFC0155_NA_PBCOR_IM.fits' ## input fits file for correct wcs definition
outfilename = 'HDFC0155_NA_PBCOR_IMPy.fits'
convertAIPStoPythonImage(filename,outfilename)
hdu = fits.open(outfilename)
wcs = WCS(hdu[0].header)
fig = plt.figure(1)
ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)
#ax1 = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)
lon = ax.coords['ra']
lat = ax.coords['dec']
#lon.grid(color='k',alpha=0.5,linestyle='dashed')
#lat.grid(color='k',alpha=0.5,linestyle='dashed')
lon.set_major_formatter('hh:mm')
lat.set_major_formatter('dd:mm')
lon.set_axislabel('Right Ascension (J2000)', minpad=1.5)
lat.set_axislabel('Declination (J2000)', minpad=1)
ax.set_xlim(-540000, 480000)
ax.set_ylim(-540000, 480000)
#ax1.set_xlim(-500000, 500000)HDFC0155_PBCOR_IM.fits
#ax1.set_ylim(-500000, 500000)X, Y, Z = grid(c.ra.degree, c.dec.degree, rms)'magma',vmin=np.amin(rms), vmax=50E-6, alpha=1)
s_mm_range = ax.scatter(s.ra.degree, s.dec.degree, transform=ax.get_transform('icrs'), s=1000,color='y',alpha=1,label=r'Sub-mm')
s_mm_pos = ax.scatter(s.ra.degree, s.dec.degree, transform=ax.get_transform('icrs'), marker='+',color='k',alpha=1,label=r'Sub-mm')
eMERGE_range = ax.scatter(c.ra.degree, c.dec.degree, transform=ax.get_transform('icrs'), s=1000,color='r',alpha=1,label=r'eMERGE')
eMERGE_pos = ax.scatter(c.ra.degree, c.dec.degree, transform=ax.get_transform('icrs'), marker='+',color='k',alpha=1,label=r'eMERGE')
legacy_range = ax.scatter(l.ra.degree, l.dec.degree, transform=ax.get_transform('icrs'), s=1000,color='g',alpha=1,label=r'Legacy')
legacy_pos = ax.scatter(l.ra.degree, l.dec.degree, transform=ax.get_transform('icrs'), marker='+',color='k',alpha=1,label=r'Legacy')
ax.set_axisbelow(True)
ax.xaxis.grid(color='k',transform=ax.get_transform('icrs'))
leg = ax.legend(((eMERGE_range,eMERGE_pos),(s_mm_range,s_mm_pos),(legacy_range,legacy_pos)),('eMERGE','Sub-mm','Legacy'),labelspacing=1.2,borderpad=0.7,loc=4,prop={'size':16})
leg.get_frame().set_facecolor('white')

#print np.amin(rms), np.amax(rms)
fig.savefig('positions_SFXC.pdf',bbox_inches='tight',dpi=fig.dpi,format="pdf")
#plt.show()
'''
c = SkyCoord(RA,DEC,unit='deg',frame='icrs')
hdu = fits.open(outfilename)
wcs = WCS(hdu[0].header)
rms = 1/(np.array(rms)/5.5e-06)
print np.amin(rms), np.amax(rms)
fig = plt.figure()
ax = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)
#ax1 = fig.add_axes([0.15, 0.1, 0.8, 0.8], projection=wcs)
lon = ax.coords['ra']
lat = ax.coords['dec']
lon.set_major_formatter('hh:mm:ss')
lat.set_major_formatter('dd:mm:ss')
lon.set_axislabel('Right Ascension (J2000)', minpad=1.5)
lat.set_axislabel('Declination (J2000)', minpad=1)
ax.set_xlim(-1100000, 1100000)
ax.set_ylim(-1100000, 1100000)
#ax1.set_xlim(-500000, 500000)HDFC0155_PBCOR_IM.fits
#ax1.set_ylim(-500000, 500000)
X, Y, Z = grid(c.ra.degree, c.dec.degree, rms)
im = ax.pcolormesh(X,Y,Z,transform=ax.get_transform('icrs'),cmap='magma',vmin=0.3, vmax=1)
ax.scatter(c.ra.degree, c.dec.degree, transform=ax.get_transform('icrs'), marker=markers,color='g',norm=matplotlib.colors.LogNorm())
#CS = ax.contour(X,Y,gaussian_filter(Z,10),levels=np.linspace(0,1,8),colors='w',transform=ax.get_transform('icrs'), interpolation='none')
divider = make_axes_locatable(ax)
cax = divider.append_axes("top", size="5%", pad=0.00,axes_class=matplotlib.axes.Axes)
cb = plt.colorbar(orientation="horizontal",mappable=im, cax=cax,format='%.1e',ticks=np.linspace(0.5,1,8),extend='min')
#cb.add_lines(CS)
cb.ax.xaxis.set_ticks_position('top')
cax.set_xlabel("1$\sigma$ r.m.s. sensitivity ($\mu$Jy/bm)", labelpad=-60)

#fig.savefig('rms_taper_cal_weights_PB.pdf',bbox_inches='tight',dpi=50,format="pdf")
plt.show()
'''
#above code tranforms AIPS into python and then sets out an axis based upon the wcs coords of the fits file
# need to get a fits file that is in the center of the field to replace HDFA....fits
# also need to set up an array of coordinates & rms of the fitsfiles to feed into an
# ax.scatter and set up the colour scale to reflect r.m.s.
# Also need to interpolate look at http://cgcooke.github.io/Scattered-Interpolation/
