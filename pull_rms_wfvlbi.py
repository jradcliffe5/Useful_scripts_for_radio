import os,sys
from astropy.io import fits
import numpy as np
try:
    path = str(sys.argv[sys.argv.index('pull_rms_wfvlbi.py')+1])
    output_numpy = str(sys.argv[sys.argv.index('pull_rms_wfvlbi.py')+2])
except IndexError:
    print 'Usage pull_rms_wfvlbi.py <fitsfile_path> <output_numpy>'
RA= []
DEC = []
rms = []
os.system('rm *Py.fits')
os.system('rm *Py.fits')
f = []
print 'Making numpy array of rms values using all fits files in path: %s' % path
for file in os.listdir(path):
    if file.endswith('.fits'):
        hdu = fits.open(path+file)
        f.append(file)
        #print 360+hdu[0].header['CRVAL1']
        RA.append(360+float(hdu[0].header['CRVAL1']))
        DEC.append(float(hdu[0].header['CRVAL2']))
        image_data = hdu[0].data
        image_data = image_data[0,0,:,:]
        rms.append(np.sqrt(np.mean(np.square(image_data[150:1200,150:1200])))*1E6)
print 'Saving array in current directory as: %s.npy' % output_numpy
np.save('%s.npy' % output_numpy,[f,c.ra.degree,c.dec.degree,rms])
