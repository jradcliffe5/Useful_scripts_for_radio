from astropy.io import fits
import numpy as np
import sys, os

def adjust_model(modelimage,thresh):
    print 'Removing negatives from %s' % modelimage
    hdu = fits.open(modelimage)
    imagedata = hdu['PRIMARY'].data
    imagedata[imagedata<thresh] = 0
    header = hdu['PRIMARY'].header
    hdu.writeto('%s_modified.fits' % modelimage.split('.fits')[0], overwrite=True)
    hdu.close()
    print('DONE.\nOutput image: %s_modified.fits' % modelimage.split('.fits')[0])

try:
    modelimage = str(sys.argv[sys.argv.index('adjust_model.py')+1])
    mode = str(sys.argv[sys.argv.index('adjust_model.py')+2])
    thresh = 0.0
    if mode == 'thresh':
        try:
            thresh = float(sys.argv[sys.argv.index('adjust_model.py')+3])
        except IndexError:
            print 'Usage of thresh needs to have value and final cmd line argument'
    adjust_model(modelimage,thresh)
except IndexError:
    print 'Usage python adjust_model.py <modelfitsfile> <mode> <thresh>\nCurrent mode(s) are del_neg/thresh: delete negatives in model/threshold'
