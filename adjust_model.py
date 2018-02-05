from astropy.io import fits
import numpy as np
import sys, os

def adjust_model(modelimage,mode):
    if mode == 'del_neg':
        print 'Removing negatives from %s' % modelimage
        hdu = fits.open(modelimage)
        imagedata = hdu['PRIMARY'].data
        imagedata[imagedata<0] = 0
        header = hdu['PRIMARY'].header
        hdu.writeto('%s_no_neg.fits' % modelimage.split('.fits')[0], overwrite=True)
        hdu.close()
        print('DONE.\nOutput image: %s_no_neg.fits' % modelimage.split('.fits')[0])

try:
    modelimage = str(sys.argv[sys.argv.index('adjust_model.py')+1])
    mode = str(sys.argv[sys.argv.index('adjust_model.py')+2])
    adjust_model(modelimage,mode)
except IndexError:
    print 'Usage python adjust_model.py <modelfitsfile> <mode>\nCurrent mode(s) are del_neg: delete negatives in model'
