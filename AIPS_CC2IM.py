import sys,re,os
import numpy as np
try:
    image_file = sys.argv[sys.argv.index('AIPS_CC2IM.py')+1]
    cutoff = float(sys.argv[sys.argv.index('AIPS_CC2IM.py')+2])
    aipsuser = int(sys.argv[sys.argv.index('AIPS_CC2IM.py')+3])
except IndexError:
    print('Usage Parseltongue AIPS_CC2IM.py <image_file> <cutoff> <aipsuser#>')
try:
    from AIPS import AIPS, AIPSDisk
    from AIPSTask import AIPSTask, AIPSList
    from AIPSData import AIPSUVData, AIPSImage, AIPSCat
    from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
    AIPS.userno = aipsuser
except ImportError:
    print 'No AIPS/Parseltongue available'

from astropy.io import fits


fitld = AIPSTask('FITLD')
fitld.datain = 'PWD:%s' % image_file
fitld.outname = 'IM'
fitld.outclass = 'IM'
fitld.go()
image_data = AIPSImage('IM','IM',1,1)
cc2im = AIPSTask('CC2IM')
cc2im.indata = image_data
cc2im.go()
image_data.zap()
image_data = AIPSImage('IM','CC2IM',1,1)
fittp = AIPSTask('FITTP')
fittp.indata = image_data
fittp.dataout = 'PWD:%s_cc2im.fits' % image_file.split('.fits')[0]
fittp.go()
image_data.zap()

modelimage = '%s_cc2im.fits' % image_file.split('.fits')[0]
hdu = fits.open(modelimage)
imagedata = hdu['PRIMARY'].data
imagedata[imagedata<cutoff] = 0
header = hdu['PRIMARY'].header
hdu.writeto('%s_adjusted.fits' % modelimage.split('.fits')[0], overwrite=True)
hdu.close()
print('DONE.\nOutput image: %s_adjusted.fits' % modelimage.split('.fits')[0])
