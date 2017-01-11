from recipes import makepb
import sys, os

msfile = str(sys.argv[3])
imagetemplate = str(sys.argv[4])

print sys.argv
print imagetemplate[-4:]

if imagetemplate[-1] == '/':
	imagetemplate = imagetemplate[:-1]

if imagetemplate[-4:]== 'fits' or imagetemplate[-4:]=='FITS':
	importfits(fitsimage=imagetemplate, imagename=imagetemplate[:-4]+'image')
	imagetemplate = imagetemplate[:-4]+'image'
	makepb.makePB(vis=msfile,field='0',imtemplate=imagetemplate,outimage=msfile[:-2]+'pb',pblimit=0.1)
	impbcor(imagename=imagetemplate, pbimage=msfile[:-2]+'pb', outfile=imagetemplate[:-6]+'_pbcor.image')
	exportfits(imagename=imagetemplate[:-6]+'_pbcor.image', fitsimage=imagetemplate[:-6]+'_pbcor.fits')
	exportfits(imagename=msfile[:-2]+'pb', fitsimage=msfile[:-2]+'pb.fits')
else:
	makepb.makePB(vis=msfile,field='0',imtemplate=imagetemplate,outimage=msfile[:-2]+'pb',pblimit=0.1)
	impbcor(imagename=imagetemplate, pbimage=msfile[:-2]+'pb', outfile=imagetemplate[:-6]+'_pbcor.image')
	exportfits(imagename=imagetemplate[:-6]+'_pbcor.image', fitsimage=imagetemplate[:-6]+'_pbcor.fits')
	exportfits(imagename=msfile[:-2]+'pb', fitsimage=msfile[:-2]+'pb.fits')
