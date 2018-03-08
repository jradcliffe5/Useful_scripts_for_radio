import os,sys,re
import numpy as np
try:
    measurement_set = str(sys.argv[sys.argv.index('restore_model.py')+1])
    nsources = int(sys.argv[sys.argv.index('restore_model.py')+2])
    nterms = int(sys.argv[sys.argv.index('restore_model.py')+3])
    image_fits = np.array(sys.argv[sys.argv.index('restore_model.py')+4:])
except IndexError or ValueError:
    print('Usage casa -c restore_model.py <measurement_set> <nsources> <nterms> <source1.tt0.fits> <source1.tt1.fits>  <source2.tt0.fits> etc.')


print('Checking existence of data sets and images\n')
if os.path.exists(measurement_set) == True:
    print 'Found %s' % measurement_set
else:
    print('%s not found... please run again\n' % i)
    sys.exit()
for i in sys.argv[sys.argv.index('restore_model.py')+4:]:
    if os.path.exists(i) == True:
        print 'Found %s' % i
    else:
        print('%s not found... please run again\n' % i)
        sys.exit()

fitsimages_for_restore = {}
print('%d source(s) with %d Taylor Term(s): %d images should be provided\n' % (nsources,nterms,nsources*nterms))
if len(image_fits) != nsources*nterms:
    print('# images provided not equal to # sources * number Taylor Terms')
    sys.exit()
else:
    print('All ok.. we may proceed\n')

print('Checking whether files are fits file via .fits extension\n')
for i in range(len(image_fits)):
    if image_fits[i].endswith('.fits'):
        print('Fitsfiles found... converting %s into CASA image\n' % image_fits[i])
        importfits(fitsimage=image_fits[i], imagename=image_fits[i].split('.fits')[0])
        image_fits[i] = image_fits[i].split('.fits')[0]

print('Complete\n Your images are %s' % ' '.join(image_fits))

print('Checking for old MODEL visibilities')
delmod(vis=measurement_set,scr=True)


image_fits = np.split(image_fits,nterms) ## split on nterms
print('Assuming Taylor term images are in order... nterms=%d' % len(image_fits))

print('Making new MODEL visibilities from TT models')
for i in range(len(image_fits)):
    print('Inputting source %d into %s' % (i+1,measurement_set))
    ft(vis=measurement_set, nterms=nterms, model=image_fits[i].tolist(), usescratch=True, incremental=True)

print('Models to restore are in the data\n')
print('Putting the models into the uvplane')
uvsub(vis=measurement_set,reverse=True)

print('NOTE THAT THE UVSUBBED DATA IS IN THE CORRECTED DATA COLUMN')
