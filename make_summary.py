import os, sys, re
import numpy as np

def make_summary(mset):
    ms.open(mset)
    ms.summary(verbose=True,listfile=mset.split('.ms')[0]+'.listobs')
    ms.close()
try:
    mset = str(sys.argv[sys.argv.index('make_summary.py')+1])
    make_summary(mset)
except IndexError:
    print 'Usage casa -c make_summary.py <measurement_set>'
