#!/usr/bin/env python

# rad2deg.py, task to convert an angle in radians
# into degrees

import sys, math, string

if len(sys.argv)==1:
    print "\n rad2deg.py written by Enno Middelberg 2002"
    print "\n Task to convert an angle in radians"
    print " into degrees. Type rad2deg.py followed by"
    print " a blank-separated list of numbers you want"
    print " to convert.\n"
    sys.exit()

for x in sys.argv[1:]:
    rad=float(x)
    deg=(180/math.pi)*rad
    hh=int(deg)
    mm=int((deg-int(deg))*60)
    ss=((deg-int(deg))*60-mm)*60
    print " %1.3e rad = %3.3e degrees" % (rad, deg),
    print "= "+string.zfill(`hh`,2)+':'+string.zfill(`mm`,2)+':'+'%10.8f' % ss
