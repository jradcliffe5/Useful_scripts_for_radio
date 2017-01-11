#!/usr/bin/env python

# deg2rad.py, task to convert an angle in degrees
# into radians

import sys, math, string

if len(sys.argv)==1:
    print "\n deg2rad.py written by Enno Middelberg 2002"
    print "\n Task to convert an angle in degrees"
    print " into radians. Type deg2rad.py followed by"
    print " a blank-separated list of numbers you want"
    print " to convert. Format may be dd:mm:ss.sssss"
    print " or a plain number.\n"
    sys.exit()

for x in sys.argv[1:]:
    if ":" in x:
	list=string.split(x, ":")
	if len(list)<3:
	    print "Too few numbers in data format."
	    sys.exit()
	deg=float(list[0])+(1.0/60.0)*float(list[1])+(1.0/3600.0)*float(list[2])
    else:
	deg=float(x)
    rad=deg*2*math.pi/360
    print x+" = %1.3e degrees = %1.3e radians" % (deg, rad)
