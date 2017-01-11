#!/usr/bin/python

# task to convert Difmap win files into AIPS-readable BOXFILE format
# written by Enno Middelberg 2002

import sys, string, fileinput, math

if len(sys.argv)==1:
    print"\n difwin2aipswin.py written by Enno Middelberg 2002"
    print"\n Task to convert Difmap .win files into AIPS-readable BOXFILE format"
    print" Usage: difwin2aipswin.py difmapwinfile difmapcellsize aipsimsize\n"
    sys.exit()

c=string.atof(sys.argv[2])
i=string.atof(sys.argv[3])

print "Found Difmap cellsize of "+`c`+" mas"
print "Found AIPS image size of "+`i`+" pixels"


c=1.0/c
i=i/2.0

# browse list, sort out comments and print new list
for line in fileinput.input(sys.argv[1]):
    if string.find(line, "!")<>0:
	newline=string.split(line)
	blcx=i-c*string.atof(newline[1])
	blcy=i+c*string.atof(newline[2])
	trcx=i-c*string.atof(newline[0])
	trcy=i+c*string.atof(newline[3])
	print "1 %4d %4d %4d %4d" % (blcx, blcy, trcx, trcy)

