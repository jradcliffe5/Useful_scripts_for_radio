#!/usr/bin/python
# deg2ra
# Program to convert a number into the format
# hh mm ss.ssssssss

print"\n deg2ra Written by Enno Middelberg 2001\n"

import math
import sys
import string

inp=sys.argv[0:]
del inp[0]
if len(inp)==0:
        print" Program to convert a float into a right ascension"
        print" of the form hh:mm:ss.ssssssss\n"
        print" Type 'deg2ra.py' followed by a list of floats"
        print" to convert them into right ascensions."

for x in inp:
        deg=float(x)

      # test whether the input numbers are sane:

        if deg < 0:
                deg=deg+360

        if deg > 360:
                print `deg`+": inputs may not exceed 360!\n"
                continue

        hh=int(deg/15)
	mm=int((deg-15*hh)*4)
	ss=(4*deg-60*hh-mm)*60
	print `deg`+":"
        print '\t'+string.zfill(`hh`,2)+':'+string.zfill(`mm`,2)+':'+'%10.8f' % ss
        print '\t'+string.zfill(`hh`,2)+' '+string.zfill(`mm`,2)+' '+'%10.8f' % ss
        print '\t'+string.zfill(`hh`,2)+'h'+string.zfill(`mm`,2)+'m'+'%10.8fs\n' % ss
