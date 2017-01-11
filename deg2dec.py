#!/usr/bin/python
#
# deg2dec 
# Program to convert a number into the format
# hh mm ss.sssss

print" deg2dec Written by Enno Middelberg 2001\n"

import sys
import string

inp=sys.argv[0:]
del inp[0]
if len(inp)==0:
	print" Program to convert a float into a declination"
	print" of the form hh:mm:ss.ssssssss\n"
	print" Type 'deg2dec.py' followed by a list of floats"
	print" to convert them into declinations."

for x in inp:
	deg=float(x)
	sign="+"

      # test whether the input numbers are sane:

      # if negative, store "-" in sign and continue calulation
      # with positive value

	if deg < 0:
		sign="-"
                deg=deg*(-1)

	if deg > 180:
		print `deg`+": inputs may not exceed 180!\n"
		continue

	if deg > 90:
		print `deg`+" exceeds 90, will convert it to negative dec\n"
		deg=deg-90
		sign="-"

	hh=int(deg)
	mm=int((deg-int(deg))*60)
	ss=((deg-int(deg))*60-mm)*60
	if sign=="-":
		print str(deg+90)+":"
	else:
		print `deg`+":"
	print '\t'+sign+string.zfill(`hh`,2)+':'+string.zfill(`mm`,2)+':'+'%10.8f' % ss
	print '\t'+sign+string.zfill(`hh`,2)+' '+string.zfill(`mm`,2)+' '+'%10.8f' % ss
	print '\t'+sign+string.zfill(`hh`,2)+'h'+string.zfill(`mm`,2)+'m'+'%10.8fs\n' % ss
