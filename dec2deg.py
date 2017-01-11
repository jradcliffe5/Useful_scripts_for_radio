#!/usr/bin/python
# dec2deg
# Program to convert a Dec of the form hh mm ss.ssssssss
# into a float 

print" dec2deg Written by Enno Middelberg 2001\n"

import sys
import string

inp=sys.argv[0:]
del inp[0]
if len(inp)==0:
        print" Program to convert a declination of the form hh mm ss.ssssssss"
        print" into degrees of a circle.\n"
        print" Type 'dec2deg.py hh mm ss.ssssssss' or hh:mm:ss.ssssssss"
        print" to calculate a float.\n"
	sys.exit()

if len(inp)==3:
        dec=inp

if len(inp)==1:
        dec=string.split(inp[0], ":")

print dec[0]+" hours, "+dec[1]+" minutes, "+dec[2]+" seconds, ok?\n"
hh=abs(float(dec[0]))
mm=float(dec[1])/60
ss=float(dec[2])/3600
if float(dec[0]) < 0:
        print" That's -"+str(hh+mm+ss)+" degrees"
else:
        print" That's "+str(hh+mm+ss)+" degrees"
