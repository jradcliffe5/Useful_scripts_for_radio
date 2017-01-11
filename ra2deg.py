#!/usr/bin/python
# ra2deg
# Program to convert an RA of the form hh mm ss.ssssssss
# into a float 

print" ra2deg Written by Enno Middelberg 2001\n"

import sys
import string

inp=sys.argv[0:]
del inp[0]
if len(inp)==0:
        print" Program to convert an RA of the form hh mm ss.ssssssss"
        print" into degrees of a circle.\n"
        print" Type 'ra2deg.py hh mm ss.ssssssss' or hh:mm:ss.ssssssss"
        print" to calculate a float."

if len(inp)==3:
        ra=inp

if len(inp)==1:
        ra=string.split(inp[0], ":")

print ra[0]+" hours, "+ra[1]+" minutes, "+ra[2]+" seconds, ok?"
hh=float(ra[0])*15
mm=(float(ra[1])/60)*15
ss=(float(ra[2])/3600)*15
print" That's "+str(hh+mm+ss)+" degrees"
