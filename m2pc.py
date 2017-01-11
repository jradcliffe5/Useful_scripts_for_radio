#!/usr/bin/python

# Program to transform meters into parsecs
# "m2pc.py 1 2 3" transforms all three numbers from 
# m to pc

import sys
import string

lines = sys.argv
del lines[0]

if len(lines)==0:
    print"\n m2pc.py written by Enno Middelberg 2001"
    print"\n Program to transform meters into parsecs"
    print" type 'm2pc.py x y z' to transform x,y,z from"
    print" m to pc"

for x in range(len(lines)):
    print lines[x],"m =",float(lines[x])/3.085677567E16,"pc"
