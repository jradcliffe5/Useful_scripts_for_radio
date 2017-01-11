#!/usr/bin/python

# Program to transform parsecs into meters
# "pc2m.py 1 2 3" transforms all three numbers from 
# pc to m

import sys
import string

lines = sys.argv
del lines[0]

if len(lines)==0:
    print"\n pc2m.py written by Enno Middelberg 2001"
    print"\n Program to transform parsecs into meters"
    print" type 'pc2m.py x y z' to transform x,y,z from"
    print" pc to m\n"

for x in range(len(lines)):
    print lines[x],"parsec =",float(lines[x])*3.085677567E16,"m"
