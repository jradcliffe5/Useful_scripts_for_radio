#!/usr/bin/env python

# rotation_measure.py written by Enno Middelberg, Aug2003
# Task to calculate the rotation measure in ionized gas

import math, sys

if len(sys.argv)<4:
    print "\n Task to calculate rotation measure in ionized gas."
    print " On the command line, specify electron density in cm^-3,"
    print " B field in mG and path length in pc.\n"
    sys.exit()

n_e=float(sys.argv[1])
B=float(sys.argv[2])
L=float(sys.argv[3])

print "\n Electron density: %1.1f / cm^3" % n_e
print " B field:          %1.1e mG" % B
print " Path length:      %1.1f pc" % L

RM=812 * n_e * B * L

print " Rotation measure: %1.3f rad/m^2\n" % RM
