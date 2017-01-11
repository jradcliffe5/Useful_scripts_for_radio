#!/usr/bin/env python

# This is a program to calculate spectral indices

import sys, math

if len(sys.argv)<5:
    print "\n Program to calculate spectral indices. Usage:"
    print "\n spix.py S(low) S(high) nu(low) nu(high)\n"
    sys.exit()

S1 =float(sys.argv[1])
S2 =float(sys.argv[2])
nu1=float(sys.argv[3])
nu2=float(sys.argv[4])

si=math.log10(S1/S2)/math.log10(nu1/nu2)

print "\nnu1 = %4.3f" % nu1
print "nu2 = %4.3f" % nu2
print "S1  = %4.3f" % S1
print "S2  = %4.3f" % S2
print "\n Spectral index is %4.3f\n" % si

