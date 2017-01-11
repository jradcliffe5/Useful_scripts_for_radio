#!/usr/bin/env python

# tau_ff.py written by Enno Middelberg, Aug2003
# Task to calculate the free-free optical depth of ionized gas

import math, sys

if len(sys.argv)<4:
    print "\n Task to calculate free-free optical depth of ionized gas."
    print " Uses Eq. 4.32 from Osterbrock (1989). On the command line,"
    print " specify temperature in K, frequency in GHz, electron "
    print " density in cm^-3 and line of sight in pc\n"
    sys.exit()

#Osterbrock equation:
#tau_ff = 8.24E-2 * T^(-1.35 ) * nu^(-2.1) * n_e^2

T  =float(sys.argv[1])
nu =float(sys.argv[2])
n_e=float(sys.argv[3])
l  =float(sys.argv[4])

print "\n Temperature:        %1.1f K" % T
print " Frequency:          %1.1f GHz" % nu
print " Electron density:   %1.1f / cm^3" % n_e
print " Line of sight:      %1.1f pc" % l

tau_ff=8.24e-2 * T**(-1.35 ) * nu**(-2.1) * n_e**2 * l

print " Optical depth:      %1.3e\n" % tau_ff
