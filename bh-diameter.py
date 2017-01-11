#!/usr/bin/env python

import sys, math

if len(sys.argv)<4:
    print "\nbh-diameter written by Enno Middelberg 2002"
    print "\nTask to read in a Distance in Mpc, black hole mass in 10^6 Msun"
    print "and beam size in mas to caluculate the beam size in units of R_s.\n"
    sys.exit()
    
distance=3.085677567e+22*float(sys.argv[1])
R_s=3e+9*float(sys.argv[2])
beam=float(sys.argv[3])*math.pi/(1000*60*60*180)

apparent_R_s=R_s/distance
nbeams=beam/apparent_R_s

print "\nDistance: \t %4.3e pc" % (1e6*float(sys.argv[1]))
print "BH mass: \t %4.3e Msun" % (1e6*float(sys.argv[2]))
print "beam size: \t %4.3f mas" % (float(sys.argv[3]))

print "\nApparent R_s: \t %4.3e rad" % (apparent_R_s)
print "beam size: \t %4.3e rad" % (beam)

print "\nApparent size of BH: "+`nbeams`+" beams.\n"

