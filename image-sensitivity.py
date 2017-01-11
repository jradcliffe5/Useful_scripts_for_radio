#!/usr/bin/env python

import math,sys

if len(sys.argv)<7:
    print "\n baseline_sensitivity.py written by Enno Middelberg 2002"
    print "\n Task to compute the image sensitivity and SNR of a number of"
    print " equal antennas. Provide Eta (usually 0.5), SEFD [Jy], # of stations,"
    print " bandwidth [MHz], integration time [s] and target source flux density [Jy]\n"
    sys.exit()

eta=float(sys.argv[1])
SEFD=float(sys.argv[2])
nst=float(sys.argv[3])
bandwidth=1e+06*float(sys.argv[4])
inttime=float(sys.argv[5])
flux=float(sys.argv[6])

baseline_time=inttime*nst*(nst-1)/2

print "\nEta              = %4.3f" % (eta)
print "SEFD             = %4.3f Jy" % (SEFD)
print "# of stations    = %4.3f" % (nst)
print "bandwidth        = %4.3f Hz" % (bandwidth)
print "integration time = %4.3f s" % (inttime)
print "baseline time    = %4.3f s" % (baseline_time)
print "source flux      = %4.3f Jy" % (flux)

sensitivity=(1/eta)*SEFD/(math.sqrt(2*bandwidth*baseline_time))
sensitivity_mJy=1000*sensitivity
SNR=flux/sensitivity

print "\nThe image sensitivity is %4.3f mJy and you will make a %4.1f sigma detection\n" % (sensitivity_mJy, SNR)
