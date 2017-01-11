#!/usr/bin/env python

import math,sys

if len(sys.argv)<6:
    print "\nbaseline_sensitivity.py written by Enno Middelberg 2002"
    print "\nTask to compute the baseline sensitivity and SNR of two equal antennas."
    print "Provide Eta (0.5 for 1-bit sampling, 0.69 for two-bit sampling), SEFD [Jy],"
    print "bandwidth [MHz], integration time [s] and target source flux density [Jy]\n"
    sys.exit()

eta=float(sys.argv[1])
SEFD=float(sys.argv[2])
bandwidth=1e+06*float(sys.argv[3])
inttime=float(sys.argv[4])
flux=float(sys.argv[5])

print "\nEta              = %4.3f" % (eta)
print "SEFD             = %4.3f Jy" % (SEFD)
print "bandwidth        = %4.3f Hz" % (bandwidth)
print "integration time = %4.3f s" % (inttime)
print "source flux      = %4.3f Jy" % (flux)

sensitivity=(1/eta)*SEFD/(math.sqrt(2*bandwidth*inttime))
sensitivity_mJy=1000*sensitivity
SNR=flux/sensitivity

print "\nThe baselinie sensitivity is %4.3f mJy and you will make a %4.1f sigma detection\n" % (sensitivity_mJy, SNR)
