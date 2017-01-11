#!/usr/bin/python

# Program to calculate proper motions

# Command line args: distance_[pc] time_interval_[yr] ang_separation_[mas]

import sys
import math

if len(sys.argv)==1:
    print "\npropmotion.py written by Enno Middelberg 2002"
    print "\nEnter distance in Mpc, time in yr and angular separation in mas"
    print "to calculate proper motion in the source:\n"
    print "$> propmotion.py Mpc time angsep\n"
    sys.exit()

d=float(sys.argv[1])
t=float(sys.argv[2])
a=float(sys.argv[3])


print "Distance: "+`d`+" Mpc, time interval: "+`t`+" yr, angular separation: "+`a`+" mas\n"

# convert angle to radians
new_a=(a/3.6E6)*2*math.pi/360

#convert Distance to pc instead of Mpc
d=d*1000000

prop_d=d*math.tan(new_a)
prop_d_per_year=prop_d/t
beta=(prop_d_per_year*3.085677567E16/31536000)/299792458

print "1 mas = "+`math.tan(4.8481368111e-09)*d`+" pc "
print "1 pc  = "+`1/(math.tan(4.8481368111e-09)*d)`+" mas\n"
print `a`+" mas = "+`new_a`+" radians\n"
print "That's "+`prop_d`+" pc or "+`prop_d_per_year`+" pc / yr or "+`beta`+" c"
