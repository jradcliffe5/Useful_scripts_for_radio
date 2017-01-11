#!/homes/emiddelb/local2/bin/python

# Task to rotate a vector specified on command line n times through
# angle p and calculate the vector sum

import sys, math

x=float(sys.argv[1])
y=float(sys.argv[2])
n=int(sys.argv[3])
angle=float(sys.argv[4])

print "\nRotate vector (%f, %f) %i times by %f degrees " % (x, y, n, angle)

# From now on, angle is in rad
angle=math.pi*angle/180

length_0=math.sqrt(x**2+y**2)

print "Length of vector before rotation: %f" % length_0

x_sum=0
y_sum=0
for i in range(n):
    x_n=x*math.cos(angle)-y*math.sin(angle)
    y_n=x*math.sin(angle)+y*math.cos(angle)
    x=x_n
    y=y_n
    x_sum=x_sum+x
    y_sum=y_sum+y
    print "%f %f %f %f" % (x,y, x_sum, y_sum)
#    print "Vector is (%f, %f) after %i rotations." % (x,y,i+1)

length=math.sqrt(x_sum**2 + y_sum**2)
angle=180*math.acos(x_sum/length)/math.pi
if y < 0.:
    angle = 360-angle

print "Final vector length is %f in orientation %f degrees." % (length, angle)
