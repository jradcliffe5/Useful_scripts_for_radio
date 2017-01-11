#!/usr/bin/python
#
# Written by Enno Middelberg on 28.6.2001
#
# reads data from a file of the form
# 20  000000000000000000000000111111111100000000000000
# 296 000000000000000000000000111111111100000000000000
# where the first number is the antenna SEFD and the 48
# digits represent the source uptime at the antenna 
# during 24 h
#
# The data file is the only argument:
#
# 'sensitivity.py datafile'

import string
import sys
#system efficiency factor:
EFF=0.5

# this line opens the file that is passed with the call
lines = open(sys.argv[1]).readlines()

SEFD=[]
up=[]
for x in lines:
    if x=="\012":
	break
    SEFD.append(float(string.split(x)[0]))
    up.append(list(string.split(x)[1]))

nstn=len(up)

print "\nNumber of stations: ",nstn,"\n"

print "SEFD"'\t'"uptime"
print "--------------------------------------------------------"
for i in range(len(SEFD)):
    print SEFD[i],'\t',
    z=""
    for j in up[i]:
	z=z+j
    print z

BW = input("Bandwidth / MHz: ")
BW = 1000000*BW
f  = input("Source duty cycle: ") 

# Calculate baseline sensitivity for 2 min fringe-fit

i  = 0
ds = 0
while i < nstn:
    j=i+1
    while j < nstn:
	ds=(1.0/EFF)*pow( SEFD[i] * SEFD[j] / (2.0 * BW * 120.0),0.5)
#	print "Sensitivity on baseline ",i,"-",j,": ",ds
	j=j+1
    i=i+1

# Calculate image sensitivity:

k  = 0
i  = 0
j  = 0
dS=[]
count = range(0,47)

while i < nstn:
    j=i+1
    while j < nstn:
	dT=0
	for x in count:
	    dT = dT + (eval(up[i][x]) * eval(up[j][x]) * 1800.0 * f)
#	print "Time on baseline ",i,"-",j,": ",dT
	dS.append((1.0 / EFF) * pow( SEFD[i] * SEFD[j] / (2.0 * BW * dT),0.5))
	j=j+1
    i=i+1

# print dS, len(dS)

# Calculate variance in image

dI = 0
C  = 0
Tk = 1 # no taper
Wk = 1 # natural weight

a=0
while a > -3:

    i=0
    while i < len(dS):
	wk = pow(dS[i], a)
	dI = dI + pow(Tk, 2.0) * pow(Wk, 2.0) * pow(wk, 2.0) * pow(dS[i], 2.0)
#	print dI
	C  = C + Tk * Wk * wk
	i=i+1

    C = 1.0 / (2.0 * C)
    dI = 2.0 * C * pow(dI, 0.5)
    print "dI(",a,") = ",dI, "Jy"
    a = a-2
