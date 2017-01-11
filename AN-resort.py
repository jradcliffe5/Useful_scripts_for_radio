#!/usr/bin/python

# task to resort AIPS AN tables into a format for plotting and averaging
# Enno Middelberg, January 2003

import string, fileinput, sys

if len(sys.argv)<4:
    print "\n Usage: AN-resort.py filename no_IFs no_antennas [-long]\n"
    sys.exit()

infile=sys.argv[1]
noifs=int(sys.argv[2])
noantennas=int(sys.argv[3])

list=[]
dowrite=0
for x in fileinput.input(infile):
    parts=string.split(x)
    if "***END*PASS***" in parts:
        dowrite=0
    if dowrite==1:
        list.append(x)
    if "***BEGIN*PASS***" in parts:
        dowrite=1

if "-long" in sys.argv:
    for x in range(noifs):
        print "\t\t\t  IF "+`x+1`+"\t\t\t ",
    print
    for x in range(noifs):
        print "\t  Rrrrriiiiggghhhtt    Lllleeeffffttt  ",
    print "\n\t",
    print
    for x in range(noifs):
        print "\t   Re        Im        Re        Im\t",
else:
    print "\t         Rrrrriiiiggghhhtt    Lllleeeffffttt  "
    print "\t            Re        Im        Re        Im"   

print

dterms=[]
for x in range(0, noantennas*noifs*2, noifs*2):
    antenna_name=list[x][10:12]
    if "-long" in sys.argv:
        print antenna_name,
    dterms.append([antenna_name])
    for y in range(0, noifs*2, 2):
        real_R=float(list[x+y][130:143])
        real_L=float(list[x+y][168:182])
        imag_R=float(list[x+y+1][130:143])
        imag_L=float(list[x+y+1][168:182])
        dterms[x/(noifs*2)].append([real_R, imag_R, real_L, imag_L])
        if "-long" in sys.argv:
            print "\t% 7.6f % 7.6f % 7.6f % 7.6f   " % (real_R, imag_R, real_L, imag_L),
        else:
            print antenna_name+"\t IF #%i  % 7.6f % 7.6f % 7.6f % 7.6f" % (y/2, real_R, imag_R, real_L, imag_L)
    print
