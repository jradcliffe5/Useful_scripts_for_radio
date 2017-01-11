#!/usr/bin/python
# avg.py
# Program to calculate the average and standard deviations
# of columns and lines in a rectangle of numbers

print "\n avg.py Written by Enno Middelberg 2002\n"
print " Calculates the average and standard errors of"
print " columns and lines in a rectangle of numbers\n"
print " avg.py [-cols] [-lines]\n"
print " sx     = sqrt( 1/(n-1) * sum ( x_i - x )^2 )"
print " sigmax = sqrt(  1/n  * sum ( x_i - x )^2 )\n"

import math
import string
import sys

# Function to work out avg, sx and sigmax from a list of numbers
def stat(numbers):
    sum=0
    # Calculate avg, sx and sigmax
    if len(numbers)>1:
	for i in numbers:
	    sum=sum+i
	    # Store total sum for later
	    total=sum
	avg=sum/len(numbers)
	sum=0
	for i in numbers:
	    sum=sum+(i-avg)**2
	sx=math.sqrt(sum/(len(numbers)-1))
	sigmax=math.sqrt(sum/len(numbers))
    else:
	avg=numbers[0]
	sx=0
	sigmax=0
	sum=numbers[0]
    # Calculate mean
    numbers.sort()
    if (len(numbers) % 2)==1:
	mean=numbers[(len(numbers)-1)/2]
    else:
	mean=(numbers[(len(numbers)/2)-1]+numbers[(len(numbers)/2)])/2
    # Return results
    return avg, sx, sigmax, total, mean

# Get numbers from stdin
list=[]
for line in sys.stdin.readlines():
    if line[0]!="#":
	list.append(string.split(line))

ncols=len(list[0])
nlines=len(list)

print "\nNumber of columns: %i, number of lines: %i" % (ncols, nlines)

# Sort columns into 1D-lists and pass them to stat()
if "-cols" in sys.argv:
    colstats=[]
    for i in range(ncols):
	colnumbers=[]
	for j in range(nlines):
	    try:
		colnumbers.append(float(list[j][i]))
	    except (ValueError):
		continue
	colstats.append(stat(colnumbers))

print "\n",

# Sort lines into 1D-lists and pass them to stat()
if "-lines" in sys.argv:
    linestats=[]
    for i in range(nlines):
	linenumbers=[]
	for j in range(ncols):
	    try:
		linenumbers.append(float(list[i][j]))
	    except (ValueError):
		continue
	linestats.append(stat(linenumbers))

# Print results
print "            ",
for i in range(ncols):
    print "%i          " %i,
print
for i in range(nlines):
    print "%i -     " %i,
    for j in range(ncols):
	try:
	    print "% 5.4e" %(float(list[i][j])),
	except (ValueError):
	    print "           ",
    print

if "-cols" in sys.argv:
    print "\nColumn statistics:"
    print "\nAvg:    ",
    for i in colstats:
	print "% 5.4e" %i[0],
    print "\nsx:     ",
    for i in colstats:
	print "% 5.4e" %i[1],
    print "\nsigmax: ",
    for i in colstats:
	print "% 5.4e" %i[2],
    print "\nSum:    ",
    for i in colstats:
	print "% 5.4e" %i[3],
    print "\nMean:   ",
    for i in colstats:
	print "% 5.4e" %i[4],


if "-lines" in sys.argv:
    print "\n\nLine statistics:    Avg         sx          sigmax      Sum         Mean\n"
    for i in range(nlines):
	print "        Line  %i - % 5.4e % 5.4e % 5.4e % 5.4e % 5.4e" % (i,linestats[i][0], linestats[i][1], linestats[i][2], linestats[i][3], linestats[i][4])


