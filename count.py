#!/usr/bin/python
# task to count from 1. nummber to 2. nummber with the specified increment

import math, sys, string

if len(sys.argv)==1:
    print"\n count.py written by Uwe Bach 2003"
    print"\n Task to count from 1st nummber to 2nd nummber with the specified increment"
    print"\n Usage: count.py [option] 1st# 2nd# increment"
    print" Options: c=output in a column, l=output in a line \n"
    sys.exit()

format=(sys.argv[1])
first=string.atof(sys.argv[2])
second=string.atof(sys.argv[3])
incr=string.atof(sys.argv[4])

if format == 'c':
	while first <= second:
		print "%7.2f \n" % first,
		first=first+incr

if format == 'l':
	while first <= second:
		print "%7.2f" % first,
		first=first+incr

