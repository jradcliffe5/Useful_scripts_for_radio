#!/usr/bin/python
# task to convert a complex number into radius/phi

print"\n cplx2rad.py Written by Enno Middelberg 2002\n"

import math
import sys

re=float(sys.argv[1])
im=float(sys.argv[2])

r=math.sqrt(re**2+im**2)
phi=(180/math.pi)*math.atan(im/re)

if im>0 and re<0:
    phi=phi+180

if im<0 and re<0:
    phi=phi-180

print r, phi
