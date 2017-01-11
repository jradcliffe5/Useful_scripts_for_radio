#!/usr/bin/python
# brightness_temperature.py
# Program to calculate the brightness temperature
# from flux density, frequency and beam size

# Uses T = h nu / k * 1/ln(1 + 2 h nu^3 omega/(B c^2))
# (directly derived from Planck's law)
#     where h = Plank's const
#             nu in Hz
#             k = Boltzmann's const
#             c = 299272000 m/s
#             omega = beam area in steradians
#             B = source brightness in W m^-2 Hz^-1 beam^-1


# set some constants
h =  6.626e-34     # /* Plank's const in Js */
k =  1.38e-23      # /* Boltzmanns const in J/K */
c =  299792000     # /* speed of light in m/s */


print"\n brightness_temperature.py written by Enno Middelberg 2003\n"

import math
import string
import sys
import fileinput

# Control whether to use file or command-line data
id=0
if "-file" in sys.argv:
    id=1

if len(sys.argv)<4 and id==0:
    print" Program to calculate the brightness temperature of a radio source.\n"
    print" Usage:\n brightness_temperature.py flux_density [mJy] frequency [GHz] beam size [mas]"
    print" Or:      brightness_temperature.py filename\n"
    sys.exit()


data=[]
if id==1:
    for x in fileinput.input(sys.argv[2]):
        list=string.split(x)
        print list
        data.append([float(list[0]), float(list[1]), float(list[2])])
else:
    data.append([float(sys.argv[1]), float(sys.argv[2]), float(sys.argv[3])])

for x in data:
    flx_mJy =float(x[0])
    freq_GHz=float(x[1])
    beam_mas=float(x[2])

    # Convert flx to W m^-2 Hz^-1
    flx=flx_mJy*1e-29

    # Convert freq to Hz
    freq=freq_GHz*1e9

    # Convert beam to sterad
    beam=0.001*beam_mas
    beam=math.pi*math.pow(math.pi*beam/(2*3600.0*180.0),2.0)

    T = (h*freq/k)/math.log(1 + 2*h*math.pow(freq,3.0)*beam/(flx*math.pow(c,2.0)))

    print "%f mJy  %f GHz  %f mas = %5.5e K" %(flx_mJy, freq_GHz, beam_mas, T)


