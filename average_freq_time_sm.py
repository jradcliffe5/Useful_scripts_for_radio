##############Packages##############################################################################################################
import numpy as np
from scipy.special import erf
from scipy.constants import c as clight
from matplotlib import pyplot as plt
import os

##########Import casa packages######################################################################################################
try:
    import casac
except ImportError, e:
    print "failed to load casa:\n", e
    exit()

#mstool = casac.homefinder.find_home_by_name('msHome')
#ms = casac.ms = mstool.create()
#tbtool = casac.homefinder.find_home_by_name('tableHome')
#tb = casac.tb #= tbtool.create()

# It's a good habit to create a new instance of the tool first, to
# avoid possible collisions:
#############################################

####################################################################################################################################
###########Define_Functions#################
############################################

###########Maximum baseline#################

def maxbaseline(msName):
    maxbaseline = 0
    tb.open(msName + '/ANTENNA')
    antab = np.ndarray.transpose(tb.getcol('POSITION'))
    antab=antab[~np.all(antab == 0, axis=1)] #removes any zero rows..eMERLIN!
    tb.close()
    for row in antab :
        for row2 in antab:
            xsep = row[0] - row2[0]
            ysep = row[1] - row2[1]
            zsep = row[2] - row2[2]
            hypxy =  np.sqrt((xsep * xsep) + (ysep * ysep))
            hypxyz = np.sqrt((zsep * zsep) + (hypxy * hypxy))
            if hypxyz > maxbaseline :
                maxbaseline = hypxyz
    return maxbaseline
############################################

###########Bandwidth Smearing###############
################################################################
# Based on Bridle & Schwab (B&S) ch 18 in Synthesis Imaging II
# (Taylor, Carilli & Perley 1999)
# Derivation information from A. Richards. 
######################################################
# 1.3 Gaussian distribution of uv coverage, square bandpass
# (denoted where required using _GS)
#     This is probably closest to a good range of hour angle coverage
#     and more than a few channels averaged e.g. an entire spw
# B&S Eq (18-24) gives
# R = sqrt(math.pi)/(1.665*beta) * erf(beta*1.665/2.)
# for observing frequency nu 
#
# This can be inverted using the first few terms of the Maclaurin
# series for erf (http://mathworld.wolfram.com/Erf.html)
#
# (2) R = (sqrt(pi)/(1.665*beta))*((2./sqrt(pi))*(z-z**3/3.+z**5/10.))
# where z = beta1*1.665/2.
# This reduces to a quadratic equation in z^2 which has real roots for
# R > 13/18 and beta is accurate to a few percent or better for R > 0.8
# 
# (3) beta_GS = (2./1.665)*sqrt((10.- sqrt(360.*R - 260.))/6.)
# Comparing (1) and (3) gives
#
# (4) dnu_GS = beta_GS/(1.E+06 * theta*pi/(180.*3600.) * B/3.E+08)
#
#######################################################
def BWSmear_FoV(max_baseline,cfreq,smearing,theta):# square bandpass,circular Gaussian tapering
    cfreq = cfreq/1.e6
    HPBW = 1.22*(clight/((cfreq*1.e6)*max_baseline))*(180./np.pi)*3600.
    beta_GS = (2./1.665)*np.sqrt((10.- np.sqrt(360.*smearing - 260.))/6.)
    #    dnu_GS = beta_GS/(1.E+06 * theta*math.pi/(180.*3600.) * B/3.E+08)
    dnu_GS = beta_GS * cfreq * (HPBW/theta)
    dnu_GS = dnu_GS
    return dnu_GS

###################Time Smearing##############
#################################################################
# Calculate time smearing
#
# Other effects e.g. primary beam, phase rate, may provide worse
# distortions at theta.
# Primary beam roughly based on 6-antenna array at L & C band (no Lo)
# and 5 antennas at K-band.
# Confusing as well as wanted sources should be considered
#
# Based on Bridle & Schwab (B&S) ch 18 in Synthesis Imaging II
# (Taylor, Carilli & Perley 1999)
# Average intensity loss for a circumpolar point source is given by
# R = 1-const (theta/thetaB)^2 dt^2
#
# where const = 1.08E-09 for uniform, circular uv coverage (possibly
# suitable for the most compact configurations) and
# const = 1.22E-09 for a Gaussian distribution of uv coverage (more
# suitable for most, especially extended configurations).
# So
# dt = sqrt[(1-R)/const] (thetaB/theta) 
###############################################################
def TimeSmear_FoV(max_baseline,cfreq,smearing,theta):
    HPBW = 1.22*(clight/((cfreq)*max_baseline))*(180./np.pi)*3600.
    #t_a = np.sqrt((1.-smearing)/1.22e-9)*(HPBW/theta)
    t_a = (np.sqrt((1.-smearing)/1.22E-09))*HPBW/theta
    return t_a

####################################################################################################################################


#############Inputs#########################

# Inputs to BW & T smearing calculations..
smearingpc = 10. #percentage smearing required
theta = 60. #radial distance to get desired smearing
msName = 'XMAS_2013_E-MERGE.ms' # name of the ms file to be averaged!
############################################


############################################
###########Change inputs####################
############################################
#some functions to get smearing variables into a more computer friendly way.
smearing = (100.-smearingpc)/100.
############################################



# You need to define msName='myMS.ms' first -- then, this will
# open the table with spectral window information
tb.open(msName + '/SPECTRAL_WINDOW')

# Read the REF_FREQUENCY and NUM_CHAN and TOTAL_BANDWIDTH into np arrays
refFreq = tb.getcol('REF_FREQUENCY')
nChan = tb.getcol('NUM_CHAN')
totBW = tb.getcol('TOTAL_BANDWIDTH')
allchanwidths = tb.getcol('CHAN_WIDTH')

# Close the table
tb.close()

# The SPW index is simply the row number in SPECTRAL_WINDOW
spwObs = len(refFreq)
# and our channel width is just the total bandwidth/no channels
widthChan = totBW/nChan
cfreq = refFreq[0]+ np.sum(totBW)/2.
print cfreq

#Check that all channels & spws are the same size.

def checkEqual2(iterator):
    return len(set(iterator)) <= 1

print 'The frequency information for BW averaging of data set %s is:' % msName
print 'No. Spectral Windows:  %d' % spwObs
print 'spw same size?         %s' % checkEqual2(totBW)
print 'Bandwidth per spw:     %f Hz' % totBW[0]
print 'No. Channels per spw:  %d' % nChan[0]
print 'and channel same size? %s' % checkEqual2(allchanwidths[0])
print 'Channel width:         %f Hz \n' % widthChan[0]

#check that channels are the same size and spectral windows.
if checkEqual2(totBW)==False or checkEqual2(allchanwidths[0])==False:
    print 'Your channels or spectral windows are not the same bandwidths'
    print 'Check you don\'t have mixed mode data.. split out the continuum please'
    print 'And run again!\n'
    exit()

del allchanwidths #delete as this array can be large to stick in memory

print 'Great! Lets continue, to average your data\n'

#################Implement bandwidth & time smearing ##########################
## We need to work out what number of channels to combine together in order
## to get closest to the value that we derive in the function
## BWSmear_FoV
## Use TimeSmear_FoV for time smearing
###############################################################################

# Use function to get channel width in MHz
max_baseline = maxbaseline(msName)
reqChan_width = BWSmear_FoV(max_baseline,cfreq,smearing,theta)*1.e6
channel_ave = int(round((reqChan_width/widthChan[0]), 0))


# And time smearing
timesmear = TimeSmear_FoV(max_baseline,cfreq,smearing,theta)

print 'We will use split to average your data to give %d percent smearing at %.1f arcsec from desired source\n' % (smearingpc,theta)
print 'We have calculated a maximum baseline of %.1f km to give:\n' %(max_baseline/1000)
print 'Bandwidth averaging:'
print 'Required channel width: %.3f MHz' % (reqChan_width/1.e6)
print 'Achieved with averaging of each %d channels per spw\n' %channel_ave
print 'Time averaging'
print 'Ths data will be averaged each %.3f seconds\n' % timesmear

print max_baseline,cfreq,smearing,theta


#dynamic range for time smearing ADD IN, thats what the timesmear/3 is for!
'''
split(vis=msName, outputvis=msName[:-3]+'averaged',width=channel_ave,timebin=str(timesmear/2)+'s')
'''







