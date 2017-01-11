#########################################
### Jack Morford - Version: 240614
#########################################

### A script to create multiple plot files intended to recreate AIPS' SPFLG plotting function for e-Merlin Radio data

### A python script writen in Parseltongue for formulation within AIPS

### This script 'should' not need altering, simply edit the input file named: SPPlot_input.py

### The core plotting routines have been adapted from PLOT_LIB routines written by Hans-Rainer Klockner.


#########################################
### Modules

import os
import sys
import cPickle as pkle
import os.path

from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData

import numpy as np

import math, time, datetime
from time import gmtime, strftime, localtime
ti = time.time()    # To time the script

###########################################
### Functions

def time2hms(seconds):
    h=int(seconds/3600)
    m=int(seconds % 3600)/60
    s=seconds-(h*3600)-(m*60)
    h=`h`
    m=`m`
    s="%4.2f" % s
    hms=h.zfill(2)+":"+m.zfill(2)+":"+s.zfill(4)
    return hms

def get_ant_dictionary(AN_version=0):
	# Returns a dictionary of the antenna 'Names':Numbers
	an_table = uvdata.table('AN', AN_version)

	antnumbers = []
	for row in an_table:
		station = row.nosta
		antnumbers.append(station)
	
	count = 0
	a_dict = {}
	for antenna in uvdata.antennas:
		a_dict[antenna] = antnumbers[count]
		count += 1

	return a_dict

def get_uvcop_bslines(datav, AN_version=0):
	# Returns the each baseline present in a UV dataset.
	import Wizardry.AIPSData as wiz

	data = wiz.AIPSUVData(datav[0],datav[1],datav[2],datav[3])

	tmpbsllist = []
	
	for vis in data:
		bsl = vis.baseline
		if bsl not in tmpbsllist and bsl[0] != bsl[1]:
			tmpbsllist.append(bsl)
		
	return tmpbsllist
		



def getbsllist(AN_version=0):
	# Creates a list of all the baslines within the dataset, format = [[1,2],[1,5],etc]
	a_dict = get_ant_dictionary() 
	bsllist = []
	
	for ant1 in a_dict.values():
		for base1 in a_dict.values():
			if ant1 != base1:
				tmp = []
				tmp.append(ant1)
				tmp.append(base1)
				tmp.sort()
				bsllist.append(tmp)

	bsllist.sort()
	i = 1
	for duplicate in xrange((len(bsllist)/2)):
		del bsllist[i]
		i += 1
	
	return bsllist


def listr():
	listr = AIPSTask('LISTR')
	listr.indata = uvdata
	listr.optype = 'SCAN'
	listr.inext = 'SU'
	listr.inver = 0
	listr.dparm = AIPSList([0,0,0,0,1,0,0,0,0])
	listr.docrt = -1
	listr.outprint = path2folder + 'listrscan.txt'
	listr.go()


def isNaN(x): # Hans!
    return x != x


def isnanv(x): # Hans!
    # Returns False of True if nan or NaN
    
    from numpy import isnan
    if list(isnan(x)).count(True) > 0:
        return(True)
    else:
        return(False)


def prtime(t): # Hans!
        # Convert time of 24hrs integer 
        # Return (day, hour, min, sec)
        
        from math import floor
        
        day  = floor(t)

        hour = 24 * (t - floor(t))
        min  = (hour - floor(hour))*60
        sec  = (min - floor(min))*60
        return (int(day),int(hour),int(min),sec)


def kdemean(x,domin=1): # Hans!
    # Use the Kernel Density Estimation (KDE) to determine the mean
    # (http://jpktd.blogspot.com/2009/03/using-gaussian-kernel-density.html )
    
    from scipy.stats import gaussian_kde
    from numpy import linspace,min,max,mean,std
    from math import sqrt,log
    accu = 50   # accuracy in determine the values of the KDE
    chrangeminmax = 1.0

    if isnanv(x) == True:
	    print 'kde isnanv'
	    return(0,0)
    if len(x) == 0:
	    print 'kde len(x) = 0'
	    return(0,0)
    if max(x) == min(x):
	    print 'kde max = min'
	    return(max(x),max(x))
    if mean(x) == std(x):
	    print 'kde mean = std'
	    return(mean(x),std(x))

    maxra = max([abs(max(x)),abs(min(x))])

    # create instance of gaussian_kde class
    gk     = gaussian_kde(x)

    fwhma = 0
    fwhmb = 0
    i = 0
    while fwhma == fwhmb == 0:
	    chrangeminmax += 1.0

	    vra    = linspace(-1*chrangeminmax*maxra,chrangeminmax*maxra,accu)
	    vraval = gk.evaluate(vra)

	    idx = 0
	    for i in range(len(vra)):
		if vraval[i] == max(vraval):          
		    mval = vra[i]
		    idx = i
		    break

	    dfwa  = -1
	    dfwb  = -2
	    for i in range(len(vra)):
		    if idx -i >= 0:
			if vraval[idx - i] < max(vraval)/2. and dfwa != 1:
			    fwhma = vra[idx] - vra[idx-i]
			    dfwa = 1

		    if len(vraval) > idx + i: 
			if vraval[idx + i] < max(vraval)/2. and dfwb != 1:
			    fwhmb = vra[idx + i] - vra[idx] 
			    dfwb = 1
		    if dfwa == dfwb:
			break
	    if chrangeminmax > 1000:
		    print 'KDE runaway',chrangeminmax
		    return(mean(x),std(x))
	    
    cfifm = delzero([fwhma,fwhmb])

    if domin == 1: 
        fwhm = min(cfifm)
    else:
        fwhm = max(cfifm)
 
    # factor 2 is because only one side is evaluated
    sigma = 2*fwhm/(2*sqrt(2*log(2)))

    return(mval,abs(sigma))


def givestokes(data): # Hans!
    # Provides you with the stokes parameter 
    
    from AIPSData import AIPSUVData as UVDataS
    uvdata  = UVDataS(data[0],data[1],data[2],data[3])
    st      = uvdata.stokes
    #
    if len(st) == 0:
        print '\nCAUTION could not determine STOKES parameter (givestokes)\n'
        sys.exit(-1)
    return(st)


def datatype(data): # Hans! - edited by Jack Morford
    # Returns the datatype of the file
    from AIPSData import AIPSCat
    from string import split,find
    
    datatype = 'None'
    catalog_in = AIPSCat(int(data[2]))
    catalog    = str(catalog_in).split()
    
    mymatch = []
    for i in xrange(len(catalog)):
    	if catalog[i] == data[0]:
    		mymatch.append(i)

    for i in range(len(catalog)-4):
        if catalog[i].find(data[0]) == 0:
        	if len(catalog[i+1]) < 7 and catalog[i+4] == 'UV':
        		datatype = str(catalog[i+4])
        		return datatype
        		#return 'UV'
        	if len(catalog[i+1]) == 8 and catalog[i+3] == 'UV':
        		datatype = str(catalog[i+3])
        		return datatype
        		#return 'UV'
    '''
    matchfname.append(i)
    matchfname.append(len(catalog))
    print matchfname
    classindex = -1
    isname  = -1
    isclass = -1
    isdisk  = -1
    isseq   = -1
    i = 0
    while (isname != 1 or isclass != 1 or isdisk !=1 or isseq !=1 or classindex != -1) and i < len(matchfname)-1:
        isname  = -1
        isclass = -1
        isdisk  = -1
        isseq   = -1
        for j in range(matchfname[i],matchfname[i+1]):

            if catalog[j].replace('.','') == data[0].replace('.',''):
                isname = 1
            if catalog[j].replace('.','') == data[1]:
                isclass = 1
            if catalog[j].replace('.','') == str(int(data[2])):
                isdisk = 1
            if catalog[j].replace('.','') == str(int(data[3])):
                isseq = 1
                        
            if catalog[j] == 'UV' and isname == 1 and isclass == 1 and isdisk == 1 and isseq ==1 and classindex == -1:
                classindex = j
            if catalog[j] == 'MA' and isname == 1 and isclass == 1 and isdisk == 1 and isseq ==1 and classindex == -1:
                classindex = j
        i += 1

    if classindex != -1:
        datatype = catalog[classindex]

    return(datatype)
    '''

def mergepdfs(tmp, bsl, amporphas='A'):

	mergedfile = outfilename+"_"+str(amporphas)+"_"+str(bsl[0])+"-"+str(bsl[1])+"_"+str(flagver)+".pdf"

	if len(tmp) > 1 and len(tmp) < 3:
		#mergedfile = outfilename+"_"+str(amporphas)+"_"+str(bsl[0])+"-"+str(bsl[1])+".pdf"
		first = str(tmp[1])
		second = str(tmp[0])
		os.system('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile='+str(path2folder+mergedfile)+' '+str(path2folder+first)+' '+str(path2folder+second))

                os.remove(path2folder+first)
                os.remove(path2folder+second)
                
		print "\nMerge Completed --> "+str(path2folder+mergedfile)
	
	if len(tmp) > 2:
		for i in xrange(len(tmp)):
			tmp[i] = str(path2folder)+tmp[i]

                tmp = tmp[::-1]
	
		#mergedfile = outfilename+"_"+str(amporphas)+"_"+str(bsl[0])+"-"+str(bsl[1])+".pdf"
		mystring2 = ' '.join(tmp)

		os.system('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile='+str(path2folder+mergedfile)+' '+str(mystring2))
	
                for i in xrange(len(tmp)):
                    os.remove(tmp[i])

		print "\nMerge Completed --> "+str(path2folder+mergedfile)
	
	if len(tmp) == 1:
		print "No merging required - only one plotfile"

        finalfile = str(path2folder+mergedfile)
        return finalfile


def mergebslpdfs(fulllist,amporphas):

	mergedfile = outfilename+"_"+str(amporphas)+"_"+"ALL"+"_"+str(flagver)+".pdf"

	if len(fulllist) > 1 and len(fulllist) < 3:
		#mergedfile = outfilename+"_"+str(amporphas)+"_"+str(bsl[0])+"-"+str(bsl[1])+".pdf"
		first = str(fulllist[0])
		second = str(fulllist[1])
		os.system('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile='+str(path2folder+mergedfile)+' '+str(first)+' '+str(second))

		print "\nMerge Completed --> "+str(path2folder+mergedfile)
	
	if len(fulllist) > 2:
		# for i in xrange(len(tmp)):
			# tmp[i] = str(path2folder)+tmp[i]
	
		#mergedfile = outfilename+"_"+str(amporphas)+"_"+str(bsl[0])+"-"+str(bsl[1])+".pdf"
		mystring2 = ' '.join(fulllist)

		os.system('gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile='+str(path2folder+mergedfile)+' '+str(mystring2))
	
 		print "\nMerge Completed --> "+str(path2folder+mergedfile)
	
	if len(fulllist) == 1:
		print "No merging required - only one plotfile"





def gensinbsl(bsl, flagver, deleteprev=True):
	### Function to generate a single baseline UV datafile within AIPS


	ant1 = bsl[0]
	base1 = bsl[1]

	### Choose ALL channels (0, 0)
	bchan = 0
	echan = 0
	### Outname
	outname = 'sbsl'+str(ant1)+'-'+str(base1)+"_"+str(flagver)+'_'+str(source)
	outclass = 'UVBSL'
	outdata = 1

	
	uvbsldata = AIPSUVData(outname, outclass, int(outdata), int(1))

	if deleteprev == True:
		if uvbsldata.exists() == False:
			print "\nNo previous Single BSL UV Datafile exists silly!"
			print "Creating baselines: "+str(ant1)+"_"+str(base1)+", FG: "+str(flagver)+" UV datafile"
			try:
				uvcop = AIPSTask('UVCOP')
				uvcop.indata = uvdata
				uvcop.outname = outname
				uvcop.outclass = outclass
				uvcop.sources = AIPSList([str(srcname)])
				uvcop.antenna[1] = ant1
				uvcop.baseline[1] = base1
				uvcop.bchan = bchan
				uvcop.echan = echan
				#if len(source_list) > 1:
				#	uvcop.sources = 
				uvcop.outdisk = int(outdata)
				uvcop.flagver = flagver
				uvcop.go()

				print "\nGenerated Single baseline UV datafile: "+str(ant1)+"-"+str(base1)+" using flagversion: "+str(flagver)

				uvbsldatav = [outname, outclass, int(outdata), int(1)]
				uvbsldata = AIPSUVData(uvbsldatav[0], uvbsldatav[1], uvbsldatav[2], uvbsldatav[3])

				return uvbsldatav, uvbsldata

			except RuntimeError:
				print "\nPlease ensure your datafile includes visibilities within your specified baselines"
				print "Could not find visibilities in baseline: "+str(ant1)+"-"+str(base1)
				#print "Aborting!"
				#sys.exit(0)

		if uvbsldata.exists() == True:
			uvbsldata.zap()
			print "\nPrevious Single baseline datafile deleted!"
			print "Creating baselines: "+str(ant1)+"_"+str(base1)+", FG: "+str(flagver)+" UV datafile"
			try:
				uvcop = AIPSTask('UVCOP')
				uvcop.indata = uvdata
				uvcop.outname = outname
				uvcop.outclass = outclass
				uvcop.sources = AIPSList([str(srcname)])
				uvcop.antenna[1] = ant1
				uvcop.baseline[1] = base1
				uvcop.bchan = bchan
				uvcop.echan = echan
				uvcop.outdisk = int(outdata)
				uvcop.flagver = flagver
				uvcop.go()

				print "\nGenerated Single baseline UV datafile: "+str(ant1)+"-"+str(base1)+" using flagversion: "+str(flagver)

				uvbsldatav = [outname, outclass, int(outdata), int(1)]
				uvbsldata = AIPSUVData(uvbsldatav[0], uvbsldatav[1], uvbsldatav[2], uvbsldatav[3])

				return uvbsldatav, uvbsldata

			except RuntimeError:
				print "\nPlease ensure your datafile includes visibilities within your specified baselines"
				print "Could not find visibilities in baseline: "+str(ant1)+"-"+str(base1)
				#print "Aborting!"
				#sys.exit(0)


	if deleteprev == False:
		if uvbsldata.exists() == True:
			uvbsldatav = [outname, outclass, int(outdata), int(1)]
			uvbsldata = AIPSUVData(uvbsldatav[0], uvbsldatav[1], uvbsldatav[2], uvbsldatav[3])
			print "\nSingle basline: "+str(ant1)+"-"+str(base1)+", FG: "+str(flagver)+" UV datafile already exists - using existing datafile"
			print "\nOtherwise, quit program and select deleteprev=True in input file"
			return uvbsldatav, uvbsldata

		if uvbsldata.exists() == False:
			try:
				uvcop = AIPSTask('UVCOP')
				uvcop.indata = uvdata
				uvcop.outname = outname
				uvcop.outclass = outclass
				uvcop.antenna[1] = ant1
				uvcop.baseline[1] = base1
				uvcop.bchan = bchan
				uvcop.echan = echan
				uvcop.outdisk = int(outdata)
				uvcop.flagver = flagver
				uvcop.go()

				print "\nGenerated Single baseline UV datafile: "+str(ant1)+"-"+str(base1)+" using flagversion: "+str(flagver)

				uvbsldatav = [outname, outclass, int(outdata), int(1)]
				uvbsldata = AIPSUVData(uvbsldatav[0], uvbsldatav[1], uvbsldatav[2], uvbsldatav[3])

				return uvbsldatav, uvbsldata

			except RuntimeError:
				print "\nPlease ensure your datafile includes visibilities within your specified baselines"
				print "Could not find visibilities in baseline: "+str(ant1)+"-"+str(base1)
				#print "Aborting!"
				#sys.exit(0)
			

		
	# This created a datafile within the AIPS disk of the single chosen baseline.

def sbslplot(outfilename, bsl, flagver, stokes=['RR'], timeperpage=1057, deleteprev=True, amporphas='A', IF_start=1, IF_end=8):
	### Format is sbslplot('outfilename', which baseline, flagversion, stokes, timeperpage, deleteprevious UV datafile?, amplitude or phase)
	# timeperpage is the amount of visibilities to plot per page corresponding to time in the y-axis.
	
	ant1 = bsl[0]
	base1 = bsl[1]
	
	try:
		uvbsldatav, uvbsldata = gensinbsl(bsl, flagver, deleteprev)


		### Generate a basline sorted python dictionary from time sorted UV Data using datatobase function
		### but if it already exists (you have ran it before) then simply load in the pickle file
		picklefilename = picklepath+"/"+str(outfilename)+"-"+str(srcname)+"_FG:"+str(flagver)+"_dic.p"


		if os.path.isfile(picklefilename) == False:
			dbasebsl,dbasebslhead  = datatobase(uvbsldatav,outfilename,flagver,True)
		else:
			print "\nVisibility data exists in pickle file - loading back in..."
			dbasebsl = pkle.load(open(picklefilename, "rb"))


		print "\n----------------- Plotting your pretty pictures -----------------\n"
		
		return plotspec(dbasebsl, bsl, outfilename, stokes, timeperpage, amporphas, IF_start, IF_end)
	
	except TypeError:
		print "No visibilities in baseline: "+str(ant1)+"-"+str(base1)+" continuing..."


def datatobase(data,outfilename,flagver,doheader=False,picklesave=True):
	# Create a dictionary containing all the IFs, visibilities, channels. (per baseline, per flagversion)
	from numpy import array
	import numpy
	from math import sqrt
	from copy import copy
	#from Wizardry.AIPSData import AIPSUVData as AIPSDataVIS

	if  datatype(data) == 'UV':
		from Wizardry.AIPSData import AIPSUVData as AIPSDataVIS
	else:
		print 'Data is not UV data',data,' abort.'

	stokes = givestokes(data)
	#ant1 = bsl[0]
	#base1 = bsl[1]

	visdata = AIPSDataVIS(data[0],data[1],data[2],data[3])
	print '\n----------------- Loading the data into memory ------------------\n'

	hi = {}
	hi['HEADER'] = copy(visdata.header)
	
	h = {}

	for visibility in visdata:

		try: 
			srcnumber = visibility.source		
		except KeyError:
			srcnumber = 1.

                #print srcnumber, source_list.keys()
		if srcnumber in source_list.keys():
			
			for st in range(len(stokes)):
                            
				for IFs in xrange(1,IF+1):
					
					if isnanv(array(visibility.visibility[IFs-1,:,st,1]).flatten()) == False:

						if h.has_key(str(visibility.baseline)) == False:
							h[str(visibility.baseline)] = {}

						if h.has_key(str(visibility.baseline)) == True:
					
							if h[str(visibility.baseline)].has_key(str('sour'+str(IFs)+str(stokes[st]))) == True:
								#if visibility.source == specifysources.keys()
								h[str(visibility.baseline)]['sour'+str(IFs)+str(stokes[st])].append(srcnumber)
							else:
								h[str(visibility.baseline)][str('sour'+str(IFs)+str(stokes[st]))] = [srcnumber]

							if h[str(visibility.baseline)].has_key(str('time'+str(IFs)+str(stokes[st]))) == True:
								h[str(visibility.baseline)]['time'+str(IFs)+str(stokes[st])].append(visibility.time)
							else:
								h[str(visibility.baseline)][str('time'+str(IFs)+str(stokes[st]))] = [visibility.time]


							if h[str(visibility.baseline)].has_key(str('real'+str(IFs)+str(stokes[st]))) == True:
								h[str(visibility.baseline)]['real'+str(IFs)+str(stokes[st])].append(array(visibility.visibility[IFs-1,:,st,0]))
							else:
								h[str(visibility.baseline)][str('real'+str(IFs)+str(stokes[st]))] = [array(visibility.visibility[IFs-1,:,st,0])]


							if h[str(visibility.baseline)].has_key(str('imag'+str(IFs)+str(stokes[st]))) == True:
								h[str(visibility.baseline)]['imag'+str(IFs)+str(stokes[st])].append(array(visibility.visibility[IFs-1,:,st,1]))
							else:
								h[str(visibility.baseline)][str('imag'+str(IFs)+str(stokes[st]))] = [array(visibility.visibility[IFs-1,:,st,1])]


							if h[str(visibility.baseline)].has_key(str('weig'+str(IFs)+str(stokes[st]))) == True:
								h[str(visibility.baseline)]['weig'+str(IFs)+str(stokes[st])].append(array(visibility.visibility[IFs-1,:,st,2]))
							else:
								h[str(visibility.baseline)][str('weig'+str(IFs)+str(stokes[st]))] = [array(visibility.visibility[IFs-1,:,st,2])]


								#if h[str(visibility.baseline)].has_key(str('uvdis'+str(IFs)+str(stokes[st]))) == True:
								#	h[str(visibility.baseline)]['uvdis'+str(IFs)+str(stokes[st])].append(sqrt(visibility.uvw[0]**2+visibility.uvw[1]**2+visibility.uvw[2]**2))
								#else:
								#	h[str(visibility.baseline)][str('uvdis'+str(IFs)+str(stokes[st]))] = [sqrt(visibility.uvw[0]**2+visibility.uvw[1]**2+visibility.uvw[2]**2)]


								#if h[str(visibility.baseline)].has_key(str('uvw'+str(IFs)+str(stokes[st]))) == True:
								#	h[str(visibility.baseline)]['uvw'+str(IFs)+str(stokes[st])].append([visibility.uvw[0],visibility.uvw[1],visibility.uvw[2]])
								#else:
								#	h[str(visibility.baseline)][str('uvw'+str(IFs)+str(stokes[st]))] = [[visibility.uvw[0],visibility.uvw[1],visibility.uvw[2]]]

	

	if picklesave == True:
		try:
			os.mkdir(picklepath)

		except OSError:
			print picklepath,"--> Directory already exists!"
			pass

		except IOError:
			print "Check you have the correct permissions to use mkdir"
			print "Abort!\n"
			sys.exit(0)

		
		#picklefilename = picklepath+"/"+str(outfilename)+"-"+str(srcname)+"_FG:"+str(flagver)+'_'+str(outclass)+"_dic.p"

		if doheader == True:
			### Save the basline sorted dictionary in picklepath to speed up later runs
			pkle.dump(h, open(picklefilename, "wb"))
			print "\nPickled the visibility data into:\n"+str(picklefilename)+'\n'
			return(h,hi)
		else:
			pkle.dump(h, open(picklefilename, "wb"))
			print "\nPickled the visibility data into:\n"+str(picklefilename)+'\n'
			return(h)

	if picklesave == False:
		if doheader == True:
			return(h,hi)

		if doheader == False:
			return(h)


def plotspec(database, bsl, plfilename='None', stokes=['RR'], timeperpage=1057, amporphas='A', IF_start=1, IF_end=8):
	# Takes the required infomation from the created visibility dictionary and sends it to the plotting function makespecplot
	from numpy import sqrt, power, arctan2, array
	from numpy import std, mean, median
	from math import pi
	
	print "\n---------------- SPPloting your pretty pictures ----------------\n"

	prtfiles = []
	drea = []
	dima = []
	weig = []
	time = []
	amps = []
	phas = []

	if database.has_key(str(bsl)) == False:
		return prtfiles

	for srcnumber in source_list.keys():
		srcname = source_list[srcnumber]
		for st in xrange(len(stokes)):
                        drea = []
                        dima = []
			amps = []
			phas = []
			weig = []
			sour = []
			count = 0
			for IFs in xrange(1,IF+1):
				if database[str(bsl)].has_key('real'+str(IFs)+stokes[st]):
					drea.append(database[str(bsl)]['real'+str(IFs)+stokes[st]])
					dima.append(database[str(bsl)]['imag'+str(IFs)+stokes[st]])
					weig.append(database[str(bsl)]['weig'+str(IFs)+stokes[st]])
					#weig.append(database[str(bsl)]['weig'+str(IFs)+stokes[st]])
					time = database[str(bsl)]['time'+str(IFs)+stokes[st]]
					sour = database[str(bsl)]['sour'+str(IFs)+stokes[st]]
					#time.append(database[str(bsl)]['time'+str(IFs)+stokes[st]])
					
					amp = sqrt(power(drea[count],2)+power(dima[count], 2))
					pha = (360./(2*pi))*arctan2(dima[count], drea[count])
					
					amps.append(amp)
					phas.append(pha)
					count += 1


			amps3 = np.array(amps)
			phas2 = np.array(phas)
			weig2 = np.array(weig)
			time2 = np.array(time)
			sours = np.array(sour)


			amps2 = amps3[:,np.where(sours[:] == srcnumber)]
			amps2 = np.squeeze(amps2)
			phas = phas2[:,np.where(sours[:] == srcnumber)]
			phas = np.squeeze(phas)
			weig = weig2[:,np.where(sours[:] == srcnumber)]
			weig = np.squeeze(weig)
			
			time = time2[np.where(sours[:] == srcnumber)]

			amps = np.where(weig[:] > 0.0, amps2[:], float('NaN'))
                        #print len(time), len(amps)

			if len(amps[0]) < timeperpage:
				plotname = makespecplot(amps, phas, weig, time, bsl, srcname, stokes[st], plfilename+'_'+str(srcname)+'_0-'+str(len(amps[0])),amporphas,IF_start,IF_end)
				prtfiles.append(plotname)
			

			else:
				for timestamp in xrange(0,len(amps[0]),timeperpage):
					if len(amps[0]) < timestamp+timeperpage:
						plotname = makespecplot(amps[:,timestamp:timestamp+timeperpage,:], phas[:,timestamp:timestamp+timeperpage,:], weig[:,timestamp:timestamp+timeperpage,:], time[timestamp:timestamp+timeperpage], bsl, srcname, stokes[st], plfilename+'_'+str(srcname)+'_'+str(timestamp)+'-'+str(len(amps[0])),amporphas,IF_start,IF_end)
						prtfiles.append(plotname)
					else:
						plotname = makespecplot(amps[:,timestamp:timestamp+timeperpage,:], phas[:,timestamp:timestamp+timeperpage,:], weig[:,timestamp:timestamp+timeperpage,:], time[timestamp:timestamp+timeperpage], bsl, srcname, stokes[st], plfilename+'_'+str(srcname)+'_'+str(timestamp)+'-'+str(timestamp+timeperpage),amporphas,IF_start,IF_end)
						prtfiles.append(plotname)
			
	#print "\n Plots created include: ", prtfiles
        prtfiles[:] = [x for x in prtfiles if x !=  None]
	print "\n Plots created include: ", prtfiles
	return prtfiles		




def makespecplot(amp,pha,weight,time,bsl,srcname='None',stokes='RR',plfilename='None',amporphas='A',IF_start=1,IF_end=8,tlabspace=15.0,sigma=3):
	# Decides whether Amplitude or Phase is plotted
	if amporphas == 'A':
		tmp = makespecplotsingle(amp,time,weight,bsl,srcname,plfilename,IF_start,IF_end,textstrg = ['AMP',bsl,stokes],fixr=[[-1E9,1E9],[0,0],sigma])
		return tmp
	
	if amporphas == 'P':
		tmp = makespecplotsingle(pha,time,weight,bsl,srcname,plfilename,IF_start,IF_end,textstrg = ['PHA',bsl,stokes],fixr=[[-99,98],[0,0],sigma])
		return tmp

	#allplfiles = []
	#allplfiles.append(makespecplotsingle(amp,time,weight,plfilename,textstrg = ['AMP',bsl,stokes],fixr=[[-1E9,1E9],[0,0],sigma]))
	#allplfiles.append(makespecplotsingle(pha,time,weight,plfilename,textstrg = ['PHA',bsl,stokes],fixr=[[-99,98],[0,0],sigma]))
	#allplfiles.append(makespecplotsingle(weight,time,weight,plfilename,textstrg = ['WEIGHTS',bsl,stokes],fixr=[[0,1E9],[0,0],sigma]))
	#return(allplfiles)
	


def makespecplotsingle(datachntime,time,weights,bsl,srcname='None',plfilename='None',IF_start=1,IF_end=8,textstrg=[],fixr=[[-1E9,1E9],[0,0],3],statstype='norm',tlabspace=15.0):
	# Plot the dynamic spectrum via Matplotlib
	
	# fixr = [[lower limit, upper limit],[mean,rms],sigma]

	# Plotting 
	
	from numpy import std,mean,median,arange
	import matplotlib
	matplotlib.use('Agg') # force the antigrain backend
	from matplotlib.backends.backend_pdf import PdfPages
	from matplotlib.pyplot import figure, close, savefig, rcParams, cm, cla, clf, text
	from matplotlib.colors import LogNorm
	

	pltype    = '.pdf'
	DPI       = 500
	fsize     = 8           # fontsize
	stlaboff  = 15

	
	### Set time (y=axis) spacing
	if len(datachntime[0]) == timeperpage:
		tlabspace = 15

	if len(datachntime[0]) < timeperpage:
		tlabspace = 15


	### set channel spacings
	chnspa = len(datachntime[0][0])/4

	# garbage collector
	import gc
	gc.collect()
        
	#for IFs in xrange(8):
	# Produce time labels 
	dstype = textstrg[0]
	bsl    = textstrg[1]
	stokes = textstrg[2]
	
	### i.e. create a list of the x-axies (time) ticks, len(tticks) = tlabspace
	tticks = []	
	if len(time) > tlabspace:
		for i in xrange(int(tlabspace)):
			tticks.append(int(i * len(time)/tlabspace))
		tticks.append(len(time)-1)
	else:
		for i in xrange(len(time)):
			tticks.append(i)
	
	### Create the y-axis (time) labels [day,hour,min,secs]
	tlab = []
	for i in xrange(len(tticks)):
		tgt = prtime(time[tticks[i]])
		tlab.append([tgt[0],tgt[1],tgt[2],int(tgt[3])])

	### Create a list for the amplitude plotting limits
	ovv   = []
	ovv2 = []
	for IFs in xrange(IF_start-1,IF_end):
	        nchan = 0
	        ovv2 = []
		for i in xrange(len(datachntime[IFs])):
			for j in xrange(len(datachntime[IFs][i])):
				if weights[IFs][i][j] > 0:
					ovv2.append(datachntime[IFs][i][j])
	                        nchan = j+1
		ovv.append(ovv2)

	if len(ovv) == 0:
		return('')


	### Create both x-axis (channels) ticks and tick labels
	xtck,xtckl = [],[]
	if nchan >= 8:
		for i in xrange(int(nchan/chnspa)):
			xtckl.append(int((i+1)*chnspa))
			xtck.append(int(((i+1)*chnspa)))

	else:
		for i in xrange(nchan):
			xtckl.append(int(i+1))
			xtck.append(int(i))

	### Define the max and min of the plots for each IF
	minim = []
	maxim = []
	dmean = []
	rms = []
        missingif = []

	for IFs in xrange(noIFs):
		try:
			minim.append(min(ovv[IFs]))
			maxim.append(max(ovv[IFs]))
			dmean.append(mean(ovv[IFs]))
			rms.append(std(ovv[IFs]))
		except ValueError:
			print "\nNo visibiilties exist in IF:%i for source: %s, baseline:" % (IFs+1, srcname), bsl
                        missingif.append(IFs)

        #print minim
        #print noIFs, len(missingif)
        if len(missingif) == noIFs:
            return
	### Average over all IFs to use in plot title!
	avmin = sum(minim)/len(minim)
	avmax = sum(maxim)/len(maxim)
	avdmean = sum(dmean)/len(dmean)
	avrms = sum(rms)/len(rms)

        #print len(minim), missingif
        #print minim

        for i in range(len(missingif)):
            minim.insert(missingif[i],avmin)
            maxim.insert(missingif[i],avmax)
            dmean.insert(missingif[i],avdmean)
            rms.insert(missingif[i],avrms)

        #print len(minim)
        #print minim

	sigma   = fixr[2]
	
	if statstype == 'kde':
		kdestat  = kdemean(ovv)
		dmean    = kdestat[0]
		rms      = kdestat[1]
		
	maxvalueallifs = max(maxim)
	minvalueallifs = min(minim)
	
	ifwminamp = minim.index(min(minim))+1
	ifwmaxamp = maxim.index(max(maxim))+1

	vvmin = []
	vvmax = []
	for IFs in xrange(noIFs):
		try:
			vmin = dmean[IFs] - (abs(3)*rms[IFs])
			vmax = dmean[IFs] + (abs(3)*rms[IFs])
			vvmin.append(vmin)
			vvmax.append(vmax)
		except IndexError:
			vvmin.append(0), vvmax.append(0)
 

	if fixr[0][0]+fixr[0][1] != 0:
		if fixr[0][0] != -1E9 and fixr[0][1] != 1E9:
			vvmin = fixr[0][0]
			vvmax = fixr[0][1]
		if fixr[0][0] != -1E9 and fixr[0][1] == 1E9:
			vvmin = fixr[0][0]
			#vvmax = maxim
		if fixr[0][0] == -1E9 and fixr[0][1] != 1E9:
			vvmin = minim
			#vvmax = fixr[0][1]
		if fixr[0][0] == -99 and fixr[0][1] == 98:
			vvmin = minim
			vvmax = maxim
	if fixr[1][0]-fixr[1][1] != 0:
		dmean    = fixr[1][0]
		rms      = fixr[1][1]
		vvmin = dmean - (abs(sigma)*rms)
		vvmax = dmean + (abs(sigma)*rms)


	### Create the figure and the IF subplots
	fig1 = figure()

        # fig1.set_size_inches(8.27,11.69) # A4 Picture size in INCHES! - Portrait
	fig1.set_size_inches(11.69,8.27) # Landscape
	fig1.subplots_adjust(wspace=0) # An option close the whitespace between each IF
	
	if colourscheme == 'jet':
            cmapchoice = cm.jet
        elif colourscheme == 'spectral':
            cmapchoice = cm.spectral
        elif colourscheme == 'gist_rainbow':
            cmapchoice = cm.gist_rainbow
        elif colourscheme == 'cool':
            cmapchoice = cm.cool
        elif colourscheme == 'summer':
            cmapchoice = cm.summer
        elif colourscheme == 'winter':
            cmapchoice = cm.winter
        else:
            cmapchoice = cm.jet

	
	subplts = []
	count = 1
	for IFs in xrange(IF_start, IF_end+1):
		tmp = fig1.add_subplot(1,noIFs,count)

                #datachntime.shape

		if scale == 'std' and scale_over_all_IFs == 'yes':
			minvalueallifs = min(vvmin)
			maxvalueallifs = max(vvmax)
			cbarimage = tmp.imshow(datachntime[IFs-1],aspect='auto',interpolation='nearest',origin='lower',cmap=cmapchoice,vmin=minvalueallifs, vmax=maxvalueallifs)
		
		if scale == 'std' and scale_over_all_IFs == 'no':
			cbarimage = tmp.imshow(datachntime[IFs-1],aspect='auto',interpolation='nearest',origin='lower',cmap=cmapchoice,vmin=vvmin[count-1], vmax=vvmax[count-1])

		if scale == 'log' and scale_over_all_IFs == 'yes':
			cbarimage = tmp.imshow(datachntime[IFs-1],aspect='auto',interpolation='nearest',origin='lower',cmap=cmapchoice,norm=LogNorm(vmin=minvalueallifs, vmax=maxvalueallifs))

		if scale == 'log' and scale_over_all_IFs == 'no':
			cbarimage = tmp.imshow(datachntime[IFs-1],aspect='auto',interpolation='nearest',origin='lower',cmap=cmapchoice,norm=LogNorm(vmin=minim[count-1], vmax=maxim[count-1]))

		if scale == 'linear' and scale_over_all_IFs == 'yes':
			cbarimage = tmp.imshow(datachntime[IFs-1],aspect='auto',interpolation='nearest',origin='lower',cmap=cmapchoice,vmin=minvalueallifs, vmax=maxvalueallifs)

		if scale == 'linear' and scale_over_all_IFs == 'no':
			cbarimage = tmp.imshow(datachntime[IFs-1],aspect='auto',interpolation='nearest',origin='lower',cmap=cmapchoice,vmin=minim[count-1], vmax=maxim[count-1])
		
		
		xlabelminmax = "Channels \n Min: %.1g \n Max: %.1g \n (Jy)" % (minim[count-1],maxim[count-1])
		tmp.set_xlabel(xlabelminmax, fontsize=fsize)
		#tmp.set_xlabel('Channels' ,fontsize=fsize)
		tmp.set_ylabel('Observation Time', fontsize=fsize)
		tmp.set_title("IF "+str(IF_start), fontsize=fsize)
		tmp.set_xticks(xtck)
		tmp.set_xticklabels(xtckl,fontsize=fsize*0.75)
		tmp.set_yticks(tticks)
		tmp.set_yticklabels(tlab,fontsize=fsize)
		if count > 1:
			tmp.get_yaxis().set_visible(False)
		IF_start += 1
		count += 1
	

	
	### Create the colourbar
	cax = fig1.add_axes([0.91,0.1,0.02,0.8]) #Specify colourbar position, width and height
	cbar = fig1.colorbar(cbarimage, cax=cax)
	cbar.ax.tick_params(labelsize=fsize)


	if str(dstype).count('amp') or str(dstype).count('AMP') > 0:
		bottomtext = dstype+' BSL ['+str(bsl[0])+':'+str(bsl[1])+'] Stoke:'+str(stokes)+' Mean: '+str('%.5f'%avdmean)+' RMS: '+str('%.5f'%avrms)+' Max: '+str('%.5f'%maxvalueallifs)+' Jy'
		#fig1.suptitle(dstype+' BSL ['+str(bsl[0])+':'+str(bsl[1])+'] Stoke:'+str(stokes)+' Mean: '+str('%.5f'%avdmean)+' RMS: '+str('%.5f'%avrms)+' Max: '+str('%.5f'%avmax)+' Jy',fontsize=fsize*2)
	if str(dstype).count('pha') or str(dstype).count('PHA') > 0:
		bottomtext = dstype+' BSL ['+str(bsl[0])+':'+str(bsl[1])+'] Stoke:'+str(stokes)+' Mean: '+str('%.5f'%avdmean)+' RMS: '+str('%.5f'%avrms)+' Max: '+str('%.5f'%maxvalueallifs)+' deg'
	#	fig1.suptitle(dstype+' BSL ['+str(bsl[0])+':'+str(bsl[1])+'] Stoke:'+str(stokes)+' Mean: '+str('%.5f'%avdmean)+' RMS: '+str('%.5f'%avrms)+' Max: '+str('%.5f'%avmax)+' deg',fontsize=fsize*2)
	if str(dstype).count('nvis') or str(dstype).count('NVIS') > 0:
		bottomtext = dstype+' BSL ['+str(bsl[0])+':'+str(bsl[1])+'] Stoke:'+str(stokes)+' Mean: '+str('%.5f'%avdmean)+' RMS: '+str('%.5f'%avrms)+' Max: '+str('%.5f'%maxvalueallifs)+' counts'
	#	fig1.suptitle(dstype+' BSL ['+str(bsl[0])+':'+str(bsl[1])+'] Stoke:'+str(stokes)+' Mean: '+str('%.5f'%avdmean)+' RMS: '+str('%.5f'%avrms)+' Max: '+str('%.5f'%avmax)+' counts',fontsize=fsize*2)


	if scale == 'std' and scale_over_all_IFs == 'yes':
		righthandtext = 'Linear within 3 std: IFs all to scale'
	if scale == 'std' and scale_over_all_IFs == 'no':
		righthandtext = 'Linear within 3 std: IFs scaled independantly'
	if scale == 'log' and scale_over_all_IFs == 'yes':
		righthandtext = 'Log Scale: all IFs to scale'
	if scale == 'log' and scale_over_all_IFs == 'no':
		righthandtext = 'Log Scale: IFs scaled independantly'
	if scale == 'linear' and scale_over_all_IFs == 'yes':
		righthandtext = 'Linear Scale: all IFs to scale'
	if scale == 'linear' and scale_over_all_IFs == 'no':
		righthandtext = 'Linear Scale: IFs scaled independantly'

	plfname = plfilename.replace('.','_')+'_'+dstype[0]+'_'+str(stokes)+'_BL:'+str(bsl[0])+'_'+str(bsl[1])+"_FG:"+str(flagver)+str(pltype)
	fig1.text(0.98,0.5, righthandtext, rotation='vertical', verticalalignment='center',fontsize=fsize)
	fig1.text(0.5, 0.935, bottomtext, fontsize=fsize*1.3,horizontalalignment='center')
	fig1.suptitle(plfname, fontsize=fsize*1.7)
	pp = PdfPages(path2folder+plfname)
	savefig(pp, format='pdf', dpi=DPI)
	

	pp.close()
	cla()
	clf()
	close('all')

	return(plfname)


def aipsbaselines(baselinelist):
	tmp = []
	for i in baselinelist:
		tmp.append(i[0])
		tmp.append(i[1])

	tmp = set(tmp)
	tmp = list(tmp)
	tmp.sort()
	
	return tmp


def uvcop(outclass, baselines, flagver, deleteprev=True):
	
	outname = str(Name)+'_'+str(flagver)

	uvbsldatav = [outname, outclass, int(outdata), int(1)]
	uvbsldata = AIPSUVData(uvbsldatav[0], uvbsldatav[1], uvbsldatav[2], uvbsldatav[3])

	if deleteprev == True and uvbsldata.exists() == True:
		uvbsldata.zap()

		try:
			uvcop = AIPSTask('UVCOP')
			uvcop.indata = uvdata
			uvcop.outname = outname
			uvcop.outclass = outclass
			#uvcop.sources = AIPSList([str(srcname)])
			if baselines != 0:
				uvcop.antenna = AIPSList(baselines)
				uvcop.baseline = AIPSList(baselines)
			uvcop.bchan = 0
			uvcop.echan = 0
			uvcop.outdisk = 1
			uvcop.flagver = flagver
			uvcop.go()

			uvbsldatav = [outname, outclass, int(outdata), int(1)]
			uvbsldata = AIPSUVData(uvbsldatav[0], uvbsldatav[1], uvbsldatav[2], uvbsldatav[3])

			return uvbsldatav
		except RuntimeError:
			print "No visibilities exist for this baseline!"
			sys.exit(0)

	if deleteprev == True and uvbsldata.exists() == False:
		print "No previous UV dataset exists for these required baselines..."

		try:
			uvcop = AIPSTask('UVCOP')
			uvcop.indata = uvdata
			uvcop.outname = outname
			uvcop.outclass = outclass
			#uvcop.sources = AIPSList([str(srcname)])
			if baselines != 0:
				uvcop.antenna = AIPSList(baselines)
				uvcop.baseline = AIPSList(baselines)
			uvcop.bchan = 0
			uvcop.echan = 0
			uvcop.outdisk = 1
			uvcop.flagver = flagver
			uvcop.go()

			uvbsldatav = [outname, outclass, int(outdata), int(1)]
			uvbsldata = AIPSUVData(uvbsldatav[0], uvbsldatav[1], uvbsldatav[2], uvbsldatav[3])

			return uvbsldatav
		except RuntimeError:
			print "No visibilities exist for this baseline!"
			sys.exit(0)

	if deleteprev == False and uvbsldata.exists() == True:
		uvbslines = get_uvcop_bslines(uvbsldatav)
		
		count = 0
		for bsl in specifybaselines:
			for bline in uvbslines:
				if bline == bsl:
					count += 1
					break

		if count == len(specifybaselines):

			print "\nUsing previous uvcopped data %s" % uvbsldatav
			return uvbsldatav

		else:
			print "\nPrevious UVCOP'd dataset found but doesn't contain your desired baselines..."
			print "It only containes visibilities for baselines: %s" % uvbslines
			print "...not visibilities for baselines: %s" % specifybaselines, "as selected!"
			print "\nPlease select deletprev == True in the input file"
			print "\nAborting!"
			sys.exit(0)

	if deleteprev == False and uvbsldata.exists() == False:
		print "No previous UV dataset exists for these required baselines..."

		try:
			uvcop = AIPSTask('UVCOP')
			uvcop.indata = uvdata
			uvcop.outname = outname
			uvcop.outclass = outclass
			#uvcop.sources = AIPSList([str(srcname)])
			if baselines != 0:
				uvcop.antenna = AIPSList(baselines)
				uvcop.baseline = AIPSList(baselines)
			uvcop.bchan = 0
			uvcop.echan = 0
			uvcop.outdisk = 1
			uvcop.flagver = flagver
			uvcop.go()

			uvbsldatav = [outname, outclass, int(outdata), int(1)]
			uvbsldata = AIPSUVData(uvbsldatav[0], uvbsldatav[1], uvbsldatav[2], uvbsldatav[3])

			return uvbsldatav
		except RuntimeError:
			print "No visibilities exist for this baseline!"
			sys.exit(0)

def plotmysplots():
	
	bsllist = specifybaselines
        fulllist = []
	
	for bsl in bsllist:
		tmp = plotspec(dbasebsl, bsl, outfilename, stokes, timeperpage, amporphas, IF_start, IF_end)
		if len(tmp) == 0:
			print "No visibilities exist for baseline: "+str(bsl)+'\n'

		plotfiles.append(tmp)
 		finalfile = mergepdfs(tmp, bsl, amporphas)
                if os.path.isfile(finalfile):
                    fulllist.append(finalfile)
        #print fulllist

        if requestfullmerge and len(fulllist) > 1:
            mergebslpdfs(fulllist,amporphas)
            if removebslfile:
                print '\nRemoving individual baseline plots\n'
                for i in xrange(len(fulllist)):
                    os.remove(fulllist[i])
            else:
                print '\nKeeping individual baseline plots\n'


#################################################################################
################################# Main script ###################################


try:
    execfile("SPPlot_input.py")
    print "\nYour AIPS USER NUMBER is: ", AIPS.userno
    
except:
    print "\n Could not find SPPlot_input.py!"
    print " Make sure input file is in the same folder as this script"
    print " Aborting!\n"
    sys.exit(0)


try:
	AIPS.userno
	Name
	Klass
	Disk
	Seq
	path2folder
	picklepath
	picklesave
	picklelookfor
	IF
	IF_start
	IF_end
	choosesources
	if choosesources == 'choose':
		specifysources
	amporphas
	outfilename
	flagver
	stokes
	timeperpage
	deleteprev
	scale_over_all_IFs
	scale
	colourscheme
	choosebaselines
        if requestfullmerge == 'yes':
            removebslfile
	if choosebaselines == 'choose':
		specifybaselines
except NameError:
	print " Please ensure ALL input variables have been specified in the input file"
	print " Aborting!\n"
	sys.exit(0)


uvdatav          = [Name, Klass, Disk, Seq]
uvdata           = AIPSUVData(uvdatav[0],uvdatav[1],uvdatav[2],uvdatav[3])


indata = 1
outdata = indata
noIFs = (IF_end-IF_start)+1

if IF_start > IF_end:
	print "IF_start is larger than IF_end - please ammend!"
	print "Aborting!"
	sys.exit(0)


newsource = [] 
source_list = {}

try:
    nsource = len(uvdata.sources)
    print "\n Number of sources (nsource): %i" % nsource
    for s in xrange(len(uvdata.sources)):
        #source_list[float(s+1)] = str(uvdata.sources[s])
        print str(uvdata.sources[s])
    print " Names of sources within your dataset (with dictionary index):"

    for tab in uvdata.tables:
        if tab[1] == 'AIPS SU' or tab[1] == 'SU':
            sutable = uvdata.table('SU',tab[0])
    for row in sutable:
        if row.source.strip() in uvdata.sources:
            newsource.append([row.id__no,row.source.strip()])
    source_list = dict(newsource)

    
    for s in source_list:
        if not s == len(source_list):
        # if not s == nsource:
            print " %i: %s," % (s, source_list[s]),
        else:
            print " %i: %s" % (s, source_list[s])
    
except:
    print "\n No SU table found... "
    print " Assuming single source..."
    nsource = 1
    source_list[1] = Name.upper()
    print " %i: %s" % (1, source_list[1])


if choosesources == 'choose':
	nsource = len(specifysources)
	for source in specifysources:
		for k, v in source_list.items():
			if v != source:
				del source_list[k]

print "\nThe following sources will be plotted:"
for s in source_list:
	print " %s" % source_list[s]


if choosebaselines == 'all':
	specifybaselines = getbsllist()

plotfiles = []


if flagver == 0:
	if len(specifybaselines) <= 5:

		outclass = 'UVBSL'
		baselines = aipsbaselines(specifybaselines)
		
		
		uvbsldatav = uvcop(outclass, baselines, flagver, deleteprev)

		#dbasebsl,dbasebslhead  = datatobase(uvbsldatav,outfilename,flagver,True,False)
		
		picklefilename = picklepath+"/"+str(Name)+"_FG:"+str(flagver)+'_'+str(outclass)+"_dic.p"
		
		if picklelookfor == True:
			if os.path.isfile(picklefilename) == False:
				dbasebsl,dbasebslhead  = datatobase(uvbsldatav,outfilename,flagver,True, picklesave)
			else:
				print "\nVisibility data exists in pickle file - loading back in...\n"
				dbasebsl = pkle.load(open(picklefilename, "rb"))

		if picklelookfor == False:
			dbasebsl, dbasebslhead = datatobase(uvbsldatav,outfilename,flagver, True, picklesave)
		
		plotmysplots()


	if len(specifybaselines) > 5:
		
		outclass = 'UVABsL'
		
		picklefilename = picklepath+"/"+str(Name)+"_FG:"+str(flagver)+'_'+str(outclass)+"_dic.p"

		if picklelookfor == True:
			if os.path.isfile(picklefilename) == False:
				dbasebsl,dbasebslhead  = datatobase(uvdatav,outfilename,flagver,True, picklesave)
			else:
				print "\nVisibility data exists in pickle file - loading back in...\n"
				dbasebsl = pkle.load(open(picklefilename, "rb"))

		if picklelookfor == False:
			dbasebsl, dbasebslhead = datatobase(uvdatav,outfilename,flagver, True, picklesave)
		
		plotmysplots()


if flagver >= 1:
	if len(specifybaselines) <= 5:

		outclass = 'UVBSL'
		baselines = aipsbaselines(specifybaselines)

		uvbsldatav = uvcop(outclass, baselines, flagver, deleteprev)

		#dbasebsl, dbasebslhead = datatobase(uvbsldatav,outfilename,flagver,True,False)
		
		picklefilename = picklepath+"/"+str(Name)+"_FG:"+str(flagver)+'_'+str(outclass)+"_dic.p"
		
		if picklelookfor == True:
			if os.path.isfile(picklefilename) == False:
				dbasebsl,dbasebslhead  = datatobase(uvbsldatav,outfilename,flagver,True, picklesave)
			else:
				print "\nVisibility data exists in pickle file - loading back in...\n"
				dbasebsl = pkle.load(open(picklefilename, "rb"))

		if picklelookfor == False:
			dbasebsl, dbasebslhead = datatobase(uvbsldatav,outfilename,flagver, True, picklesave)
		
		plotmysplots()


	if len(specifybaselines) > 5:

		outclass = 'UVABsL'
		baselines = aipsbaselines(specifybaselines)

		uvbsldatav = uvcop(outclass, baselines, flagver, deleteprev)

		picklefilename = picklepath+"/"+str(Name)+"_FG:"+str(flagver)+'_'+str(outclass)+"_dic.p"

		if picklelookfor == True:
			if os.path.isfile(picklefilename) == False:
				dbasebsl,dbasebslhead  = datatobase(uvbsldatav,outfilename,flagver,True, picklesave)
			else:
				print "\nVisibility data exists in pickle file - loading back in...\n"
				dbasebsl = pkle.load(open(picklefilename, "rb"))

		if picklelookfor == False:
			dbasebsl, dbasebslhead = datatobase(uvbsldatav,outfilename,flagver, True, picklesave)
		
		plotmysplots()


plotfiles[:] = [x for x in plotfiles if x != []]
print "\nTotal plotfiles produced:\n",plotfiles

print "\nTime taken (hh:mm:ss):", time2hms(time.time()-ti)
print 'Finished on %s' % strftime("%d-%b-%y"), 'at:', strftime("%H:%M:%S", localtime()),'\n'

