# comment out this area surrounded by ~~ once installed into the eMERGE pipeline!
# For this to work you NEED to source your lofar stack & then aips
# i.e. for JBCA you would (in the terminal):
# 1) source /usr/local/lofar-2_15/init_env_release.csh; source /usr/local/lofar-2_15/lofar/release/lofarinit.csh
# 2) source /aips/LOGIN.CSH
# Adjust as according to your shell!

# Changelog
# v1.0 basic functionality
# v1.1 added function to remove intermediate files
# v1.2 added version number tracking
# v1.3 removed IO overhead by removing intermediate files in process & fixed single source problem
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

import os, re, time, datetime, sys, math, fnmatch, numpy
from os.path import join, getsize
from datetime import date
from collections import deque
import Utilities
from AIPS import AIPS, AIPSDisk
from AIPSTask import AIPSTask, AIPSList
from AIPSData import AIPSUVData, AIPSImage, AIPSCat
from Wizardry.AIPSData import AIPSUVData as WizAIPSUVData
import math, time, datetime
from numpy import *
import itertools
from time import gmtime, strftime, localtime
ti = time.time()    # To time the script
from subprocess import Popen, PIPE
from os import environ

def source(script, update=True):
    """
    http://pythonwise.blogspot.fr/2010/04/sourcing-shell-script.html (Miki Tebeka)
    http://stackoverflow.com/questions/3503719/#comment28061110_3505826 (ahal)
    """
    import subprocess
    import os
    proc = subprocess.Popen(
        ['bash', '-c', 'set -a && source {} && env -0'.format(script)], 
        stdout=subprocess.PIPE, shell=False)
    output, err = proc.communicate()
    output = output.decode('utf8')
    env = dict((line.split("=", 1) for line in output.split('\x00') if line))
    if update:
        os.environ.update(env)
    return env

def AOFlag_AIPS(uvdata,sensitivity,profile,casa,ref_freq,AOsources,NonAOsources,intermediate,rfistrategy):
	# Part 1: Split the data
	print sys.argv
	if len(sys.argv) < 2:
		print 'ERROR: please use as ParselTongue AOFlag_module_standalone.py CASA_convert_standalone.py'
		sys.exit()

	original = uvdata
	version = str(1.5)
	splat = AIPSTask('SPLAT')
	splat.indata = uvdata
	splat.docalib = 100
        splat.aparm[6]=1
        splat.aparm[7]=1
	splat.go()
	
	# Part 2: Export the data to be flagged
	for i in range(len(AOsources)):
		
		fittp = AIPSTask('FITTP')		
		fittp.indata = AIPSUVData(AOsources[i],'SPLAT',original.disk,1)
		fittp.dataout = 'PWD:'+AOsources[i]+'_pre_AO.fits'
		fittp.go()
        
        	#Change here to add TASAV file of each split so that they can be deleted
		tasav = AIPSTask('TASAV')
		tasav.indata = AIPSUVData(AOsources[i],'SPLAT',original.disk,1)
		tasav.outdata = AIPSUVData(AOsources[i],'TASAV',original.disk,1)
        	tasav.go()
                if intermediate == True:
			AIPSUVData(AOsources[i],'SPLAT',original.disk,1).zap()
        
        	# Part 3: Execute CASA script to convert data sets from fits to ms
		os.system(casa+' -c '+str(sys.argv[1])+' fits '+AOsources[i]+'_pre_AO.fits')
		if intermediate == True:
		    os.system('rm '+AOsources[i]+'_pre_AO.fits')
		    os.system('rm *log')
        
        	# Part 4: Create the rfistrategy defined
		if rfistrategy == '':
			os.system('rfistrategy -s '+str(sensitivity[i])+' '+str(profile)+' rfistrat_'+AOsources[i]+'_'+str(sensitivity[i])+'.rfis')
			rfistrat = 'rfistrat_'+AOsources[i]+'_'+str(sensitivity[i])+'.rfis'
		else:
			rfistrat = rfistrategy
		# Part 5: AOFlag the data to make it sparkly clean (add some functionality, plots? SPPlot?)        
		os.system('aoflagger -strategy '+rfistrat+' '+AOsources[i]+'_pre_AO.ms') #Aoflag the data
		os.system('mv '+AOsources[i]+'_pre_AO.ms '+AOsources[i]+'_post_AO.ms') #Change the name to avoid confusion
				
		# Part 6: Convert AOflagged data back to fits
		os.system(casa+' -c '+str(sys.argv[1])+' ms '+AOsources[i]+'_post_AO.ms')
		if intermediate == True:
		    os.system('rm -r '+AOsources[i]+'_post_AO.ms')
		    os.system('rm *log')
        
		# Part 7: Load AOFlagged data back into AIPS
		fitld = AIPSTask('FITLD')
		fitld.digicor = -1
		fitld.datain = 'PWD:'+AOsources[i]+'_post_AO.fits'
		fitld.outnam = AOsources[i]
		fitld.outclass = 'AO'
		fitld.go()
		if intermediate == True:
		    os.system('rm -r '+AOsources[i]+'_post_AO.fits')
                
		# Part 8: Move the Antenna table off the old TASAV'ed data to the new
		data = AIPSUVData(AOsources[i],'AO',original.disk,1)
		data.zap_table('AN',1)
		tacop = AIPSTask('TACOP')
		tacop.indata = AIPSUVData(AOsources[i],'TASAV',original.disk,1)
		tacop.outdata = data
		tacop.inext = 'AN'
		tacop.invers = 1
		tacop.go()
		if intermediate == True:
		    AIPSUVData(AOsources[i],'TASAV',original.disk,1).zap()
        
        	# Part 9: AXDEFINE the AOflagged data to reference pixel 1 and proper frequency
		data = WizAIPSUVData(AOsources[i],'AO',original.disk,1)
		data.header['crval'][2] = ref_freq
		data.header['crpix'][2] = 1
		data.header.update()
		
		
		# Part 10 Uvfix the AO data so the frequencies match
		uvfix = AIPSTask('UVFIX')
		uvfix.indata = AIPSUVData(AOsources[i],'AO',original.disk,1)
		uvfix.uvfixprm[14] = ref_freq
		uvfix.fqcenter = -1
		uvfix.go()
		if intermediate == True:
		    AIPSUVData(AOsources[i],'AO',original.disk,1).zap()
		print 'Next source'
		
	# Part 11: Uvfix the nonAOflagged data so the frequencies match
	if NonAOsources > 0:
		for i in range(len(NonAOsources)):
			uvfix = AIPSTask('UVFIX')
			uvfix.indata = AIPSUVData(NonAOsources[i],'SPLAT',original.disk,1)
			uvfix.uvfixprm[14] = ref_freq
			uvfix.fqcenter = -1
			uvfix.go()
			if intermediate == True:
			    AIPSUVData(NonAOsources[i],'SPLAT',original.disk,1).zap()
	
	# Part 12: Combine the sources, uvsort the data sets and index first
	sources = AOsources + NonAOsources	
	for i in range(len(sources)):	
		
		AIPSUVData(sources[i],'UVFIX',original.disk,1).zap_table('NX',1)
		uvsrt = AIPSTask('UVSRT')
		uvsrt.indata = AIPSUVData(sources[i],'UVFIX',original.disk,1)
		uvsrt.sort = 'TB'
		uvsrt.go()
		if intermediate == True:
		    AIPSUVData(sources[i],'UVFIX',original.disk,1).zap()
		
		indxr = AIPSTask('INDXR')
		indxr.indata = AIPSUVData(sources[i],'UVSRT',original.disk,1)
		indxr.go()

		# Part 13: Multi the data so DBCON cannot interfere
		multi = AIPSTask('MULTI')
		multi.indata = AIPSUVData(sources[i],'UVSRT',original.disk,1)
		multi.go()
		if intermediate == True:
		    AIPSUVData(sources[i],'UVSRT',original.disk,1).zap()
		
		uvdata = AIPSUVData(sources[i],'MULTI',original.disk,1)

		# Part 14: Re-index so multi cannot interfere either		
		uvdata.zap_table('CL',1)
		uvdata.zap_table('NX',1)
		indxr = AIPSTask('INDXR')
		indxr.indata = uvdata
		indxr.go()
		
	
	# Part 15: DBCON the data together one sources at a time
	uvdata = AIPSUVData(sources[0],'MULTI',original.disk,1)
	uvdata.rename(name='ALL_AO',klass='UV',seq=1)
        if len(sources) == 1:
		print 'Find your flagged file in AIPS no.:'+str(AIPS.userno)+' As ALL_AO.UV'
		sys.exit()
	for i in range(len(sources)-1):
		dbcon = AIPSTask('DBCON')
		dbcon.indata = AIPSUVData('ALL_AO','UV',original.disk,i+1)
		dbcon.in2data = AIPSUVData(sources[i+1],'MULTI',original.disk,1)
		dbcon.outdata = AIPSUVData('ALL_AO','UV',original.disk,i+2)
		dbcon.fqcenter = -1
		dbcon.go()
		if intermediate == True:
			AIPSUVData(sources[i+1],'MULTI',original.disk,1).zap()
	# Part 15: Sort the dbconned data set & index again and all should be good followed by a final renaming to the original file with increased sequence number
	uvdata = AIPSUVData('ALL_AO','UV',original.disk,len(sources))
	uvdata.zap_table('CL',1)
	uvdata.zap_table('NX',1)
	uvsrt = AIPSTask('UVSRT')
	uvsrt.indata = uvdata
	uvsrt.sort = 'TB'
	uvsrt.go()
    	if intermediate == True:
		for i in range(len(sources)):
        		AIPSUVData('ALL_AO','UV',original.disk,i+1).zap()
	
	indxr = AIPSTask('INDXR')
	indxr.indata = AIPSUVData('ALL_AO','UVSRT',original.disk,1)
	indxr.go()
	
	AIPSUVData('ALL_AO','UVSRT',original.disk,1).rename(name=original.name,klass=original.klass,seq=(original.seq +1))

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

#Inputs for AOmodule
AIPS.userno = 1000                                  #AIPS usernumber the file you want to flag is located
uvdata = AIPSUVData('ALL','DQUAL',1,1)            #UV data file you want to flag in Parseltongue format i.e. AIPSUVData('NAME','CLASS',DISK,SEQ) 
AOsources = ['1241+6020','1236+6212','1331+305','1407+284','0319+415']                            #Python list of the sources that you wish to flag with AOFlagger
NonAOsources = []              #Python list of sources you don't want to flag (AOFLagger isnt great for bright sources)
sensitivity = [3,2]                                 #Sensitivity profile (2.5 is normally a good bet) used in flagging corresponding element-wise to AOsources input e.g. 1241+602 will be flagged with 3 sensitivity  
profile = 'best'                                    #Profile to use in flagging, it is best to leave it as best!
casa = '/usr/local/casa-release-4.6.0-el6/bin/casa' #Path to the casa executable!
ref_freq = 1.2545250E+09                            #Frequency of the observations at channel 1, IF 1!
intermediate_files = False                           #Delete intermediate working files if set to True
rfistrategy = 'strategy3.rfis'

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

AOFlag_AIPS(uvdata,sensitivity,profile,casa,ref_freq,AOsources,NonAOsources,intermediate_files,rfistrategy)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


