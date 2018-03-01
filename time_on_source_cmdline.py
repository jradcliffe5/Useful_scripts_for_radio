import os, sys, re
import numpy as np
import time as timeUtilities
import casadef

casaVersionWithMSMD = '4.1.0'
casaVersionWithMSMDFieldnames = '4.5'
casaVersionWithUndefinedFrame = '4.3.0'
casaAvailable = True
h=6.6260755e-27
k=1.380658e-16
c=2.99792458e10
c_mks=2.99792458e8
jy2SI=1.0e-26
jy2cgs=1.0e-23
pc2cgs=3.0857e18
au2cgs=1.4960e13
solmass2g=1.989e33
earthmass2g=5.974e27
radiusEarthMeters=6371000
solLum2cgs = 3.826e33
mH = 1.673534e-24
G  = 6.67259e-8
Tcmb = 2.725
TROPICAL = 1
MID_LATITUDE_WINTER = 2
MID_LATITUDE_SUMMER = 3
ALMA_LONGITUDE=-67.754748 # ICT-4143,  -67.754694=JPL-Horizons,  -67.754929=CASA observatories
ALMA_LATITUDE=-23.029211  # ICT-4143,  -23.029167=JPL-Horizons,  -23.022886=CASA observatories
ARCSEC_PER_RAD=206264.80624709636
JPL_HORIZONS_ID = {'ALMA': '-7',
                   'VLA': '-5',
                   'GBT': '-9',
                   'MAUNAKEA': '-80',
                   'OVRO': '-81',
                   'geocentric': '500'
}
def ComputeLST(mjdsec=None, longitude=ALMA_LONGITUDE, ut=None, hms=False,
               observatory=None, date=None, verbose=False, prec=0):
    """
    Computes the LST (in hours) for a specified time and longitude.
    The input longitude is in degrees, where east of Greenwich is > 0.
    Two options to specify the time:
    mjdsec: MJD seconds
    ut: ut time on the current day as a HH:MM:SS string or floating point hours
    hms: if True, return the time as a HH:MM:SS string
    observatory: if specified, ignore the longitude argument
    prec: digits of fractional seconds to show
    date: a date/time string
    Either of these formats is valid: 2011/10/15 05:00:00
                                      2011/10/15-05:00:00
                                      2011-10-15 05:00:00
                                      2011-10-15T05:00:00
                                      2011-Oct-15T05:00:00
    If mjdsec, ut and date are all None, then it uses the current date+time.
    -- Todd Hunter
    """
    if (mjdsec == None and ut==None and date==None):
        mjdsec = getMJDSec()
    elif (ut is not None):
        if (date == None):
            datestring = getCurrentDate()
        else:
            datestring = date
        if (type(ut) == str):
            datestring += ' ' + ut
        else:
            minutes = int(60*(ut-int(ut)))
            datestring += ' %02d:%02d:%02.0f' % (int(ut),minutes,3600*(ut-int(ut)-minutes/60.))
        mjdsec = dateStringToMJDSec(datestring, verbose=False)
    elif (date is not None):
        mjdsec = dateStringToMJDSec(date, verbose=False)
    if verbose: print "MJD seconds = ", mjdsec
    JD = mjdToJD(mjdsec/86400.)
    T = (JD - 2451545.0) / 36525.0
    sidereal = 280.46061837 + 360.98564736629*(JD - 2451545.0) + 0.000387933*T*T - T*T*T/38710000.

    # now we have LST in Greenwich, need to scale back to site
    if (observatory is not None):
        longitude = getObservatoryLongitude(observatory)
        if (longitude == None): return
    sidereal += longitude
    sidereal /= 360.
    sidereal -= np.floor(sidereal)
    sidereal *= 24.0
    if (sidereal < 0):
        sidereal += 24
    if (sidereal >= 24):
        sidereal -= 24
    if (hms):
        return(hoursToHMS(sidereal, prec=prec))
    else:
        return(sidereal)

def lstRange(vis, verbose=True, vm='',intent='OBSERVE_TARGET#ON_SOURCE'):
  """
  Compute the LST of start and end of the specified ms.
  For further help and examples, see http://casaguides.nrao.edu/index.php?title=Lstrange
  -- Todd Hunter
  """
  return(lstrange(vis,verbose=verbose,vm=vm,intent=intent))

def lstrange(vis, verbose=True, vm='', intent='OBSERVE_TARGET#ON_SOURCE', fieldID=-1):
  """
  Compute the LST of start and end of the specified ms.
  intent: use only the scans observed with the specified intent
  fieldID: for computing the hour angle, use this field ID rather than the first with the
           specified intent
  -- Todd Hunter
  """
  if (os.path.exists('%s/table.dat'%vis)==False):
      print "Could not find %s/table.dat, are you sure this is an ms?" % (vis)
      return
  observatory = getObservatoryName(vis)
  if (observatory == ''):
      observatory = 'ALMA'
      print "Assuming ", observatory

  if (casadef.casa_version >= casaVersionWithMSMD):
      mymsmd = createCasaTool(msmdtool)
      mymsmd.open(vis)
      uniqueScans = mymsmd.scannumbers()
  else:
      if (vm == ''):
          print "Running ValueMapping... (this may take a minute)"
          vm = ValueMapping(vis)
      uniqueScans = np.unique(vm.scans)

# This way is really slow:
#  onsourceScans = np.unique(vm.getScansForIntent(intent))
  if (verbose):
      print "Checking intents for each scan..."
  if (casadef.casa_version >= casaVersionWithMSMD):
      intents = mymsmd.intents()
      foundIntent = False
      for i in intents:
          if (i.find(intent) >= 0):
              foundIntent = True
              break
      if (foundIntent):
          onsourceScans = mymsmd.scansforintent(intent)
      else:
          onsourceScans = []
  else:
      onsourceScans = getScansForIntentFast(vm,uniqueScans,intent)
  if (len(onsourceScans) < 1):
      # try the "old school" format
      intent = intent.replace('#','.')
      if (casadef.casa_version >= casaVersionWithMSMD):
          intents = mymsmd.intents()
          foundIntent = False
          for i in intents:
              if (i.find(intent) >= 0):
                  foundIntent = True
                  break
          if (foundIntent):
              onsourceScans = mymsmd.scansforintent(intent)
          else:
              onsourceScans = []
      else:
          onsourceScans = getScansForIntentFast(vm,uniqueScans,intent)
  if (verbose):
        print "%s scans = " % (intent), onsourceScans
  i = 0
  wikiline = ''
  wikiline2 = ''
  wikiline3 = 'no onsource time | '
  # First, examine the list of all scans (where i will be 0), then
  # the list of the on-target scans (where i will be 1).
  for scans in [uniqueScans,onsourceScans]:
      if (len(scans) > 0):
          [latitude, longitude, obs] = getObservatoryLatLong(observatory)
          if (casadef.casa_version >= casaVersionWithMSMD):
              times = mymsmd.timesforscans(scans)
              mjdsecmin = np.min(times)
              mjdsecmax = np.max(times)
              mjds = times/86400.
              intMjds = np.int32(mjds)
              uniqueMjds = np.unique(intMjds)
              LSTs = []
              for mjd in uniqueMjds:
                  idx = np.where(mjd == intMjds)[0]
                  mjdsecmin = np.min(times[idx])
                  mjdsecmax = np.max(times[idx])
                  LSTs.append(ComputeLST(mjdsecmin, longitude))
                  LSTs.append(ComputeLST(mjdsecmax, longitude))
          else:
              times = vm.getTimesForScans(scans)
              mjdsecmin = 1e12
              mjdsecmax = 0
              mjds = []
              LSTs = []
              for t in times:
    #  #  #  This is too slow:
    #  #  #        mjdsecmin = np.amin([np.amin(t),mjdsecmin])
    #  #  #        mjdsecmax = np.amax([np.amax(t),mjdsecmax])
    #  #  #  Assume the times are in ascending order:
                  mjdsecmin = np.amin([t[0],mjdsecmin])
                  mjdsecmax = np.amax([t[-1],mjdsecmax])
                  mjds += np.unique(np.int32(t/86400.))
                  LSTs.append(ComputeLST(mjdsecmin, longitude))
                  LSTs.append(ComputeLST(mjdsecmax, longitude))
              uniqueMjds = np.unique(mjds)
#          print "unique MJDs = ", uniqueMjds
          LST = np.zeros(2)
          LST[0] = np.min(LSTs)
          LST[1] = np.max(LSTs)
          if (i == 1):
              style = "on source"
          else:
              style = "of whole SB"
          duration = LST[1]-LST[0]
          duration2 = LST[1]+24-LST[0]
          if (duration2 < duration):
              duration = duration2
          if (verbose):
              print "LST range %s = %.2f to %.2f = %02d:%02d to %02d:%02d (%.1f minutes)" % (style,LST[0],LST[1],
                  np.floor(LST[0]),np.floor(60*(LST[0]-np.floor(LST[0]))), np.floor(LST[1]),
                  np.floor(60*(LST[1]-np.floor(LST[1]))), duration*60)
              if (i==1):
                  # get RA of science field
                  if (fieldID == -1):
                      if (casadef.casa_version >= casaVersionWithMSMD):
                          fieldID = mymsmd.fieldsforintent(intent)[0]
                      else:
                          fieldname = getFieldsForTime(vm.fieldsForTimes,times[0][0])
                          fieldID = vm.getFieldIdsForFieldName(fieldname)
                  RA, DEC = getRADecForField(vis, fieldID, forcePositiveRA=True, usemstool=True)
                  RA = RA*12/np.pi
                  HA = np.zeros(2)
                  HA[0] = LST[0]-RA
                  HA[1] = LST[1]-RA
                  HA0sign = ("%+f" % HA[0])[:1]
                  HA1sign = ("%+f" % HA[1])[:1]
                  aHA = np.zeros(2)
                  aHA[0] = abs(HA[0])
                  aHA[1] = abs(HA[1])
                  print "Field %d RA = %fh" % (fieldID, RA)
                  print " HA range %s = %+.2f to %+.2f = %c%02d:%02d to %c%02d:%02d" % (style,HA[0],HA[1],
                      HA0sign, np.floor(aHA[0]), np.floor(60*(aHA[0]-np.floor(aHA[0]))),
                      HA1sign, np.floor(aHA[1]), np.floor(60*(aHA[1]-np.floor(aHA[1]))))
          [mjdmin,utmin] = mjdSecondsToMJDandUT(mjdsecmin)
          [mjdmax,utmax] = mjdSecondsToMJDandUT(mjdsecmax)
          if (i==0):
              if (casadef.casa_version >= casaVersionWithMSMD):
                  alltimes = times
              else:
                  alltimes = np.array([val for subl in times for val in subl])
              deltaT = 0.5*(abs(np.min(alltimes[np.where(alltimes > mjdsecmin)]) - mjdsecmin) +
                            abs(np.max(alltimes[np.where(alltimes < mjdsecmax)]) - mjdsecmax))
              clockTimeMinutes = (mjdmax - mjdmin + deltaT/86400.)*1440.
          if (verbose):
              print "MJD range %s = %.4f to %.4f" % (style, mjdmin, mjdmax)
              meanJD = mjdToJD(0.5*(mjdmin+mjdmax))
              print "Mean JD %s = %f" % (style, meanJD)
              print " UT range %s = %s to %s" % (style, utmin, utmax)
          tb.open(vis+'/OBSERVATION')
          sbname = 'unknown'
          exec_uid = 'unknown'
          if ('SCHEDULE' in tb.colnames()):
              if (tb.iscelldefined('SCHEDULE',0)):
                  sched = tb.getcol('SCHEDULE')
                  sbname = '%s' % (sched[0][0].split()[1])  # This is the SB UID.
                  exec_uid = '%s' % (sched[1][0].split()[1])
          tb.close()
          if (i==0):
              wikiline2 += "| %s | %s | %s | %s-%s | %02d:%02d-%02d:%02d | %.1f | " % (utmin[0:10],sbname,exec_uid,utmin[10:-6],utmax[11:-6],np.floor(LST[0]),np.floor(60*(LST[0]-np.floor(LST[0]))), np.floor(LST[1]), np.floor(60*(LST[1]-np.floor(LST[1]))), clockTimeMinutes)
              csvline = "%s,%s,%.2f" % (sbname, exec_uid, clockTimeMinutes)
              wikiline += "%s-%s | %02d:%02d-%02d:%02d | %.1f |" % (utmin[0:-6],utmax[11:-6],np.floor(LST[0]),np.floor(60*(LST[0]-np.floor(LST[0]))), np.floor(LST[1]), np.floor(60*(LST[1]-np.floor(LST[1]))),clockTimeMinutes)
          else:
              # print out elevation range for mjdsecmin to mjdsecmax
              # could use TsysExplorer(vis) but it fails at line 1144
              tb.open("%s/POINTING" % vis)
              azel = 1
              try:
                  elevation = np.transpose(tb.getcol("DIRECTION")[azel])
                  elevTime  = tb.getcol("TIME")
                  tb.close()
                  if (casadef.casa_version >= casaVersionWithMSMD):
                      t = mymsmd.timesforscan(scans[0])[0]
                  else:
                      t = vm.getTimesForScans(scans[0])[0]
                  matches1 = np.where(elevTime > np.min(t[0]))[0]
                  matches2 = np.where(elevTime < np.max(t[-1]))[0]
                  matches = np.intersect1d(matches1,matches2)
                  startElev = elevation[matches[0]]*180/math.pi
                  if (casadef.casa_version >= casaVersionWithMSMD):
                      t = mymsmd.timesforscan(scans[0])[0]
                  else:
                      t = vm.getTimesForScans(scans[-1])[0]
                  matches1 = np.where(elevTime > np.min(t[0]))[0]
                  matches2 = np.where(elevTime < np.max(t[-1]))[0]
                  matches = np.intersect1d(matches1,matches2)
                  stopElev = elevation[matches[-1]]*180/math.pi
                  if (verbose):
                      print "Elevation range on %s scans = %.1f-%.1f" % (intent, startElev,stopElev)
                  wikiline += "%02d:%02d-%02d:%02d | %.0f-%.0f | " % (np.floor(LST[0]),np.floor(60*(LST[0]-np.floor(LST[0]))), np.floor(LST[1]), np.floor(60*(LST[1]-np.floor(LST[1]))), startElev, stopElev)
                  wikiline3 = "%.0f-%.0f | " % (startElev,stopElev)
              except:
                  wikiline3 = "pointing table empty | "
                  if (verbose):
                      print "The pointing table appears to be empty.  Was it deleted because this is a mosaic?"
                  wikiline += "%02d:%02d-%02d:%02d | " % (np.floor(LST[0]),np.floor(60*(LST[0]-np.floor(LST[0]))), np.floor(LST[1]), np.floor(60*(LST[1]-np.floor(LST[1]))))
          i += 1
  if (verbose):
      print "wikiline = %s" % (wikiline)
  csvline = 'csvline = %s' % (csvline)
  if (casadef.casa_version >= casaVersionWithMSMD):
      mymsmd.close()
  return (wikiline2, wikiline3,clockTimeMinutes,csvline)
  # end of lstrange()
def getObservatoryLatLong(observatory='',verbose=False, radians=False):
     """
     Opens the casa table of known observatories and returns the latitude and longitude
     in degrees for the specified observatory name string.
     observatory: string name, JPL integer, or integer string  (e.g. 'ALMA' == -7)
     -- Todd Hunter
     """
     if (observatory == ''):
        listAvailableObservatories()
        return
     repotable = os.getenv("CASAPATH").split()[0]+"/data/geodetic/Observatories"
     try:
        tb.open(repotable)
     except:
        print "Could not open table = %s, returning ALMA coordinates in au instead" % (repotable)
        longitude = ALMA_LONGITUDE
        latitude = ALMA_LATITUDE
        observatory = 'ALMA'
        if radians:
            latitude = np.radians(latitude)
            longitude = np.radians(longitude)
        return([latitude,longitude,observatory])
     if (type(observatory) == 'int' or str(observatory) in JPL_HORIZONS_ID.values()):
         if (str(observatory) in JPL_HORIZONS_ID.values()):
             observatory = JPL_HORIZONS_ID.keys()[JPL_HORIZONS_ID.values().index(str(observatory))]
             if (verbose):
                 print "Recognized observatory = %s" % observatory
             if (observatory == 'MAUNAKEA'): observatory = 'SMA'
             if (observatory == 'OVRO'): observatory = 'OVRO_MMA'
         else:
             print "Did not recognize observatory='%s' in %s, using ALMA instead." % (observatory,str(JPL_HORIZONS_ID.values()))
             observatory = 'ALMA'
#     else:
#         print "%s is not in %s" % (observatory, str(JPL_HORIZONS_ID.values()))

     Name = tb.getcol('Name')
     matches = np.where(np.array(Name)==observatory)
     if (len(matches) < 1 and str(observatory).find('500') < 0):
         print "Names = ", Name
         print "Did not find observatory='%s', using ALMA value in au instead." % (observatory)
         for n in Name:
             if (n.find(observatory) >= 0):
                 print "Partial match: ", n
         observatory = 'ALMA'
         longitude = ALMA_LONGITUDE
         latitude = ALMA_LATITUDE
     elif (str(observatory).find('500') >= 0 or
           str(observatory).lower().find('geocentric') >= 0):
         observatory = 'Geocentric'
         longitude = 0
         latitude = 0
     else:
         longitude = tb.getcol('Long')[matches[0]]
         latitude = tb.getcol('Lat')[matches[0]]
     tb.close()
     if radians:
         latitude = np.radians(latitude)
         longitude = np.radians(longitude)
     return([latitude,longitude,observatory])

def call_qa_time(arg, form='', prec=0, showform=False):
    """
    This is a wrapper for qa.time(), which in casa 4.0.0 returns a list
    of strings instead of just a scalar string.
    arg: a time quantity
    - Todd Hunter
    """
    if (type(arg) == dict):
        if (type(arg['value']) == list or
            type(arg['value']) == np.ndarray):
            if (len(arg['value']) > 1):
                print "WARNING: call_qa_time() received a dictionary containing a list of length=%d rather than a scalar. Using first value." % (len(arg['value']))
            arg['value'] = arg['value'][0]
    result = qa.time(arg, form=form, prec=prec, showform=showform)
    if (type(result) == list or type(result) == np.ndarray):
        return(result[0])
    else:
        return(result)

def computeDurationOfScan(scanNumber,t=None, vis=None, returnSubscanTimes=False,
                          verbose=False, gapFactor=None, includeDate=False,
                          mymsmd=None, scienceSpwsOnly=False):
    """
    This function is used by timeOnSource() to empirically determine the number
    of subscans and the total duration of all the subscans of a particular scan.
    Inputs:
    scanNumber: the scan number, simply for generating a file list of timestamps
    t = a sequence of integration timestamps (optional for casa >= 4.1.0)
    vis = the measurement set (not necessary if t is given)
    mymsmd = an msmd instance (as an alternative to vis)
    gapFactor: default=2.75 for integrations<1sec, 2.0 otherwise
    Returns:
    1) duration in seconds
    2) the number of subscans
    and if returnSubscanTimes==True
    3) the begin/end timestamps of each subscan
    4) the begin/end timestampStrings of each subscan
    -- Todd Hunter
    """
    if (t is None and vis is None and (mymsmd is None or mymsmd=='')):
        print "You must specify either vis, mymsmd or t."
        return
    keepmymsmd = False
    if (t is None):
        if (casadef.casa_version < casaVersionWithMSMD):
            print "For this version of casa, you must specify t rather than vis."
            return
        if (mymsmd is None or mymsmd == ''):
            mymsmd = createCasaTool(msmdtool)
            mymsmd.open(vis)
        else:
            keepmymsmd = True
    elif (mymsmd is None or mymsmd==''):  # fix for CSV-2903 on 29-Apr-2016
        mymsmd = createCasaTool(msmdtool)
        mymsmd.open(vis)
    else:
        keepmymsmd = True
    if (t is None):
        t = pickTimesForScan(mymsmd, scanNumber, scienceSpwsOnly=scienceSpwsOnly)
    else:
        t = np.unique(t)
    if (len(t) <= 1):
        if (casadef.casa_version < casaVersionWithMSMD):
            print "This version of CASA is too old for this function to handle single-dump integrations."
            return(0,0)
        else:
            if (scanNumber==1):
                duration = np.min(pickTimesForScan(mymsmd, scanNumber+1, scienceSpwsOnly=scienceSpwsOnly)) - t[0]
            else:
                duration = t[0] - np.max(pickTimesForScan(mymsmd, scanNumber-1, scienceSpwsOnly=scienceSpwsOnly))
            if (not keepmymsmd): mymsmd.close()
            if (returnSubscanTimes):
                timestampsString = mjdsecToTimerange(t[0]-0.5*duration,t[0]+0.5*duration,
                                                     decimalDigits=1,includeDate=includeDate)
                return(duration,1,{0:t},{0:timestampsString})
            else:
                return(duration,1)
    else:
        if (mymsmd is not None and mymsmd != '' and not keepmymsmd):
            mymsmd.close()
        d = np.max(t) - np.min(t)
        # initial estimate for interval
        diffs = np.diff(t)
        avgInterval = np.median(diffs)
        startTime = previousTime = t[0]
        subscans = 1
        if (gapFactor is None):
            if (avgInterval > 1):
                gapFactor = 2
            else:
                gapFactor = 2.75 # it was 3 for a long time, but failed on 2013-01-24 dataset (2014-09-23)
        startTime = previousTime = t[0]
        duration = 0
        tdiffs = []
        s = ''
        timestamps = {}
        timestampsString = {}
        droppedTimeTotal = 0
        daygaps = 0
        for i in range(1,len(t)):
            s += "%.2f " % (t[i]-t[0])
            tdiff = t[i]-previousTime
            tdiffs.append(tdiff)
            if (tdiff > gapFactor*avgInterval):
                droppedTime = t[i]-previousTime+avgInterval
                droppedTimeTotal += droppedTime
                if droppedTime > 12*3600:
                    daygaps += droppedTime
                if (verbose):
                    print "    ***********************"
                    print "    i=%d, Dropped %.1f seconds" % (i,droppedTime)
                subscanLength = t[i-1] - startTime + avgInterval
                duration += subscanLength
                if (subscanLength > 1.5*avgInterval):
                    # Don't count single point subscans because they are probably not real.
                    timestamps[subscans] = [startTime,t[i-1]]
                    timestampsString[subscans] = mjdsecToTimerange(startTime,t[i-1],decimalDigits=1,includeDate=includeDate)
                    if (verbose):
                        print "Scan %d: Subscan %d: %s, duration=%f" % (scanNumber,subscans,s,subscanLength)
                    s = ''
                    subscans += 1
                elif (verbose):
                    print "Scan %d: dropped a dump of length %f because it was between subscans" % (scanNumber,subscanLength)
                startTime = t[i]
            previousTime = t[i]
        if (droppedTimeTotal > 0 and verbose):
            print "+++++ Scan %d: total dropped time = %.1f seconds = %.1f minutes" % (scanNumber,droppedTimeTotal, droppedTimeTotal/60.)
        if daygaps > 0:
            print "+++++ Scan %d: large time gaps: %.1f sec = %.1f hrs = %.1f days" % (scanNumber,daygaps,daygaps/3600.,daygaps/86400.)
        duration += t[len(t)-1] - startTime
        timestamps[subscans] = [startTime,t[len(t)-1]]
        timestampsString[subscans] = mjdsecToTimerange(startTime,t[len(t)-1],decimalDigits=1,includeDate=includeDate)
    if (returnSubscanTimes):
        return(duration, subscans, timestamps, timestampsString)
    else:
        return(duration, subscans)
def mjdToJD(MJD=None):
    """
    Converts an MJD value to JD.  Default = now.
    """
    if (MJD==None): MJD = getMJD()
    JD = MJD + 2400000.5
    return(JD)
def mjdSecondsToMJDandUT(mjdsec, use_metool=True, debug=False, prec=6, delimiter='-'):
    """
    Converts a value of MJD seconds into MJD, and into a UT date/time string.
    prec: 6 means HH:MM:SS,  7 means HH:MM:SS.S
    example: (56000.0, '2012-03-14 00:00:00 UT')
    Caveat: only works for a scalar input value
    Todd Hunter
    """
    if (not casaAvailable or use_metool==False):
        mjd = mjdsec / 86400.
        jd = mjdToJD(mjd)
        trialUnixTime = 1200000000
        diff  = ComputeJulianDayFromUnixTime(trialUnixTime) - jd
        if (debug): print "first difference = %f days" % (diff)
        trialUnixTime -= diff*86400
        diff  = ComputeJulianDayFromUnixTime(trialUnixTime) - jd
        if (debug): print "second difference = %f seconds" % (diff*86400)
        trialUnixTime -= diff*86400
        diff  = ComputeJulianDayFromUnixTime(trialUnixTime) - jd
        if (debug): print "third difference = %f seconds" % (diff*86400)
        # Convert unixtime to date string
        utstring = timeUtilities.strftime('%Y'+delimiter+'%m'+delimiter+'%d %H:%M:%S UT',
                                          timeUtilities.gmtime(trialUnixTime))
    else:
        me = createCasaTool(metool)
        today = me.epoch('utc','today')
        mjd = np.array(mjdsec) / 86400.
        today['m0']['value'] =  mjd
        hhmmss = call_qa_time(today['m0'], prec=prec)
        date = qa.splitdate(today['m0'])
        utstring = "%s%s%02d%s%02d %s UT" % (date['year'],delimiter,date['month'],delimiter,
                                             date['monthday'],hhmmss)
    return(mjd, utstring)
def pickScansForField(mymsmd, field):
    """
    Chooses whether to run msmd.scansforfield() or the workaround for the
    bug in CASA 4.4 (CAS-7622).
    field: should be an integer
    -Todd Hunter
    """
    if (casadef.casa_version >= '4.4' and casadef.casa_version < '4.5'):
        s = getScansForField(mymsmd, field)
    else:
        s = mymsmd.scansforfield(field)
    return(s)

def pickTimesForScan(mymsmd, scan, scienceSpwsOnly=True, useTimesForSpwsIfAvailable=False):
    """
    Chooses whether to run msmd.timesforscan() or the aU workaround for the
    bug in CASA 4.4 (CAS-7622).  Called only by computeDurationOfScan.
    -Todd Hunter
    """
    if (casadef.casa_version >= '4.4' and casadef.casa_version < '4.5'):
        t = getTimesForScan(mymsmd, scan)
    else:
        t = mymsmd.timesforscan(scan)
        if scienceSpwsOnly:
            scienceSpws = getScienceSpws(mymsmd.name(), mymsmd=mymsmd, returnString=False)
            if len(scienceSpws) == 0:
                print "There are no science spws, so considering all spws instead."
                scienceSpwsOnly = False
        if scienceSpwsOnly:
            if useTimesForSpwsIfAvailable:
                try:
                    # available starting in 5.1.0-46
                    print "Calling msmd.timesforspws(%s)" % (str(scienceSpws[0]))
                    times = mymsmd.timesforspws(scienceSpws[0])
                    t = np.intersect1d(times, t)
                except:
                    mytb = createCasaTool(tbtool)
                    mytb.open(mymsmd.name())
                    print "Restricting times to first spw: %d" % (scienceSpws[0])
                    myt = mytb.query('DATA_DESC_ID == %d && SCAN_NUMBER == %d' % (scienceSpws[0], scan))
                    times = myt.getcol('TIME')
                    myt.close()
                    mytb.close()
                    idx = np.where((t >= np.min(times)) & (t <= np.max(times)))
                    print "Dropped %d points" % (len(t)-len(idx))
                    t = t[idx]
            else:
                    mytb = createCasaTool(tbtool)
                    mytb.open(mymsmd.name())
                    print "Restricting times to spw %d" % (scienceSpws[0])
                    myt = mytb.query('DATA_DESC_ID == %d && SCAN_NUMBER == %d' % (scienceSpws[0], scan))
                    times = myt.getcol('TIME')
                    myt.close()
                    mytb.close()
                    if False:
                        # This will include SQLD, and was the old version of the code
                        idx = np.where((t >= np.min(times)) & (t <= np.max(times)))
                        print "Dropped %d points" % (len(t)-len(idx))
                        t = t[idx]
                    else:
                        # This will not include SQLD.
                        t = np.unique(times)
    return (t)
def mjdsecToTimerangeComponent(mjdsec, decimalDigits=2, includeDate=True, use_metool=True, debug=False):
    """
    Converts a value of MJD seconds into a UT date/time string suitable for one
    member of the timerange argument in plotms.
    example: '2012/03/14/00:00:00.00'
    input options:
       decimalDigits: how many digits to display after the decimal point
       use_metool: True=use casa tool to convert to UT, False: use formula in aU
    Todd Hunter
    """
    if (not casaAvailable or use_metool==False):
        mjd = mjdsec / 86400.
        jd = mjdToJD(mjd)
        trialUnixTime = 1200000000
        diff  = ComputeJulianDayFromUnixTime(trialUnixTime) - jd
        if (debug): print "first difference = %f days" % (diff)
        trialUnixTime -= diff*86400
        diff  = ComputeJulianDayFromUnixTime(trialUnixTime) - jd
        if (debug): print "second difference = %f seconds" % (diff*86400)
        trialUnixTime -= diff*86400
        diff  = ComputeJulianDayFromUnixTime(trialUnixTime) - jd
        if (debug): print "third difference = %f seconds" % (diff*86400)
        # Convert unixtime to date string
        if (includeDate):
            utstring = timeUtilities.strftime('%Y/%m/%d/%H:%M:%S',
                                              timeUtilities.gmtime(trialUnixTime))
        else:
            utstring = timeUtilities.strftime('%H:%M:%S',
                                              timeUtilities.gmtime(trialUnixTime))
        utstring += '.%0*d'  % (decimalDigits, np.round(10**decimalDigits*(trialUnixTime % 1)))
    else:
        me = createCasaTool(metool)
        today = me.epoch('utc','today')
        mjd = np.array(mjdsec) / 86400.
        today['m0']['value'] =  mjd
        hhmmss = call_qa_time(today['m0'],prec=6+decimalDigits)
        date = qa.splitdate(today['m0'])
        if (includeDate):
            utstring = "%s/%02d/%02d/%s" % (date['year'],date['month'],date['monthday'],hhmmss)
        else:
            utstring = hhmmss
    return(utstring)

def mjdsecToTimerange(mjdsec1, mjdsec2=None, decimalDigits=2, includeDate=True,
                      use_metool=True, debug=False):
    """
    Converts two value of MJD seconds into a UT date/time string suitable for
    the timerange argument in plotms.  They can be entered as two separated values,
    or a single tuple.
    Example output: '2012/03/14/00:00:00.00~2012/03/14/00:10:00.00'
    input options:
       decimalDigits: how many digits to display after the decimal point
       use_metool: True=use casa tool to convert to UT, False: use formula in aU
    -Todd Hunter
    """
    if (type(mjdsec1) == list or type(mjdsec1)==np.ndarray):
        mjdsec2 = mjdsec1[1]
        mjdsec1 = mjdsec1[0]
    return(mjdsecToTimerangeComponent(mjdsec1, decimalDigits, includeDate, use_metool, debug) + '~' \
           + mjdsecToTimerangeComponent(mjdsec2, decimalDigits, includeDate, use_metool, debug))

def createCasaTool(mytool):
    """
    A wrapper to handle the changing ways in which casa tools are invoked.
    Relies on "from taskinit import *" in the preamble above.
    Todd Hunter
    """
    if (type(casac.Quantity) != type):  # casa 4.x
        myt = mytool()
    else:  # casa 3.x
        myt = mytool.create()
    return(myt)

def getObservatoryName(ms):
    """
    Returns the observatory name in the specified ms.
    -- Todd Hunter
    """
    antTable = ms+'/OBSERVATION'
    try:
        mytb = createCasaTool(tbtool)
        mytb.open(antTable)
        myName = mytb.getcell('TELESCOPE_NAME')
        mytb.close()
    except:
        print "Could not open OBSERVATION table to get the telescope name: %s" % (antTable)
        myName = ''
    return(myName)


class ValueMapping:
    """
    Input: The name of an MS dataset as a string.
    Purpose: This class provides details on the mapping of various parameters to each other.  For example, if you would like to
             know which scans observed a given field source or over which time interval a scan was executed, this is the place to look.
             Included in that are functions which map antenna name to antenna id and field name to field id.  This is useful in building
             other routines that allow you to not require the user to input one or the other type of information.  It also gets unique
             lists of items, like antennas, field sources, scans, intents, etc.

    Responsible: S. Corder and other contributors
    Example: vm = aU.ValueMapping('myInputMS.ms')
    Suggested Improvements:
          (done, 06-04-2011, scorder)1) Change some of the get methods to do methods because they aren't really returning anything
          2) Add spectral window mapping information, spectral windows, spectral windows to fields, spectral windows to scans,
             spectral windows to frequency (central and channelized): Basically make a dictionary of all the spectral line info stuff,
             per spectral window.  Rework combination of SensitivityCalculator and VM....
          3) Add integration time calculator per source
          4) Do sensitivity calculator (maybe needs to be separate function/class that inhereits this)
    """
    def __init__(self,inputMs):
        """
        Instantiation of this class calls this, i.e., vm = aU.ValueMapping('myInputMS.ms').  The dataset name is the only allowed
        input.  It generates the mappings are part of this instantiation.
        """
        self.inputMs = inputMs
        self.setInputMs(self.inputMs)

    def getInputMs(self):
        """
        Input: None
        Output: The active measurement set in the class as a string
        Responsible: S. Corder
        Purpose: Return the name of the active dataset
        """
        return self.inputMs

    def setInputMs(self,inputMs):
        """
        Input: New measurement set that you wish to become the active one, as a string and activate that change to all other parameters
        Output: None
        Responsible: S. Corder
        Purpose: This changes the active dataset and remakes the relevant mappings.  The order of the functions is very
                 important.
        """
        self.inputMs = inputMs
        self.doScanTimeMapping()
        self.doFieldsAndSources()   ;  self.doStatesAndIntents()
        self.doAntennasAndNames()
        self.doScanStateMapping()   ;  self.doFieldTimeMapping() ;
        self.doAntennaTimeMapping() ;  self.doAntennaStateMapping()
        self.doDataDescId()         ; self.doSpwAndDataDescId()
        self.doPolarizations()
        self.doSpwAndFrequency()
        self.doSpwScanMapping()     ;  self.doSpwFieldMapping()
        self.doSpwIntentMapping()

    def doSpwAndFrequency(self,ignoreWVR=True) :
        """
        Input: None
        Output: None
        Responsible: S. Corder
        Purpose: Creates a dictionary (spwInfo) of spectral windows, with keys of the spectral window number.
                 For each spectral window, another dictionary is formed that has keys of bandwidth, sideband,
                 chanFreqs, chanWidth, numChannels, and meanFreq.
        """
        self.spwInfo = {}
        tb.open("%s/SPECTRAL_WINDOW" % self.inputMs)
        specWinIds = range(tb.nrows())
        junk = []
        for i in specWinIds :
            self.spwInfo[i] = {}
            self.spwInfo[i]["bandwidth"] = tb.getcell("TOTAL_BANDWIDTH",i)
            self.spwInfo[i]["chanFreqs"] = tb.getcell("CHAN_FREQ",i)
            self.spwInfo[i]["chanWidth"] = tb.getcell("CHAN_WIDTH",i)[0]
            self.spwInfo[i]["edgeChannels"] = [min(self.spwInfo[i]["chanFreqs"]),max(self.spwInfo[i]["chanFreqs"])]
            netSideband  = tb.getcell("NET_SIDEBAND",i)
            if netSideband == 2 : self.spwInfo[i]["sideband"] = 1
            else : self.spwInfo[i]["sideband"] = -1
            self.spwInfo[i]["meanFreq"]  = self.spwInfo[i]["chanFreqs"].mean()
            self.spwInfo[i]["numChannels"] = self.spwInfo[i]["chanFreqs"].shape[0]
            if ((ignoreWVR) and (self.spwInfo[i]['numChannels'] == 4)) :
                junk.append(i)
                self.spwInfo.pop(i)
        tb.close()
        if ignoreWVR:
            if (len(junk) > 0):
                print "Ignoring spectral window %s because it is WVR related" % junk

    def doSpwAndDataDescId(self) :
        tb.open("%s/DATA_DESCRIPTION" % self.inputMs)
        self.spwForDataDescId = tb.getcol('SPECTRAL_WINDOW_ID')
        tb.close()

    def doSpwFieldMapping(self) :
        tb.open("%s" % self.inputMs)
        self.fieldsForSpw = {}
        for i in self.spwForDataDescId :
            spw = self.spwForDataDescId[i]
            indices = np.where(self.dataDescId == i)
            self.fieldsForSpw[spw] = np.unique(self.fields[indices])
        tb.close()
        return

    def doDataDescId(self) :
        tb.open("%s" % self.inputMs)
        self.dataDescId = tb.getcol('DATA_DESC_ID')
        tb.close()

    def doPolarizations(self) :
        # Determine the number of polarizations for the first OBSERVE_TARGET intent.
        # Used by plotbandpass for BPOLY plots since the number of pols cannot be inferred
        # correctly from the caltable alone.  You cannot not simply use the first row, because
        # it may be a pointing scan which may have different number of polarizations than what
        # the TARGET and BANDPASS calibrator will have.
        # -- T. Hunter
        myscan = -1
        starttime = timeUtilities.time()
        for s in self.uniqueScans:
            intents = self.getIntentsForScan(s)
            for i in intents:
                if (i.find('OBSERVE_TARGET')>=0):
                    myscan = s
#                    print "First OBSERVE_TARGET scan = ", myscan
                    break
            if (myscan >= 0):
                break
        if (myscan == -1):
            # if there is no OBSERVE_TARGET, then just use the first scan
            myscan = 0
        self.getDataColumnNames()
        tb.open("%s" % self.inputMs)
        if (myscan == 0):
            # assume the first row in the table is for the first scan, to save time
            self.nPolarizations = np.shape(tb.getcell(self.dataColumnName,0))[0]
        else:
            scans = tb.getcol('SCAN_NUMBER')
            self.nPolarizations = 0
            for s in range(len(scans)):
                if (scans[s]==myscan):
                    self.nPolarizations = np.shape(tb.getcell(self.dataColumnName,s))[0]
                    break
        tb.close()
        donetime = timeUtilities.time()
#        print "doPolarizations took %.1f sec" % (donetime-starttime)

    def getDataColumnNames(self):
        tb.open(self.inputMs)
        colnames = tb.colnames()
        self.correctedDataColumnName = ''
        self.modelDataColumnName = ''
        if 'FLOAT_DATA' in colnames:
            self.dataColumnName = 'FLOAT_DATA'
            self.correctedDataColumnName = 'FLOAT_DATA'
        elif 'DATA' in colnames:
            self.dataColumnName = 'DATA'
        if 'CORRECTED_DATA' in colnames:
            self.correctedDataColumnName = 'CORRECTED_DATA'
        if 'MODEL_DATA' in colnames:
            self.modelDataColumnName = 'MODEL_DATA'
        tb.close()
        return

    def doSpwScanMapping(self) :
        tb.open("%s" % self.inputMs)
        self.scansForSpw = {}
        for i in self.spwForDataDescId :
            spw = self.spwForDataDescId[i]
            indices = np.where(self.dataDescId == i)
            self.scansForSpw[spw] = np.unique(self.scans[indices])
        tb.close()
        return

    def doSpwIntentMapping(self) :
        tb.open("%s" % self.inputMs)
        self.intentsForSpw = {}
        for i in self.spwForDataDescId :
            spw = self.spwForDataDescId[i]
            indices = np.where(self.dataDescId == i)
            statesForSpw = np.unique(self.states[indices])
            _intent = []
            for i in statesForSpw :
                __intent = []
# The 'if' statement is needed to support telescopes w/o intents. -T. Hunter
                if (len(self.intentsForStates) > 0):
                  for j in self.intentsForStates[i] :
#                    __map = j.split('#')[0]
                    __map = j
                    __intent.append(__map)
                  _intent += __intent
            self.intentsForSpw[spw] = np.unique(np.array(_intent))
        tb.close()


    def getSpwsForIntent(self,intent) :
        spwsForIntent = []
        for i in self.intentsForSpw.keys() :
            if (intent in self.intentsForSpw[i]) : spwsForIntent.append(i)
        return spwsForIntent

    def getIntentsForSpw(self,spw) :
        return self.intentsForSpw[spw]

    def getSpwsForField(self,field) :
        if not str(field).isdigit() : field = self.getFieldIdsForFieldName(field)
        spwsForField = []
        for i in self.fieldsForSpw.keys() :
            if (field in self.fieldsForSpw[i]) : spwsForField.append(i)
        return spwsForField

    def getFieldsForSpw(self,spw,returnName = True) :
        if returnName :
            return self.getFieldNamesForFieldId(np.unique(np.array(self.fieldsForSpw[spw])))
        else :
            return np.unique(np.array(self.fieldsForSpw[spw]))

    def getSpwsForScan(self,scan):
        spwsForScan = []
        for i in self.scansForSpw.keys() :
            if (scan in self.scansForSpw[i]) : spwsForScan.append(i)
        return spwsForScan

    def getScansForSpw(self,spw) :
        return self.scansForSpw[spw]

    def getAntennaNamesForAntennaId(self,id):
        """
        Input: Antenna id as an integer or string
        Output: Antenna name as a string
        Responsible: S. Corder
        Purpose: Allows translation between antenna id and antenna name.
        """

        return self.antennaNamesForAntennaIds[int(id)]

    def getAntennaIdsForAntennaName(self,antennaName):
        """
        Input: Antenna names as a string
        Output: Antenna index as an integer.
        Responsible: S. Corder
        Purpose: This allows translation between antenna name and antenna id
        """

        return np.where(self.antennaNamesForAntennaIds == antennaName)[0][0]

    def doStatesAndIntents(self):
        """
        Input: None
        Output: None
        Responsible: S. Corder
        Purpose: This function defines two attributes, uniqueStates, which is a python list of the different intents present in the dataset,
                 and intentsForStates is another python list which give the intents for state id as a nested list.  The first index of
                 the intentsForStates is the state id.  If you choose a state id, then the result is a list of intents for that state.
        """
        tb.open("%s/STATE" % self.inputMs)
        intents = tb.getcol("OBS_MODE")
        tb.close()
        _intents = []
        for i in intents : _intents.append(i.split(','))
        self.intentsForStates = _intents
        self.uniqueIntents = []
        for i in self.intentsForStates : self.uniqueIntents.extend(i)
        self.uniqueIntents = np.unique(np.array(self.uniqueIntents))

    def doFieldsAndSources(self):
        """
        Input: None
        Output: None
        Responsible: S. Corder
        Purpose: This function defines two attributes, uniqueField and fieldNamesForFieldIds.  For the time being these are identical.
                 fieldNamesForFieldIds is simply a numpy array where the index is the field id and the value is the name of the field source.
        """
        tb.open("%s/FIELD" % self.inputMs)
        self.fieldNamesForFieldIds = tb.getcol('NAME')
#        print '%d field names = '%len(self.fieldNamesForFieldIds), self.fieldNamesForFieldIds
        self.uniqueFields = self.fieldNamesForFieldIds
        tb.close()

    def doAntennasAndNames(self) :
        """
        Input: None
        Output: None
        Responsible: S. Corder
        Purpose: This function defines two attributes, uniqueAntennas (which is a little excessive) and antennaNamesForAntennaIds.
                 antennaNamesForAntennaIds is a numpy array and has indices that are the antenna ids and values that are the antenna names.
        """
        tb.open("%s/ANTENNA" % self.inputMs)
        self.antennaNamesForAntennaIds = tb.getcol('NAME')
        self.uniqueAntennas = np.unique(self.antennaNamesForAntennaIds)
        self.numAntennas = len(self.uniqueAntennas)
        tb.close()

    def doScanStateMapping(self):
        """
        Input: None
        Output: None
        Responsible: S. Corder
        Purpose: This function defines an attribute, statesForScans, that is the mapping between states and scan numbers.  It is
                 a python dictionary that has keys of the scan number and values, in a list, of the states used in that scan.
        """
        tb.open("%s" % self.inputMs)
        self.states = tb.getcol("STATE_ID")
        tb.close()
        self.statesForScans = {}
        for i in self.uniqueScans :
            indices = np.where(self.scans == i)
            self.statesForScans[i] = np.unique(self.states[indices])

    def doScanTimeMapping(self):
        """
        Input: None
        Output: None
        Responsible: S. Corder
        Purpose: This function defines four attributes, scans, time, uniqueScans and scansForTiems.  scans and time are simply
                 the scan number and time table from the main data table as python arrays.  uniqueScans is a numpy array of the independent scans
                 in the table.  scansForTimes is a python dictionary with keys as scans and values as times in which data was taken
                 for that scan.
        """
        tb.open(self.inputMs)
        self.scans = tb.getcol('SCAN_NUMBER')
        self.time = tb.getcol('TIME')
        tb.close()
        self.uniqueScans = np.unique(self.scans)
        self.scansForTimes = {}
        for i in self.uniqueScans :
            indices = np.where(self.scans == i)
            self.scansForTimes[i] = self.time[indices]

    def doFieldTimeMapping(self):
        """
        Input: None
        Output: None
        Responsible: S. Corder
        Purpose: This function defines two attributes, fields, a numpy array, and fieldsForTimes, a dictionary with keys of field name.
                 The fields is just the field id from the data table.  The values of fieldsForTimes are the times during which data was
                 collected for that field source.
        """

        tb.open(self.inputMs)
        self.fields = tb.getcol('FIELD_ID')
        tb.close()
        self.fieldsForTimes = {}
        for i in range(len(self.fieldNamesForFieldIds)) :
            indices = np.where(self.fields == i)
            self.fieldsForTimes[self.fieldNamesForFieldIds[i]] = self.time[indices]

    def doAntennaTimeMapping(self):
        """
        Input: None
        Output: None
        Responsible: S. Corder
        Purpose: This function defines three attributes. antenna1 and antenna2 are numpy arrays containing the antenna1 and antenna2 columns
                 from the data table.  antennasForTimes defines the times over which data was collected for that antenna.  It is a python
                 dictionary with keys of baseline (using antenna names) and values of numpy array of times.
        """

        tb.open(self.inputMs)
        self.antennas1 = tb.getcol('ANTENNA1')
        self.antennas2 = tb.getcol('ANTENNA2')
        tb.close()
        self.antennasForTimes = {}
        for i in range(len(self.uniqueAntennas)) :
            for j in range(len(self.uniqueAntennas)) :
                if i <= j :
                    antennaKey = "%s-%s" % (str(self.uniqueAntennas[j]),str(self.uniqueAntennas[i]))
                    indices = np.where((self.antennas1 == list(self.antennaNamesForAntennaIds).index(self.uniqueAntennas[i])) *
                                       (self.antennas2 == list(self.antennaNamesForAntennaIds).index(self.uniqueAntennas[j])))
                    self.antennasForTimes[antennaKey] = self.time[indices]

    def doAntennaStateMapping(self):
        """
        Input: None
        Output: None
        Responsible: S. Corder
        Purpose: This function defines one attribute, antennasForStates, a python dictionary.  The keys are baselines (using
                 antenna names) and values are staes used for that baseline.  Usually the autocorrelations are most useful.
        """

        self.antennasForStates = {}
        for i in range(len(self.uniqueAntennas)) :
            for j in range(len(self.uniqueAntennas)) :
                if i <= j :
                    antennaKey = "%s-%s" % (str(self.uniqueAntennas[j]),str(self.uniqueAntennas[i]))
                    indices = np.where((self.antennas1 == list(self.antennaNamesForAntennaIds).index(self.uniqueAntennas[i])) *
                                       (self.antennas2 == list(self.antennaNamesForAntennaIds).index(self.uniqueAntennas[j])))
                    self.antennasForStates[antennaKey] = self.states[indices]


    def getScansForTime(self,time,fudge=0.0):
        """
        Input: Time stamp in CASA native units as a string or float
        Output: Scan number associated with that time stamp.
        Responsible: S. Corder
        Purpose: This function returns the scan number for a specific timestamp.  It allows translation between time and scan.
        """
        for i in self.scansForTimes.keys() :
            if ((float(time) >= self.scansForTimes[i][0]-fudge) and (float(time) <= self.scansForTimes[i][-1]+fudge)) :
                return i

    def getTimesForScans(self,scans):
        """
        Input: Scan number as an integer or string or list of integers
        Output: A list of time ranges over which data exists for that scan (or those scans), each as a numpy array.
        Responsible: S. Corder, copied from getTimesForScan and modified by T. Hunter
        Purpose: Return the times associated with a given timestamp.  This allows translation between scan and time.
        """
        times = []
        if (type(scans) == int or type(scans) == np.int32):
            scans = [scans]
        for scan in scans:
            times.append(self.scansForTimes[int(scan)])
        return (times)

    def getTimesForScan(self,scan):
        """
        Input: Scan number as an integer or string
        Output: Time range over which data exists for that scan as a numpy array.
        Responsible: S. Corder
        Purpose: Return the times associated with a given timestamp.  This allows translation between scan and time.
        """
        return self.scansForTimes[int(scan)]

    def getScansForState(self,state):
        """
        Input: State id as an integer or string
        Output: The scan numbers, as a list, that use that specific state.
        Responsible: S. Corder
        Purpose: Return the scans that used a specific state.  This allos translation between state and scan.
        """

        scansForState = []
        for i in self.uniqueScans :
            if int(state) in self.statesForScans[i] : scansForState.append(i)
        return scansForState

    def getStatesForScan(self,scan):
        """
        Input: Scan number as a string or integer
        Output: States used during that scan.
        Responsible: S. Corder
        Purpose: Returns the states used during a given scan.  This allows translation between scan and state
        """

        return self.statesForScans[int(scan)]

    def getIntentsForScan(self,scan) :
        """
        Input: Scan as an integer or string.
        Output: Intent as a an array of strings with the names of the intents as values.
        Responsible: S. Corder
        Purpose: This returns the intents used in a specific scan allowing translation between scan and intent.
        """

        intentsForScan = []
        for i in range(len(self.intentsForStates)) :
            subIntents = self.intentsForStates[i]
            if int(scan) in self.getScansForState(i) : intentsForScan.extend(subIntents)
        return np.unique(intentsForScan)

    def getScansForIntent(self,intent) :
        """
        Input: Intent (as a string)
        Output: A numpy array of scans using the input intent.
        Responsible: S. Corder
        Purpose: This returns the scans using a specific intent.  This allows flagging based on intent and translation
                 between intent and scan.
        """

        scansForIntent = []
        for i in range(len(self.states)) :
            if intent in self.intentsForStates[self.states[i]] :
                scansForIntent.extend(self.getScansForState(self.states[i]))
        return np.unique(scansForIntent)

    def getScansForFieldID(self,field):
        """
        Input: Field, as an id.
        Output: Scans using that field
        Responsible: T. Hunter
        Purpose: This takes a field ID and tells you what scans it was used in.  It was
                 created to avoid a strange behavior of getScansForField for integer inputs.
        """
        indices = np.where(self.fields == field)
        return np.unique(self.scans[indices])

    def getScansForField(self,field):
        """
        Input: Field, as a name or id.
        Output: Scans using that field
        Responsible: S. Corder
        Purpose: This takes a field source and tells you what scans it was used in.
        """

        if not str(field).isdigit() : field = self.getFieldIdsForFieldName(field)
        indices = np.where(self.fields == field)
        return np.unique(self.scans[indices])

    def getFieldsForScans(self,scans,returnName=True):
        slist = []
        for scan in scans:
            slist.append(self.getFieldsForScan(scan))
        return([item for sublist in slist for item in sublist])

    def getFieldsForScan(self,scan,returnName=True):
        """
        Input: Scan as an integer or string
        Output: Field ids observed during that scan.
        Responsible: S. Corder
        Purpose: This takes a scan number and returns a field observed during that scan.  This allows translation between
                 scan and field.
        """

        indices = np.where(self.scans == int(scan))
        if returnName : return self.getFieldNamesForFieldId(np.unique(self.fields[indices]))
        else : return np.unique(self.fields[indices])

    def getFieldsForIntent(self,intent,returnName=True):
        """
        Input: intent as a string
        Output: field id as integer or array of names
        Responsible: S. Corder
        Purpose: This retrieves all of the fields that have been assigned a given intent during an observation.
        """
        _fields = []
        scans = self.getScansForIntent(intent)
        for i in scans :
            _field = self.getFieldsForScan(i)
            if _field not in _fields : _fields.append(_field)
        if returnName :
            return _fields
        else :
            return self.getFieldIdsForFieldName(_fields)

    def getFieldIdsForFieldName(self,sourceName):
        """
        Input: source name as string
        Output: field id as integer (actually it is returning the source id-Todd)
        Responsible: S. Corder
        Purpose: This translates between source/field name and field id.
        """
# The following fails because the case varies when .title() is applied: QSO vs. Qso, and TW Hya vs Tw Hya
#        return np.where(upper == sourceName.title())[0][0]
        if (type(sourceName) == list):
            sourceName = sourceName[0]
#        print "looking for %s in " % (sourceName), self.fieldNamesForFieldIds
        return np.where(self.fieldNamesForFieldIds == sourceName)[0]

    def getFieldNamesForFieldId(self,sourceId):
        """
        Input: field id (as string or integer)
        Output: field name
        Responsible: S. Corder
        Purpose: This translates between field id and field/source name.
        """
        if (type(sourceId) == int or type(sourceId) == np.int32 or type(sourceId) == np.int64 or type(sourceId) == str):
            if (len(self.fieldNamesForFieldIds) > int(sourceId)):
                # prevent "index out of bounds" error if field is not present
                return self.fieldNamesForFieldIds[int(sourceId)]
            else:
                return (None)
        else:
            # Todd added this check which was necessary for the Antennae Band 7 mosaic
            return [self.fieldNamesForFieldIds[s] for s in sourceId]

    def getFieldsForTime(self,time,returnName=True):
        """
        Input: Time in casa native units, returnName (boolean).  If returnName is true, the name is returned, else the id is returned.
               Default is returnName=True
        Output: Field name or id (depending on value of returnName).
        Responsible: S. Corder
        Purpose: Allows the field id/name to be returned for a specific observation time.
        """

        for i in self.fieldsForTimes.keys() :
            if (time in self.fieldsForTimes[i]) :
                if returnName : return i
                else : return self.getFieldNamesForFieldId(i)

    def getTimesForField(self,field):
        """
        Input: Field name or id (as a string or integer)
        Output: Times as a numpy array over which that field was observed in casa native units.
        Responsible: S. Corder
        Purpose: This allows you to determine the data time stamps for observations of a specific field source.
        """

        if str(field).isdigit() : field = self.fieldNamesForFieldIds(int(field))
        return self.fieldsForTimes[field]


def timeOnSource(ms='', field='', verbose=True,
                 asdm='', help=False, vm='', gapFactor=None, debug=False,
                 verboseComputeDuration=False, scienceSpwsOnly=False):
    """
    Uses ValueMapping to get the integration timestamps and computes a
    list of durations on the fields specified, attempting to detect and
    account for the inter-subscan latency.
    Inputs:
    ms: measurement set
    field: integer ID or name
    vm: a pre-existing ValueMapping structure for this measurement set
    scienceSpwsOnly: by default (False) WVR is ignored, but not SQLD,
          set to True to extract times using only the first science spw

    Returns a dictionary indexed by the source ID integer:
    {0: {'num_of_subscans': 2,
         'scans': [5],
         'field_id': [0],
         'source_name': '3c279',
         'minutes_on_source': 58.716000556945801,
         'minutes_on_source_with_latency': 65.664000511169434
         },
     ...
     'clock_time': 120
     'minutes_on_science': 58.716
     'minutes_on_science_with_latency': 65.66
     'percentage_time_on_science': 20.40
    }

It also prints a string convenient to a wikitable of format:
| date | SB name | exec UID | UT start-end | LST start-end | Total time | Time on source | el range | med pwv | antennas |

Notes on this string:
1) If the pointing table has been deleted it will print "pointing table empty" in the elevation range column.
2) If there is no on-source time, it will print "no onsource time" in the elevation range column.
3) If neither the CalWVR.xml nor ASDM_CALWVR files are present, it will print "unknown" in the pwv column.
4) The median pwv is over the entire dataset, not just the on-source scans.
5) For the number of antennas, it does not detect whether an antenna has been totally flagged!

    For further help and examples, see http://casaguides.nrao.edu/index.php?title=TimeOnSource
    -- Todd Hunter
    """
    if (os.path.exists(ms) == False):
        print "Could not open ms = %s" % (ms)
        return({})
#    if (usemsmd):
#        return(timeOnSourceMSMD(vis,field,verbose,asdm,gapFactor)
    if (vm==''):
#        Someday I should implement this:  au.timeOnSourceMSMD
        if (verbose):
            print "Running ValueMapping... (this may take a minute)"
        vm = ValueMapping(ms)
    tb.open(ms+'/FIELD')
    sourceIDs = tb.getcol('SOURCE_ID')
#    print "sourceIDs: ", sourceIDs
    tb.close()
    telescopeName = getObservatoryName(ms)
    mydict = {}
    if (field == ''):
        field = range(len(sourceIDs))
    elif (type(field) == int):
        field = [field]
    elif (type(field) == str):
        if (field.find(',') >= 0):
            field = field.split(',')
            fd = []
            for f in field:
                fd.append(int(f))
            field = fd
        else:
            try:
                field = int(field)
                field = [field]
            except:
                field = vm.getFieldIdsForFieldName(field)
                if (len(field) < 1):
                    print "No match for field name = " % (field)
                    return({})
    elif (type(field) == list):
        if (type(field[0]) == str):
            nf = []
            for f in field:
                nf.append(vm.getFieldIdsForFieldName(f))
            field = nf
    if (verbose):
        print "Considering non-CalAtmosphere, non-CalPointing, non-CalSideband scans."
        print "Total time on scans (including inter-subscan latency) in %s:" % (os.path.basename(ms))
    durations = []
    subscansAll = []
    fieldId = -1
    previousField = -1
    legend = "Field "
    multiField = False
    totalMinutes = 0
    scienceMinutesWithLatency = 0
    mydict = {}
    legend = ""
    multiFieldScansObserved = []
    if (casadef.casa_version >= casaVersionWithMSMD):
        mymsmd = createCasaTool(msmdtool)
        mymsmd.open(ms)
    else:
        mymsmd = ''
    if debug: print "field = ", field
    times = {}
    for findex in range(len(field)):
        times[findex] = {}
        f = field[findex]
        durationWithLatency = 0
        scienceDurationWithLatency = 0
        scans = vm.getScansForField(f)
        if debug:
            print scans, "= scans on field ", f
        scansObserved = []
        totalSubscans = 0
        multiField = False
        if (findex < len(field)-1):
          if (sourceIDs[f] == sourceIDs[field[findex+1]]): # and scans[0] == vm.getScansForField(field[findex+1])[0]):
              multiField = True
              legend += '%2d,' % (f)
#              print multiFieldScansObserved, scans, len(scans)
              for s in scans:
                  multiFieldScansObserved.append(s)
              continue
          else:
              legend += '%2d' % (f)
        else:
          if (sourceIDs[f] == sourceIDs[field[findex-1]]): #  and scans[0] == vm.getScansForField(field[findex-1])[0]):
              multiField = True
              for s in scans:
                  multiFieldScansObserved.append(s)
              legend += '%2d' % (f)
          else:
              legend += '%2d' % (f)
        subscansPerScan = []
        if (multiField):
            scans = np.unique(multiFieldScansObserved)
        if (debug):
            print "scans to check = ", scans
        for s in scans:
            intents = vm.getIntentsForScan(s)
            if debug: print "Scan %d intents = " % (s), intents
            skip = False
            for i in intents:
                if (i.find('CALIBRATE_ATMOSPHERE') >= 0 or
                    i.find('CALIBRATE_POINTING') >= 0 or
                    i.find('CALIBRATE_SIDEBAND') >= 0):
                    if debug: print "Skipping scan %d because it is calibration" % (s)
                    skip = True
                if (casadef.casa_version >= casaVersionWithMSMD):
                    if (s not in mymsmd.scansforfield(f)):
                        skip = False
            if (skip):
                if debug:
                    if (s not in mymsmd.scansforfield(f)):
                        print "Skipping scan %d, since msmd says it is not on field %d" % (s,d)
                    else:
                        print "Skipping calibration scan %d" % (s)
                continue
            scansObserved.append(s)
            times[findex][s] = vm.getTimesForScans(s)
            # times will be a list of length 1 (because s is a single integer scan)
            for t in times[findex][s]:
                (d,subscans) = computeDurationOfScan(s,None,verbose=verboseComputeDuration,gapFactor=gapFactor,
                                                     vis=ms, mymsmd=mymsmd, scienceSpwsOnly=scienceSpwsOnly)
                scanLength = np.max(t) - np.min(t)
                if (scanLength == 0):
                    # This happens with simulated data as it has no subscans
                    scanLength = d
                durationWithLatency += scanLength
                totalSubscans += subscans
                subscansPerScan.append(subscans)
                if ('OBSERVE_TARGET#ON_SOURCE' in intents or 'OBSERVE_TARGET.ON_SOURCE' in intents
                    or 'OBSERVE_TARGET#UNSPECIFIED' in intents or telescopeName == 'SMA'):
                    scienceDurationWithLatency += scanLength
        totalMinutes += durationWithLatency/60.
        scienceMinutesWithLatency += scienceDurationWithLatency/60.
        if (verbose):
            if (totalSubscans > 1 and multiField==False):
                print "Source %2d = Field %s = %s: %.1f sec = %.2f min (%d scan%s: %s, %d subscans)" % (sourceIDs[f],
                   legend, vm.getFieldNamesForFieldId(f),
                   durationWithLatency, durationWithLatency/60.,
                   len(scansObserved), 's' if len(scansObserved)>1 else '', scansObserved, totalSubscans)
            else:
                print "Source %2d = Field %s = %s: %.1f sec = %.2f min (%d scan%s: %s)" % (sourceIDs[f],
                   legend, vm.getFieldNamesForFieldId(f),
                   durationWithLatency, durationWithLatency/60.,
                   len(scansObserved), 's' if len(scansObserved)>1 else '', scansObserved)
        if (len(legend.split(',')) > 1):
            myfield = [int(x) for x in legend.split(',')]
        else:
            myfield = [f]
        mydict[sourceIDs[f]] = {'field_ids': myfield,
                     'source_name':vm.getFieldNamesForFieldId(f),
                     'scans': scansObserved,
                     'num_of_scans': len(scansObserved),
                     'num_of_fields': len(legend.split(',')),
                     'num_of_subscans': subscansPerScan,
                     'minutes_on_source_with_latency': durationWithLatency/60.,
                     'minutes_on_source': 0  # fill this in later
                     }
        legend = ""
    fullreport = True
    if (fullreport):
        [wikiline2,wikiline3,clockTimeMinutes,csvline] = lstrange(ms,verbose=False,vm=vm)
    else:
        clockTimeMinutes = computeClockTimeOfMS(vm=vm)
    print "Clock time = %.2f min, Total time = %.2f min,  science time = %.2f min" % (clockTimeMinutes, totalMinutes, scienceMinutesWithLatency)
    csvline += ',%.2f,%.2f' % (totalMinutes, scienceMinutesWithLatency)
    ##################################################################################3
    if (verbose):
        print "\nMy attempt to detect and account for inter-subscan latency follows:"
    ##################################################################################3
    legend = ""
    multiField = False
    totalMinutes = 0
    scienceMinutes = 0
    fid = 0
    multiFieldScansObserved = []
    if (debug): print "for findex in ", range(len(field))
    for findex in range(len(field)):
        f = field[findex]
        # f is the field ID
        duration = 0
        scienceDuration = 0
        scans = vm.getScansForField(f)
        scansObserved = []
        totalSubscans = 0
        multiField = False
        if (findex < len(field)-1):
            if debug: print "on field = ", findex, field[findex]
            if (sourceIDs[f] == sourceIDs[field[findex+1]]): #  and scans[0] == vm.getScansForField(field[findex+1])[0]):
                # non-final fields
                multiField = True
                legend += '%2d,' % (f)
                for s in scans:
                    multiFieldScansObserved.append(s)
                # By continuing here, we effectively use the final field in the Field table as the
                # effective field.  Sometimes, the final field is never observed, so we need to
                # avoid that possibility, but only if the prior field is not entirely a CAL field.
                if (s in pickScansForField(mymsmd,f)):
                    # use fieldToUse unless it has no science scans
                    myintents = mymsmd.intentsforfield(f)
                    if ('OBSERVE_TARGET#ON_SOURCE' in myintents or
                        'OBSERVE_TARGET.ON_SOURCE' in myintents or
                        'OBSERVE_TARGET#UNSPECIFIED' in myintents):
                        fieldToUse = f
                if debug: print "On field %d, continuing before checking scans" % (f)
                continue
            else:
                legend += '%2d' % (f)
                fieldToUse = f
        else:
            # final field
            if debug: print "on final field = ", f
            if (findex > 0):
                # This is the same logic as the non-final fields, but in reverse sense
                if (sourceIDs[f] != sourceIDs[field[findex-1]]):
                    fieldToUse = f
                    if debug: print "final 2 source IDs not equal: Setting fieldToUse = ", f
                else:
                    myintents = mymsmd.intentsforfield(fieldToUse)
                    if ('OBSERVE_TARGET#ON_SOURCE' not in myintents and
                        'OBSERVE_TARGET.ON_SOURCE' not in myintents and
                        'OBSERVE_TARGET#UNSPECIFIED' not in myintents):
                        fieldToUse = f
                        if debug: print "Setting fieldToUse = ", f
            elif (s in pickScansForField(mymsmd,f)):
                if debug: print "final elif: set fieldToUse = ", f
                fieldToUse = f
            else:
                if debug: print "Not setting fieldToUse"
            if (sourceIDs[f] == sourceIDs[field[findex-1]]):
                multiField = True
                if debug: print "multiField=True"
                legend += '%2d' % (f)
                for s in scans:
                    multiFieldScansObserved.append(s)
                if debug: print "multiFieldScansObserved=", multiFieldScansObserved
            else:
                legend += '%2d' % (f)

        if (multiField):
            scans = np.unique(multiFieldScansObserved)
        if debug: print "Checking fieldToUse=%d in scans = " % (fieldToUse), scans
        scienceMinutesPerScan = {}
        for s in scans:
            intents = vm.getIntentsForScan(s)
            skip = False
            for i in intents:
                if (i.find('CALIBRATE_ATMOSPHERE') >= 0 or
                    i.find('CALIBRATE_POINTING') >= 0 or
                    i.find('CALIBRATE_SIDEBAND') >= 0):
                    skip = True
                    if debug:
                        print "Skipping calibration scan %d" % (s)
                if (casadef.casa_version >= casaVersionWithMSMD):
                    if (s not in pickScansForField(mymsmd,fieldToUse)):  # change f to fieldToUse Jun 29, 2015
                        skip = True
                        if debug:
                            print "Skipping scan %d because not in scans for field %d: %s" % (s, fieldToUse, pickScansForField(mymsmd,fieldToUse))
            if (skip):
                continue
            scansObserved.append(s)
            scienceDurationThisScan = 0
            # times[findex][s] will be a list of length 1 (because s is a single integer scan)
            for t in times[findex][s]:
#                print "Running au.computeDurationOfScan(%d,%s,verbose=False,gapFactor=%g,vis='%s')" % (s,str(t),gapFactor,ms)
                (d,subscans) = computeDurationOfScan(s,None,verbose=False,gapFactor=gapFactor,vis=ms,
                                                     mymsmd=mymsmd, scienceSpwsOnly=scienceSpwsOnly)
                duration += d
                totalSubscans += subscans
                if ('OBSERVE_TARGET#ON_SOURCE' in intents or 'OBSERVE_TARGET.ON_SOURCE' in intents
                     or 'OBSERVE_TARGET#UNSPECIFIED' in intents or telescopeName == 'SMA'):
                    scienceDurationThisScan += d
            scienceMinutesPerScan[s] = scienceDurationThisScan/60.
            scienceDuration += scienceDurationThisScan
        totalMinutes += duration/60.
        scienceMinutes += scienceDuration/60.
        if (verbose):
            fieldWarning = ""
            if (totalSubscans > 1 and multiField==False):
                print "Source %2d = Field %s = %s: %.1f sec = %.2f min (%d scan%s: %s, %d subscans%s)" % (sourceIDs[f],
                 legend, vm.getFieldNamesForFieldId(f), duration, duration/60.,
                 len(scansObserved), 's' if len(scansObserved)>1 else '', scansObserved, totalSubscans, fieldWarning)
            else:
                print "Source %2d = Field %s = %s: %.1f sec = %.2f min (%d scan%s: %s)" % (sourceIDs[f],
                    legend, vm.getFieldNamesForFieldId(f), duration, duration/60.,
                    len(scansObserved), 's' if len(scansObserved)>1 else '', scansObserved)
        legend = ""
        mydict[sourceIDs[f]]['minutes_on_source'] = duration/60.
        mydict['minutes_on_science'] = scienceMinutes
        mydict['minutes_on_science_per_scan'] = scienceMinutesPerScan
        mydict['minutes_on_science_with_latency'] = scienceMinutesWithLatency
        mydict['percentage_time_on_science'] = 100*scienceMinutes/clockTimeMinutes
        mydict['clock_time'] = clockTimeMinutes
        mydict['source_ids'] = list(sourceIDs)
        mydict['num_of_sources'] = len(sourceIDs)
# might add this someday
#        elevs = csvline.split(',')
#        startElev = float(elevs[-2])
#        stopElev = float(elevs[-1])
#        mydict['elevation_range'] = [startElev,stopElev]
    if (fullreport):
        wikiline2 += '%.1f | %.1f | ' % (totalMinutes, scienceMinutes)
        wikiline2 += wikiline3
    # Now get the PWV if possible
    pwvmean = -1
    if (casadef.casa_version >= casaVersionWithMSMD):
        mymsmd.close()
    if (telescopeName.find('ALMA')>=0):
        if (os.path.exists('ASDM_CALWVR') or os.path.exists(ms+'/ASDM_CALWVR')):
            if (os.path.exists(ms+'/ASDM_CALWVR')):
                [pwvmean, pwvstd]  = getMedianPWV(ms)
            else:
                [pwvmean, pwvstd]  = getMedianPWV('.')
        elif (os.path.exists('ASDM_CALATMOSPHERE') or os.path.exists(ms+'/ASDM_CALATMOSPHERE')):
            if (os.path.exists(ms+'/ASDM_CALATMOSPHERE')):
                [pwvmean, pwvstd]  = getMedianPWV(ms)
            else:
                [pwvmean, pwvstd]  = getMedianPWV('.')
        elif (os.path.exists(ms+'/CalWVR.xml')):
            [pwvtime, pwv, antenna] = readpwv(ms)
            pwvmean = np.mean(pwv)
        elif (os.path.exists('CalWVR.xml')):
            [pwvtime, pwv, antenna] = readpwv('.')
            pwvmean = np.mean(pwv)
        else:
            print "No ASDM_CALWVR, ASDM_CALATMOSPHERE or CalWVR.xml table found.  You should importasdm with asis='*' or copy the CalWVR.xml file from your ASDM to your working directory (or your ms directory)."
    tb.open(ms+'/ANTENNA')
    nAntennas = len(tb.getcol('NAME'))
    tb.close()
    if (fullreport):
        if (pwvmean < 0):
            wikiline2 += ' unknown_PWV | %d ' % (nAntennas)
        else:
            wikiline2 += ' %.2f | %d |' % (pwvmean, nAntennas)
        wikiline2 += '   |'   # Requested by Andreas Lundgren on 2012-05-23
    print "Latency removed: Total time = %.2f min,   science time = %.2f min" % (totalMinutes, scienceMinutes)
    csvline += ',%.1f,%.1f' % (totalMinutes, scienceMinutes)
    if (fullreport):
        print "wikiline = %s" % (wikiline2)
        print csvline
    print "WARNING: This task does not account for any flagging."
    return(mydict)
    # end of timeOnSource  vm.

try:
    mset = str(sys.argv[sys.argv.index('time_on_source_cmdline.py')+1])
    timeOnSource(ms=mset)
except IndexError:
    print 'Usage casa -c time_on_source_cmdline.py <measurement_set>'
