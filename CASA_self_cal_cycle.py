msfile = str(sys.argv[3]) 
modelfile = str(sys.argv[4])
cycle = str(sys.argv[5])
caltype = str(sys.argv[6])
soli = str(sys.argv[7])
mysteps = map(int,sys.argv[8:])
print sys.argv
print mysteps


if len(sys.argv)<3:
	print "Usage casapy -c self_cal_cycle.py msfile modelfile cycleno ap/p solint mysteps"
	sys.exit()


thesteps = [0]
step_title = {0: 'Model in the data',
	      1: 'gaincal',
              2: 'applycal'}

try:
 	print 'List of steps to be executed ...', mysteps
	thesteps = mysteps
except:
	print 'global variable mysteps not set.'
if (thesteps==[]):
	thesteps = range(0,len(step_title))
	print 'Executing all steps: ', thesteps

mystep=0
if(mystep in thesteps):
	casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
	print 'Step ', mystep, step_title[mystep]
	
	ft(vis=msfile, model=modelfile, usescratch=True)
mystep=1
if(mystep in thesteps):
	casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
	print 'Step ', mystep, step_title[mystep]
	minbl = 3	
	if caltype == 'ap':
		minbl = 4

	gaincal(vis=msfile, caltable=msfile[:-3]+'_SC_cycle'+cycle+'.'+caltype, solint=soli, calmode=caltype, minsnr=5.0, minblperant=minbl)

mystep=2
if(mystep in thesteps):
	casalog.post('Step '+str(mystep)+' '+step_title[mystep],'INFO')
	print 'Step ', mystep, step_title[mystep]
	sctable = msfile[:-3]+'_SC_cycle'+cycle+'.'+caltype

	applycal(vis=msfile, docallib=False, gaintable=sctable)
