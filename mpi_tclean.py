import time
start_time = time.time()

bands = ['L','C','X','K','Ka']
spw_bands = ['0~15','16~31','32~47','49~63','64~79']
cellsizes = ['0.3arcsec','0.1arcsec','0.06arcsec','0.02arcsec','0.01arcsec']
image = [[512,512],[1560,1000],[2560,2000],[5560,4256],[10000,8000]]

for i in range(len(bands)):
	msfile = 'All_bands.ms'
	imname = '2012_M82_'+bands[i]+'_IM'
	fields = 'M82'
	spws = spw_bands[i]
	datacol = 'corrected'
	imagesize = image[i]
	cellsize = cellsizes[i]
	grid = 'widefield'
	facet = 1
	wproj = -1
	deconvolve = 'clarkstokes'
	scale = []
	weights = 'natural'
	niterations = 50000

	tclean(vis=msfile,selectdata=True,field=fields,spw=spws,timerange="",uvrange="",antenna="",scan="",observation="",intent="",datacolumn=datacol,imagename=imname,imsize=imagesize,cell=cellsize,gridder=grid,facets=facet,wprojplanes=wproj,deconvolver=deconvolve,scales=[],nterms=2,restoringbeam=[],outlierfile="",weighting=weights,robust=0.5,npixels=0,uvtaper=[],niter=niterations,gain=0.1,threshold=0.0,cycleniter=-1,cyclefactor=1.0,calcres=True,calcpsf=True,parallel=True)

print '--- %s seconds ---' % (time.time() - start_time)
