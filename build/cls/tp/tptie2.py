# Multiple 2D well ties
# @author Andrew Munoz, Colorado School of Mines
# @version 04.10.2013

from dtw import *
from tputils import *
from wtutils import *
from imports import *

_pres = '_cwpp13'

#printWellLengths()
# Params:
fp = 35
q  = 100
nlr = 1.0
# Time warping
vrmin = 0.98 # v*/v min 0.5 #0.99,1.01,.08,195 
vrmax = 1.02 # v*/v max inf
dr    = 0.02 # 1D smoothness vertically
# Image warping
dr1   = 0.02 # 2D smoothness vertically
dr2   = 0.50 # 2D smoothness horizontally
dr3   = 0.50 # 3D smoothness horizontally
rmin1,rmax1 = -0.08,0.08 # 2D vertical constraint
rmin2,rmax2 = -1.00,1.00 # 2D horizontal constraint
rmin3,rmax3 = -0.50,0.50 # 3D horizontal constraint
smin,smax = -0.20,0.20 # min and max time shifts for 2D/3D warping
#print "# of vertical slopes = ",1+2*rmax1/dr1 # base numbers on the vmin and vmax
#print "# of horiz y  slopes = ",1+2*rmax2/dr2 # base number of slopes on misalignment of wells 
#print "# of horiz x  slopes ,dw= ",rmax3/dr3

def main(args):
	goMultiWell2D()

def goMultiWell2D():
	# 1. Get the 3D seismic image
	global s1,s2,s3
	global wells
	tpimg,s3,s2,s1 = getImage(normalize=True,smooth=True)
	#tpimg,s3,s2,s1 = getImage(normalize=True,cut1=[0.10,1.25],cut2=[2.5,8.0],cut3=[0.5,3.5])
	#tpimgs,s3,s2,s1 = getImage(normalize=True,smooth=True)
	dt = s1.delta
	global smin,smax
	smin,smax = infl(smin/dt),ince(smax/dt)

	# 2. Get the wells, coordinates, and tops
	################
	# TESTING FOR 3D PREP:
	# wells less than 127: 8,6,7,13,2,17
	###wells = [0,1,2,3,4,6,7,8,9,10,11,12,13,15,16,17] # 
	#wells =    [0,1,  3,4,      9,10,11,12,   15,16   ] #>=100 all 
	#wells =    [      3,          10,         15,16   ] #>=100 filt1
	#wells =    [      3,        9,10,11               ] #>=100 filt2 #throw 11,1?
	#wells =    [0,    3,4,       9,     12        ] #>=100 filt3
	#wells =    [0,    3,                12,17,2, 13] #>=100 filt3
	#wells = [0,3,4,9,15,16,12,17,2,13,11]
	#wells = [0,3,  9,15,16,10           ]
	#wells = [0,3,4,9,10,                ]
	#wells = [3,10,15,16                ]
	wells = [3,15,16,11,4]

	#wells = [0,1,3,4,9,10,11,12,15,16]
	## potential slices:
	#wells =  [3,6,10,11,12] 
	#wells = [3,16,15,4,9] # All long wells (>200)
	######wells = [3,16,15,4,9,12,17] # 3D
	## test lines all with 3:
	#wells = [8,10,3,6,0,12,4,11,9] # 11,10,8
	#wells = [3,12,4,9] # 8,11,10,0,6
	#wells = [3,1,11,9] # 1,11 no tie
	#wells = [3,16,15,7,13,2,17] 
	################
	# 2D slices for display
	##wells = [3,16,15]
	#wells = [3,12,4,17]
	#wells = [3,4,9]
	#wells = [3,16,4,9] # for 2D figs: 1.02,.98,.02,.02,.05,-.10,.10,-1.,1.,-.05,.05,ph153,fp35,q80
	#wells = [3,4,9]
	#phases = [172]*len(wells)
	phases = [71]*len(wells)
	#phases = [100]*len(wells) #4%
	#phases = [60]*len(wells) # 2%
	#print 'Average phase is ',phases[0]
	wells = sortWells(wells)
	tops,fmsp,datums = getWellTops(wells)
	ail,vl,dl,gl,sz,x2w,x3w = getMultipleLogs(wells,datums)
	#showMaxWellDist(x2w,x3w)
	#showMinWellDist(x2w,x3w)

	# 3. Extract the image between the wells
	x2s,x3s,sn,x2wn = getCurveThruWells(x2w,x3w,s2.delta)
	f = getImageAlongCurve(tpimg,x2s,x3s)
	imgh = WTutils().getSlice(470,tpimg)
	plotCoordSlice(s2,s3,x2s,x3s,x2w,x3w,imgh,slides=_pres+"hcoordslice")

	# Flatten Image
	#goSlopes2(f,sn)
	#goFlat2(f,sn)

	# 4. Get horizons 
	gh = get2DHorizons(x2s,x3s,sz)

	# 5. Get the synthetics
	ssy,sy,tzl,ztl,sztl = getSeismograms(sz,vl,dl,q,phase=phases)

	#lims = [sn.first,0.2,sn.last,1.8]
	#lims = [sn.first,0.2,sn.last,1.6]
	#lims = [sn.first,s1.first,sn.last,1.8]
	lims = [sn.first,0.1,sn.last,1.2]
	plotSlice(sn,s1,f,lims=lims,title='seismic image')
	#plot2DLogs(s1,sn,f,ssy,sy,x2wn,lims=lims,title='initial synthetics')
	#plot2DLogs(s1,sn,f,ssy,sy,x2wn,lims=lims,gh=gh,tops=tops,fmsp=fmsp,sz=sz,tz=tzl,title='wells untied')
	#plotLogPanelsTops(s1,sz,ail,gl,wells,fmsp,tops,ztl,tzl,gh,x2w)
	
	# 6. Warp the synthetics individually
	#ncl = 20
	#wp = constraintsForTops(tops,gh,tzl,x2w,s1,sz,ncl)
	u,h,sh,tz,vint,zt,szt = warpIndividually(\
		tpimg,ssy,sy,sz,x2w,x3w,tzl,vl,vrmin,vrmax,dr,vb=[0.25,0.35])#,ph2=True)
	#plot2DLogs(s1,sn,f,sh,h,x2wn,lims=lims,title='single tie synthetics')
	#plot2DLogs(s1,sn,f,sh,h,x2w,lims=lims,title='wells tied individually',\
	#					gh=gh,tops=tops,fmsp=fmsp,sz=sz,tz=tz)
	#plotLogPanelsTops(s1,sz,ail,gl,wells,fmsp,tops,zt,tz,gh,x2w)

	print '2D warp refinement'
	# 7. Interpolate and warp the synthetics
	bgs = 2.0 
	tensors = makeTensors(f)
	cs = None#[-2.5,2.5]
	#plotSlice(sn,s1,f,tens=tensors,lims=lims)

	vlogs=None
	vinterps=True
	dinterps=None
	ginterps=None	
	tdepth=True

	fplot  = False
	slides = True

	wonly = True
	wsngl = None

	#########################################################
	### PLOTS

	if slides:

		lims = [sn.first,0.1,sn.last,1.2]
	
		plotSlice(sn,s1,f,lims=lims,paper=_pres+'seisextract',psize=[712,524])
		plot2DLogs(s1,sn,f,ssy,sy,x2wn,lims=lims,paper='2dwuntindv2',rmb=True)#,psize=[712,394])
		#plot2DLogs(s1,sn,f,ssy,sy,x2wn,lims=lims,\
		#	gh=gh,tops=tops,fmsp=fmsp,sz=sz,tz=tzl,paper='2dwuntindvwh')#,psize=[712,394])

		plot2DLogs(s1,sn,f,sh,h,x2wn,lims=lims,title='wells tied individually',paper='2dwtindv2',rmb=True)#,psize=[712,394])
		#plot2DLogs(s1,sn,f,sh,h,x2wn,lims=lims,title='wells tied individually',\
		#						gh=gh,tops=tops,fmsp=fmsp,sz=sz,tz=tz,paper='2dwtindvwh')#,psize=[712,394])

		if wsngl:
			#dr1 = 0.02
			#dr2 = 0.50 
			#rmin1,rmax1 = -0.12,0.12 
			#rmin2,rmax2 = -1.00,1.00 
			#smin,smax = -0.10,0.10 
			gd = interpolateSynthetics2(sy,ssy,sn,x2wn,bgs,ts=tensors)
			plotSlice(sn,s1,gd,lims=lims,title='interpolated initial synthetics ',clips=cs,paper='intinitsy')
			ga = interpolateSynthetics2(h,sh,sn,x2wn,bgs,ts=tensors)
			plotSlice(sn,s1,ga,lims=lims,title='interpolated single tie synthetics',clips=cs,paper='intsgltiesy')
			h2d,u2d,hn,shn,tzn2,vintn2,ztn2,sztn2,gs = warpAllin2D(f,ga,sn,x2wn,h,sh,vint,tz,sz,ssy)
			gc = interpolateSynthetics2(hn,shn,sn,x2wn,bgs,ts=tensors)
			plotSlice(sn,s1,gc,lims=lims,title='interpolated multi-tie synthetics',clips=cs,paper='intmltsy')
			plotSlice(sn,s1,h2d,lims=lims,title='warped single tied synthetics',clips=cs,paper='warpmltsy')
			plot2DLogs(s1,sn,f,shn,hn,x2wn,lims=lims,title='multi-tie synthetics',paper='2dwtsglmult2',rmb=True)
			votl2 = makeVot(sz,vintn2,tzn2,shn)
			plot2DLogs(s1,sn,f,shn,votl2,x2wn,lims=lims,title='velocity logs and seismic',cbar='Velocity (km/s)',lcmap=jet,lclips=[2.5,5.5],paper='vnsmultsgl',velc=True)

			print 'nrms of single ties is ',nrms(ga,f)
			print 'nrms of refined single ties is ',nrms(gc,f)

		# Iterate over multiple image warpings 
		if wonly: 
			gd = interpolateSynthetics2(sy,ssy,sn,x2wn,bgs,ts=tensors)
			#plotSlice(sn,s1,gd,lims=lims,title='interpolated initial synthetics')
			h2d,u2d,hn,shn,tzn,vintn,ztn,sztn,gs = warpAllin2Donly(tpimg,sn,f,gd,x2wn,x3w,sy,ssy,vl,tzl,sz)#,ph=True)
			gb = interpolateSynthetics2(hn,shn,sn,x2wn,bgs,ts=tensors)
			plotSlice(sn,s1,gb,lims=lims,title='interpolated multi tied synthetics '+str(1),clips=cs,paper='intmltonlysy')
			plotSlice(sn,s1,gd,lims=lims,title='interpolated inital synthetics '+str(1),clips=cs,paper='intinitsy')
			print 'nrms of initial ties is ',nrms(gd,f)
			if wsngl:
				print 'nrms of single ties is ',nrms(ga,f)
				print 'nrms of refined single ties is ',nrms(gc,f)
			print 'nrms of multiple ties is ',nrms(gb,f)
			votl = makeVot(sz,vintn,tzn,shn)
			plot2DLogs(s1,sn,f,shn,votl,x2wn,lims=lims,title='velocity logs and seismic',cbar='Velocity (km/s)',lcmap=jet,lclips=[2.5,5.5],paper='vnsmulto',velc=True)
		
		
		plotSlice(sn,s1,f,lims=lims,title='seismic image',paper='seismic')
		plotSlice(sn,s1,f,lims=lims,title='seismic image',tens=tensors,paper='seismicwtens')
		#plotSlice(sn,s1,dist2,pts=gs,lims=lims,title="distances from wells",cmap=jet,cbar=True)
		plotSlice(sn,s1,h2d,lims=lims,title='warped multi-tie synthetics',paper='warpmltonly')
		#plotSlice(sn,s1,u2d,pts=gs,cbar=True,title='Interpolated tied wells 2D shifts')
		plot2DLogs(s1,sn,f,shn,hn,x2wn,lims=lims,title='multi-tie synthetics',paper='2dwtmult2',rmb=True)
		plot2DLogs(s1,sn,f,shn,hn,x2wn,lims=lims,cbar='Amplitude',title='multi-tie synthetics',paper='seiswcb',rml=True)
		#plot2DLogs(s1,sn,f,shn,hn,x2wn,lims=lims,title='wells tied with 2D refinement',\
		#						gh=gh,tops=tops,fmsp=fmsp,sz=sz,tz=tzn)
	
		if dinterps:
			print 'convert time to depth...'
			if tdepth:
				bgs = 1.0
				stz,tzs = computetz2(tzn,sz,sn,x2wn,bgs)
				plotSlice(sn,stz,tzs,cmap=jet,cbar='Time (s)',title='Times',paper='tzmulo',vlabel='Depth (km)')
				plotSlice(sn,stz,tzs,cmap=jet,cbar='Time (s)',title='Times',paper='tzmulocont',vlabel='Depth (km)',cont=True)
				
				stz,tsc = computetz2(tzn,sz,sn,x2wn,bgs,scat=True)
				plotSlice(sn,stz,tsc,cmap=jet,cbar='Time (s)',title='Times',paper='tzscat',vlabel='Depth (km)',linear=True)

				dimg = imageTimeToDepth2(f,tzs,sn)
				plotSlice(sn,stz,dimg,title='seismic depth image',paper='seismicdepth',vlabel='Depth (km)')
				plotSlice(sn,s1,f,lims=[sn.first,0.0,sn.last,max(tzs)],title='seismic image',paper='seismic0')

			print 'interpolate depths...'
			bgs = 0.5
			tensors = makeTensors(zerofloat(s1.count,s2.count))#,p1=0.6)
			ztlmg = interpolateSynthetics2(ztl,sztl,sn,x2wn,bgs,ts=tensors)
			plotSlice(sn,s1,ztlmg,title='log depths',cbar='Depth (km)',cmap=jet,lims=lims,paper='ztinit')
			if wonly:
				ztimg = interpolateSynthetics2(zt,szt,sn,x2wn,bgs,ts=tensors)
				plotSlice(sn,s1,ztimg,cmap=jet,lims=lims,cbar='Depth (km)',title='2D multi depths',paper='ztmulo')
				plotSlice(sn,s1,ztimg,cmap=jet,lims=lims,cbar='Depth (km)',title='2D multi depths',paper='ztmulocont',cont=True)
			if wsngl:
				ztomg = interpolateSynthetics2(zt,szt,sn,x2wn,bgs,ts=tensors)
				plotSlice(sn,s1,ztomg,title='1D tie depths',cbar='Depth (km)',cmap=jet,lims=lims,paper='ztsgl')
				ztimg2 = interpolateSynthetics2(ztn2,sztn2,sn,x2wn,bgs,ts=tensors)
				plotSlice(sn,s1,ztimg2,cmap=jet,lims=lims,cbar='Depth (km)',title='2D refined depths',paper='ztmult')

			#tptz = getTptz()
			#tpdepths = getImageAlongCurve(tptz,x2s,x3s,s1,s2,s3)
			
		if vinterps:
			print 'interpolate velocities...'
			bgs = 1.0
			cv = [2.5,5.5]
			tensors = makeTensors(f,p1=1.0)
			vlmg = interpolateLogs2(vl,sz,sn,x2wn,tzl,bgs,ts=tensors)
			plotSlice(sn,s1,vlmg,title='initial log velocities',cbar='Velocity (km/s)',cmap=jet,lims=lims,clips=cv,paper='velinit')
			plotSlice(sn,s1,vlmg,cmap=jet,lims=lims,cbar='Velocity (km/s)',title='log velocities',clips=cv,paper=_pres+"velinittens",tens=tensors)
			if wonly:
				vimg = interpolateLogs2(vintn,sz,sn,x2wn,tzn,bgs,ts=tensors)
				plotSlice(sn,s1,vimg,cmap=jet,lims=lims,cbar='Velocity (km/s)',title='1 step multi-tie velocities',clips=cv,paper='velmultonly')
			if wsngl:
				vintmg = interpolateLogs2(vint,sz,sn,x2wn,tz,bgs,ts=tensors)
				plotSlice(sn,s1,vintmg,title='single-tie updated velocities',cbar='Velocity (km/s)',cmap=jet,lims=lims,clips=cv,paper='velsgltie')
				vimg2 = interpolateLogs2(vintn2,sz,sn,x2wn,tzn2,bgs,ts=tensors)
				plotSlice(sn,s1,vimg2,cmap=jet,lims=lims,cbar='Velocity (km/s)',title='multi-tie updated velocities',clips=cv,paper='velmult')


def warpAllin2D(f,g,s2n,x2w,syh,ssyh,vl,tzl,sz,ssy):
	dt = s1.delta
	dwc = DynamicWarpingWTM(f,g,smin,smax)
	u = dwc.findShifts(rmin1,rmax1,dr1,rmin2,rmax2,dr2)
	h = dwc.applyShifts(u)
	u = mul(u,-1.0*dt)
	g1=dwc.getG1();g2=dwc.getG2();
	gs= getGS(dt,s2n.delta,s2n.first,g1,g2)
	pred = 0.0
	wh,wsh,tz,vint,zt,szt=[],[],[],[],[],[]
	nw = len(x2w)
	for i in range(nw):
		x2i = int((x2w[i]-s2n.first)/s2n.delta)
		fsh = ssyh[i].first; fshi=int(round(fsh))
		h0,u2 = h[x2i],u[x2i]
		#plotCurve(u2,s1)
		wsh0 = Sampling(s1.count,dt,0.0)#fsh+u2[fshi])
		tz0,vint0 = updateTDandVelocity(sz[i],u2,tzl[i])
		zt0,szt0 = getzt(sz[i],tz0)
		wh.append(h0); wsh.append(wsh0); szt.append(szt0);
		tz.append(tz0); vint.append(vint0); zt.append(zt0);
		#pred+=printError(s1,wsh[i],h0,f[x2i]); 
		print 'R = ',correlation(h0,f[x2i],wsh[i],s1)
	#print 'Average Predictability =',pred/nw
	return h,u,wh,wsh,tz,vint,zt,szt,gs

def warpAllin2Donly(tpimg,s2n,f,g,x2w,x3w,sy,ssy,vl,tzl,sz,nr=1,ph=None):
	dt = s1.delta
	nw = len(x2w)
	if ph:
		traceMulti = []
		for iw in range(nw):
			traceMulti.append(getTrace(tpimg,s1,s2,s3,x2w[iw],x3w[iw],tonly=True,nr=nr))
		phase = rotatePhaseMulti(sy,traceMulti,rmin1+1,rmax1+1,dr1,dt)
		for iw in range(nw):
			sy[iw] = applyPhaseRot(phase,sy[iw])
	dwc = DynamicWarpingWTM(f,g,smin,smax)
	u = dwc.findShifts(rmin1,rmax1,dr1,rmin2,rmax2,dr2)
	h = dwc.applyShifts(u)
	u = mul(u,-1.0*dt)
	#plotSlice(s2n,s1,u,cbar=True)
	g1=dwc.getG1();g2=dwc.getG2();
	gs= getGS(dt,s2n.delta,s2n.first,g1,g2)
	pred = 0.0; rsum=0
	wh,wsh,tz,vint,zt,szt=[],[],[],[],[],[]
	for i in range(nw):
		x2i = inro((x2w[i]-s2n.first)/s2n.delta)
		h2,u2 = h[x2i],u[x2i]
		h0,wsh0 = getNewSyntheticSampling(h2,u2,ssy[i])
		tz0,vint0 = updateTDandVelocity(sz[i],u2,tzl[i])
		zt0,szt0 = getzt(sz[i],tz0)
		wh.append(h0); szt.append(szt0);
		tz.append(tz0); vint.append(vint0); zt.append(zt0);
		wsh.append(wsh0)
		trc = f[x2i]
		r = correlation(h0,trc,wsh0,s1)
		rsum += r
		wn = str(wells[i])
		print 'R = ',r,'for well ',wn
		#plotCurve(f[x2i],s1,h0,s1,title=str(i)+" h vs g")
		#pred+=printError(s1,wsh[i],h0,f[x2i]); 
		print 'R = ',correlation(h0,f[x2i],wsh[i],s1)
	#print 'Average Predictability =',pred/nw
	return h,u,wh,wsh,tz,vint,zt,szt,gs

def getNewSyntheticSampling(h,u,ssy):
	ni = ssy.count
	fi = ssy.first
	li = ssy.last
	n1 = s1.count
	d1 = s1.delta
	f1 = s1.first
	fii = infl((fi-f1)/d1)
	lii = ince((li-f1)/d1)
	tm = add(floats(s1.getValues()),u)
	fo = tm[fii]
	lo = tm[lii]
	np = infl((fi-tm[0])/d1)
	fio = infl((fo-f1)/d1)
	no = infl((lo-fo)/d1)
	if no+fio>n1: no -= no+fio-n1
	sh = Sampling(no,d1,fo)
	ho = copy(no,fio,h)
	return ho,sh

#------------------------------------------------------------------------------#

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

