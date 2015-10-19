"""
"""
from tp import * 
from wt import *
from dtw import DynamicWarpingWT
from dtw import SyntheticSeismogram
from dtw import SeismicModel1D
from tputils import *
from wtutils import *
from imports import *
from tpwelltie import *
from edu.mines.jtk.dsp import DynamicWarping as DynamicWarpingO

_pres = 'cwpr13_'
#_pres = 'cwpp13_'

# Run types for a single well
Comparisonsf = None	
Warpsf = True

def main(args):
	print '\n'
	#goSingleWellf()
	#goMultiWellf()
	#printWellIDList()

def goSingleWellf():
	#Params
	well = [3]
	phase = 49 #	well 3
	_ID,_set = getWellFromSet(well[0]) 
	ail,yv,xv,zv,vl,dl = getLogs(_ID,_set)
	tops,fmsp,datum = getWellTops(well)
	nz = len(vl)
	fz = zv#+datum[0]
	sz = Sampling(nz,dz,zv) # in km
	fp = 35 #(Hz) peak frequency
	tpimg,sx,sy,st = getImage(normalize=True)
	traces,sx,smy,st,xslice,yc,xc = getTrace(tpimg,st,sy,sx,yv,xv,nr=2) #(x3,x2,x1)=(x,y,z)
	tr = traces[0]
	dt = st.delta
	q1 = 10000000
	q2 = 80
	#vrmin = 0.85 # v*/v min 0.5 
	#vrmax = 1.25 # v*/v max inf
	#dr    = 0.02 # resolution
	vrmin = 0.98 # v*/v min 0.5 
	vrmax = 1.02 # v*/v max inf
	dr    = 0.02 # resolution
	#dr = computedr(vrmin,vrmax,dr) 
	#plotLogPanel(sz,vl,dl,reflectivity(ail),paper='plogpanel')
	
	def warpAndUpdatef(ssy,sy,tzl,title):
		# Rotate the phase of the seismogram to find the optimum phase
		#sy = rotatePhase(sy,traces,vrmin,vrmax,dr)
		#sy = rotateAllPhase(sy,traces,vrmin,vrmax,dr)
		#sy = rotatePhaseHv(st,sz,sx,smy,xc,yc,ssy,sy,tzl,traces,vrmin,vrmax,dr,vl)
		#vrstep = 0.02; drstep = 0.002
		#sy,vrmint,vrmaxt,drt = reduceAll(sy,traces,vrmin,vrmax,dr,vrstep,drstep)
		
		# Apply Dynamic time warping to get shifts and shifted seismogram
 		e,u,h,sh = getWarpingShifts(ssy,st,sy,traces,vrmin,vrmax,dr,vl,smax=[0.0,0.2])
		print 'nrms '+title+'=',nrms(h,copy(sh.count,inro(sh.first/sh.delta),traces[0]))
		velpanel=None
		if velpanel:
 			e1,u1,h1,sh1 = getWarpingShifts(ssy,st,sy,traces,0.5,1.5,0.5,vl)
 			e2,u2,h2,sh2 = getWarpingShifts(ssy,st,sy,traces,0.9,1.1,0.1,vl)
 			e3,u3,h3,sh3 = getWarpingShifts(ssy,st,sy,traces,0.9,1.1,0.02,vl)
 			tz1,vint1 = updateTDandVelocity(sz,st,u1,sub(tzl,tzl[0]))
 			tz2,vint2 = updateTDandVelocity(sz,st,u2,sub(tzl,tzl[0]))
 			tz3,vint3 = updateTDandVelocity(sz,st,u3,sub(tzl,tzl[0]))
			plotVelPanel(sz,vint1,vint2,vint3,paper=_pres+'velcomppanel')
 		# Get updated time-depth curve and velocity curve
 		tz,vint = updateTDandVelocity(sz,st,u,sub(tzl,tzl[0]))
 		sl = Sampling(len(tr)-len(sy)+1,dt,0.0)

 		# Get horizons and well tops for the trace and well
 		gh = getHorizonsSlice(sx,smy,sz,st,xc,yc,tz)

		ulim=[min(u)*0.5,max(u)*1.5]
		#plotSequences(tr,sy,st,ssy,title='f & g')
		plotSequences(tr,h,st,sh,title='f & h')
		#plotCurve(u,ssy,d="t",hlabel="u")
		#plotMatrix(etran(e),ssy,sl,lim=ulim)
		#plotMatrix(etran(e),ssy,sl,u=u,lim=ulim)
		#plotTDCurves(sz,tzl,tz,paper=_pres+title+'tdcurves')
		#plotCurve(div(1.0,vl),sz,c2=vint,s2=sz,d="z",hlabel="slowness (s/km)",paper=_pres+title+'slcurves')
		plotCurvePanel2(ssy,st,sh,sy,tr,h,lims=[0.0,1.4],paper=_pres+'sigtie')

		# TD AND VELOCITY PLOTS
		plotCurve(vl,sz,c2=vint,s2=sz,d="z",hlabel="velocity (km/s)",paper=_pres+title+'vcurves')
		plotTDandV(sz,vl,vint,tzl,tz,h2label='velocity (km/s)',paper=_pres+title+'vtdplot')
		plotTDandV(sz,div(1.0,vl),div(1.0,vint),tzl,tz,h2label='slowness (s/km)',paper=_pres+title+'stdplot')
		plotCurveHorz(u,ssy,d="t",lim=[0.35,0.55],paper=_pres+title+'shifts')
		
		#lims = [yv-0.65,ssy.first-0.25,yv+0.65,ssy.last+0.25]
		lims = [yv-0.65,0.1,yv+0.65,1.35]
		plot2DLogs(st,smy,xslice,[ssy],[sy],[yv],tops=tops,fmsp=fmsp,lims=lims,\
								sz=[sz],tz=[tzl],gh=gh,paper=_pres+title+'1duntiedwtnh',psize=[723,578])
		plot2DLogs(st,smy,xslice,[ssy],[sy],[yv],lims=lims,paper=_pres+title+'1duntied',psize=[723,578])
		plot2DLogs(st,smy,xslice,[sh],[h],[yv],tops=tops,fmsp=fmsp,lims=lims,\
								sz=[sz],tz=[tz],gh=gh,paper=_pres+title+'1dtiedwtnh',psize=[723,578])
		plot2DLogs(st,smy,xslice,[sh],[h],[yv],lims=lims,paper=_pres+title+'1dtied',psize=[723,578])
		plotCurvePanel2(ssy,st,sh,sy,tr,h,lims=[0.0,1.4])
		plotSlice(smy,st,xslice,lims=lims,paper=_pres+title+'1dseismic',psize=[723,578])

	if Warpsf:
		#title='Simple method- no multiples or Q'
		#print '\n'+title
		#title='smpl'
		#ssyc,syc,tzlc = getSimpleSeismogram(sz,st,vl,ail,fp,phase=phase,normalize=True)
		#warpAndUpdatef(ssyc,syc,tzlc,title)
		#title='Propagator method 1- multiples added'
		#print '\n'+title
		#ssyp,syp,tzlp = getPropagatorSeismogram(sz,dt,vl,dl,fp,q1,phase=0,normalize=True)
		#warpAndUpdatef(ssyp,syp,tzlp,title)
		#title='Propagator method 2- no multiples, attenuation added'
		#print '\n'+title
		#ssyp,syp,tzlp = getPropagatorSeismogram(sz,dt,vl,dl,fp,q2,phase=0,normalize=True,nomult=True)
		#warpAndUpdatef(ssyp,syp,tzlp,title)
		title='Propagator method 3- multiples and attenuation added'
		print '\n'+title
		title = 'prop'
		ssyp,syp,tzlp = getPropagatorSeismogram(sz,dt,vl,dl,fp,q2,phase=phase,normalize=True)
		warpAndUpdatef(ssyp,syp,tzlp,title)

	if Comparisonsf:
		#ssyp,syp,tzlp = getPropagatorSeismogram(sz,dt,vl,dl,fp,q1,phase=0,normalize=False,nomult=True,nocut=True)
		#ssyc,syc,tzlc = getSimpleSeismogram(sz,st,vl,ail,fp,phase=0,normalize=False)
		#plotCompare(ssyc,ssyp,syc,syp,'S(blk) & P(red) nm')

		#ssyp,syp,tzlp = getPropagatorSeismogram(sz,dt,vl,dl,fp,q1,phase=0,normalize=False,nocut=True)
		#syp = div(syp,2.0)
		#plotCompare(ssyc,ssyp,syc,syp,'S(blk) & P(red) m')

		### Seismogram Panel	
		ssyp0,syp0,tzlp0 = getSimpleSeismogram(sz,st,vl,ail,fp,phase=0,normalize=False)
		ssyp1,syp1,tzlp1 = getPropagatorSeismogram(sz,dt,vl,dl,fp,q1,phase=0,normalize=False,nomult=True)
		ssyp2,syp2,tzlp2 = getPropagatorSeismogram(sz,dt,vl,dl,fp,q1,phase=0,normalize=False)
		ssyp3,syp3,tzlp3 = getPropagatorSeismogram(sz,dt,vl,dl,fp,q2,phase=0,normalize=False)
		plotSeismogramPanel(ssyp1,ssyp2,ssyp3,syp1,div(syp2,2.0),syp3,paper=_pres+'seismopanel')
		plotCompare(ssyp2,ssyp3,div(syp2,5.0),syp3,'S(blk) & P(red) m')
		plotCompare(ssyp1,ssyp3,div(syp1,5.0),syp3,'S(blk) & P(red) m')

		#ssyp,syp,tzlp = getPropagatorSeismogram(sz,dt,vl,dl,fp,q2,phase=0,normalize=True,nomult=True,nocut=True)
		#plotCompare(ssyc,ssyp,syc,syp,'S(blk) & P(red) nm, q=100')

		#ssyp,syp,tzlp = getPropagatorSeismogram(sz,dt,vl,dl,fp,q2,phase=0,normalize=True,nomult=True,nocut=True)
		#plotCompare(ssyc,ssyp,syc,syp,'S(blk) & P(red) m, q=100')

		#ssyp,syp,tzlp = getPropagatorSeismogram(sz,dt,vl,dl,fp,q,phase=0,normalize=True,geophone=True,nomult=True,nocut=True)
		#plotCompare(ssyc,ssyp,syc,syp,'S(blk) & P(red) geophone, nm, norm')

		#ssyp,syp,tzlp = getPropagatorSeismogram(sz,dt,vl,dl,fp,q,phase=0,normalize=False,nomult=True,anisotropic=True,nocut=True)
		#ssyc,syc,tzlc = getSimpleSeismogram(sz,st,vl,ail,fp,phase=0,normalize=False)
		#plotCompare(ssyc,ssyp,syc,syp,'Simple vs. Prop anis, no mult')

################################################################################
## Multiple well tie methods
vlogs=None
vinterps=True
dinterps=True
ginterps=None
def goMultiWellf():
	#printWellLengths()
	# Params:
	fp = 35
	q  = 100
	nlr = 1.0
	# Time warping
	vrmin = 0.98 # v*/v min 0.5 #0.99,1.01,.08,195 
	vrmax = 1.02 # v*/v max inf
	dr    = 0.01 # 1D smoothness vertically
	# Image warping
	dr1   = 0.02 # 2D smoothness vertically
	dr2   = 0.10 # 2D smoothness horizontally
	dr3   = 0.10 # 3D smoothness horizontally
	rmin1,rmax1 = -0.08,0.08 # 2D vertical constraint
	rmin2,rmax2 = -0.30,0.30 # 2D horizontal constraint
	rmin3,rmax3 = -0.30,0.30 # 3D horizontal constraint
	smin,smax = -0.12,0.02 # min and max time shifts for 2D/3D warping
	#dr = computedr(vrmin,vrmax,dr)

	# 1. Get the 3D seismic image
	tpimg,s3,s2,s1 = getImage(normalize=True)
	dt = s1.delta
	smin,smax = int(floor(smin/dt)),int(ceil(smax/dt))

	# 2. Get the wells, coordinates, and tops
	################
	# TESTING FOR 3D PREP:
	# wells less than 127: 8,6,7,13,2,17
	####wells = [0,1,2,3,4,6,7,8,9,10,11,12,13,15,16,17]
	#wells = [0,1,3,4,9,10,11,12,15,16]
	## potential slices:
	#wells =  [3,6,10,11,12] 
	#wells = [3,16,15,4,9] # All long wells (>200)
	######wells = [3,16,15,4,9,12,17] # 3D
	## test lines all with 3:
	#wells = [8,10,3,6,0,12,4,11,9] # 11,10,0,6
	#wells = [3,12,4,9] # 8,11,10,0,6
	#wells = [3,1,11,9] # 1,11 no tie
	#wells = [3,16,15,17]#7,13,2,17] 
	################
	# 2D slices for display
	##wells = [3,16,15]
	##wells = [3,12,4,17]
	#wells = [3,4,9]
	wells = [3,15,4,9] # for 2D figs: 1.02,.98,.02,.02,.05,-.10,.10,-1.,1.,-.05,.05,ph153,fp35,q80
	#wells = [3,4,9]
	#phases = [195]*len(wells)
	##phases = [306]*len(wells)
	phases = [49]*len(wells)
	#phases = [0]*len(wells)
	wells = sortWells(wells)
	tops,fmsp,datums = getWellTops(wells)
	ail,vl,dl,gl,sz,x2w,x3w = getMultWells(wells,datums)

	# 3. Get the wavelets (if desired) and synthetics
	ssy,sy,tzl,ztl,sztl = getSeismograms(sz,s1,vl,dl,fp,q,sz,phase=phases)

	# 4. Extract the image between the wells
	x2s,x3s,sn = getCurveThruWells(x2w,x3w,s2.delta)
	f = getImageAlongCurve(tpimg,x2s,x3s,s1,s2,s3)
	imgh = wtu.getSlice(470,tpimg)
	plotCoordSlice(s2,s3,x2s,x3s,x2w,x3w,imgh,paper=_pres+'hwellslice')
	lims = [sn.first,0.10,sn.last,1.35]
	
	# 5. Get horizons 
	gh = getMultiHorizons(x2s,x3s,s2,sz,s1)
	
	plotSlice(sn,s1,f,lims=lims,paper=_pres+'seisextract',psize=[712,524])
	plot2DLogs(s1,sn,f,ssy,sy,x2w,lims=lims,paper='2dwuntindv')#,psize=[712,394])
	plot2DLogs(s1,sn,f,ssy,sy,x2w,lims=lims,\
		gh=gh,tops=tops,fmsp=fmsp,sz=sz,tz=tzl,paper='2dwuntindvwh')#,psize=[712,394])

	# 6. Warp the synthetics individually
	u,h,sh,tz,vint,zt,szt = warpIndividually(tpimg,s2,s3,\
												 	 ssy,s1,sz,fp,x2w,x3w,sy,tzl,vl,vrmin,vrmax,dr,smax=[0.00,0.15])	
	plot2DLogs(s1,sn,f,sh,h,x2w,lims=lims,title='wells tied individually',paper='2dwtindv')#,psize=[712,394])
	plot2DLogs(s1,sn,f,sh,h,x2w,lims=lims,title='wells tied individually',\
							gh=gh,tops=tops,fmsp=fmsp,sz=sz,tz=tz,paper='2dwtindvwh')#,psize=[712,394])

	print '2D Warp refinement'
	# 7. Interpolate and warp the synthetics
	bgs = 0.8 
	tensors = makeTensors(f)

	gb = interpolateWells(s1,f,sy,ssy,sn,x2w,bgs,ts=tensors)
	plotSlice(sn,s1,gb,lims=lims,title='Interpolated untied wells',paper=_pres+'interpsyuntied')
	ga = interpolateWells(s1,f,h,sh,sn,x2w,bgs,ts=tensors)
	plotSlice(sn,s1,ga,lims=lims,title='interpolated single tie synthetics',paper=_pres+'interp1sytied')
	plotSlice(sn,s1,f,tens=tensors,lims=lims,paper=_pres+'seistens')
	gs = interpolateWells(s1,f,h,sh,sn,x2w,bgs,ts=tensors,scat=True)
	plotSlice(sn,s1,gs,lims=lims,paper=_pres+'tie1scatwells2')
	gs = interpolateWells(s1,f,sy,ssy,sn,x2w,bgs,ts=tensors,scat=True)
	plotSlice(sn,s1,gs,lims=lims,paper=_pres+'untiescatwells2')

	h2d,u2d,hn,shn,tzn,vintn,ztn,sztn,gs = warpAllin2D(\
				s1,sn,f,ga,rmin1,rmax1,rmin2,rmax2,smin,smax,dr1,dr2,x2w,h,sh,vl,tz,sz,ssy)
	
	# Single well warp and then multi-well warp
	h2d,u2d,hn,shn,tzn2,vintn2,ztn2,sztn2,gs = warpAllin2D(\
			s1,sn,f,ga,rmin1,rmax1,rmin2,rmax2,smin,smax,dr1,dr2,x2w,h,sh,vint,tz,sz,ssy)
	gc = interpolateWells(s1,f,hn,shn,sn,x2w,bgs,ts=tensors)
	plotSlice(sn,s1,gc,lims=lims,title='interpolated single then multi tied synthetics',paper=_pres+'interp1and2fin')
	plot2DLogs(s1,sn,f,shn,hn,x2w,lims=lims,title='wells tied with 2D refinement',paper=_pres+'2dp1dt')
	
	# Warp all at once
	h2d,u2d,hn,shn,tzn,vintn,ztn,sztn,gs = warpAllin2Donly(tpimg,s1,s2,s3,sn,f,gb,\
			rmin1,rmax1,rmin2,rmax2,smin,smax,dr1,dr2,x2w,x3w,sy,ssy,vl,tzl,sz,fp)
	gb = interpolateWells(s1,f,hn,shn,sn,x2w,bgs,ts=tensors)
	plotSlice(sn,s1,gb,lims=lims,title='interpolated multi tied synthetics',paper=_pres+'interp2fin')
	plot2DLogs(s1,sn,f,shn,hn,x2w,lims=lims,title='wells tied with 2D only',paper=_pres+'2dto')

	#plotSlice(sn,s1,h2d,lims=lims,title='interpolated 2D tie synthetics',paper=_pres+'interp2sytied')
	#plotSlice(sn,s1,u2d,cbar=True,title='Interpolated tied wells 2D shifts',paper=_pres+'interpu')
	#plot2DLogs(s1,sn,f,shn,hn,x2w,lims=lims,title='wells tied with 2D refinement',psize=[712,394],paper='2dt2dv')
	#plot2DLogs(s1,sn,f,shn,hn,x2w,lims=lims,title='wells tied with 2D refinement',\
	#						gh=gh,tops=tops,fmsp=fmsp,sz=sz,tz=tzn,psize=[712,394],paper='2dt2dvwh')
	bgs = 0.5
	if dinterps:
		ztlmg = interpolateWells(s1,f,ztl,sztl,sn,x2w,bgs,ts=tensors)
		plotSlice(sn,s1,ztlmg,title='log depths',cbar="Depth (km)",cmap=jet,lims=lims,paper=_pres+'ztn0')
		ztomg = interpolateWells(s1,f,zt,szt,sn,x2w,bgs,ts=tensors)
		plotSlice(sn,s1,ztomg,title='interpolated 1D tie depths',cbar="Depth (km)",cmap=jet,lims=lims,paper=_pres+'ztn1')
		ztimg = interpolateWells(s1,f,ztn,sztn,sn,x2w,bgs,ts=tensors)
		plotSlice(sn,s1,ztimg,cmap=jet,lims=lims,cbar="Depth (km)",title='2D warp depths',paper=_pres+'ztn2')
		ztimg2 = interpolateWells(s1,f,ztn2,sztn2,sn,x2w,bgs,ts=tensors)
		plotSlice(sn,s1,ztimg2,cmap=jet,lims=lims,cbar="Depth (km)",title='2D refined depths',paper=_pres+'ztn2r')
		#plotSlice(sn,s1,sub(ztimg,ztomg),cmap=jet,lims=lims,cbar=True,title='difference of refined depths')
		
	if vinterps:
		bgs = 0.8
		#tensors = makeTensors(f,p0=0.6,p1=1.0)
		vlmg = interpolateWells(s1,f,vl,ssy,sn,x2w,bgs,tzl,ts=tensors)
		plotSlice(sn,s1,vlmg,title='interpolated original log velocities',cbar='Velocity (km/s)',cmap=jet,lims=lims,clips=[2.5,6.8],paper='vint0tie')
		vintmg = interpolateWells(s1,f,vint,sh,sn,x2w,bgs,tzl,ts=tensors)
		plotSlice(sn,s1,vintmg,title='interpolated 1D tie velocities',cbar='Velocity (km/s)',cmap=jet,lims=lims,paper='vint1tie',clips=[2.5,6.8])
		vimg = interpolateWells(s1,f,vintn2,shn,sn,x2w,bgs,tzn2,ts=tensors)
		plotSlice(sn,s1,vimg,cmap=jet,lims=lims,cbar='Velocity (km/s)',title='interpolated 2D refined velocities',paper='vint2rtie',clips=[2.5,6.8])
		vimg2 = interpolateWells(s1,f,vintn,shn,sn,x2w,bgs,tzn,ts=tensors)
		plotSlice(sn,s1,vimg2,cmap=jet,lims=lims,cbar='Velocity (km/s)',title='interpolated 2D warp velocities',paper='vint2tie',clips=[2.5,6.8])
		#plotSlice(sn,s1,vimg,cmap=jet,lims=lims,cbar=True,title='interpolated 2D refined velocities',clips=[2.5,7.0],tens=tensors)

	if ginterps:
		gimg = interpolateWells(s1,f,gl,shn,sn,x2w,bgs,tzn,ts=tensors)
		plotSlice(sn,s1,gimg,cmap=jet,lims=lims,cbar=True,title='interpolated gamma logs')
		plot2DLogs(s1,sn,gimg,shn,hn,x2w,lims=lims,title='interpolated gamma logs'\
			,cbar=True,cmap=jet,gh=gh,tops=tops,fmsp=fmsp,sz=sz,tz=tzn,wwidth=1,nolog=True)
	
	if vlogs:
		for i in range(len(wells)):
			plotCurve(vl[i],sz[i],c2=vint[i],s2=sz[i],d="z",hlabel="velocity (km/s)")
			plotTDCurves(sz[i],tzl[i],tz[i])	
		for i in range(len(wells)):
			plotCurve(vl[i],sz[i],c2=vintn[i],s2=sz[i],d="z",hlabel="velocity (km/s)")
			plotTDCurves(sz[i],tzl[i],tzn[i])
			plotCurve(ztn[i],sztn[i],d="t",hlabel="zt")





	# reinterpolate wells with nearest neighbor to check consitency
	#ga = interpolateWells(s1,f,hn,shn,sn,x2w,bgs=0.0,ts=tensors)
	#plotSlice(sn,s1,ga,tens=tensors,lims=lims,title='interpolated tied synthetics')#,slides=_pres+'2dinterptens-4.0')
	#plotSlice(sn,s1,ga,lims=lims,title='interpolated tied synthetics')#,slides=_pres+'2dinterpw1ds')


# TODO Make methods for calling interpolated logs/seismic and plotting

#------------------------------------------------------------------------------#

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())



