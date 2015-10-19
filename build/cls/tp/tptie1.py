# Single well ties
# @author Andrew Munoz, Colorado School of Mines
# @version 04.10.2013

from dtw import *
from tputils import *
from wtutils import *
from imports import *

well = [3]
_pres = '_cwps13'

# Global variables set-up
_ID,_set = getWellFromSet(well[0]) 
ail,yv,xv,zv,vl,dl,sz = getLogs(_ID,_set)
tops,fmsp,datum = getWellTops(well)
tpimg,s3,s2,s1  = getImage(normalize=True)
traces,yc,xc = getTraces(tpimg,yv,xv,nr=3) 
gh = getXSliceHorizons(xc,sz)
xslice = tpimg[xc]
dt = s1.delta

def main(args):
	goWarps()
	#goComparisons()

def goWarps():
	phase = 44 #	well 3 prop
	q1 = 10000000
	q2 = 100
	vrmin = 0.98 # v*/v min 0.5 
	vrmax = 1.02 # v*/v max inf
	dr    = 0.02 # resolution

	simple = True
	monly  = None
	qonly  = None
	prop   = True

	if simple:
		title='Simple method- no multiples or Q'
		print '\n'+title
		title='smpl '
		ssyc,syc,tzlc = getSimpleSeismogram(sz,vl,ail,phase=phase,normalize=True)
		warpAndUpdate(ssyc,syc,tzlc,title,vrmin,vrmax,dr)
	if monly:
		title='Propagator method 1- multiples added'
		print '\n'+title
		ssyp,syp,tzlp = getPropagatorSeismogram(sz,vl,dl,q1,phase=0,normalize=True)
		warpAndUpdate(ssyp,syp,tzlp,title,vrmin,vrmax,dr)
	if qonly:
		title='Propagator method 2- no multiples, attenuation added'
		print '\n'+title
		ssyp,syp,tzlp = getPropagatorSeismogram(sz,vl,dl,q2,phase=0,normalize=True,nomult=True)
		warpAndUpdate(ssyp,syp,tzlp,title,vrmin,vrmax,dr)
	if prop:
		title='Propagator method 3- multiples and attenuation added'
		print '\n'+title
		title = 'prop'
		ssyp,syp,tzlp = getPropagatorSeismogram(sz,vl,dl,q2,phase=phase,normalize=True)
		warpAndUpdate(ssyp,syp,tzlp,title,vrmin,vrmax,dr)
	#plotCurve(syp,ssyp,syp2,ssyp2,d="t",title="compare norm")


def warpAndUpdate(ssy,sy,tzl,title,vrmin,vrmax,dr,phaser=None):

	# Rotate the phase of the seismogram to find the optimum phase
	#sy = rotatePhase(sy,traces,vrmin,vrmax,dr)
	#sy = rotateAllPhase(sy,traces,vrmin,vrmax,dr)
			
	tr = traces[0]
	#tzl = sub(tzl,tzl[0])
	# Apply Dynamic time warping to get shifts and shifted seismogram
 	e,u,h,sh = getWarpingShifts(ssy,s1,sy,traces,vrmin,vrmax,dr,vl,sz=sz)

	# Get updated time-depth curve and velocity curve
	#tz,vint,pv = updateTDandVelocity(sz,u,sub(tzl,tzl[0]))
	tz,vint,pv = updateTDandVelocity(sz,u,tzl)

	#goWarpPlots(tr,tzl,tz,vint,u,sy,ssy,h,sh,title)
	goPaper(tr,tzl,tz,vint,u,sy,ssy,h,sh,title)
	#goSlides(tr,tzl,tz,vint,u,sy,ssy,h,sh,title)


def goComparisons():
		
		ssyp,syp,tzlp = getPropagatorSeismogram(sz,vl,dl,q1,phase=0,normalize=False,nomult=True,nocut=True)
		ssyc,syc,tzlc = getSimpleSeismogram(sz,vl,ail,phase=0,normalize=False)
		plotCompare(ssyc,ssyp,syc,syp,'S(blk) & P(red) nm')

		ssyp,syp,tzlp = getPropagatorSeismogram(sz,vl,dl,q1,phase=0,normalize=False,nocut=True)
		syp = div(syp,2.0)
		plotCompare(ssyc,ssyp,syc,syp,'S(blk) & P(red) m')

		#ssyp,syp,tzlp = getPropagatorSeismogram(sz,vl,dl,q2,phase=0,normalize=True,nomult=True,nocut=True)
		#plotCompare(ssyc,ssyp,syc,syp,'S(blk) & P(red) nm, q=100')

		#ssyp,syp,tzlp = getPropagatorSeismogram(sz,vl,dl,q2,phase=0,normalize=True,nomult=True,nocut=True)
		#plotCompare(ssyc,ssyp,syc,syp,'S(blk) & P(red) m, q=100')

		#ssyp,syp,tzlp = getPropagatorSeismogram(sz,vl,dl,q,phase=0,normalize=True,geophone=True,nomult=True,nocut=True)
		#plotCompare(ssyc,ssyp,syc,syp,'S(blk) & P(red) geophone, nm, norm')

		#ssyp,syp,tzlp = getPropagatorSeismogram(sz,vl,dl,q,phase=0,normalize=False,nomult=True,anisotropic=True,nocut=True)
		#ssyc,syc,tzlc = getSimpleSeismogram(sz,vl,ail,phase=0,normalize=False)
		#plotCompare(ssyc,ssyp,syc,syp,'Simple vs. Prop anis, no mult')

def goWarpPlots(tr,tzl,tz,vint,u,sy,ssy,h,sh,title):
	#plotSequences(tr,sy,s1,ssy,title='f & g')
	plotSequences(tr,h,s1,sh,title='f & h')
	#plotCurve(u,ssy,d="t",hlabel="u")
	#plotMatrix(etran(e),ssy,sl,lim=ulim)
	#plotMatrix(etran(e),ssy,sl,u=u,lim=ulim)
	plotTDCurves(sz,tzl,tz)
	plotCurve(vl,sz,c2=vint,s2=sz,d="z",hlabel="v (km/s)")
	
	lims = [3.0,0.1,4.4,1.8]
	#plot2DLogs(s1,s2,xslice,[ssy],[sy],[yv],tops=tops,fmsp=fmsp,lims=lims,\
	#					sz=[sz],tz=[tzl],gh=gh,title=title+'untied')
	plot2DLogs(s1,s2,xslice,[ssy],[sy],[yv],lims=lims,title=title+'untied',wwidth=1)
	#plot2DLogs(s1,s2,xslice,[sh],[h],[yv],tops=tops,fmsp=fmsp,lims=lims,\
	#					sz=[sz],tz=[tz],gh=gh,title=title+'tied')
	plot2DLogs(s1,s2,xslice,[sh],[h],[yv],lims=lims,title=title+'tied',wwidth=1)
	plotSlice(s2,s1,xslice,lims=lims)
	plotCurvePanel2(ssy,s1,sh,sy,tr,h,lims=[0.0,1.4])
	#ssyr,syr,tzlr = getPropagatorSeismogram(sz,vint,dl,q2,phase=phase,normalize=True)
	#plotCurve(syr,ssyr,h,sh,d="t")

def goPaper(tr,tzl,tz,vint,u,sy,ssy,h,sh,title):
	velpanel=True
	if velpanel:
 		e1,u1,h1,sh1 = getWarpingShifts(ssy,s1,sy,traces,0.5,1.5,0.5,vl) # Classic
 		e2,u2,h2,sh2 = getWarpingShifts(ssy,s1,sy,traces,0.9,1.1,0.1,vl) # Constrained
 		e3,u3,h3,sh3 = getWarpingShifts(ssy,s1,sy,traces,0.9,1.1,0.02,vl) # Smoothly Const
 		tz1,vint1,pv1 = updateTDandVelocity(sz,u1,tzl)
 		tz2,vint2,pv2 = updateTDandVelocity(sz,u2,tzl)
 		tz3,vint3,pv3 = updateTDandVelocity(sz,u3,tzl)
		plotVelPanel(sz,vl,vl,vl,vint1,vint2,vint3,avel=True,paper=title+' velcomppanel')

	plotLogPanel(sz,vl,dl,reflectivity(ail),paper='plogpanel')

	##plotSequences(tr,sy,s1,ssy,title='f & g')
	#plotSequences(tr,h,s1,sh,title='f & h')
	##plotCurve(u,ssy,d="t",hlabel="u")
	##plotMatrix(etran(e),ssy,sl,lim=ulim)
	##plotMatrix(etran(e),ssy,sl,u=u,lim=ulim)
	##plotTDCurves(sz,tzl,tz,paper=_pres+title+'tdcurves')
	##plotCurve(div(1.0,vl),sz,c2=vint,s2=sz,d="z",hlabel="slowness (s/km)",paper=_pres+title+'slcurves')
	#plotCurvePanel2(ssy,s1,sh,sy,tr,h,lims=[0.0,1.4],paper=_pres+title+'sigtie')

	## TD AND VELOCITY PLOTS
	#plotCurve(vl,sz,c2=vint,s2=sz,d="z",hlabel="Velocity (km/s)",paper=_pres+title+'vcurves')
	#plotTDandV(sz,vl,vint,tzl,tz,h2label='Velocity (km/s)',paper=_pres+title+'vtdplot')
	#plotTDandV(sz,div(1.0,vl),div(1.0,vint),tzl,tz,h2label='slowness (s/km)',paper=_pres+title+'stdplot')
	#plotCurveHorz(u,ssy,d="t",lim=[0.35,0.55],paper=_pres+title+'shifts')
	#
	## Warp and tie plots
	##lims = [yv-0.65,ssy.first-0.25,yv+0.65,ssy.last+0.25]
	#lims = [yv-0.65,0.1,yv+0.65,1.35]
	##plot2DLogs(s1,s2,xslice,[ssy],[sy],[yv],tops=tops,fmsp=fmsp,lims=lims,\
	##						sz=[sz],tz=[tzl],gh=gh,paper=_pres+title+'1duntiedwtnh',psize=[723,578])
	#plot2DLogs(s1,s2,xslice,[ssy],[sy],[yv],lims=lims,paper=_pres+title+'1duntied',psize=[723,578],wwidth=1)
	##plot2DLogs(s1,s2,xslice,[sh],[h],[yv],tops=tops,fmsp=fmsp,lims=lims,\
	##						sz=[sz],tz=[tz],gh=gh,paper=_pres+title+'1dtiedwtnh',psize=[723,578])
	#plot2DLogs(s1,s2,xslice,[sh],[h],[yv],lims=lims,paper=_pres+title+'1dtied',psize=[723,578],wwidth=1)
	#plotCurvePanel2(ssy,s1,sh,sy,tr,h,lims=[0.0,1.4])
	#plotSlice(s2,s1,xslice,lims=lims,paper=_pres+title+'1dseismic',psize=[723,578])

def goSlides(tr,tzl,tz,vint,u,sy,ssy,h,sh,title):
	velpanel=None
	if velpanel:
 		e1,u1,h1,sh1 = getWarpingShifts(ssy,s1,sy,traces,0.5,1.5,0.5,vl)
 		e2,u2,h2,sh2 = getWarpingShifts(ssy,s1,sy,traces,0.9,1.1,0.1,vl)
 		e3,u3,h3,sh3 = getWarpingShifts(ssy,s1,sy,traces,0.9,1.1,0.02,vl)
 		tz1,vint1,pv1 = updateTDandVelocity(sz,u1,sub(tzl,tzl[0]))
 		tz2,vint2,pv2 = updateTDandVelocity(sz,u2,sub(tzl,tzl[0]))
 		tz3,vint3,pv3 = updateTDandVelocity(sz,u3,sub(tzl,tzl[0]))
		plotVelPanel(sz,vint1,vint2,vint3,slides=_pres+'velcomppanel')

	plotLogPanel(sz,vl,dl,reflectivity(ail),slides='plogpanel')

	#plotSequences(tr,sy,s1,ssy,title='f & g')
	plotSequences(tr,h,s1,sh,title='f & h')
	#plotCurve(u,ssy,d="t",hlabel="u")
	#plotMatrix(etran(e),ssy,sl,lim=ulim)
	#plotMatrix(etran(e),ssy,sl,u=u,lim=ulim)
	plotTDCurves(sz,tzl,tz,slides=_pres+title+'tdcurves')
	#plotCurve(div(1.0,vl),sz,c2=vint,s2=sz,d="z",hlabel="slowness (s/km)",slides=_pres+title+'slcurves')
	plotCurvePanel2(ssy,s1,sh,sy,tr,h,lims=[0.0,1.4],slides=_pres+'sigtie')

	# TD AND VELOCITY PLOTS
	plotCurve(vl,sz,c2=vint,s2=sz,d="z",hlabel="Velocity (km/s)",slides=_pres+title+'vcurves')
	plotTDandV(sz,vl,vint,tzl,tz,h2label='Velocity (km/s)',slides=_pres+title+'vtdplot')
	plotTDandV(sz,div(1.0,vl),div(1.0,vint),tzl,tz,h2label='Slowness (s/km)',slides=_pres+title+'stdplot')
	plotCurveHorz(u,ssy,d="t",lim=[0.35,0.55],slides=_pres+title+'shifts')
	
	# Warp and tie plots
	#lims = [yv-0.65,ssy.first-0.25,yv+0.65,ssy.last+0.25]
	lims = [yv-0.65,0.1,yv+0.65,1.35]
	#plot2DLogs(s1,s2,xslice,[ssy],[sy],[yv],tops=tops,fmsp=fmsp,lims=lims,\
	#						sz=[sz],tz=[tzl],gh=gh,slides=_pres+title+'1duntiedwtnh',psize=[723,578])
	plot2DLogs(s1,s2,xslice,[ssy],[sy],[yv],lims=lims,slides=_pres+title+'1duntied',psize=[723,578],wwidth=1)
	#plot2DLogs(s1,s2,xslice,[sh],[h],[yv],tops=tops,fmsp=fmsp,lims=lims,\
	#						sz=[sz],tz=[tz],gh=gh,slides=_pres+title+'1dtiedwtnh',psize=[723,578])
	plot2DLogs(s1,s2,xslice,[sh],[h],[yv],lims=lims,slides=_pres+title+'1dtied',psize=[723,578],wwidth=1)
	plotCurvePanel2(ssy,s1,sh,sy,tr,h,lims=[0.0,1.4])
	plotSlice(s2,s1,xslice,lims=lims,slides=_pres+title+'1dseismic',psize=[723,578])



#------------------------------------------------------------------------------#

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

