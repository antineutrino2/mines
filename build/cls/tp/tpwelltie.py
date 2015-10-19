# Dynamic Time Warping plotter for teapot dome data
# This runs all Dtw methods and plots
# @author Andrew Munoz, Colorado School of Mines
# @version 09.11.2012

from tp import * 
from wt import *
from dtw import DynamicWarpingWT
from dtw import DynamicWarpingWTM
from dtw import SyntheticSeismogram
from dtw import SeismicModel1D
from tputils import *
from wtutils import *
from imports import *
from edu.mines.jtk.dsp import DynamicWarping as DynamicWarpingO

# My Mac
_tpDir = "/Users/amunoz/Home/data/tp/"
# Backus
#_tpDir = "/data/amunoz/tp/"

_pngDir = './tppics/2013/'
#_pres = 'testing_'
_pres = 'sega13_'
#_pres = 'cwpr13_'
#_pres = 'ctmsp13_'
#olddatum = 0.4041648
#seismicDatum = 1.98120-0.003048 # 6500-10 ft 
seismicDatum = 1.6764-0.003048 # 5500-10 ft
dataDir = _tpDir+"wellties/"
topsDir = _tpDir+"csm/welllogs/welltops/"
wtu = WTutils()
we = WaveletExtractor()
dz = 0.0001524
sigmanorm = 100
fp = 35

# Run types for a single well
Comparisons =	False
Warps = True

#fms = ["F2WC","DKOT","LKOT","CRMT","RDPK","A Sand","C1 Dolo"]
fms = ["F2WC","DKOT","A Sand","C1 Dolo"]

def main(args):
	print '\n'
	#goSingleWell()
	goMultiWell()
	#printWellIDList()

def goSingleWell():
	#Params
	well = [3]
	phase = 112 #	well 3 prop
	_ID,_set = getWellFromSet(well[0]) 
	ail,yv,xv,zv,vl,dl = getLogs(_ID,_set)
	tops,fmsp,datum = getWellTops(well)
	nz = len(vl)
	fz = zv#+datum[0]
	szd = Sampling(nz,dz,fz) # in km
	sz = Sampling(nz,dz,zv) # in km
	fp = 35 #(Hz) peak frequency
	tpimg,s3,s2,s1 = getImage(normalize=True)
	traces,sx,smy,st,xslice,yc,xc = getTrace(tpimg,s1,s2,s3,yv,xv,nr=3) #(x3,x2,x1)=(x,y,z)
	tr = traces[0]
	dt = st.delta
	q1 = 10000000
	q2 = 80
	vrmin = 0.98 # v*/v min 0.5 
	vrmax = 1.02 # v*/v max inf
	dr    = 0.02 # resolution
	
	def warpAndUpdate(ssy,sy,tzl,title):
		# Rotate the phase of the seismogram to find the optimum phase
		#sy = rotatePhase(sy,traces,vrmin,vrmax,dr)
		#sy = rotateAllPhase(sy,traces,vrmin,vrmax,dr)
		#sy = rotatePhaseHv(st,sz,sx,smy,xc,yc,ssy,sy,tzl,traces,vrmin,vrmax,dr,vl)
		#vrstep = 0.02; drstep = 0.002
		#sy,vrmint,vrmaxt,drt = reduceAll(sy,traces,vrmin,vrmax,dr,vrstep,drstep)
		
		# Apply Dynamic time warping to get shifts and shifted seismogram
 		e,u,h,sh = getWarpingShifts(ssy,st,sy,traces,vrmin,vrmax,dr,vl,sz=sz,lv=2.7)
		print 'nrms '+title+'=',nrms(h,copy(sh.count,inro(sh.first/sh.delta),traces[0]))

 		# Get updated time-depth curve and velocity curve
 		tz,vint = updateTDandVelocity(sz,st,u,sub(tzl,tzl[0]))
 		sl = Sampling(len(tr)-len(sy)+1,dt,0.0)

 		# Get horizons and well tops for the trace and well
		gh = getXSliceHorizons(xc,sx,sz,st)

		pred = printError(st,sh,h,traces[0])	

		ulim=[min(u)*0.5,max(u)*1.5]
		#plotSequences(tr,sy,st,ssy,title='f & g')
		plotSequences(tr,h,st,sh,title='f & h')
		#plotCurve(u,ssy,d="t",hlabel="u")
		#plotMatrix(etran(e),ssy,sl,lim=ulim)
		#plotMatrix(etran(e),ssy,sl,u=u,lim=ulim)
		plotTDCurves(sz,tzl,tz)
		plotCurve(vl,sz,c2=vint,s2=sz,d="z",hlabel="v (km/s)")
		
		lims = [3.0,0.1,4.4,1.8]
		plot2DLogs(st,smy,xslice,[ssy],[sy],[yv],tops=tops,fmsp=fmsp,lims=lims,\
							sz=[sz],tz=[tzl],gh=gh,title=title+'untied')
		plot2DLogs(st,smy,xslice,[ssy],[sy],[yv],lims=lims,title=title+'untied')
		plot2DLogs(st,smy,xslice,[sh],[h],[yv],tops=tops,fmsp=fmsp,lims=lims,\
							sz=[sz],tz=[tz],gh=gh,title=title+'tied')
		plot2DLogs(st,smy,xslice,[sh],[h],[yv],lims=lims,title=title+'tied')
		plotSlice(smy,st,xslice,lims=lims)
		plotCurvePanel2(ssy,st,sh,sy,tr,h,lims=[0.0,1.4])
		#ssyr,syr,tzlr = getPropagatorSeismogram(sz,dt,vint,dl,fp,q2,phase=phase,normalize=True)
		#plotCurve(syr,ssyr,h,sh,d="t")

	if Warps:
		title='Simple method- no multiples or Q'
		print '\n'+title
		title='smpl '
		ssyc,syc,tzlc = getSimpleSeismogram(szd,st,vl,ail,fp,phase=phase,normalize=True)
		warpAndUpdate(ssyc,syc,tzlc,title)
		#title='Propagator method 1- multiples added'
		#print '\n'+title
		#ssyp,syp,tzlp = getPropagatorSeismogram(sz,dt,vl,dl,fp,q1,phase=0,normalize=True)
		#warpAndUpdate(ssyp,syp,tzlp,title)
		#title='Propagator method 2- no multiples, attenuation added'
		#print '\n'+title
		#ssyp,syp,tzlp = getPropagatorSeismogram(sz,dt,vl,dl,fp,q2,phase=0,normalize=True,nomult=True)
		#warpAndUpdate(ssyp,syp,tzlp,title)
		title='Propagator method 3- multiples and attenuation added'
		print '\n'+title
		title = 'prop'
		ssyp,syp,tzlp = getPropagatorSeismogram(szd,dt,vl,dl,fp,q2,phase=phase,normalize=True)
		warpAndUpdate(ssyp,syp,tzlp,title)
		ssyp,syp,tzlp = getPropagatorSeismogram(szd,dt,vl,dl,fp,q2,phase=phase,normalize=True,nocut=True)
		plotCurve(syc,ssyc,syp,ssyp,d="t",title="compare norm")

	if Comparisons:
		ssyp,syp,tzlp = getPropagatorSeismogram(sz,dt,vl,dl,fp,q1,phase=0,normalize=False,nomult=True,nocut=True)
		ssyc,syc,tzlc = getSimpleSeismogram(sz,st,vl,ail,fp,phase=0,normalize=False)
		plotCompare(ssyc,ssyp,syc,syp,'S(blk) & P(red) nm')

		ssyp,syp,tzlp = getPropagatorSeismogram(sz,dt,vl,dl,fp,q1,phase=0,normalize=False,nocut=True)
		syp = div(syp,2.0)
		plotCompare(ssyc,ssyp,syc,syp,'S(blk) & P(red) m')

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

def goMultiWell():
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
	dr1   = 0.01 # 2D smoothness vertically
	dr2   = 0.25 # 2D smoothness horizontally
	dr3   = 0.10 # 3D smoothness horizontally
	rmin1,rmax1 = -0.00,0.00 # 2D vertical constraint
	rmin2,rmax2 = -0.50,0.50 # 2D horizontal constraint
	rmin3,rmax3 = -0.30,0.30 # 3D horizontal constraint
	smin,smax = -0.15,0.15 # min and max time shifts for 2D/3D warping
	#print "# of vertical slopes = ",1+2*rmax1/dr1 # base numbers on the vmin and vmax
	#print "# of horiz y  slopes = ",1+2*rmax2/dr2 # base number of slopes on misalignment of wells 
	#print "# of horiz x  slopes = ",rmax3/dr3
	goMultiWell2D(fp,q,nlr,vrmin,vrmax,dr,dr1,dr2,rmin1,rmin2,rmax1,rmax2,smin,smax)
	#goMultiWell3D(fp,q,nlr,vrmin,vrmax,dr,dr1,dr2,dr3,rmin1,rmin2,rmin3,rmax1,rmax2,rmax3,smin,smax)



################################################################################
################################################################################
## Warp methods

################################################################################
## Seismogram and 1D Warp methods



#------------------------------------------------------------------------------#

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

