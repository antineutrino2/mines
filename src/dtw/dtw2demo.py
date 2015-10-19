#############################################################################
# Dynamic warping 2D images

import sys
from java.awt import *
from java.awt.image import *
from java.io import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.dsp import DynamicWarping as DynamicWarpingO
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.interp import *
from edu.mines.jtk.util.ArrayMath import *

from dtw import DynamicWarpingR
from dtw import DynamicWarpingWTM
from wt import WTutils

#############################################################################

pngDir = "./dtwpics/demo2d/" 
wtu = WTutils()

# 2D modeling
nt = 501
nx = 301
dt = 1.0
dx = 1.0
nr = 2 # number of evenly spaced model layers
fpeak = 3.0/35
fv = 250.0
dv = 10.0/nr
fd = 2.5
dd = 0.5/nr
sv = Sampling(nr,fd,fv) # velocity modeling
sd = Sampling(nr,dd,fd) # density modeling
theta = 0.0  # layer tilt in degrees
st = Sampling(nt,dt,0.0)
sx = Sampling(nx,dx,0.0)
hx = (nx-1)/2 + 1 #center

setSinStrainMax1 = 0.1
setSinStrainMax2 = 0.1
seed1 = 32421
noise1= 0.01
hgrad=None

# Warp the wells to the seismic
smin,smax = -50,50
r1min,r1max = -0.5,0.5
r2min,r2max = -0.5,0.5
d1r,d2r 		= 0.10,0.10
dwr = DynamicWarpingR(smin,smax,st,sx)
dwr.setStrainLimits(r1min,r1max,r2min,r2max)
dwr.setSmoothness(d1r,d2r)

dwo = DynamicWarpingO(int(smin/dt),int(smax/dt))
dwo.setStrainMax(r1max,r2max)
dwo.setErrorSmoothing(2)
dwo.setShiftSmoothing(1/d1r,1/d2r)

def main(args):
	v,d = makeFlat2Dmodel()
	f = getImage(v,d)
	n1 = nt; n2 = nt
	w,sw,xc = extract2Wells(f,n1,n2)
	w[0] = warpSin(w[0],setSinStrainMax1)
	w[1] = warpSin(w[1],setSinStrainMax2)

	# Do a first pass tie attempt in 1D
	#lot2D(sx,st,f,sw,w,xc,title="1st pass wells and seismic")
	#w = warp1Dpass(f,w,xc)
	#lot2D(sx,st,f,sw,w,xc,title="1st pass shifted wells and seismic")
	
	# Interpolate the Wells
	bgs = 0.0
	g = interpolateWells(w,sw,xc,bgs)

	#wxc = xc[0]
	#warp1Dtest(f[wxc],g[wxc])

	#warpR2Dtest(f,g,w,sw,xc)
	
	#warpO2Dtest(f,g,w,sw,xc)
	
	warpC2Dtest(f,g,w,sw,xc)


def warp1Dpass(f,w,xc):
	nnw = len(w)
	wn = w
	for i in range(nnw):
		g = w[i]
		u = dwo.findShifts(f[xc[i]],g)
		wn[i] = dwo.applyShifts(u,g)
	return wn

def	warpO2Dtest(f,g,w,sw,xc):
	e = dwo.computeErrors(f,g)
	u = dwo.findShifts(f,g)
	h = dwo.applyShifts(u,g)
	wn1 = h[xc[0]]
	wn2 = h[xc[1]]
	wn = [wn1,wn2]
	swn = [Sampling(len(wn1),dt,0.0),#u[xc[0]][0]),
				 Sampling(len(wn2),dt,0.0)]#u[xc[1]][0])]
	plot2D(sx,st,f,title='seismic',cmap=gray)
	plot2D(sx,st,f,swn,wn,xc,title="shifted wells and seismic")
	plot2D(sx,st,h,title="shifted wells")
	plot2D(sx,st,u,title="shifts",cmap=jet)
	plot2D(sx,st,f,sw,w,xc,title="wells and seismic")
	#plot(sx,st,sm,sy=wh,sst=sh,title='seismic',cmap=gray,perc=99,vaxis="time (s)")
	#plot1(w[0],sw[0])
	#plot1(w[1],sw[1])

def	warpC2Dtest(f,g,w,sw,xc):
	dwc = DynamicWarpingWTM(f,g)
	u = dwc.findShifts(r1min,r1max,d1r,r2min,r2max,d2r)
	h = dwc.applyShifts(u)
	wn1 = h[xc[0]]
	wn2 = h[xc[1]]
	wn = [wn1,wn2]
	swn = [Sampling(len(wn1),dt,0.0),#u[xc[0]][0]),
				 Sampling(len(wn2),dt,0.0)]#u[xc[1]][0])]
	plot2D(sx,st,f,title='seismic',cmap=gray)
	plot2D(sx,st,f,swn,wn,xc,title="shifted wells and seismic")
	plot2D(sx,st,h,title="shifted wells")
	plot2D(sx,st,u,title="shifts",cmap=jet)
	plot2D(sx,st,f,sw,w,xc,title="wells and seismic")
	#plot(sx,st,sm,sy=wh,sst=sh,title='seismic',cmap=gray,perc=99,vaxis="time (s)")
	#plot1(w[0],sw[0])
	#plot1(w[1],sw[1])

def	warpR2Dtest(f,g,w,sw,xc):
	u = dwr.findShifts(st,f,st,g)
	h = zerofloat(nt,nx)
	for ix in range(nx):
		h[ix] = dwr.applyShifts(st,g[ix],u[ix])
	wn1 = h[xc[0]]
	wn2 = h[xc[1]]
	wn = [wn1,wn2]
	swn = [Sampling(len(wn1),dt,0.0),#u[xc[0]][0]),
				 Sampling(len(wn2),dt,0.0)]#u[xc[1]][0])]
	plot2D(sx,st,f,title='seismic',cmap=gray)
	plot2D(sx,st,f,swn,wn,xc,title="shifted wells and seismic")
	plot2D(sx,st,h,title="shifted wells")
	plot2D(sx,st,u,title="shifts",cmap=jet)
	plot2D(sx,st,f,sw,w,xc,title="wells and seismic")
	#plot(sx,st,sm,sy=wh,sst=sh,title='seismic',cmap=gray,perc=99,vaxis="time (s)")
	#plot1(w[0],sw[0])
	#plot1(w[1],sw[1])


def warp1Dtest(f,g):
	e = dwr.computeErrors(st,f,st,g)
	u = dwr.findShifts(st,f,st,g)
	h = dwr.applyShifts(st,g,u)
	plot1(u,haxis="u")
	plot1(h,haxis="h")
	plot1(f,haxis="f")
	plot1(g,haxis="g")
	SimplePlot().addPixels(e)


##################################################################
## Utils

def interpolateWells(w,sw,xc,bgs=0.0):
	pnull = -999.0
	g = zerofloat(nt,nx)
	gv = fillfloat(pnull,nt,nx)
	for i in range(len(w)):
		gv[xc[i]] = w[i]
	#plot2D(sx,st,g)
	xa1,xa2,wa = getXcoords(sw,xc,w)
	bg = BlendedGridder2(wa,xa1,xa2)
	bg.setSmoothness(bgs)
	gt = bg.gridNearest(pnull,gv)
	bg.gridBlended(gt,gv,g)
	plot2D(sx,st,g,title="interpolated wells",cmap=jet)
	return g

def getXcoords(sw,xc,w):
	sw1,sw2 = sw[0],sw[1]
	xc1,xc2 = xc[0],xc[1]
	w1,w2 	=	w[0],w[1]
	n1,n2 	=	sw1.count,sw2.count
	n = n1+n2
	xa1,xa2,wa = zerofloat(n),zerofloat(n),zerofloat(n)
	for i in range(n1):
		xa1[i] = sw1.delta*i+sw1.first
		xa2[i] = xc1
		wa[i]  = w1[i]
	for i in range(n2):
		xa1[n1+i] = sw2.delta*i+sw2.first
		xa2[n1+i] = xc2
		wa[n1+i]  = w2[i]
	return xa1,xa2,wa
			
def warpSin(x,strain):
	n = len(x)
	start=0
	amp = strain*n/(2*PI)
	w = float(2*PI/n)
	y = zerofloat(n)
	t = zerofloat(n)
	for i in range(0,n):
		t[i]=start+i-10#amp*sin(i*w)   
	si = SincInterp()
	si.setExtrapolation(SincInterp.Extrapolation.CONSTANT)
	si.interpolate(n,1.0,0.0,x,n,t,y)    
	if noise1:
		y = addNoise(noise1,fpeak,y,seed1)
	return y

def extract2Wells(f,n1,n2):
	wi = nx/3
	x1 = wi
	x2 = wi*2
	w1 = f[x1]
	w2 = f[x2]
	sw1 = st
	sw2 = st
	w = [w1,w2]
	sw = [sw1,sw2]
	xc = [x1,x2]
	return w,sw,xc

def getImage(v,d):
	r = make2Drefl(v,d)
	s = make2Dseismic(r)
	plot2D(sx,st,v,cmap=jet,title='velocity m/s')
	plot2D(sx,st,d,cmap=jet,title='density g/cc')
	#plot(sx,st,r,title='reflectivity',cmap=gray)
	return s

def makeFlat2Dmodel():
	v = zerofloat(nt,nx)
	d = zerofloat(nt,nx)
	# Create velocity and density model
	for ix in range(nx):
		dr=nt/nr
		ir=0
		for it in range(nt):
			v[ix][it] = fv+dv*ir
			d[ix][it] = fd+dd*ir
			if it==dr:
				ir+=1
				dr+=nr
	return v,d

def make2Drefl(v,d):
	r = zerofloat(nt,nx)
	for ix in range(nx):
		imp0 = v[ix][0]*d[ix][0]
		for it in range(1,nt):
			imp1 = v[ix][it]*d[ix][it]	
			r[ix][it] = (imp1-imp0)/(imp1+imp0)
			imp0 = imp1
	return r

def make2Dseismic(r):
	s = zerofloat(nt,nx)
	for ix in range(nx):
		f = addRickerWavelet(fpeak,r[ix])
		f = mul(1.0/max(abs(f)),f)
		s[ix] = addNoise(noise1,fpeak,f,seed1+ix)
	return s


def make2Dmodels():
	theta = 0.0000001
	nl = nt/nr
	xl = nl/tan((FLT_PI/180.0)*theta)  
	gr = 0.0
	if hgrad:	
		gr = hgrad/nx
	if (xl>=nx):
		xl = float(nx)
	 	theta = atan(nl/xl)*180/FLT_PI
		print 'xl>=nx'
		print 'max theta= ',theta
	dz = nl
	cx = (nx-int(xl))/2
	v = zerofloat(nt,nx)
	d = zerofloat(nt,nx)
	# Create array of z values
	zv = zeroint(nx)
	nxmc = nx-cx
	for i in range(cx):
		zv[i] = nl
	for i in range(nxmc,nx):
		zv[i] = 2*nl
	j=0
	for i in range(cx,nxmc):
		zv[i] = int((dt/xl)*j+nl)
		j+=1
	# Create velocity and density model
	for ix in range(nx):
		for izl in range(0,nr):
			zinp = zv[ix]+(izl  )*dt+1
			zinm = zv[ix]+(izl-1)*dt
			if (izl==0): zinm=0;
			if (zinp>nt): zinp = nt
			for iz in range(zinm,zinp):
				v[ix][iz] = fv+dv*izl-gr*ix*fv
				d[ix][iz] = fd+dd*izl
	return v,d

def normalize(x):
  xmin = min(x)
  xmax = max(x)
  return mul(sub(x,xmin),1.0/(xmax-xmin))

def makeRandomEvents(n,seed=0):
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  return pow(mul(2.0,sub(randfloat(r,n),0.5)),15.0)

def addRickerWavelet(fpeak,f):
  n = len(f)
  ih = int(3.0/fpeak)
  nh = 1+2*ih
  h = zerofloat(nh)
  for jh in range(nh):
    h[jh] = ricker(fpeak,jh-ih)
  g = zerofloat(n)
  Conv.conv(nh,-ih,h,n,0,f,n,0,g)
  return g

def ricker(fpeak,time):
  x = PI*fpeak*time
  return (1.0-2.0*x*x)*exp(-x*x)

def addNoise(nrms,fpeak,f,seed=0):
  n = len(f)
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  nrms *= max(abs(f))
  g = mul(2.0,sub(randfloat(r,n),0.5))
  g = addRickerWavelet(fpeak,g)
  #rgf = RecursiveGaussianFilter(3.0)
  #rgf.apply1(g,g)
  frms = sqrt(sum(mul(f,f))/n)
  grms = sqrt(sum(mul(g,g))/n)
  g = mul(g,nrms*frms/grms)
  return add(f,g)

def floats(x):
  n1 = len(x)
  y = zerofloat(n1)
  for i1 in range(n1):
    y[i1] = float(x[i1])
  return y

##################################################################

##################################################################
## plots

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(sx,sz,x,cmap=jet,perc=100,title=None,png=None,vaxis=None,sy=None,sst=None,ssx=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.addColorBar()
  sp.setSize(1200,500)
  pv = sp.addPixels(sz,sx,x)
  pv.setColorModel(cmap)
  sp.setHLabel("Distance (m)")
  sp.setVLabel("Depth (m)")
  if sy and sst and not ssx:
		asy = zerofloat(len(sy),1)
		copy(sy,asy[0])
		ssx = Sampling(1,sx.delta,sx.delta*(sx.count/2))
		sp.addPixels(sst,ssx,asy)
  if sy and sst and ssx:
		asy = zerofloat(len(sy),1)
		copy(sy,asy[0])
		sp.addPixels(sst,ssx,asy)
  if vaxis:
  	sp.setVLabel(vaxis)
  if title: 
		sp.setTitle(title)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if perc<100:
    pv.setPercentiles(100-perc,perc)
  if png:
    sp.paintToPng(360,3.33,pngDir+png+'.png')

def plot2D(sx,st,f,sw=None,w=None,xc=None,cmap=gray,title=None):
	pp = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
	pp.addColorBar()
	pp.setColorBarWidthMinimum(100)
	img = pp.addPixels(st,sx,f)
	img.setColorModel(cmap)
	pp.setHLabel("Distance (m)")
	pp.setVLabel("Time (s)")
	if title:
		pp.setTitle(title)
	if sw and w and xc:
		nw = len(w)
		for iw in range(nw):
			swi = sw[iw]
			sxw = Sampling(1,sx.delta,xc[iw]*sx.delta)
			nw = swi.count
			w1 = zerofloat(nw,1)
			copy(nw,w[iw],w1[0])
			log = pp.addPixels(swi,sxw,w1)
			log.setColorModel(cmap)
			log.setInterpolation(PixelsView.Interpolation.NEAREST)
			log.setClips(min(f),max(f))
	frame = PlotFrame(pp)
	frame.setSize(1200,500)
	frame.setVisible(True)
	frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);

def plot1(x,s=None,png=None,haxis=None,vaxis=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if s:
		sp.addPoints(s,x)
  else:
		sp.addPoints(x)
  sp.setSize(500,1000)
  if haxis: sp.setHLabel(haxis)
  sp.setVLabel("Depth (m)")
  if vaxis: sp.setVLabel(vaxis)
  if png:
    sp.paintToPng(360,3.33,pngDir+png+'.png')

#############################################################################
# Do everything on Swing thread.

import sys,time
class RunMain(Runnable):
  def run(self):
    start = time.time()
    main(sys.argv)
    s = time.time()-start
    h = int(s/3600); s -= h*3600
    m = int(s/60); s -= m*60
    print '%02d:%02d:%02ds'%(h,m,s)
SwingUtilities.invokeLater(RunMain())

