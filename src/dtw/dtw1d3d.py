#############################################################################
# Dynamic time warping for 1D to 2D or 1D to 3D sequences

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
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from dtw import DynamicWarping1d3d
from wt import WTutils

#############################################################################

pngDir = "./dtwpics/1d2d/" 
#pngDir = "./dtwpics/1d3d/" 

# 1D to 2D modeling
nz = 501
nx = 301
dz = 8.0
dx = 8.0
nr = 10 # number of model reflectors
fpeak = 0.031 
fp = .7797/fpeak 
dt = 0.002
sv = Sampling(nr,2500/nr,2200) # velocity modeling
sd = Sampling(nr,0.7/nr,2.50) # density modeling
theta = 35 #degrees
sz = Sampling(nz,dz,0.0)
sx = Sampling(nx,dx,0.0)
hx = (nx-1)/2 + 1 #centered seismogram


# 1D to 3D modeling
#
#
#

def main(args):
	#2D
	st,wst,sm,ws = go1D2Dmodels()
	go1D2Dwarp(st,wst,sm,ws)
	#3D
	#go1D3D()

def go1D2Dwarp(st,wst,sm,ws):
	nt1 = st.count
	nt2 = wst.count
	ods1 = 10
	ods2 = 4 #even
	nl1 = (nt1-nt2)*ods1+1
	nl2 = 5 #odd
	nl2i = nl2*ods2+1
	sms = copy(nt1,nl2i,0,hx-(nl2i-1)/2,sm)
	hnl2 = (nl2i-1)/2+1
	
	dw = DynamicWarping1d3d(nl1,nl2i,ods1,ods2)
	# Set constraints
	k1min,k1max = -10,10 
	j1min,j1max = -5,5
	k2min,k2max = -1,1
	j2min,j2max = -1,1
	dw.setl1(k1min,k1max,j1min,j1max)
	dw.setl2(k2min,k2max,j2min,j2max)
	# Apply warping
	e = dw.computeErrors1d2d(ws,sms)
	d = dw.accumulate(e)
	u1,u2 = dw.findShifts(d)
	wh = dw.applyShifts(u1,ws)
	u1 = mul(u1,dt)
	u2 = div(sub(hnl2,u2),ods2)
	sh = Sampling(len(wh),dt,u1[0])

	#sf2 = SimpleFrame()
	#sf2.addImagePanels(d)
	plot1(u1,s=wst,haxis="u1",vaxis="time (s)")
	plot1(u2,s=wst,haxis="u2",vaxis="time (s)")
	plot(sx,st,sm,sy=wh,sst=sh,title='seismic',cmap=gray,perc=99,vaxis="time (s)")

def go1D2Dmodels():
	# Get 2D models
	v,d = make2Dmodels(nz,nx,sv,sd,theta,hgrad=0.10)
	toz,st = make2Dtoz(v,sz,dt)
	r = make2Drefl(v,d)
	s,st = make2Dseismic(st,sz,toz,r,fp)
	plot(sx,sz,v,title='velocity m/s')
	plot(sx,sz,d,title='density g/cc')
	#plot(sx,sz,r,title='reflectivity',cmap=gray)
	#plot(sx,sz,toz,title="t(z)")
	plot(sx,st,s,title='seismic',cmap=gray,perc=99,vaxis="time (s)")
	
	# Get 1D model
	wtheta = theta 
	wv,wd = make1Dmodel(nz,nx,sv,sd,wtheta,hx)
	ws,wst = make1Dseismic(wv,wd,sz,dt,fp)
	plot1(ws,s=wst,haxis="seismic theta= "+str(wtheta))
	plot(sx,st,s,sy=ws,sst=wst,title='seismic',cmap=gray,perc=99,vaxis="time (s)")
	return st,wst,s,ws

##################################################################
## Utils

def make2Dmodels(nz,nx,sv,sd,theta,noise=None,hgrad=None):
	nr = sv.count
	dv = sv.delta
	fv = sv.first
	dd = sd.delta
	fd = sd.first
	nl = nz/nr
	xl = nl/tan((FLT_PI/180.0)*theta)  
	if hgrad:	
		gr = hgrad/nx
	if (xl>=nx):
		xl = float(nx)
	 	theta = atan(nl/xl)*180/FLT_PI
		print 'xl>=nx'
		print 'max theta= ',theta
	dz = nl
	cx = (nx-int(xl))/2
	v = zerofloat(nz,nx)
	d = zerofloat(nz,nx)
	# Create array of z values
	zv = zeroint(nx)
	nxmc = nx-cx
	for i in range(cx):
		zv[i] = nl
	for i in range(nxmc,nx):
		zv[i] = 2*nl
	j=0
	for i in range(cx,nxmc):
		zv[i] = int((dz/xl)*j+nl)
		j+=1
	# Create velocity and density model
	for ix in range(nx):
		for izl in range(0,nr):
			zinp = zv[ix]+(izl  )*dz+1
			zinm = zv[ix]+(izl-1)*dz
			if (izl==0): zinm=0;
			if (zinp>nz): zinp = nz
			for iz in range(zinm,zinp):
				v[ix][iz] = fv+dv*izl-gr*ix*fv
				d[ix][iz] = fd+dd*izl
	return v,d

def make2Dtoz(v,sz,dt):
	wtu = WTutils()
	nx = len(v)
	nz = len(v[0])
	toz = zerofloat(nz,nx)
	for ix in range(nx):
		toz[ix] = wtu.tzVlog(v[ix],sz)
	nt = int(max(toz)/dt)
	st = Sampling(nt,dt,0.0)
	return toz,st

def make2Drefl(v,d):
	nx = len(v)
	nz = len(v[0])
	r = zerofloat(nz,nx)
	for ix in range(nx):
		imp0 = v[ix][0]*d[ix][0]
		for iz in range(1,nz):
			imp1 = v[ix][iz]*d[ix][iz]	
			r[ix][iz] = (imp1-imp0)/(imp1+imp0)
			imp0 = imp1
	return r

def make2Dseismic(st,sz,toz,r,fpeak):
	wtu = WTutils()
	nx = len(r)
	s = zerofloat(st.count,nx)
	ml = 0
	for ix in range(nx):
		s[ix] = wtu.syntheticSeismogram(fpeak,r[ix],toz[ix],st,sz,"ricker")
		if (len(s[ix])>ml): 
			ml = len(s[ix])
	nt = ml
	s2 = zerofloat(nt,nx)
	for ix in range(nx):
		for it in range(len(s[ix])):
			s2[ix][it] = s[ix][it]
	st = Sampling(nt,st.delta,st.first)
	#dz = sz.delta
	#nx = len(r)
	#nz = len(r[0])
	#t = zerofloat(nz)
	#s = zerofloat(nz,nx)
	#for ix in range(nx):
	#	s[ix] = addRickerWavelet(fpeak,r[ix])
  #	s[ix] = mul(1.0/max(abs(r[ix])),r[ix])
	return s2,st

def make1Dmodel(nz,nx,sv,sd,theta,hx,noise=None):
	nr = sv.count
	dv = sv.delta
	fv = sv.first
	dd = sd.delta
	fd = sd.first
	nl = nz/nr
	xl = nl/tan((FLT_PI/180.0)*theta)  
	if (xl>=nx):
		xl = float(nx)
	 	theta = atan(nl/xl)*180/FLT_PI
		print 'xl>=nx'
		print 'max theta= ',theta
	dz = nl
	cx = (nx-int(xl))/2
	v = zerofloat(nz)
	d = zerofloat(nz)
	# Create array of z values
	zv = zeroint(nx)
	nxmc = nx-cx
	for i in range(cx):
		zv[i] = nl
	for i in range(nxmc,nx):
		zv[i] = 2*nl
	j=0
	for i in range(cx,nxmc):
		zv[i] = int((dz/xl)*j+nl)
		j+=1
	# Create velocity and density model
	ix = hx
	for izl in range(0,nr):
		zinp = zv[ix]+(izl  )*dz+1
		zinm = zv[ix]+(izl-1)*dz
		if (izl==0): zinm=0;
		if (zinp>nz): zinp = nz
		for iz in range(zinm,zinp):
			v[iz] = fv+dv*izl
			d[iz] = fd+dd*izl
	return v,d

def make1Dseismic(v,d,sz,dt,fpeak):
	wtu = WTutils()
	nz = len(v)
	toz = wtu.tzVlog(v,sz)
	r = zerofloat(nz)
	s = zerofloat(nz)
	imp0 = v[0]*d[0]
	for iz in range(1,nz):
		imp1 = v[iz]*d[iz]	
		r[iz] = (imp1-imp0)/(imp1+imp0)
		imp0 = imp1
	st = Sampling(int(max(toz)/dt),dt,0.0)
	s = wtu.syntheticSeismogram(fpeak,r,toz,st,sz,"ricker")
	st = Sampling(len(s),dt,0.0)
	#s = addRickerWavelet(fpeak,r)
 	#s = mul(1.0/max(abs(s)),s)
	return s,st

def sig1(n1,rmsNoise, fpeak):
  seed   = 23232
  seedN1 = 31415 
  f = makeRandomEvents(n1,seed) 
  f = addRickerWavelet(fpeak,f)
  f = mul(1.0/max(abs(f)),f)
  f = addNoise(rmsNoise,fpeak,f,seedN1)
  return f

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

