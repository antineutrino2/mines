"""
Jython utilities for Teapot Dome.
Author: Andrew Munoz, Colorado School of Mines
Version: 2012.09.11
"""
from dtw import *
from imports import *

#############################################################################
#############################################################################
# Phase tools

"""
Rotates the phase of the input synthetic seismogram by a constant specified 
in degrees.
"""
def applyPhaseRot(ph,xr):
  n = len(xr)
  xim = zerofloat(n)
  y = zerofloat(n)
  hb = HilbertTransformFilter()
  hb.apply(n,xr,xim)
  phr = (FLT_PI/180)*ph 
  y = add(mul(cos(phr),xr),mul(sin(phr),xim))
  return y

"""
Finds constant phase rotation of synthetic more efficiently.
"""
def rotatePhase(tr,sy,fr):
  n1 = len(sy)
  n2 = len(tr)
  nl = 1+fr*(n2-n1)
  lmin = 0
  lmax = nl-fr+1
  stretchMax = 1.0
  dw = DynamicWarpingF(lmin,lmax,fr)
  dw.setStrainMax(stretchMax)
  nh = 360
  mset = [30,5,1]
  mphr = zerofloat(nh); 
  pmrms=1; sa=0; en=nh; ihm=0;
  for mu in mset:
    for ih in range(sa,en,mu):
      print ih
      syr = applyPhaseRot(ih,sy)
      e = dw.computeErrorsFrac2(syr,tr)
      d = dw.accumulateForward(e)
      u = dw.findShiftsFrac(syr,tr,e,d)
      mphr[ih] = u[0]
      if ih==0: pmrms=mphr[ih]
      if (mphr[ih]<=pmrms): 
        ihm=ih
        pmrms = mphr[ih]
    sa = ihm-mu
    en = sa+2*mu  
    if sa<0: sa=0
    if en>nh: en=nh
  print 'Optimum Phase is',ihm
  return applyPhaseRot(ihm,sy)

"""
Finds constant phase rotation of synthetic by rotating through and checking
360 degrees using dtw.
"""
def rotateAllPhase(tr,sy,fr):
  n1 = len(sy)
  n2 = len(tr)
  nl = 1+fr*(n2-n1)
  lmin = 0
  lmax = nl-fr+1
  stretchMax = 1.0
  dw = DynamicWarpingF(lmin,lmax,fr)
  dw.setStrainMax(stretchMax)
  nh = 361
  mphr = zerofloat(nh); 
  pmrms = 1; ihm=0;
  for ih in range(nh):
    print ih
    syr = applyPhaseRot(ih,sy)
    e = dw.computeErrorsFrac2(syr,tr)
    d = dw.accumulateForward(e)
    u = dw.findShiftsFrac(syr,tr,e,d)
    mphr[ih] = u[0]    
    if ih==0: pmrms=mphr[ih]
    if (mphr[ih]<=pmrms): 
      ihm=ih
      pmrms = mphr[ih]
  print 'Optimum Phase is',ihm
  #plotPhaseRMS(mphr,nh,1,slides='sphaseplot',paper='phaseplot')
  return applyPhaseRot(ihm,sy)


#############################################################################
#############################################################################
# Normalizations 

"""
Local RMS normalization.
See java method documentation. 
"""
def normalizeRMSlocal(x,sig):
  return WTutils().localRMSnorm(x,sig)

"""
Global RMS normalization
"""
def normalizeRMS(x):
  n = len(x)
  rms = sqrt(sum(mul(x,x))/n)
  div(x,rms)
  return x

"""
Global normalization
"""
def globNorm(a):
  amax = max(a)
  div(a,amax)
  return a

#############################################################################
#############################################################################
# General utility methods

"""
Plots the amplitude spectrum of the input sequence using
the sampling, the values, and the title of the plot (s,x,l)
"""
def plotSpectrum(s,x,l):
  Spectrum().apply(s,x,l)

"""
Used to over sample a signal
@param x1 the signal
@param xs1 the uniform Sampling
@param xs2 the new uniform Sampling
"""
def overSample(x1,x2,xs1,xs2):
  n1 = len(x1)
  n2 = xs2.getCount() 
  si = SincInterpolator()
  si.setUniform(n1,xs1.getDelta(),xs1.getFirst(),x1)
  si.interpolate(n2,xs2.getDelta(),xs2.getFirst(),x2)
  return x2

"""
Pad 2nd array if
len(f1)>len(f2) 
"""
def addPad(f1,f2,z2,y2,x2):
  n1,n2 = len(f1),len(f2)
  ft,zt,yt,xt = zerofloat(n1),zerofloat(n1),zerofloat(n1),zerofloat(n1)
  copy(n2,f2,ft)
  copy(n2,x2,xt)
  copy(n2,y2,yt)
  copy(n2,z2,zt)
  return ft,xt,yt,zt

"""
Converts x,y,z coords from kilometers to meters
"""
def km2m(x,y,z):
  x = mul(x[i],1000.0)
  y = mul(y[i],1000.0)
  z = mul(z[i],1000.0)
  return x,y,z

"""
Converts 1D float array to 1D double array
"""
def double(x):
  xd = zerodouble(len(x))
  xd = add(x,0.0)
  return xd

"""
Converts 1D double array to 1D float array
"""
def afloat(x):
  xf = zerofloat(len(x))
  xf = add(x,0.)
  return xf

"""
Clips a signal based on the average to remove the outliers
"""
def Clip(a):
  return WTutils().clip(a)

def normalize(e):
  emin = min(e)
  emax = max(e)
  return mul(sub(e,emin),1.0/(emax-emin))

def etran(e):
  e = normalize(e)
  return transpose(pow(e,0.25))


#############################################################################
#############################################################################
# Test signals
# @author Luming Liang

def sig1L(n1):
  seed   = 23232
  seedN1 = 31415 
  fpeak  = 0.125
  rsmNoise = 0.1
  wtu = WTutils()
  X = wtu.makeEvents(n1, seed) 
  X = wtu.addRickerWavelet(fpeak, X)
  X = mul(1.0/max(abs(X)),X)
  X = wtu.addNoise(rsmNoise, seedN1, X)
  return X

def sig1S(n1,n2):
  _si = SincInterpolator()
  _si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT)
  seedN2 = 10000
  rsmNoise = 0.5
  amp = 30
  Y = zerofloat(n2)
  t = zerofloat(n2)
  X = sig1L(n1);
  w = float((2*PI/n2))
  start = 300
  for i in range (0,n2):
    t[i]=start+i-amp*sin(i*w) 
  _si.setUniform(n1,1.0,0.0,X)
  _si.interpolate(n2,t,Y)   
  Y = WTutils().addNoise(rsmNoise, seedN2, Y)
  return X,Y

#############################################################################
#############################################################################
# Plotting

def plotMatrix(c,sf,slag,u=None,png=None):
  n1,nlag = len(c[0]),len(c)
  #slag = Sampling(nlag,1.0,-(nlag-1)/2)
  panel = PlotPanel(PlotPanel.Orientation.X1RIGHT_X2UP)
  cv = panel.addPixels(sf,slag,c)
  cv.setInterpolation(PixelsView.Interpolation.NEAREST)
  cv.setColorModel(ColorMap.JET)
  if u:
    uv = panel.addPoints(sf,u)
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
  panel.setVLimits(slag.first,slag.last)
  panel.setHLimits(sf.first,sf.last)
  panel.setHLabel("sample")
  panel.setVLabel("lag")
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(18)
  frame.setSize(1000,800)
  frame.setVisible(True)
  if png:
    frame.paintToPng(720,3.33,_pngDir+png+'.png')
  
def plotSequences(f,g,s1=None,s2=None):
  n1 = len(f)
  n2 = len(g)
  if s1 is None:
    s1 = Sampling(n1,1.0,0.0)
  if s2 is None:
    s2 = Sampling(n2,1.0,0.0)
  panel = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  fv = panel.addPoints(s1,f)
  gv = panel.addPoints(s2,g)
  gv.setLineColor(Color.RED)
  panel.setHLabel("f & g")
  panel.setVLabel("sample")
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(18)
  frame.setSize(300,800)
  frame.setVisible(True)

"""
Plots a panel of well logs
"""
def plotLogPanel(sz,v,d,rf,paper=None,slides=None):
  p1 = PlotPanel(1,3,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  p1.setVLabel("Depth (km)")
  vp = p1.addPoints(0,0,sz,v)
  dp = p1.addPoints(0,1,sz,d)
  rf = p1.addPoints(0,2,sz,rf)
  p1.setHLabel(0,"Velocity (km/s)")
  p1.setHLabel(1,"Density (g/cc)")
  p1.setHLabel(2,"Reflectivity")
  p1.setHLimits(0,min(v),max(v))
  p1.setHLimits(1,1.51,max(d))
  p1.setHInterval(1,0.5)
  p1.setHLimits(2,-0.13,0.13)
  p1.setHInterval(2,0.1)
  #dp.setLineColor(BLUE)
  #gp.setLineColor(GREEN)
  frame = PlotFrame(p1)
  frame.setSize(800,1500)
  frame.setFontSize(28)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  if paper:
    frame.setSize(548,738)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    frame.setSize(690,561)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)

"""
Plots comparison of reflectvity to its resulting synthetic seismogram
"""
def plotCurveW(sz,rf,ss,sy,paper=None):
  p1 = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  p2 = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  p1.addPoints(0,0,sz,rf)
  p2.addPoints(0,0,ss,sy)
  p1.setVLabel(0,"Depth (km)")
  p2.setVLabel(0,"t (s)")
  p1.setHLimits(0,-0.13,0.13)
  p2.setHLimits(0,-0.5,0.5)
  frame = PlotFrame(p1,p2,PlotFrame.Split.HORIZONTAL)
  frame.setSize(300,1500)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True)
  if paper:
    frame.setSize(584,757)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')

"""
Plots comparison of synthetic seismogram before and after normalization
"""
def plotCurveN(st,tr,trn,ss,sy,syn,lim=None,paper=None):
  p1 = PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  p2 = PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  p1.addPoints(0,0,st,tr)
  p1.addPoints(0,1,st,trn)
  p2.addPoints(0,0,ss,sy)
  p2.addPoints(0,1,ss,syn)
  p1.setVLabel(0,"Time (s)")
  p2.setVLabel(0,"t (s)")
  p2.setHLimits(0,-1.1,1.1)
  frame = PlotFrame(p1,p2,PlotFrame.Split.HORIZONTAL)
  frame.setSize(300,1500)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True)
  if paper:
    frame.setSize(584,757)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')

"""
Plots time-depth curves before and after warping
"""
def plotTDCurves(sz,tauz,tz,paper=None,slides=None):
  p1 = PlotPanel()
  l1 = p1.addPoints(sz,tz)
  l2 = p1.addPoints(sz,tauz)
  l1.setLineColor(RED)
  p1.setVLabel('Time (s)')
  p1.setHLabel('Depth (km)')
  frame = PlotFrame(p1)
  frame.setSize(700,500)
  frame.setFontSizeForPrint(8.0,240.0)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  if paper:
    frame.setSize(513,379) #700,500)
    frame.setFontSizeForPrint(8.0,192.0)
    frame.paintToPng(720,2.165,_pngDir+paper+'.png')
  if slides:
    frame.setSize(513,379)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)

"""
Plots teaser figure for CWP report 2012 of synthetic and trace before and after alignment
"""
def plotTeaser(tr,st,y,ys,sst,ss,tm=None,lims=None,paper=None,slides=None,paper2=None,paper3=None):
  if tm:
    p1 = PlotPanel(1,3,PlotPanel.Orientation.X1DOWN_X2RIGHT)
    p1.setHLabel(2,"tau")
    stau = Sampling(ss.getCount(),ss.getDelta(),tm[0])
    p1.addPoints(0,2,stau,tm)
  else: 
    p1 = PlotPanel(2,1)#,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  l1 = p1.addPoints(0,0,st,tr)
  l1.setLineWidth(2)
  p1.setHLabel("Time (s)")
  #ss2 = Sampling(ss.getCount(),ss.getDelta(),sst.getFirst())
  l2 = p1.addPoints(0,0,ss,y)
  l2.setStyle('r-')
  l2.setLineWidth(2)
  l3 = p1.addPoints(1,0,st,tr)
  l3.setLineWidth(2)
  l4 = p1.addPoints(1,0,sst,ys)
  l4.setStyle('r-')
  l4.setLineWidth(2)
  if lims:
    p1.setVLimits(0,-lims,lims)
    p1.setVLimits(1,-lims,lims)
    if tm: p1.setVLimits(2,-lims,lims)
  p1.setHLimits(0,sst.getLast()+0.5)
  frame = PlotFrame(p1)
  frame.setSize(900,1500)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True)
  if paper:
    frame.setSize(1360,354)
    frame.setFontSizeForPrint(8.0,504.0)
    frame.paintToPng(720,7.00,_pngDir+paper+'.png')
  if paper2:
    frame.setSize(904,379)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper2+'.png')
  if paper3:
    frame.setSize(904,379)
    frame.setFontSizeForPrint(8.0,312.0)
    frame.paintToPng(720,4.33,_pngDir+paper3+'.png')
  if slides:
    frame.setSize(500,900)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)

"""
Plots a curve vertically and optionally another as a red in the same panel
"""
def plotCurve(c1,s1,c2=None,s2=None,title=None,lim=None,png=None,d=None):
  p1 = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  l1 = p1.addPoints(s1,c1)
  l1.setLineColor(BLACK)
  if c2 and s2:
    l2 = p1.addPoints(s2,c2)
    l2.setLineColor(RED)
  p1.setVLabel("samples")
  if d=="z":
    p1.setVLabel("depth (km)")
  if d=="t":
    p1.setVLabel("time (s)")
  if lim:
    p1.setHLimits(-lim,lim)
  if title: 
    p1.setTitle(title)
  frame = PlotFrame(p1)
  frame.setSize(300,1500)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True)
  if png:
    frame.paintToPng(720,3.33,_pngDir+png+'.png') 

"""
Plots a 2-D seismic slice, a well log through that slice, and seismic horizons on the slice
"""
def plot2DLog(img,x,yc,st,sy,zot,tops=None,tz=None,ghz=None,ghy=None,nolog=None,paper=None,slides=None):
  nx = len(x)
  ny = len(img)
  x2 = zerofloat(nx,3)
  copy(x,x2[0])         
  copy(x,x2[1])
  copy(x,x2[2])         
  sy2 = Sampling(3,sy.getDelta(),(yc-1)*sy.getDelta())
  sp = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT) 
  im1 = sp.addPixels(st,sy,img)
  im1.setColorModel(ColorMap.GRAY)#RED_WHITE_BLUE)#GRAY)
  if ghz:
    # Add Horizons
    pv1 = sp.addPoints(ghz[0],ghy[0])
    pv1.setStyle("g-")
    pv1.setLineWidth(2)
    pv2 = sp.addPoints(ghz[1],ghy[1])
    pv2.setStyle("c-")
    pv2.setLineWidth(2)
    pv3 = sp.addPoints(ghz[2],ghy[2])
    pv3.setStyle("y-")
    pv3.setLineWidth(2)
    pv4 = sp.addPoints(ghz[3],ghy[3])
    pv4.setStyle("m-")
    pv4.setLineWidth(2)
  if tops:
    for i in range(len(tops)):
      if (tops[i]==1):
        it = int((tz[i]-tz[0])/st.getDelta())
	x2[0][it] = 100; x2[1][it] = 100; x2[2][it] = 100;
	x2[0][it-1] = 100; x2[1][it-1] = 100; x2[2][it-1] = 100;
	x2[0][it+1] = 100; x2[1][it+1] = 100; x2[2][it+1] = 100;
  log = sp.addPixels(zot,sy2,x2)
  log.setColorModel(ColorMap.GRAY)#RED_WHITE_BLUE)#GRAY)
  log.setInterpolation(PixelsView.Interpolation.NEAREST)
  log.setClips(min(img),max(img))
  if nolog: sp.remove(log)
  if ((zot.getFirst()-0.2)<0): zotgf = 0
  else: zotgf = zot.getFirst()-0.2
  sp.setLimits(yc*sy.getDelta()-1.2,zotgf,yc*sy.getDelta()+1.2,zot.getLast()+0.2)
  #sp.setLimits(yc*sy.getDelta()-1.2,0.45,yc*sy.getDelta()+1.2,1.4)
  sp.setHLabel("Distance (km)")
  sp.setVLabel("Time (s)")
  frame = PlotFrame(sp)
  frame.setSize(1500,1000)
  frame.setFontSize(24)
  if paper:
    frame.setSize(1200,800)
    frame.setFontSizeForPrint(8.0,504.0)
    frame.paintToPng(720,7.00,_pngDir+paper+'.png')
  if slides:
    frame.setSize(1000,700)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);

"""
Plots a 2-D seismic slice
"""
def plotSlice(sy,st,img,png=None,slides=None):
  p1 = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  p1.addPixels(st,sy,img)
  p1.setHLabel("Distance (km)")
  p1.setVLabel("Time (s)")
  frame = PlotFrame(p1)
  frame.setSize(1000,1000)
  frame.setFontSize(24)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True)
  if png:
    frame.paintToPng(1000,5.00,_pngDir+png+'.png')  
  if slides:
    frame.setSize(900,600)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)

"""
Plots phase rotation errors
"""
def plotPhaseRMS(mrms,nph,phf,png=None,paper=None,slides=None):
  sp = SimplePlot()
  sph = Sampling(nph/phf,phf,0.0)
  sp.addPoints(mrms)
  sp.setHLabel("Phase (degrees)")
  sp.setVLabel("Minimum Dn")
  sp.setTitle("Minumum RMS vs. Phase")
  sp.setSize(1500,1000)
  if png:
    sp.paintToPng(1200,6.00,_pngDir+png+'.png')
  if paper:
    sp.removeTitle()
    sp.setSize(1000,800)
    sp.setFontSizeForPrint(8.0,240.0)
    sp.paintToPng(720,3.33,_pngDir+paper+'.png')  
  if slides:
    sp.removeTitle()
    sp.setSize(950,800)
    sp.setFontSizeForSlide(1.0,1.0)
    sp.paintToPng(720,7.00,_pngDir+slides+'.png')

"""
Plots just a gaussian
"""
def gaussianPlot(sig):
  nf = 801
  f1 = zerofloat(nf)
  f2 = zerofloat(nf)
  f1[int(nf/2)] = 100000.0
  RecursiveGaussianFilter(sig).apply2(f1,f2)
  f2 = mul(-1,f2)
  si = SimplePlot()
  si.addPoints(f2)
  si.setSize(700,500)
  si.paintToPng(900,3.33,_pngDir+"RICKER"+".png")

"""
Histogram plot 
@author Farhad Bazarghani
"""
def goHistogram(data,xlabel,png=None):
  hg = Histogram(data)
  h = hg.getCounts()
  hf = zerofloat(len(h))
  for i in range (1,len(h)):
    hf[i]=float(h[i])
  sam = hg.getBinSampling()
  sp = SimplePlot()
  sp1 = sp.addSequence(sam,hf)
  sp.setTitle("histogram"+" "+xlabel)
  sp.setVLabel("# of points")
  sp.setHLabel(xlabel)
  if png:
    sp.paintToPng(360,3.33,_pngDir+png+'.png') 


#############################################################################
# Run the function main on the Swing thread
import sys
class _RunMain(Runnable):
  def __init__(self,main):
    self.main = main
  def run(self):
    self.main(sys.argv)
def run(main):
  SwingUtilities.invokeLater(_RunMain(main)) 
