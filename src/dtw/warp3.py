#############################################################################
# Dynamic time warping for 3D images

from imports import *
from dtw import Displacement3
from util import FakeData

#############################################################################

#pngDir = "./png"
pngDir = None

#seed = abs(Random().nextInt()/1000000)
seed = 588
print "seed =",seed

nrms = 0.0 # rms noise/signal ratio
r1min,r1max = -1.0,1.0 
r2min,r2max = -1.0,1.0 
r3min,r3max = -1.0,1.0 
dr1,dr2,dr3 = 1.00,1.00,1.00
n1,n2,n3 = 101,51,101 # numbers of samples
#d1,d2,d3 = 0.004,0.025,0.025 # sampling intervals
d1,d2,d3 = 1.000,1.000,1.000 # sampling intervals
#d1,d2,d3 = 0.100,1.000,1.000 # sampling intervals
f1,f2,f3 = 0.000,0.000,0.000 # first samples
s1,s2,s3 = Sampling(n1,d1,f1),Sampling(n2,d2,f2),Sampling(n3,d3,f3)
label1,label2,label3 = "Time (s)","Inline (km)","Crossline (km)"
dataDir = "/Users/amunoz/Home/data/seis/fake/"
#dataDir = "/data/amunoz/fake/"
sic = 50
smin,smax = -12*d1,12*d1 # shifts for fake images are in samples

def main(args):
  #goFakeImages()
  #goFakeShifts3()
  goFakeShifts2()
  goFakeShifts2O()
  #goFakeShiftsBench()

def goFakeImages():
  f,g,s = makeFakeImages(smax/d1,nrms)
  print "s: min =",min(s),"max =",max(s)
  writeImage(dataDir+"fakef.dat",f)
  writeImage(dataDir+"fakeg.dat",g)
  writeImage(dataDir+"fakes.dat",s)
  show(f)
  show(g)
  show(s)

def goFakeShifts3():
  f = readImage(dataDir+"fakef.dat",n1,n2,n3)
  g = readImage(dataDir+"fakeg.dat",n1,n2,n3)
  s = readImage(dataDir+"fakes.dat",n1,n2,n3); s = mul(s,d1)
  #f,g,s = makeFakeImages(smax,nrms)
  mlag = 4+smax
  #dw = DynamicWarping(-mlag,mlag)
  #dw.setStrainMax(strainMax1,strainMax2,strainMax3)
  #dw.setErrorSmoothing(esmooth)
  #dw.setShiftSmoothing(usmooth)
  #dw.setWindowSizeAndOverlap(51,51,0.5,0.5)
  #u = dw.findShifts(f,g)
  #h = dw.applyShifts(u,g)
  dw = DynamicWarpingWT(-mlag,mlag,s1,s2,s3)
  dw.setStrainLimits(r1min,r1max,r2min,r2max,r3min,r3max)
  dw.setSmoothing(dr1,dr2,dr3)
  u = dw.findShifts(s1,s1,f,g)
  h = dw.applyShifts(f,u)
  print "s: min =",min(s),"max =",max(s)
  print "u: min =",min(u),"max =",max(u)
  usmin,usmax = min(s),max(s)
  show(s,usmin,usmax)
  show(u,usmin,usmax)
  print 'difference=',sum(sub(u,s))/(len(u[0][0])*len(u[0])*len(u))
  show(g)
  show(f)
  print type(h)
  show(h)

def goFakeShifts2():
  f = readImage(dataDir+"fakef.dat",n1,n2,n3)
  g = readImage(dataDir+"fakeg.dat",n1,n2,n3)
  s = readImage(dataDir+"fakes.dat",n1,n2,n3)
  f = f[sic]
  g = g[sic]
  s = s[sic]; s = mul(s,d1)
  pad = 0.1
  dw = DynamicWarpingWT(smin-pad,smax+pad,s1,s2)
  dw.setStrainLimits(r1min,r1max,r2min,r2max)
  dw.setSmoothing(dr1,dr2)
  dw.setInterpolationSpline()
  #dw.setInterpolationLinear()
  u = dw.findShifts(s1,s1,g,f)
  h = dw.applyShifts(g,u)
  dwb = DynamicWarpingR(smin-pad,smax+pad,s1,s2)
  dwb.setStrainLimits(r1min,r1max,r2min,r2max)
  dwb.setSmoothness(dr1,dr2)
  ub = dwb.findShifts(s1,g,s1,f)
  hb = dw.applyShifts(g,ub)
  print "smin =",smin-pad," smax =",smax+pad
  print "s: min =",min(s),"max =",max(s)
  print "u: min =",min(u),"max =",max(u)
  usmin,usmax = min(s),max(s)
  show2(ub,usmin,usmax,cmap=jet,title="ub")
  show2(hb,title="hb")
  show2(s,usmin,usmax,cmap=jet,title="s")
  show2(u,usmin,usmax,cmap=jet,title="u")
  #print 'difference=',sum(sub(us,s))/(len(us[0])*len(us))
  show2(f,title="f")
  show2(g,title="g")
  show2(h,title="h")

def goFakeShiftsBench():
  f = readImage(dataDir+"fakef.dat",n1,n2,n3)
  g = readImage(dataDir+"fakeg.dat",n1,n2,n3)
  s = readImage(dataDir+"fakes.dat",n1,n2,n3)
  f = f[sic]
  g = g[sic]
  s = s[sic]
  pad = 0.1
  dw = DynamicWarpingR(smin-pad,smax+pad,s1,s2)
  dw.setStrainLimits(r1min,r1max,r2min,r2max)
  dw.setSmoothness(dr1,dr2)
  u = dw.findShifts(s1,f,s1,g)
  h = dw.applyShifts(s1,g,u)
  print "smin =",smin-pad," smax =",smax+pad
  print "s: min =",min(s),"max =",max(s)
  print "u: min =",min(u),"max =",max(u)
  usmin,usmax = min(s),max(s)
  show2(s,usmin,usmax,cmap=jet,title="sb")
  show2(u,usmin,usmax,cmap=jet,title="ub")
  print 'difference=',sum(sub(u,s))/(len(u[0])*len(u))
  show2(g,title="gb")
  show2(f,title="fb")
  show2(h,title="hb")

def goFakeShifts2O():
  f = readImage(dataDir+"fakef.dat",n1,n2,n3)
  g = readImage(dataDir+"fakeg.dat",n1,n2,n3)
  s = readImage(dataDir+"fakes.dat",n1,n2,n3)
  f = f[sic]
  g = g[sic]
  s = s[sic]; s = mul(s,d1)
  pad = 0.1
  dw = DynamicWarpingWTM(f,g,round((smin-pad)/d1),round((smax+pad)/d1))
  u = dw.findShifts(r1min,r1max,dr1,r2min,r2max,dr2)
  #h,us = dw.applyShifts(u)
  h = dw.applyShiftsF(u)
  u = mul(u,d1)
  #us = mul(us,d1)
  print "smin =",smin-pad," smax =",smax+pad
  print "s: min =",min(s),"max =",max(s)
  print "u: min =",min(u),"max =",max(u)
  usmin,usmax = min(s),max(s)
  show2(s,usmin,usmax,cmap=jet,title="so")
  show2(u,usmin,usmax,cmap=jet,title="uo")
  #print 'difference=',sum(sub(us,s))/(len(us[0])*len(us))
  show2(f,title="fo")
  show2(g,title="go")
  show2(h,title="ho")



def makeFakeImages(smax,nrms):
  f = FakeData.seismic3d2010A(n1,n2,n3,20.0,10.0,30.0,0.5,0.0);
  w = Displacement3.sinusoid(0.5*smax,0.0,0.0,n1,n2,n3)
  #w = Warp3.constant(smax,0.0,0.0,n1,n2,n3)
  g = w.warp(f)
  f = addNoise(nrms,f,seed=10*seed+1)
  g = addNoise(nrms,g,seed=10*seed+2)
  s = w.u1x()
  return f,g,s

#############################################################################
# utilities

def readImage(fileName,n1,n2,n3):
  x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x

def writeImage(fileName,x):
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(x)
  aos.close()

def addRickerWavelet(fpeak,f):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  ih = int(3.0/fpeak)
  nh = 1+2*ih
  h = zerofloat(nh)
  for jh in range(nh):
    h[jh] = ricker(fpeak,jh-ih)
  g = zerofloat(n1,n2,n3)
  for i3 in range(n3):
    for i2 in range(n2):
      Conv.conv(nh,-ih,h,n1,0,f[i3][i2],n1,0,g[i3][i2])
  return g

def ricker(fpeak,time):
  x = PI*fpeak*time
  return (1.0-2.0*x*x)*exp(-x*x)

def addNoise(nrms,f,seed=0):
  n1,n2,n3 = len(f[0][0]),len(f[0]),len(f)
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  g = mul(2.0,sub(randfloat(r,n1,n2,n3),0.5))
  g = addRickerWavelet(0.125,g) # bandlimited in time
  rgf = RecursiveGaussianFilter(1.0)
  rgf.applyX0X(g,g)
  rgf.applyXX0(g,g)
  frms = sqrt(sum(mul(f,f))/n1/n2/n3)
  grms = sqrt(sum(mul(g,g))/n1/n2/n3)
  g = mul(g,nrms*frms/grms)
  return add(f,g)
 
#############################################################################
# plotting

jet = ColorMap.JET
gray = ColorMap.GRAY

def show(f,cmin=0.0,cmax=0.0):
  frame = SimpleFrame()
  ip = frame.addImagePanels(f)
  ip.setSlices(n1-1,n2/2,n3/2)
  if cmin<cmax:
    ip.setClips(cmin,cmax)
  frame.getOrbitView().setScale(2.0)
  frame.setSize(900,900)

def show2(f,cmin=0.0,cmax=0.0,cmap=gray,title=""):
  pp = PlotPanel()
  px = pp.addPixels(s1,s2,f)
  px.setColorModel(cmap)
  if cmin<cmax:
    px.setClips(cmin,cmax)
  pp.addColorBar()
  pp.setColorBarWidthMinimum(100)
  pp.setTitle(title)
  frame = PlotFrame(pp)
  frame.setSize(900,900)
  frame.setVisible(True)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
