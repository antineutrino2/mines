#############################################################################
# Dynamic time warping for 2D images

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
from edu.mines.jtk.sgl import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

from dtw import *
#from dtw.Util import *

#############################################################################

#pngDir = "./png"
pngDir = './migpics/Layered/'
#pngDir = './migpics/Marmousi/'
saveDir = './tmpdata/'

seed = abs(Random().nextInt()/1000000)
seed = 580
seed = 588
#print "seed =",seed

nrms = 2.00
npass = 5
stretchMax = 0.25
showLsf = False
smoothShifts = False
nv = 51 # number of datasets
_s = 1  # number of shifts
J1 = zerofloat(nv)
J2 = zerofloat(nv)
J3 = zerofloat(nv)
J4 = zerofloat(nv)

J41 = zerofloat(nv)
y0 = 150

def main(args):
  #go1Dwrite()
  #goTestImages()
  #goFaultImages()
  #goTestShifts() #smax = 0.20, nrms = 2.0
  #goFaultShifts() #smax = 0.25, npass = 3
  readJs()
  goImageShift2()
  writeJs()

def goImageShift2():
  for s in range(1,_s+1):
    for it in range(5,10):
      #f = readMMImage(25)
      #g = readMMImage(it)
      f = readLayImage(50)
      g = readLayImage(it)
      _shift = s
      h1 = goImageShift(f,g,it,_shift)
      if (s==2 or s>2):
        _shift = s
        h2 = goImageShift(f,h1,it,_shift)
        if (s==3 or s>3):
          _shift = s
          h3 = goImageShift(f,h2,it,_shift)
          if (s==4 or s>4):
            _shift = s
            h4 = goImageShift(f,h3,it,_shift)
            if s==5:
              _shift = s
              h5 = goImageShift(f,h4,it,_shift)
    goPlotJ(s)
  """
  g = zerofloat(len(f[0]),len(f))
  readImage("h5",g)
  h1 = goImageShift(f,g,it)
  h2 = goImageShift(f,h1,it)
  h3 = goImageShift(f,h2,it)
  h4 = goImageShift(f,h3,it)
  h5 = goImageShift(f,h4,it)
  #writeImage("h5",h5)
  writeImage("h10",h5)
  """
def goPlotJ(s):
  plotJ(J1,s,png='F-H')
  plotJ(J2,s,png='U*H')
  plotJ(J3,s,png='U')
  plotJ(J4,s,png='F-G')
  #plotJ(J41,s,png='F-G_1D')

def writeJs():
  writeImage('J1',J1)
  writeImage('J2',J2)
  writeImage('J3',J3)
  writeImage('J4',J4)

def readJs():
  readJ('J1',J1)
  readJ('J2',J2)
  readJ('J3',J3)
  readJ('J4',J4)

def readJ(name,J):
  fileName = saveDir+name+'.dat'
  ais = ArrayInputStream(fileName)
  ais.readFloats(J)
  ais.close()

def goImageShift(f,g,it,_shift):
  nm = ":v"+str(it)
  shift = 25 
  sigma = shift
  ml = 2*shift
  uclips = (-shift,shift)
  dw = DynamicWarping(-ml,ml)
  dw.setStretchMax(stretchMax)
  #fclips = (min(f),max(f))
  fclips = None
  plot(f,"f"+nm,fclips,png='f'+str(it))#+'_shift='+str(_shift))
  plot(g,"g"+nm,fclips,png='g'+str(it))#+'_shift='+str(_shift))
  plot(sub(f,g),"dc-do"+nm,fclips,png='f-g'+str(it))#+'_shift='+str(_shift))
  #plot1(f[y0],"f"+nm)
  #plot1(g[y0],"g"+nm)
  #plot1(sub(f[y0],g[y0]),"dc-do"+nm)
  #objectiveFunctionFG1(f,g,it)
  #h = 1
  e = dw.computeErrors(f,g)
  if npass>=5:
    u = shifts12121(dw,e)
    print "errors  u12121: sum =",dw.sumErrors(e,u)
    plot(u,"u12121"+nm,uclips,png='u'+str(it))#+'_shift='+str(_shift))
    objectiveFunctionU(u,it)
  h = align(u,f,g)
  objectiveFunctionFH(f,h,it)
  objectiveFunctionFG(f,g,it)
  objectiveFunctionUH(u,h,it)
  plot(h,"h"+nm,fclips,png='h'+str(it))#+'_shift='+str(_shift))
  plot(sub(f,h),"dc-h"+nm,fclips,png='f-h'+str(it))#+'_shift='+str(_shift))
  plot(mul(u,g),"u*g"+nm,fclips,png='u*g'+str(it))#+'_shift='+str(_shift))
  return h

def readMMImage(t):
  ddir = "../../../../../data/migclass/"
  num = str(t)
  n1,n2 = 421,5101
  v = zerofloat(n2,n1)
  fileName = ddir+"marmousiTest-v"+num+".dat"
  ais = ArrayInputStream(fileName)
  ais.readFloats(v)
  ais.close()
  return v
  
def readLayImage(t):
  ddir = "../../../../../data/migclass/"
  num = str(t)
  n1,n2 = 301,3501
  v = zerofloat(n2,n1)
  fileName = ddir+"layered-sineSourceTapered2-"+num+".dat"
  ais = ArrayInputStream(fileName)
  ais.readFloats(v)
  ais.close()
  return v

def objectiveFunctionFH(f,h,it):
  n1,n2 = len(f[0]),len(f)
  sumu = 0
  for i2 in range(n2):
    for i1 in range(n1):
      u = f[i2][i1]-h[i2][i1]
      uu = u*u 
      sumu += uu
  J1[it] = 0.5*sumu

def objectiveFunctionUH(u,h,it):
  n1,n2 = len(u[0]),len(u)
  sumu = 0
  for i2 in range(n2):
    for i1 in range(n1):
      uu = u[i2][i1]*u[i2][i1]
      hh = h[i2][i1]*h[i2][i1]
      sumu += uu*hh
  J2[it] = 0.5*sumu

def objectiveFunctionU(u,it):
  n1,n2 = len(u[0]),len(u)
  sumu = 0
  for i2 in range(n2):
    for i1 in range(n1):
      uu = u[i2][i1]*u[i2][i1]
      sumu += uu
  J3[it] = 0.5*sumu

def objectiveFunctionFG(f,g,it):
  n1,n2 = len(f[0]),len(f)
  sumu = 0
  for i2 in range(n2):
    for i1 in range(n1):
      u = f[i2][i1]-g[i2][i1]
      uu = u*u 
      sumu += uu
  J4[it] = 0.5*sumu

def objectiveFunctionFG1(f,g,it):
  n = len(f[0])
  sumu = 0
  for i in range(n):
    u = f[y0][i]-g[y0][i]
    uu = u*u 
    sumu += uu
  J41[it] = 0.5*sumu

def goTestShifts():
  shift = 16
  sigma = shift
  ml = 2*shift
  uclips = (-shift,shift)
  dw = DynamicWarping(-ml,ml)
  dw.setStretchMax(stretchMax)
  f,g,s = makeTestImages()
  fclips = (min(f),max(f))
  plot(f,"f",fclips)
  plot(g,"g",fclips)
  e = dw.computeErrors(f,g)
  if npass>=1:
    u = shifts1(dw,e);
    print "errors      u1: sum =",dw.sumErrors(e,u)
    plot(u,"u1",uclips)
  if npass>=2:
    u = shifts12(dw,e)
    print "errors     u12: sum =",dw.sumErrors(e,u)
    plot(u,"u12",uclips)
  if npass>=3:
    u = shifts121(dw,e)
    print "errors    u121: sum =",dw.sumErrors(e,u)
    plot(u,"u121",uclips)
  if npass>=4:
    u = shifts1212(dw,e)
    print "errors   u1212: sum =",dw.sumErrors(e,u)
    plot(u,"u1212",uclips)
  if npass>=5:
    u = shifts12121(dw,e)
    print "errors  u12121: sum =",dw.sumErrors(e,u)
    plot(u,"u12121",uclips)
  if showLsf:
    v = copy(u)
    LocalShiftFinder(ml,sigma).find1(-ml,ml,f,g,v)
    print "errors   u lsf: sum =",dw.sumErrors(e,v)
    plot(v,"u lsf",uclips)
  if s:
    print "errors       s: sum =",dw.sumErrors(e,s)
    plot(s,"s")
  h = align(u,f,g)
  plot(h,"h",fclips)

def goTestImages():
  f,g,s = makeTestImages()
  plot(f,"f")
  plot(g,"g")
  plot(s,"s")


def makeTestImages():
  dip = 30.0
  shift = 16
  n1,n2 = 501,501; f = FakeData.seismic2d2011A(n1,n2,dip)
  #n1,n2 = 462,951; f = readImage("/data/seis/f3d/f3d75.dat",n1,n2)
  f = sub(f,sum(f)/n1/n2)
  #w = Warp2.constant(shift,0.0,n1,n2)
  w = Warp2.sinusoid(shift,0.0,n1,n2)
  g = w.warp(f)
  f = addNoise(nrms,f,seed=10*seed+1)
  g = addNoise(nrms,g,seed=10*seed+2)
  s = zerofloat(n1,n2)
  for i2 in range(n2):
    for i1 in range(n1):
      s[i2][i1] = w.u1x(i1,i2)
  return f,g,s

def go1Dwrite():
  y = 150
  f = readLayImage(5)
  g = readLayImage(50)
  fng = zerofloat(len(f[y]),2)
  for i in range(len(f[y])):
    fng[0][i] = f[y][i]
    fng[1][i] = g[y][i]
  writeImage("f=5_g=50_y=150_len="+str(len(f[y])),fng)

#############################################################################
# shifts

def smooth(u):
  v = copy(u)
  if smoothShifts:
    rgf = RecursiveGaussianFilter(8); rgf.apply00(u,v)
  return v

def normalize(x):
  xmin = min(x)
  xmax = max(x)
  return mul(sub(x,xmin),1.0/(xmax-xmin))

def align(u,f,g):
  n1,n2 = len(u[0]),len(u)
  si = SincInterpolator()
  si.setUniformSampling(n1,1.0,0.0)
  h = copy(g)
  r = rampfloat(0.0,1.0,n1)
  for i2 in range(n2):
    t = add(r,u[i2])
    si.setUniformSamples(g[i2])
    si.interpolate(n1,t,h[i2])
  return h

def shifts1(dw,e):
  d = dw.accumulateForward1(e)
  u = dw.findShiftsReverse1(d,e)
  return smooth(u)

def shifts12(dw,e):
  e = dw.accumulate1(e)
  e = normalize(e)
  d = dw.accumulateForward2(e)
  u = dw.findShiftsReverse2(d,e)
  return smooth(u)

def shifts121(dw,e):
  e = dw.accumulate1(e)
  e = normalize(e)
  e = dw.accumulate2(e)
  e = normalize(e)
  d = dw.accumulateForward1(e)
  u = dw.findShiftsReverse1(d,e)
  return smooth(u)

def shifts1212(dw,e):
  e = dw.accumulate1(e)
  e = normalize(e)
  e = dw.accumulate2(e)
  e = normalize(e)
  e = dw.accumulate1(e)
  e = normalize(e)
  d = dw.accumulateForward2(e)
  u = dw.findShiftsReverse2(d,e)
  return smooth(u)

def shifts12121(dw,e):
  e = dw.accumulate1(e)
  e = normalize(e)
  e = dw.accumulate2(e)
  e = normalize(e)
  e = dw.accumulate1(e)
  e = normalize(e)
  e = dw.accumulate2(e)
  e = normalize(e)
  d = dw.accumulateForward1(e)
  u = dw.findShiftsReverse1(d,e)
  return smooth(u)

#############################################################################
# utilities

def readImage(fileName,n1,n2):
  x = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()
  return x

def readImage(name,x):
  fileName = saveDir+name+'.dat'
  ais = ArrayInputStream(fileName)
  ais.readFloats(x)
  ais.close()

def writeImage(name,image):
  fileName = saveDir+name+'.dat'
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()

def addRickerWavelet(fpeak,f):
  n1,n2 = len(f[0]),len(f)
  ih = int(3.0/fpeak)
  nh = 1+2*ih
  h = zerofloat(nh)
  for jh in range(nh):
    h[jh] = ricker(fpeak,jh-ih)
  g = zerofloat(n1,n2)
  for i2 in range(n2):
    Conv.conv(nh,-ih,h,n1,0,f[i2],n1,0,g[i2])
  return g

def ricker(fpeak,time):
  x = PI*fpeak*time
  return (1.0-2.0*x*x)*exp(-x*x)

def addNoise(nrms,f,seed=0):
  n1,n2 = len(f[0]),len(f)
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  g = mul(2.0,sub(randfloat(r,n1,n2),0.5))
  g = addRickerWavelet(0.125,g) # same wavelet used in signal
  rgf = RecursiveGaussianFilter(1.0)
  rgf.applyX0(g,g)
  frms = sqrt(sum(mul(f,f))/n1/n2)
  #frms = max(abs(f))
  grms = sqrt(sum(mul(g,g))/n1/n2)
  g = mul(g,nrms*frms/grms)
  return add(f,g)
 
#############################################################################
# plotting

def plot(f,title=None,clips=None,png=None):
  #sp = SimplePlot.asPixels(f)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(f)
  if clips:
    pv.setClips(clips[0],clips[1])
  if title:
    sp.setTitle(title)
  sp.addColorBar()
  sp.setSize(900,900)
  if png:
    sp.paintToPng(360,3.33,pngDir+png+'.png') 

def plot1(f,title=None):
  sp = SimplePlot()
  sp.addPoints(f)
  sp.setSize(900,500)
  if title:
    sp.setTitle(title)

def plotJ(J,it,png=None):
  it = str(it)
  sp = SimplePlot()
  pv = sp.addPoints(J)
  sp.setHLabel("Velocity Models")
  #sp.setHInterval(1.0)
  sp.setVLabel("J- Objective Function")
  sp.setTitle("Data misfit "+png)
  sp.setSize(1200,500)
  if png:
    sp.paintToPng(360,3.33,pngDir+png+'-'+it+'.png')

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

    #main(sys.argv)
SwingUtilities.invokeLater(RunMain())
