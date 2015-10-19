"""
Jython utilities for BP 2D TTI dataset
Author: Andrew Munoz, Colorado School of Mines
Version: 2014.01.30
"""
from imports import *

#############################################################################
# Internal constants

_home = "/Users/amunoz/Home/data/bp/"
_sgydir = _home+"sgy/"
_datdir = _home+"dat/"

sigmanorm = 100

#############################################################################
# get data

def getGathers():
  """
  Reads in gathers with the correct samplings 
    samplings s1,s2,s3
  Example: getModel("vp")
  """
  global s1,so,s2
  n1,no,n2 = 1151,100,3599#12901
  d1,do,d2 = 0.008,0.1,0.00625 # s,km,km
  f1,fo,f2 = 0.0,0.0,0.0 # s,km,km
  s1,so,s2 = Sampling(n1,d1,f1),Sampling(no,do,fo),Sampling(n2,d2,f2)
  anidir = _datdir+"ani/"
  image = zerofloat(n1,no,n2)
  g = readImage(image,anidir+"cmps")
  return g

def getModel(which):
  """
  Reads in a model with the correct samplings 
    samplings s1,s2,s3
  Example: getModel("vp")
  """
  global s1,s2
  n1,n2 = 1801,12596
  d1,d2 = 0.00625,0.00625 # km
  f1,f2 = 0.0,0.0 # km
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  modeldir = _datdir+"model/"
  image = zerofloat(n1,n2)
  g = readImage(image,modeldir+which)
  return g

def getModelt(which):
  """
  Reads in a model with the correct samplings 
    samplings s1,s2,s3
  Example: getModel("vp")
  """
  global s1,s2
  n1,n2 = 1151,12596
  d1,d2 = 0.008,0.00625 # s,km
  f1,f2 = 0.0,0.0 # s,km
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  modeldir = _datdir+"model/"
  image = zerofloat(n1,n2)
  g = readImage(image,modeldir+which+"t")
  return g


def getStack():
  """
  Reads in a seismic stack with the correct samplings 
    samplings s1,s2,s3
  Example: getModel("vp")
  """
  global s1,s2
  n1,n2 = 1151,12596
  d1,d2 = 0.008,0.00625 # s, km
  f1,f2 = 0.0,0.0 # s, km
  s1,s2 = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  isodir = _datdir+"iso/"
  image = zerofloat(n1,n2)
  g = readImage(image,isodir+"bp_pstm_iso_stk")
  return g

def getSamplings():
  return s1,s2

def makeSubset(ix2):
  filePath = _datdir+"sub/"
  fileName = filePath+"sub_"+str(ix2)+".dat"
  if path.exists(fileName) and path.isfile(fileName): 
    print "Subset exists at",ix2,"..."
    return
  gat = getGathers()[ix2]
  vp = getModelt("vp")[ix2]
  de = getModelt("delta")[ix2]
  ep = getModelt("epsilon")[ix2]
  th = getModelt("theta")[ix2]
  print "Writing subset at",ix2,"..."
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(gat)
  aos.writeFloats(vp)
  aos.writeFloats(de)
  aos.writeFloats(ep)
  aos.writeFloats(th)
  aos.close()

def readSubset(ix2,xmax):
  n1,no = 1151,100
  d1,do = 0.008,0.1
  f1,fo = 0.0,-10.0
  s1,so = Sampling(n1,d1,f1),Sampling(no,do,fo)
  gat = zerofloat(n1,no)
  vp = zerofloat(n1)
  de = zerofloat(n1)
  ep = zerofloat(n1)
  th = zerofloat(n1)
  print "Reading subset at",ix2,"..."
  filePath = _datdir+"sub/"
  fileName = filePath+"sub_"+str(ix2)+".dat"
  ais = ArrayInputStream(fileName)
  ais.readFloats(gat)
  ais.readFloats(vp)
  ais.readFloats(de)
  ais.readFloats(ep)
  ais.readFloats(th)
  ais.close()
  #xod = makexod(so,s1,vp)
  gat = mutefa(s1,so,gat)
  wb = getwb(so,s1,gat)
  so,gat = fixOffset(so,gat,xmax)
  gat = normalizeRMSlocal(gat)
  return s1,so,gat,vp,de,ep,th,wb
  
def makexod(so,s1,vp):
  """
  computes max offset to depth ratio for gather
  """
  n1 = s1.count
  d1 = s1.delta
  mo = so.last*1000
  zt = zerofloat(n1)
  hd1 = 0.5*d1
  for i1 in range(1,n1):
    zt[i1] = zt[i1-1]+hd1*vp[i1]
  return div(zt,mo)

def getwb(so,s1,g):
  g0 = g[-1]
  n1 = len(g0)
  a = 0; wb=0;
  for i in range(n1):
    a = abs(g0[i])
    if a>0.1:
      wb = s1.getValue(i)
      break
  print 'water bottom is',wb,'s'
  return wb

def fixOffset(so,gat,xmax):
  no = len(gat)
  n1 = len(gat[0])
  li = 0
  lx = so.getValue(li)
  for io in range(no-1,0,-1):
    sumg = abs(sum(gat[io]))
    if sumg>0:
      continue
    else:
      lx = so.getValue(io);
      li = io
      print 'new xmax is',lx 
      break
  if xmax: 
      lx=xmax
      li =  so.indexOfNearest(lx)
      print 'new xmax is',lx 
  so = Sampling(so.count-int(((lx-so.first)/so.delta)),so.delta,lx)
  gat = copy(n1,so.count,0,li,gat)
  return so,gat

def mutefa(s1,so,g):
  nx = so.count
  nt = s1.count
  w = 50
  for ix in range(nx):
    t = -so.getValue(ix)/1.492
    ia = s1.indexOfNearest(t)
    ib = 0#max(0,ia-w)
    ie = min(nt,ia+w)
    for it in range(ib,ie):
      g[ix][it] = 0
  return g

def normalizeRMSlocal(x):
  return Semblance.localRMSnorm(x,sigmanorm)

#############################################################################
# read/write files

def readImage(image,fname):
  """ 
  Reads an image from a file with specified name.
  name: base name of image file; e.g., "tpsz"
  """
  print "Reading",fname,"..."
  fileName = fname+".dat"
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  print "min=",min(image),"max=",max(image)
  return image

def writeImage(image,fname):
  """ 
  Writes an image to a file with specified name.
  name: base name of image file; e.g., "tpgp"
  image: the image
  """
  print "Writing",fname,"..."
  print "min=",min(image),"max=",max(image)
  fileName = fname+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
  return image


#############################################################################
# graphics

def seeModelt(which):
  g = getModelt(which)
  plotModel(g,which)

def seeModelz(which):
  g = getModel(which)
  plotModel(g,which,d="Depth (km)")

jet = ColorMap.JET
gray = ColorMap.GRAY
def plotModel(g,which,sd=None,d="Time (s)",cmap=jet,perc=100):
  si = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  if sd:
    pi = si.addPixels(sd,s2,g)
  else:
    pi = si.addPixels(s1,s2,g)
  pi.setPercentiles(100-perc,perc)
  pi.setColorModel(cmap)
  si.setHLabel("Distance (km)")
  si.setVLabel(d)
  si.addColorBar(which)
  si.setSize(1200,500)

#############################################################################
# other

def makeTimeModels():
  global s1,s2
  n1z = 1801
  n1t,n2 = 1151,12596
  d1t,d1z,d2 = 0.008,0.00625,0.00625 # km
  f1t,f1z,f2 = 0.0,0.0,0.0 # km
  s1z,s1t,s2 = Sampling(n1z,d1z,f1z),Sampling(n1t,d1t,f1t),Sampling(n2,d2,f2)
  modeldir = _datdir+"model/"
  zt,tz = makezt(s1z,s1t,s2) # z(t)
  which = ["vp","delta","epsilon","theta"]
  for img in which:
    image = zerofloat(n1z,n2)
    fz = readImage(image,modeldir+img)
    ft = convertModelDepthToTime(zt,fz,s1z,s1t,s2,tz)
    writeImage(ft,modeldir+img+"t")

def makezt(s1z,s1t,s2):
  modeldir = _datdir+"model/"
  n2 = s2.count
  n1z = s1z.count; n1t = s1t.count
  d1z = s1z.delta; f1t = s1t.first
  image = zerofloat(n1z,n2)
  vpz = readImage(image,modeldir+"vp")
  tz = maketz(vpz,d1z,f1t)
  zt = zerofloat(n1t,n2)
  ii = InverseInterpolator(s1z,s1t)
  ii.invert(tz,zt)
  z = floats(s1z.getValues())
  t = floats(s1t.getValues())
  plotModel(zt,"zt",sd=s1t)
  plotModel(tz,"tz",sd=s1z,d="Depth (km)")
  return zt,tz

def maketz(vp,d1z,f1t):
  n2 = len(vp)
  n1 = len(vp[0])
  tz = zerofloat(n1,n2)
  d1z *= 2000
  for i2 in range(n2):
    tz[i2][0] = f1t
    for i1 in range(1,n1):
      tz[i2][i1] = tz[i2][i1-1]+d1z/vp[i2][i1]
  return tz

def convertModelDepthToTime(zt,fz,s1z,s1t,s2,tz):
  n2 = s2.count
  n1 = s1t.count
  nz = s1z.count 
  ft = zerofloat(n1,n2)
  #z = floats(s1z.getValues())
  #t = floats(s1t.getValues())
  #si = SincInterp()
  for i2 in range(n2):
    #ci = CubicInterpolator(CubicInterpolator.Method.LINEAR,nz,z,fz[i2])
    #ci.interpolate(n1,zt[i2],ft[i2])
    for i1 in range(n1):
      #ft[i2][i1] = si.interpolate(s1z,s2,fz,zt[i2][i1],s2.getValue(i2))
      #ft[i2][i1] = ci.interpolate(zt[i2][i1])
      ft[i2][i1] = fz[i2][s1z.indexOfNearest(zt[i2][i1])]
  return ft

def makeEffectiveParamsAll(s1,vp,eta,de):
  n2 = s2.getCount()
  n1 = s1.getCount()
  vnmo2 = zerofloat(n1,n2)
  eeta2 = zerofloat(n1,n2)
  for i2 in range(n2):
    vnmo2[i2],eeta2[i2] = makeEffectiveParams(s1,vp[i2],eta[i2],de[i2])
  plotModel(vnmo2,"Effective Vnmo")
  plotModel(eeta2,"Effective Eta")

  
def makeEffectiveParams(s1,vp,eta,de):
  n = len(vp)
  dt = s1.delta
  vrms = zerofloat(n)
  eeta = zerofloat(n)
  vpa = mul(vp,sqrt(add(1,mul(2,de))))
  vps = mul(vpa,vpa)
  vpq = mul(vps,vps)
  to = floats(s1.getValues())
  vrms[0] = vps[0]*dt
  for i in range(1,n):
    vrms[i] = vrms[i-1]+vps[i]*dt
  dts = dt
  for i in range(n):
    vrms[i] /= dts
    dts += dt
  eeta[0] = vpq[0]*(1+8*eta[0])*dt
  for i in range(1,n):
    eeta[i] = eeta[i-1]+vpq[i]*(1+8*eta[i])*dt
  vrms4 = zerofloat(n); dts = dt
  for i in range(n):
    vrms4[i] = vrms[i]*vrms[i]*dts
    dts += dt
  eeta = mul(sub(div(eeta,vrms4),1),0.125)
  return mul(sqrt(vrms),0.001),eeta

def addNoise(s,x):
  """ adds s percent random noise to reduce problems with zero traces """
  xdif = max(x)-min(x)
  n1,n2,n3 = s1.count,s2.count,s3.count
  r = randfloat(n1,n2,n3)
  x = add(mul(sub(r,0.5),s*xdif),x)
  return x

def floats(x):
  """
  Converts 1D array to 1D float array
  """
  n = len(x)
  xd = zerofloat(n)
  for i in range(n):
    xd[i] = float(x[i])
  return xd

#############################################################################
# testing 

def main(args):
  #makeTimeModels()
  #seeModelt("vp")
  #seeModelz("vp")
  #seeModelt("delta")
  #seeModelz("delta")
  #seeModelt("epsilon")
  #seeModelz("epsilon")
  #seeModelt("theta")
  #seeModelz("theta")

  #vp = getModelt("vp")
  #de = getModelt("delta")
  #ep = getModelt("epsilon")
  #eta = div(sub(ep,de),add(1,mul(2,de)))
  #makeEffectiveParamsAll(s1,vp,eta,de)

  #g = getStack()
  #plotModel(g,"amplitude",cmap=gray,perc=98)

  

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

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
