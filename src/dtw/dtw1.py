#############################################################################
# Dynamic time warping for 1D sequences

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

from dtw import DynamicWarpingWT
from wt import WTutils

#############################################################################

_pngDir = "./dtwpics/" 

def main(args):
	n1 = 301 
	n2 = 1001
	f0 = 200 #true f0
	##dt = 0.01
	rmsNoise = 0.0
	maxStrain = 0.5
	g,f = makeSin(n2,n1,f0,maxStrain,rmsNoise,C=1.0) #true g0
	#g,f = makeLin(n2,n1,f0,maxStrain,rmsNoise) #true g0
	
	#goDTW1(n1,n2,f0,dt,rmsNoise,maxStrain,g,f)
	#goDTW2(n1,n2,f0,dt,rmsNoise,maxStrain,g,f)
	#goDTW3()
	#goDTWR(n1,n2,f0,f,g)
	goTest(n1,n2,f0,f,g)



def goTest(nf,ng,f0,f,g):
	dt = 0.001
	smin = 0
	smax = ng-nf
	strainlim = 1.0
	dw = DynamicWarpingWells(smin,smax)
	dw.setStrainMax(strainlim)
	u = dw.findShifts(f,g)
	h = dw.applyShifts(u,f)
	u = mul(u,dt)

	sh = Sampling(len(h),dt,u[0])	
	sf = Sampling(nf,dt,0.0)
	sg = Sampling(ng,dt,0.0) 

	plotSequences(g,f,sg,sf)
	plotSequences(g,h,sg,sh)


def goDTWR(nf,ng,f0,f,g):
	nl = ng-nf+1
	dr = 0.5
	dt = 0.01
	dst = dt
	sf = Sampling(nf,dt,0.0)
	sg = Sampling(ng,dt,0.0) #guess g0
	sl = Sampling(nl,dt,0.0)
	rmin = -1.0
	rmax =  1.0
	print 'kmax=',int(ceil(rmax/dr))
	print 'kmin=',int(floor(rmin/dr))
	
	dw = DynamicWarpingWT(nl,rmin,rmax,dr)
	e = dw.computeErrors(f,g)
	u = dw.findShiftsR(e)
	h = dw.applyShifts(u,f)
	u = mul(u,dt)
	sh = Sampling(len(h),dt,u[0])
  
	plotu(u)
	plotSequences(g,f,sg,sf)
	plotSequences(g,h,sg,sh)
	plotMatrix(etran(e),sf,sl)
	plotMatrix(etran(e),sf,sl,u=u)


def goDTW1(n1,n2,f0,dt,rmsNoise,maxStrain,g,f):
  vrmin  = 0.001# minimum velocity ratio v*/v (>=0.5) 
  vrmax  = 2.00 # maximum velocity ratio v*/v (<infinity)
  dvrmin = 0.001# minimum change of velocity (>=0.0)
  dvrmax = 2.00 # maximum change of velocity (<=2.0)
  dvr    = 0.50 # percent velocity change
  
  fr,ds,kmin,kmax,jmin,jmax = paramDTW(vrmin,vrmax,dvrmin,dvrmax,dvr) 
  
  nl = 1+fr*(n2-n1)
  dst = dt/fr
  sf = Sampling(n1,dt,0.0)
  sg = Sampling(n2,dt,0.0) #guess g0
  sl = Sampling(nl,dst,0.0)
  
  #jmin,jmax = 0,0
  #kmin,kmax = 0,0
  sw = Stopwatch()
  dw = DynamicWarpingWT(nl,fr,kmin,kmax,jmin,jmax)
  sw.start()
  e = dw.computeErrors(f,g)
  d = dw.accumulate(e)
  u = dw.findShifts(d)
  #u = dw.findShiftsFast(f,g)
  sw.stop()
  print 'dtw done at '+str(sw.time())+' seconds'
  h = dw.applyShifts(u,f)
  u = mul(u,dt)
  sh = Sampling(len(h),dt,u[0])
  #m = floats(dw.getm())


  #plotSequences(g,f,sg,sf)
  #plotSequences(g,h,sg,sh)
  #plotMatrix(etran(e),sf,sl,slides='sse')
  #plotMatrix(etran(d),sf,sl,cp=99.9,slides='ssd')
  ##plotMatrix(mtran(m),sf,sl)
  #plotMatrix(etran(e),sf,sl,u=u,slides='sseu')
  #plotMatrix(etran(d),sf,sl,u=u,cp=99.9,slides='sseu')
  #plotMatrix(etran(d),sf,sl,u=u,lim=1,cp=99.9)
  #plotu(u,sf)

  # For presenations/papers
  #plotCurvePanel3(sg,sf,g,f,slides='slagmin')
  #plotCurvePanel3(sg,Sampling(n1,dt,(n2-n1)*dt),g,f,slides='slagmax')
  #plotCurvePanel3(sg,sf,g,f,sh=sh,h=h,slides='spaneltie')

def goDTW2(n1,n2,f0,dt,rmsNoise,maxStrain,g,f):
 
  #for uo
  fr = 2
  kmin = -2
  kmax =  2

  rmin = -1.5
  rmax =  1.5

  nl = 1+fr*(n2-n1)
  dst = dt/fr
  sf = Sampling(n1,dt,0.0)
  sl = Sampling(nl,dst,0.0)
  dw = DynamicWarpingWT(nl,fr,kmin,kmax)
  uoi,udmin = dw.finduo(f,g)
  uo = uoi[0]
  #uo = 350 
  print 'uo= ',uo,'\n'
  e = dw.computeErrors(f,g)
  d = dw.accumulate(e)
  for ir in range(4):
    dr = 0.001*pow(10,ir)
    nr = int(floor((rmax-rmin)/dr)) + 1
    sr = Sampling(nr,dr,rmin)

    u = dw.findShiftsX(e,sr,uo)
    u = mul(u,dt)
    dmin = dw.findShiftsX(e,sr,uo,True)
     
    print 'dr= ',dr
    print 'nr= ',nr
    print 'esum= ',dmin[0],'\n'

    #plotMatrix(etran(d),sf,sl)
    plotMatrix(etran(d),sf,sl,u=u,title='dr= '+str(dr))

def goDTW3():
  # errors
  nt = 501
  ns = 101
  lw = 1
  uo = 20
  bg = 10
  #e = makeErrorTest1(nt,ns,uo,bg,lw)
  e = makeErrorTest2(nt,ns,uo,bg,5,0)
  st = Sampling(nt,1,0)
  ss = Sampling(ns,1,0)

  # warping
  rmax =  1.5
  rmin = -1.5
  dw = DynamicWarpingWT(ns,1)
  for ir in range(3):
    dr = 0.01*pow(10,ir)
    nr = int(floor((rmax-rmin)/dr)) + 1
    sr = Sampling(nr,dr,rmin)

    u = dw.findShiftsX(e,sr,uo)
    #es = esum(e,u)
     
    print 'dr= ',dr
    print 'nr= ',nr
    #print 'esum= ',es,'\n'

    plotMatrix(etran(e),st,ss,u=u,title='dr= '+str(dr)+', rmin='+str(rmin)+', rmax='+str(rmax),haxis="time samples",vaxis="shift samples")#,slides='etest'+str(ir+1))

  # plotting
  plotMatrix(etran(e),st,ss,haxis="time samples",vaxis="shift samples",title="errors")#,slides='etest0')



def makeErrorTest1(nt,ns,uo,bg,lw):
  e = popnoise(ns,nt,old='etest1')
  n2 = len(e)
  uo = uo-lw+1
  for i2 in range(bg,n2):
    for i1 in range(lw):
      e[i2][uo+i1] = 0
  return e

def makeErrorTest2(nt,ns,uo,bg,mh,mv):
  e = popnoise(ns,nt,old='etest1')
  i2 = bg
  i1 = uo
  ntm = nt-1
  nsm = ns-1
  mh+=1
  mv+=1
  while i2<ntm:
    for ih in range(mh):   
      i2+=1 
      if i2>=ntm: 
        i2=ntm
      e[i2][i1] = 0
    for iv in range(mv):   
      i1+=1 
      if i1>=nsm: 
        i1=nsm
      e[i2][i1] = 0
  return e

def popnoise(n1,n2,new=None,old=None):
  if new:
    r = Random(235324)
    e = randfloat(r,n1,n2)
    rgf = RecursiveGaussianFilter(1)
    rgf.apply00(e,e)
    aos = ArrayOutputStream("./tmpdata/"+new+".dat")
    aos.writeFloats(e)
    aos.close()
  if old:
    e = zerofloat(n1,n2)
    ais = ArrayInputStream("./tmpdata/"+old+".dat")
    ais.readFloats(e)
    ais.close()
  return e


"""
Parameterize dtw
"""
def paramDTW(vrmin,vrmax,dvrmin,dvrmax,dvr):
  dumax = vrmax - 1
  dumin = vrmin - 1
  durmax = dvrmax - 1
  durmin = dvrmin - 1
  #ds1 = 1/vrmax - 1/(vrmax+dvr)
  #ds2 = 1/(vrmin-dvr) - 1/(vrmin)
  # for small dvr:
  ds1 = dvr/(vrmax*vrmax)
  ds2 = dvr/(vrmin*vrmin)
  ds = min(ds1,ds2)
  fr = int(ceil(1/ds))
  kmax= int(floor(dumax/ds))
  kmin= int(ceil(dumin/ds))
  jmax= min(kmax,int(floor(durmax/ds)))
  jmin= max(kmin,int(ceil(durmin/ds)))
  print 'fr=',fr
  print 'ds=',ds
  print 'kmax=',kmax
  print 'kmin=',kmin
  print 'jmax=',jmax
  print 'jmin=',jmin
  return fr,ds,kmin,kmax,jmin,jmax

def smooth(u):
  v = copy(u)
  rgf = RecursiveGaussianFilter(4)
  rgf.apply0(u,v)
  return v

def normalize(e):
  emin = min(e)
  emax = max(e)
  return mul(sub(e,emin),1.0/(emax-emin))

def etran(e):
  e = normalize(e)
  return transpose(pow(e,0.25))
def mtran(e):
  return transpose(e)

def makeSequences():
  n1 = 1501
  n2 = 501
  fpeak = 0.0125
  fpeak2 = 0.0110
  fg = 100
  shift = 2.0/fpeak
  #w = Warp1.constant(shift,n)
  #w = Warp1.sinusoid(shift,n)
  f = makeCosine(fpeak,n1)
  #f = makeRandomEvents(n1,seed=seed) 
  g = makeCosine(fpeak2,n2,fg)
  #g = w.warp(f)
  #f = addRickerWavelet(fpeak,f)
  #g = addRickerWavelet(fpeak,g)
  #f = addNoise(nrms,fpeak,f,seed=10*seed+1)
  #g = addNoise(nrms,fpeak,g,seed=10*seed+2)
  #s = zerofloat(n)
  #for i in range(n):
  #  s[i] = w.ux(i)
  sf = Sampling(n1,1.0,0.0)
  sg = Sampling(n2,1.0,fg)
  return f,g,sf,sg#,s

def getSequences():
  name = "./tmpdata/f=5_g=50_y=150_len=3501.dat"
  nf = 3501
  fng = zerofloat(nf,2)
  ais = ArrayInputStream(name)
  ais.readFloats(fng)
  ais.close()
  f = zerofloat(nf)
  g = zerofloat(nf)
  for i in range(nf):
    f[i] = fng[0][i]
    g[i] = fng[1][i]
  return f,g

def makeSin(n1,n2,start,maxStrain,rmsNoise,C=None):
  si = SincInterp()
  si.setExtrapolation(SincInterp.Extrapolation.CONSTANT)
  seedN2 = 10000
  fpeak  = 0.125
  amp = maxStrain*(n2/(2*PI))
  g = zerofloat(n2)
  t = zerofloat(n2)
  f = sig1(n1,rmsNoise,fpeak)
  w = float((2*PI/n2))
  if C:
    w = float((2*PI/n2))*C
  for i in range (0,n2):
    t[i]=start+i-amp*sin(i*w)   
  si.interpolate(n1,1.0,0.0,f,n2,t,g)    
  g = addNoise(rmsNoise,fpeak,g,seedN2)
  #g = add(g,50) # To test DDTW
  return f,g

def makeLin(n1,n2,start,maxStrain,rmsNoise):
  si = SincInterp()
  si.setExtrapolation(SincInterp.Extrapolation.CONSTANT)
  seedN2 = 10000
  fpeak  = 0.125
  amp = maxStrain
  g = zerofloat(n2)
  t = zerofloat(n2)
  f = sig1(n1,rmsNoise,fpeak)
  for i in range (0,n2):
    t[i]=start+i-amp*i
  si.interpolate(n1,1.0,0.0,f,n2,t,g)    
  g = addNoise(rmsNoise,fpeak,g,seedN2)
  #g = add(g,50) # To test DDTW
  return f,g
 
def sig1(n1,rmsNoise, fpeak):
  seed   = 23232
  seedN1 = 31415 
  f = makeRandomEvents(n1,seed) 
  f = addRickerWavelet(fpeak,f)
  f = mul(1.0/max(abs(f)),f)
  f = addNoise(rmsNoise,fpeak,f,seedN1)
  return f

def makeCosine(freq,n,i=None):
  if i:
    return cos(mul(2.0*PI*freq,rampfloat(i,1.0,n)))
  return cos(mul(2.0*PI*freq,rampfloat(0.0,1.0,n)))

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
  n2 = len(x[0])
  y = zerofloat(n2,n1)
  for i1 in range(n1):
    for i2 in range(n2):
      y[i1][i2] = float(x[i1][i2])
  return y

 
#############################################################################
# plotting

def plotWithMatrix(f,g,c,uf,ug):
  n = len(f)
  panel = PlotPanel(2,2,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  panel.mosaic.setWidthElastic(0,25)
  panel.mosaic.setHeightElastic(0,25)
  panel.mosaic.setWidthElastic(1,100)
  panel.mosaic.setHeightElastic(1,100)
  fv = panel.addPoints(0,1,f)
  gv = panel.addPoints(0,1,g)
  gv.setLineColor(Color.RED)
  fv.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
  gv.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
  gv = panel.addPoints(1,0,g)
  fv = panel.addPoints(1,0,f)
  fv.setLineColor(Color.RED)
  cv = panel.addPixels(1,1,c)
  cv.setInterpolation(PixelsView.Interpolation.NEAREST)
  #cv.setClips(-1.0,1.0)
  cv.setColorModel(ColorMap.GRAY)
  uv = panel.addPoints(1,1,ug,uf)
  uv.setLineColor(Color.RED)
  dv = panel.addPoints(1,1,(0.0,n-1.0),(0.0,n-1.0))
  dv.setLineColor(Color.WHITE)
  panel.setHLabel(1,"sample")
  panel.setVLabel(1,"sample")
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(18)
  frame.setSize(1000,1000)
  frame.setVisible(True)

def plot(f,g,c,s1=None,s2=None,h=None,sh=None,u=None,clip=None):
  n2,nlag = len(c[0]),len(c)
  n1 = len(f)
  if s1 is None:
    s1 = Sampling(n1,1.0,0.0)
  if s2 is None:
    s2 = Sampling(n2,1.0,0.0)
  fs2 = s2.getFirst()
  slag = Sampling(nlag,1.0,-(nlag-1)/2-fs2)
  panel = PlotPanel(2,1,PlotPanel.Orientation.X1RIGHT_X2UP)
  panel.mosaic.setHeightElastic(0,50)
  panel.mosaic.setHeightElastic(1,50)
  panel.setVLimits(1,slag.first,slag.last)
  fv = panel.addPoints(0,0,s1,f)
  gv = panel.addPoints(0,0,s2,g)
  gv.setLineColor(Color.RED)
  if h:
    hv = panel.addPoints(0,0,sh,h)
    hv.setLineColor(Color.BLUE)
  cv = panel.addPixels(1,0,s2,slag,c)
  cv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if clip:
    cv.setClips(0.0,clip)
  cv.setColorModel(ColorMap.JET)
  if u:
    uv = panel.addPoints(1,0,s2,u)
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
  panel.setHLabel("sample")
  panel.setVLabel(0,"f & g")
  panel.setVLabel(1,"lag")
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(18)
  frame.setSize(1000,800)
  frame.setVisible(True)
  
def plotMatrix(c,sf,sl,u=None,lim=None,vaxis=None,haxis=None,cp=None,title=None,png=None,slides=None):
  n1,nlag = len(c[0]),len(c)
  #slag = Sampling(nlag,1.0,-(nlag-1)/2)
  panel = PlotPanel(PlotPanel.Orientation.X1RIGHT_X2UP)
  cv = panel.addPixels(sf,sl,c)
  cv.setInterpolation(PixelsView.Interpolation.NEAREST)
  cv.setColorModel(ColorMap.JET)
  if u:
    uv = panel.addPoints(sf,u)
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
  panel.setVLimits(sl.first,sl.last)
  panel.setHLimits(sf.first,sf.last)
  panel.setHLabel("time (s)")
  panel.setVLabel("time lag (s)")
  panel.addColorBar()
  if title:
    panel.setTitle(title)
  if vaxis:
    panel.setVLabel(vaxis)
  if haxis:
    panel.setHLabel(haxis)
  if lim:
    panel.setVLimits(lim*min(u),lim*max(u))
  if cp:
    cv.setPercentiles(0,cp)
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(18)
  frame.setSize(1000,800)
  frame.setVisible(True)
  if png:
    frame.paintToPng(720,3.33,_pngDir+png+'.png')
  if slides:
    frame.setSize(550,500)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,3.50,_pngDir+slides+'.png')
  
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
  panel.setVLabel("time (s)")
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(18)
  frame.setSize(300,800)
  frame.setVisible(True)

def plotu(u,s=None):
  si = SimplePlot()
  if s:
    si.addPoints(s,u)
  else:
    si.addPoints(u)

"""
Plots 3 curves on top of eachother
"""
def plotCurvePanel3(sf,sg,f,g,sh=None,h=None,tlim=None,paper=None,slides=None):
  if h:
    p1 = PlotPanel(3,1,PlotPanel.Orientation.X1RIGHT_X2UP,PlotPanel.AxesPlacement.LEFT_TOP)
    pf = p1.addPoints(0,0,sf,f)
    pg = p1.addPoints(1,0,sg,g)
    pg2f = p1.addPoints(0,0,sh,h)
    pg2g = p1.addPoints(2,0,sh,h)
    pg2f.setLineColor(Color.RED)
    pg2g.setLineColor(Color.RED)
    #p1.setVLabel(0,"f & gw")
    #p1.setVLabel(1,"g & gw")
  else: 
    p1 = PlotPanel(2,1,PlotPanel.Orientation.X1RIGHT_X2UP,PlotPanel.AxesPlacement.LEFT_TOP)
    pf = p1.addPoints(0,0,sf,f)
    pg = p1.addPoints(1,0,sg,g)
  p1.setHLabel(0,"Time (s)")
  #p1.setVLabel(0,"f")
  #p1.setVLabel(1,"g")
  if tlim:
    p1.setHLimits(0,tlim)
  else: p1.setHLimits(0,sf.getLast())
  frame = PlotFrame(p1)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(28)
  frame.setSize(1500,350)
  if paper:
    frame.setSize(1100,413)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    frame.setSize(1000,350)
    #frame.setSize(1000,500)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
    #if g2:
    #  frame.setSize(965,455)
    #else: frame.setSize(965,349)
    #frame.setFontSizeForSlide(1.0,1.0)
    #frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)

  


#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
