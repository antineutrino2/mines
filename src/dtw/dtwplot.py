# Dynamic Time Warping plotter for tests
# This runs all Dtw methods and plots
# @author Andrew Munoz, Colorado School of Mines
# @version 10.20.2011

from imports import *

from dtw import * 
from tp import TDutils

_pngDir = '../tp/tppics/dtwmodels/'

def main(args):
  goDtw1dif()
  #goDtw1lag()
  #goDtw2pseudo()
  #goDtw2Tree()
  #goDtw2ModTree


# DTW for 1D signals of two different lengths

def goDtw1dif():
  n1,n2 = 1001, 301     # n1 > n2 always... a
  g0 = 350
  tui=0; tli=n1-1;
  r = 4 
  x = sig1L(n1) 
  y = sig1S(n1,n2,250) 
  sx = Sampling(n1)
  sy = Sampling(n2,1,g0)
  cm = Dtw1().cm(x,y) 
  acm = Dtw1().acm(cm) 
  #cm = Dtw1n().cm(x,y) 
  for sl in [2]:
    #acm = Dtw1n().acm(cm,sl) 
    #pe,tpi,zpi = Dtw1n().paths(cm,sl)
    #pe,tpi,zpi = Dtw1().paths(cm)
    #pe1,tpi1,zpi1 = Dtw1().pathRemoval(pe,tpi,zpi,tui,tli)
    #pe1,tpi1,zpi1 = Dtw1n().pathRemoval(pe,tpi,zpi,tui,tli)
    #pe2,tpi2,zpi2 = Dtw1().pathRemoval(pe,tpi,zpi,acm[len(acm)-1],r)
    #pe3,tpi3,zpi3 = Dtw1n().pathRemoval(pe,tpi,zpi,acm[len(acm)-1],1)
    #pe3,tpi3,zpi3 = Dtw1().pathRemoval(pe,tpi,zpi,acm[len(acm)-1],1)

    #tm,tn,tz = TDutils().tdcurve(Sampling(1),sx,sy,tpi3[0],zpi3[0],zerofloat(1))
    #y2 = TDutils().stretch(tm,y,sy)
    #sy2 = Sampling(len(y2),sy.getDelta(),tm[0])
    
    #plotSubsq(sx,sy,x,y,st3=sy2,g2=y2,slides='sdtwmcs')#paper='pdtwmc2')#,slides='sdtwmcs')
    #plotSubsq(sx,sy,x,y,slides='sdtwmc')#paper='pdtwmc2')#,slides='sdtwmcs')
    #plotSubsq(sx,sy,x,y,paper='pdtwmc')#,slides='sdtwmc')
    #plotImgB     (sx,sy,acm,pe=pe                   )#,slides='sdtwmdb')#,paper='dtwmd')
    #plotImgB     (sx,sy,cm                         )#,slides='sdtwmdb')#,paper='dtwmd')
    #plotImgB     (sx,sy,acm,cn=len(pe),pt=tpi,pd=zpi)#,slides='sdtwmdbca')#,paper='dtwmd')
    #plotImgB     (sx,sy,acm,pe=pe,cn=len(pe1),pt=tpi1,pd=zpi1)#,slides='sdtwmdbct')#,paper='dtwmd')
    #plotImgB     (sx,sy,acm,pe=pe,cn=len(pe2),pt=tpi2,pd=zpi2)#,slides='sdtwmdbcr')#,paper='dtwmd')
    #plotImgB     (sx,sy,acm,pe=pe,cn=len(pe3),pt=tpi3,pd=zpi3)#,slides='sdtwmdbcm')#,paper='dtwmd')

    #plotImgCurves(sx,sy, cm,x,y,cn=len(pe3),pt=tpi3,pd=zpi3)#,slides='sdtwmecm')#,paper='dtwmec')
    #plotImgCurves(sx,sy,acm,x,y,cn=len(pe3),pt=tpi3,pd=zpi3)#,slides='sdtwmdcm')#,paper='dtwmdc')
    #plotImgCurves(sx,sy,acm,x,y,cn=len(pe2),pt=tpi2,pd=zpi2)#,slides='sdtwmdcr')#,paper='dtwmdc')
    #plotImgCurves(sx,sy,acm,x,y,cn=len(pe1),pt=tpi1,pd=zpi1)#,slides='sdtwmdct')#,paper='dtwmdc')
    #plotImgCurves(sx,sy,acm,x,y,cn=len(pe),pt=tpi,pd=zpi)#,slides='sdtwmdca')#,paper='dtwmdc')
    plotImgCurves(sx,sy,acm,x,y                         )#,slides='sdtwmd')#,paper='dtwmd')
    #plotImgCurves(sx,sy, cm,x,y                         )#,slides='sdtwme')#,paper='dtwme')
    
    #plotb(sx,acm,pe,paper='pdtwmbn')
    #plot2Dp(sx,sy, cm,cbwmin=70,slides='sdtwme2')#,paper='pdtwme')
    #plot2Dp(sx,sy,acm,cbintv=20,slides='sdtwmd2')#,paper='pdtwmd')
    #plot2Dp(sx,sy,acm,cn=len(pe),pt=tpi,pd=zpi,cbintv=20,paper='pdtwmdca')
    #plot2Dp(sx,sy,acm,cn=len(pe1),pt=tpi1,pd=zpi1,cbintv=20)#,paper='pdtwmdcs')
    #plot2Dp(sx,sy,acm,cn=len(pe2),pt=tpi2,pd=zpi2,cbintv=20,paper='pdtwmdcr')
    #plot2Dp(sx,sy,acm,cn=len(pe3),pt=tpi3,pd=zpi3,cbintv=20,slides='sdtwmdcm2')#,paper='pdtwmdcm')

"""
DTW for 1D signals of the same length and in lag space
Using Dave's modified dtw code
"""
def goDtw1lag():
  #TODO
  pass


# DTW for 2D using Tree-Serial Dynamic Programming

def goDtw2tree():
  #TODO
  pass


# DTW for 2D using pseudo2D dynamic programming

def goDtw2pseudo():
  n1, n2 = 2001, 301
  #x = readImage(n1,n2,'/Users/drumdude42/Home/data/migclass/layered1.dat')
  #y = readImage(n1,n2,'/Users/drumdude42/Home/data/migclass/layered2.dat')
  x = readImage(n1,n2,'/Users/amunoz/Home/data/migclass/layeredVert1.dat')
  #y = readImage(n1,n2,'/Users/amunoz/Home/data/migclass/layered2.dat')
  y = readImage(n1,n2,'/Users/amunoz/Home/data/migclass/layeredVert2.dat')
  path = Dtw2.dtw2PS(x,y)
  plot2D(x,"Image 1")
  plot2D(y,"Image 2")
  plot2D(path,"1D Warp paths")


##########################################################################################
## Thanks to Luming Liang for this signal format:  

def sig1L(n1):
  seed   = 23232
  seedN1 = 31415 
  fpeak  = 0.125
  rsmNoise = 0.1
  X = Dtw1().makeEvents(n1, seed) 
  X = Dtw1().addRickerWavelet(fpeak, X)
  X = mul(1.0/max(abs(X)),X)
  X = Dtw1().addNoise(rsmNoise, seedN1, X)
  return X

def sig1S(n1,n2,start):
  _si = SincInterpolator()
  _si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT)
  seedN2 = 10000
  rsmNoise = 0.1
  amp = 30
  Y = zerofloat(n2)
  t = zerofloat(n2)
  X = sig1L(n1);
  w = float((2*PI/n2))
  for i in range (0,n2):
    t[i]=start+i-amp*sin(i*w)   
  _si.setUniform(n1,1.0,0.0,X) 
  _si.interpolate(n2,t,Y)    
  Y = Dtw1().addNoise(rsmNoise, seedN2, Y)
  #Y = add(Y,50) # To test DDTW
  return Y

# Test signals:

def sig2S(n2):
  seed   = 20000
  seedN1 = 31415 
  #seed   = 23232
  #seedN1 = 31415 
  fpeak  = 0.125
  rsmNoise = 0.1
  Y = Dtw1().makeEvents(n2, seed) 
  Y = Dtw1().addRickerWavelet(fpeak, Y)
  Y = mul(1.0/max(abs(Y)),Y)
  Y = Dtw1().addNoise(rsmNoise, seedN1, Y)
  return Y

def sig2L(n1,n2):
  _si = SincInterpolator()
  _si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT)
  seedN2 = 10000
  rsmNoise = 0.1
  amp = 30
  X = zerofloat(n1)
  t = zerofloat(n1)
  w = float((2*PI/n1))
  Y = sig2S(n2);
  for i in range (0,n1):
    t[i]=i-amp*sin(i*w)
  _si.setUniform(n2,1.0,0.0,Y)
  _si.interpolate(n1,t,X)
  X = Dtw1().addNoise(rsmNoise, seedN2, X)
  return X

def sig3S(n2,X):
  Y = zerofloat(n2)
  for i in range (0,n2):
    Y[i] = X[i] + 1.0 
  return Y
  
###########################################################################

def plotImgCurves(st,sz,img,t,z,c1=None,cn=None,pt=None,pd=None,paper=None,slides=None):
  p1 = PlotPanel(2,2,PlotPanel.Orientation.X1RIGHT_X2UP,PlotPanel.AxesPlacement.LEFT_TOP)
  p1.mosaic.setWidthElastic(0,25)
  p1.mosaic.setHeightElastic(0,25)
  p1.mosaic.setWidthElastic(1,100)
  p1.mosaic.setHeightElastic(1,100)
  pp1 = p1.addPoints(0,1,st,t)
  pp1.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
  p1.setVLabel(0,"f")
  p1.setVLabel(1,"Tau (s)")
  p1.setHLabel(0,"g")
  p1.setHLabel(1,"Time (s)")
  pp2 = p1.addPoints(1,0,sz,z)
  pp2.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
  im = p1.addPixels(1,1,st,sz,img)
  im.setInterpolation(PixelsView.Interpolation.NEAREST)
  im.setColorModel(ColorMap.JET)
  im.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
  #p1.addColorBar()
  #p1.setColorBarWidthMinimum(100)
  # ADDED PATHS AND CURVES
  dz=sz.getDelta(); dt=st.getDelta();
  fz=sz.getFirst(); ft=st.getDelta();
  if pd and pt:
    cp = [0]
    if c1: 
      cp = [c1]
    if cn: 
      cp = range(cn)
    for i in cp:
      pt2 = mul(pt[i],dt)
      pd2 = mul(pd[i],dz)
      pd2 = add(pd2,fz)
      ln = p1.addPoints(1,1,pt2,pd2)
      ln.setLineColor(WHITE)
      ln.setLineWidth(3)
      ln.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
  p1.setHLimits(1,st.getFirst(),st.getLast())
  p1.setVLimits(1,sz.getFirst(),sz.getLast())
  p1.addColorBar()
  frame = PlotFrame(p1)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setSize(1500,1000)
  if paper:
    frame.setSize(1000,700)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    frame.setSize(1050,610)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)

def plotImgB(st,sz,img,pe=None,c1=None,cn=None,pt=None,pd=None,paper=None,slides=None):
  p1 = PlotPanel(2,1,PlotPanel.Orientation.X1RIGHT_X2UP,PlotPanel.AxesPlacement.LEFT_TOP)
  p1.mosaic.setHeightElastic(0,25)
  p1.mosaic.setHeightElastic(1,100)
  p1.setVLabel(1,"Tau (s)")
  p1.setHLabel(0,"Time (s)")
  p1.setVLabel(0,"D")
  b = add(img[len(img)-1],0)
  if pe:
    b2 = zerofloat(len(b))
    for i in range(len(b)):
      b2[i] = b[i]/len(pe[i])
    p1.setVLabel(0,"Dn")
    p1.addPoints(0,0,st,b2)
    p1.setVLimits(0,min(b2)*0.5,min(b2)*3.0)
    p1.setVFormat(0,"%.1f")
  else: 
    p1.addPoints(0,0,st,b)
    p1.setVInterval(0,40)
    p1.setVLimits(0,0,max(b))
    p1.setVFormat(0,"%1.0f")
  im = p1.addPixels(1,0,st,sz,img)
  im.setInterpolation(PixelsView.Interpolation.NEAREST)
  im.setColorModel(ColorMap.JET)
  im.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
  #p1.addColorBar()
  #p1.setColorBarWidthMinimum(100)
  # ADDED PATHS AND CURVES
  dz=sz.getDelta(); dt=st.getDelta();
  fz=sz.getFirst(); ft=st.getDelta();
  if pd and pt:
    cp = [0]
    if c1: 
      cp = [c1]
    if cn: 
      cp = range(cn)
    for i in cp:
      pt2 = mul(pt[i],dt)
      pd2 = mul(pd[i],dz)
      pd2 = add(pd2,fz)
      ln = p1.addPoints(1,0,pt2,pd2)
      ln.setLineColor(WHITE)
      ln.setLineWidth(3)
      ln.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
  p1.setHLimits(0,st.getFirst(),st.getLast())
  p1.setVLimits(1,sz.getFirst(),sz.getLast())
  p1.addColorBar()
  frame = PlotFrame(p1)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setSize(1500,1000)
  if paper:
    frame.setSize(1000,700)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    frame.setSize(900,560)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)


def plotb(sx,img,pe,paper=None,slides=None):
  n1 = len(img[0]) #long
  n2 = len(img)    #short
  
  p2 = PlotPanel(2,1,PlotPanel.Orientation.X1RIGHT_X2UP,PlotPanel.AxesPlacement.LEFT_TOP)
  #img2 = transpose(img)
  ba = add(img[n2-1],0.)
  b2 = zerofloat(len(ba))
  for i in range(len(ba)):
    b2[i] = ba[i]/len(pe[i])
  p2.addPoints(0,0,sx,b2)
  p2.setVLimits(0,min(b2)*0.5,min(b2)*3.0)
  p2.setVFormat(0,"%.1f")
  p2.mosaic.setHeightElastic(0,25)
  p2.mosaic.setWidthElastic(0,100)
  p2.mosaic.setHeightElastic(1,100)
  #p2.setVLabel(1,"f(x)")
  #p2.setHLabel(1,"g(x)")
  p2.addPixels(1,0,img)
  p2.addColorBar()
  frame = PlotFrame(p2)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setSize(1500,300)
  if paper:
    frame.setSize(772,500)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    frame.setSize(1500,1000)
    frame.setFontSizeForSlide(8.0,504.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)

def plot(x,y,img,title):
  p1 = PlotPanel() 
  ln1 = p1.addPoints(x)
  ln2 = p1.addPoints(y)
  ln1.setLineColor(BLACK)
  ln2.setLineColor(BLUE)
  p1.setHLabel("Samples")
  p1.setVLabel("Amplitude")
  p1.setTitle(title)
  
  p2 = PlotPanel(PlotPanel.Orientation.X1RIGHT_X2UP)
  pv = p2.addPixels(img)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.JET)
  #p2.addColorBar()
  #p2.setColorBarWidthMinimum(100)
  p2.setHLabel("Long Signal")
  p2.setVLabel("Short Signal")
  frame = PlotFrame(p1,p2,PlotFrame.Split.VERTICAL)
  frame.setSize(1500,1000)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True)

def plotSubsq(st,st2,f,g,st3=None,g2=None,paper=None,slides=None):
  if g2:
    p1 = PlotPanel(3,1,PlotPanel.Orientation.X1RIGHT_X2UP,PlotPanel.AxesPlacement.LEFT_TOP)
    pf = p1.addPoints(0,0,st,f)
    pg = p1.addPoints(1,0,st2,g)
    pg2f = p1.addPoints(0,0,st3,g2)
    pg2g = p1.addPoints(2,0,st3,g2)
    pg2f.setLineColor(RED)
    pg2g.setLineColor(RED)
    #p1.setVLabel(0,"f & gw")
    #p1.setVLabel(1,"g & gw")
  else: 
    p1 = PlotPanel(2,1,PlotPanel.Orientation.X1RIGHT_X2UP,PlotPanel.AxesPlacement.LEFT_TOP)
    pf = p1.addPoints(0,0,st,f)
    pg = p1.addPoints(1,0,st2,g)
  p1.setHLabel(0,"Time (s)")
  #p1.setVLabel(0,"f")
  #p1.setVLabel(1,"g")
  p1.setHLimits(0,st.getLast())
  frame = PlotFrame(p1)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(28)
  frame.setSize(1500,350)
  if paper:
    frame.setSize(1100,413)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    if g2:
      frame.setSize(965,455)
    else: frame.setSize(965,349)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)


def plot2Dp(st,sz,img,title=None,pt=None,pd=None,cn=None,c1=None,cbintv=None,cbm=None,cbwmin=None,paper=None,slides=None):
  # For v1 and v2
  dt = st.getDelta(); dz = sz.getDelta();
  ft = st.getFirst(); fz = sz.getFirst();
  lt = st.getLast(); lz = sz.getLast();
  pp = PlotPanel()
  pv = pp.addPixels(st,sz,img)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.JET)
  if pd and pt:
    cp = [0]
    if c1: 
      cp = [c1]
    if cn: 
      cp = range(cn)
    for i in cp:
        pt2 = mul(pt[i],dt)
        pd2 = mul(pd[i],dz)
        pd2 = add(pd2,fz)
        ln = pp.addPoints(pt2,pd2)
        ln.setLineColor(WHITE)
        ln.setLineWidth(2)
  #pp.setHLabel("Time (s)")
  #pp.setVLabel("Tau (s)")
  if title: pp.setTitle(title)
  pp.setLimits(ft,fz,lt,lz)
  cb = pp.addColorBar()#"Error")
  if cbintv:
    cb.setInterval(cbintv)
  if cbwmin: pp.setColorBarWidthMinimum(cbwmin)
  frame = PlotFrame(pp)
  frame.setSize(1169,344) # Perfect value to frame ratio
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  if paper:
    frame.setSize(772,253)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    frame.setSize(700,350)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)


def plotbn(sx,b,pe,paper=None):
  pp = PlotPanel()
  b2 = zerofloat(len(b))
  for i in range(len(b)):
    b2[i] = b[i]/len(pe[i])
  pp.addPoints(0,0,sx,b2)
  pp.setVLimits(0,min(b2)*0.5,min(b2)*3.0)
  pp.setVFormat(0,"%.1f")
  frame = PlotFrame(pp)
  frame.setSize(1169,344) # Perfect value to frame ratio
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  if paper:
    frame.setSize(772,253)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  frame.setVisible(True)


def plot1D(x,title):
  sp = SimplePlot()
  pv = PointsView(x)
  sp.setTitle(title)
  sp.setSize(1500,500)
  sp.add(pv)

def plot2D(img,title):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  im1 = sp.addPixels(img)
  im1.setPercentiles(2,98)
  sp.setTitle(title)
  sp.setSize(1200,500)
  sp.addColorBar()
  

def readImage(n1,n2,filename):
  ais = ArrayInputStream(filename)
  x = zerofloat(n1,n2)
  ais.readFloats(x)
  ais.close()
  return x

def double(x):
  y = zerodouble(len(x))
  for i in range(len(x)):
    y[i] = x[i]
  return y


#---------------------------------------------------------#

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

 


