"""
Converts sgy image files (SEG-Y format) to dat files (3D arrays of floats)
Removes all SEG-Y headers and, if necessary, converts the data format to
IEEE floats. The byte order for floats in dat files is BIG_ENDIAN. 

Author: Dave Hale, Colorado School of Mines
Edited by: Andrew Munoz, CSM, (put in BP dataset params)
Version: 2014.01.31
"""
from imports import *
global n1,n2,no

basedir = "/Users/amunoz/Home/data/bp/"

#############################################################################
def main(args):
  #goModels()
  #goIsoStack()
  goGathers()
  #goVSP() #TODO

def goModels():
  goModel("vp")
  goModel("ep")
  goModel("de")
  goModel("th")

def goGathers():
  x1 = goGats("1",makeCmp=True)
  #x2 = goGats("2",makeCmp=True)
  #x3 = goGats("3",makeCmp=True)
  #x4 = goGats("4",makeCmp=True)
  fname = basedir+"dat/ani/cmps1.dat"
  makeCmps(fname,x1)


def goIsoStack():
  """
  ****** beginning of SEG-Y file info ******
  file name = /Users/amunoz/Home/data/bp/sgy/iso/FD_Model_PSTM_STK.sgy
  byte order = BIG_ENDIAN
  number of bytes = 61018624
  number of traces = 12596
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: m (will be converted to km)
  indices and coordinates from trace headers:
    i2min =   803, i2max = 13398 (inline indices)
    i3min =     0, i3max =     0 (crossline indices)
    xmin =    0.000000, xmax =    0.787187 (x coordinates, in km)
    ymin =    0.000000, ymax =    0.000000 (y coordinates, in km)
  grid sampling:
    n1 =  1151 (number of samples per trace)
    n2 = 12596 (number of traces in inline direction)
    n3 =     1 (number of traces in crossline direction)
    d1 = 0.008000 (time sampling interval, in s)
    d2 = 0.000063 (inline sampling interval, in km)
    d3 = 0.000000 (crossline sampling interval, in km)
  grid corner points:
    i2min =   803, i3min =     0, x =    0.000000, y =    0.000000
    i2max = 13398, i3min =     0, x =    0.787188, y =    0.000000
    i2min =   803, i3max =     0, x =    0.000000, y =    0.000000
    i2max = 13398, i3max =     0, x =    0.787188, y =    0.000000
  grid azimuth: 90.00 degrees
  ****** end of SEG-Y file info ******
  """
  firstLook = False # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = True # reads all traces, writes an image
  showImage = True # displays the image
  sgyfile = basedir+"sgy/iso/FD_Model_PSTM_STK.sgy"
  datfile = basedir+"dat/iso/bp_pstm_iso_stk.dat"
  i1min,i1max = 0,1150
  i2min,i2max = 803,13398
  n1,n2 = 1+i1max-i1min,1+i2max-i2min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
    plotIbmIeeeFloats(si)
  if secondLook:
    si.printAllInfo()
    plot23(si)
    plotXY(si)
  if writeImage:
    scale = 1.0
    si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,0,0)
  si.close()
  if showImage:
    f = readImage(datfile,n1,n2)
    show2d(f,cbar="Amplitude")

def goModel(which):
  """
  ****** beginning of SEG-Y file info ******
  file name = /Users/amunoz/Home/data/bp/sgy/model/Vp_Model.sgy
  file name = /Users/amunoz/Home/data/bp/sgy/model/Epsilon_Model.sgy
  file name = /Users/amunoz/Home/data/bp/sgy/model/Delta_Model.sgy
  file name = /Users/amunoz/Home/data/bp/sgy/model/Theta_Model.sgy
  byte order = BIG_ENDIAN
  number of bytes = 93768224
  number of traces = 12596
  format = 1 (4-byte IBM floating point)
  units for spatial coordinates: m (will be converted to km)
  indices and coordinates from trace headers:
    i2min =   803, i2max = 13398 (inline indices)
    i3min =     0, i3max =     0 (crossline indices)
    xmin =    0.000000, xmax =    0.000000 (x coordinates, in km)
    ymin =    0.000000, ymax =    0.000000 (y coordinates, in km)
  grid sampling:
    n1 =  1801 (number of samples per trace)
    n2 = 12596 (number of traces in inline direction)
    n3 =     1 (number of traces in crossline direction)
    d1 = 0.006250 (time sampling interval, in s)
    d2 = 0.000000 (inline sampling interval, in km)
    d3 = 0.000000 (crossline sampling interval, in km)
  grid corner points:
    i2min =   803, i3min =     0, x =    0.000000, y =    0.000000
    i2max = 13398, i3min =     0, x =    0.000000, y =    0.000000
    i2min =   803, i3max =     0, x =    0.000000, y =    0.000000
    i2max = 13398, i3max =     0, x =    0.000000, y =    0.000000
  grid azimuth: 90.00 degrees
  ****** end of SEG-Y file info ******
  """
  firstLook = False # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = False # reads all traces, writes an image
  showImage = True # displays the image
  if which=="vp":
    sgyfile = basedir+"sgy/model/Vp_Model.sgy"
    datfile = basedir+"dat/model/vp.dat"
  elif which=="ep":
    sgyfile = basedir+"sgy/model/Epsilon_Model.sgy"
    datfile = basedir+"dat/model/epsilon.dat"
  elif which=="de":
    sgyfile = basedir+"sgy/model/Delta_Model.sgy"
    datfile = basedir+"dat/model/delta.dat"
  elif which=="th":
    sgyfile = basedir+"sgy/model/Theta_Model.sgy"
    datfile = basedir+"dat/model/theta.dat"
  i1min,i1max = 0,1800
  i2min,i2max = 803,13398
  n1,n2 = 1+i1max-i1min,1+i2max-i2min
  si = SegyImage(sgyfile)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
    plotIbmIeeeFloats(si)
  if secondLook:
    si.printAllInfo()
    #plot23(si)
    #plotXY(si)
  if writeImage:
    scale = 1.0
    si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,0,0)
  si.close()
  if showImage:
    f = readImage(datfile,n1,n2)
    show2d(f,cmap=jet,cbar=which)

def goGats(which,makeCmp=False):
  """
  %%%%%%%%%% From data info file %%%%%%%%%%
  1641 Shots @ 50m SP interval
  Source depth = 6m
  First SP = 1:  X = 0m; Y = 0m: = CDP No. 803
  Azimuth 90 deg
  
  Number of channels = 800
  Channel Separation = 12.5m
  Near offset = 37.5m
  Far offset = 10025m
  Receiver Depth = 8m
  
  Sample rate = 8ms
  Record length = 9200ms
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  """
  firstLook = False # fast, does not read all trace headers
  secondLook = False # slow, must read all trace headers
  writeImage = False # reads all traces, writes an image
  showImage = False # displays the image
  i1min,i1max = 0,1150
  i2min,i2max = 1,800
  if which=="1":
    sgyfile = basedir+"sgy/ani/Anisotropic_FD_Model_Shots_part1.sgy"
    datfile = basedir+"dat/ani/shots1.dat"
    itmin,itmax = 1,500
  elif which=="2":
    sgyfile = basedir+"sgy/ani/Anisotropic_FD_Model_Shots_part2.sgy"
    datfile = basedir+"dat/ani/shots2.dat"
    itmin,itmax = 501,1000
  elif which=="3":
    sgyfile = basedir+"sgy/ani/Anisotropic_FD_Model_Shots_part3.sgy"
    datfile = basedir+"dat/ani/shots3.dat"
    itmin,itmax = 1001,1500
  elif which=="4":
    sgyfile = basedir+"sgy/ani/Anisotropic_FD_Model_Shots_part4.sgy"
    datfile = basedir+"dat/ani/shots4.dat"
    itmin,itmax = 1501,1641
  i3min,i3max = itmin,itmax
  # n3 is number of shots
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  si = SegyImage(sgyfile)
  si.setInlineXlineBytes(9,13)
  if firstLook:
    si.printSummaryInfo();
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
    si.printTraceHeader(3)
    si.printTraceHeader(4)
    si.printTraceHeader(901)
    si.printTraceHeader(902)
    #plotIbmIeeeFloats(si)
  if secondLook:
    si.printAllInfo()
    #plot23(si)
    #plotXY(si)
  if writeImage:
    scale = 1.0
    si.writeFloats(datfile,scale,i1min,i1max,i2min,i2max,i3min,i3max,True)
  si.close()
  if showImage:
    x = readImage(datfile,n1,n2,n3)
    show3d(x,clip=1.0)
    show2d(x[100],clip=1.0)
  if makeCmp:
    x = readImage(datfile,n1,n2,n3)
    return x


def makeCmps(fname,x1,x2=None,x3=None,x4=None):
  ds = 50
  ns = 1641
  fr = 37.5
  dr = 12.5
  lr = 10025
  nr = 800
  ss = Sampling(ns,ds,lr)
  sr = Sampling(nr,dr,fr)
  ns1 = 501 
  ss1 = Sampling(ns1,ds,lr)
  SegyImage.cmpSort(fname,ss1,sr,x1);
  nm = int((sr.last-sr1.first+sr.delta)/(ds*2.0));
  ng = ns1*nm
  nt = 1151
  x = readImage(datfile,nt,nm,ng)
  show3d(x,clip=1.0)

  



#############################################################################

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE

def show2d(f,clip=None,title=None,cmap=gray,cbar=None):
  print "show2d: f min =",min(f)," max =",max(f)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  pv = sp.addPixels(f)
  pv.setColorModel(cmap)
  if clip:
    pv.setClips(-clip,clip)
  else:
    pv.setPercentiles(2,98)
  if title:
    sp.setTitle(title)
  if cbar:
    sp.addColorBar(cbar)
  sp.setSize(1100,600)
  
def show3d(f,clip=None):
  print "show3d: f min =",min(f)," max =",max(f)
  frame = SimpleFrame()
  ipg = frame.addImagePanels(f)
  if clip:
    ipg.setClips(-clip,clip)
  frame.orbitView.setScale(2.0)
  frame.setSize(1000,1000)

def plot23(si):
  i2 = si.getI2sAsFloats()
  i3 = si.getI3sAsFloats()
  sp = SimplePlot()
  sp.setHLabel("inline sample index i2")
  sp.setVLabel("crossline sample index i3")
  pv = sp.addPoints(i2,i3)
  pv.setMarkStyle(PointsView.Mark.POINT);
  pv.setLineStyle(PointsView.Line.NONE);
  w,h = goodWidthHeight(i2,i3)
  sp.setSize(w,h)

def plotXY(si):
  x = si.getXs()
  y = si.getYs()
  sp = SimplePlot()
  sp.setHLabel("x coordinate (km)")
  sp.setVLabel("y coordinate (km)")
  pv = sp.addPoints(x,y)
  pv.setMarkStyle(PointsView.Mark.POINT);
  pv.setLineStyle(PointsView.Line.NONE);
  w,h = goodWidthHeight(x,y)
  sp.setSize(w,h)

def goodWidthHeight(x,y):
  xmin,xmax = min(x),max(x)
  ymin,ymax = min(y),max(y)
  w,h = 1000,1000
  if (xmax-xmin)>(ymax-ymin):
    h = int(h*(ymax-ymin)/(xmax-xmin))
  else:
    w = int(w*(xmax-xmin)/(ymax-ymin))
  return w,h

def readImage(datfile,n1,n2,n3=1):
  if n3==1:
    x = zerofloat(n1,n2)
  else:
    x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(datfile)
  ais.readFloats(x)
  ais.close()
  return x

def writeImage(datfile,x):
  aos = ArrayOutputStream(datfile)
  aos.writeFloats(x)
  aos.close()

def plotIbmIeeeFloats(si):
  ntrace = si.countTraces()
  itrace = ntrace/2
  fmt = si.getFormat()
  si.setFormat(1) # IBM floats
  fibm = si.getTrace(itrace)
  si.setFormat(5) # IEEE floats
  fieee = si.getTrace(itrace)
  si.setFormat(fmt)
  pp = PlotPanel(2,1)
  pp.setTitle("IBM (top) versus IEEE (bottom)")
  pp.setHLabel(0,"Sample index")
  pp.setVLabel(0,"IBM amplitude")
  pp.setVLabel(1,"IEEE amplitude")
  pp.addPoints(0,0,fibm)
  pp.addPoints(1,0,fieee)
  pf = PlotFrame(pp)
  pf.setSize(1000,800)
  pf.setVisible(True)

def lowpass(f3db,f):
  bf = ButterworthFilter(f3db,6,ButterworthFilter.Type.LOW_PASS)
  bf.apply1ForwardReverse(f,f)

def gain(hw,f):
  g = mul(f,f)
  RecursiveExponentialFilter(hw).apply1(g,g)
  div(f,sqrt(g),f)

def stretch(c,f):
  n1,n2 = len(f[0]),len(f)
  t = rampfloat(0.0,1.0/c,n1)
  si = SincInterp()
  g = zerofloat(n1)
  for i2 in range(n2):
    si.interpolate(n1,1.0,0.0,f[i2],n1,1.0/c,0.0,g)
    copy(g,f[i2])

def spow(p,f):
  return mul(sgn(f),pow(abs(f),p))

def slog(f):
  return mul(sgn(f),log(add(1.0,abs(f))))

def sexp(f):
  return mul(sgn(f),sub(exp(abs(f)),1.0))



#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain()) 
