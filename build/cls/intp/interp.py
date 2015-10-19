"""
This is an example jython script for spatially interpolating segy volumes. The 
interpolation schemes are basic and include: linear, monotonic, and 
cubic spline. This retains the same orientation, corners, origin point, and 
starting inline and crossline numbers. To interpolate the segy, you must choose
a maximum spacing for X and Y.

When first using a SEGY, run firstlook to print the details, then copy them 
above your information for reference later. This will help you substantially as
some future interpolation inputs are needed from this method. 

THIS IS NOT FOR SEISMIC. This is simply a trilinear/tricubic interpolation for 
SEGY volumes (ie, use for smooth velocity models, ect.). 

Author: Andrew Munoz, Newfield Exploration and Colorado School of Mines
Version: 7/9/2013
"""

"""
*****************************************************************
SEGY setups (just copy, paste, and replace to make things easier)
*****************************************************************
"""

######## for testing only, remove later:
from segy import *
###########################################
from imports import *

def main(args):
  goSGYDemo()
  #insert new method here or just replace current one
  
def goSGYDemo(): 
  """
  SEGY Demo
  """
  global sgyInFile,sgyOutFile
  #sgyInFile  = "../data/sgyDemoIn.sgy"
  #sgyOutFile = "../data/sgyDemoOut.sgy"
  sgyInFile  = "E:/data/demo/sgyDemoIn.sgy"
  sgyOutFile = "E:/data/demo/sgyDemoOut.sgy"

  firstlook = True
  interp = True
  writesgy = True
  display = True

  # Set Byte Locations for SEGY
  inlineByte = 17
  xlineByte  = 25
  xByte      = 73
  yByte      = 77
  si = SegyImage(sgyInFile)
  si.setInlineXlineBytes(inlineByte,xlineByte)
  si.setXYBytes(xByte,yByte)

  """
****** beginning of SEG-Y file info ******
file name = E:/data/demo/sgyDemoIn.sgy
byte order = BIG_ENDIAN
number of bytes = 4577688
number of traces = 10302
format = 1 (4-byte IBM floating point)
WARNING  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
WARNING: format may actually be 5 (IEEE float)
WARNING  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
indices and coordinates from trace headers:
  i2min =   500, i2max =   600 (crossline indices)
  i3min =   500, i3max =   601 (inline indices)
  xmin = 1896995.000000, xmax = 1909098.000000 (x coordinates)
  ymin = 14577969.000000, ymax = 14586316.000000 (y coordinates)
grid sampling:
  n1 =    51 (number of samples per trace)
  n2 =   101 (number of traces in crossline direction)
  n3 =   102 (number of traces in inline direction)
  d1 = 0.010000 (time sampling interval, in s)
  d2 = 82.482803 (crossline sampling interval)
  d3 = 119.172348 (inline sampling interval)
grid corner points:
  i2min =   500, i3min =   500, x = 1896995.000000, y = 14586217.000000
  i2max =   600, i3min =   500, x = 1897063.000000, y = 14577969.000000
  i2min =   500, i3max =   601, x = 1909031.000000, y = 14586316.000000
  i2max =   600, i3max =   601, x = 1909099.000000, y = 14578068.000000
grid azimuth: 179.53 degrees
****** end of SEG-Y file info ******
  """
  # Parameters obtained from firstlook
  global i1min,i2min,i3min
  i1min,i1max = 0,50    # Time samples
  i2min,i2max = 500,600 # Crossline samples
  i3min,i3max = 500,601 # Inline
  global n1,n2,n3
  n1,n2,n3 = 1+i1max-i1min,1+i2max-i2min,1+i3max-i3min
  
  # Interpolations parameters
  interpType = 'linear'
  #interpType = 'monotonic'
  #interpType = 'spline'
    
  # Manually set to desired maximum x/y sampling output 
  d2Out = 82.48
  d3Out = 82.48
    
  if firstlook:
    si.printAllInfo()
    #si.printSummaryInfo()
    si.printBinaryHeader()
    si.printTraceHeader(0)
    si.printTraceHeader(1)
  if interp:
    g = goInterp(si,d2Out,d3Out,interpType)
    if writesgy:
      goWrite(si,d2Out,d3Out,g,sgyOutFile)
  if display:
    if interp:
      goDisplay(g)
    else:
      goDisplay(si.getImage())

def goInterp(si,d2o,d3o,inty):
  print 'Interpolating...'
  d1 = si.getD1(); d2 = si.getD2(); d3 = si.getD3();
  f1 = 0.0; f2 = 0.0; f3 = 0.0;
  n2o = int(ceil(n2*d2/d2o))
  n3o = int(ceil(n3*d3/d3o))
  print 'old n2= ',n2 ,' old n3= ',n3 
  print 'new n2= ',n2o,' new n3= ',n3o
  s1o = Sampling(n1,d1,f1)
  s2o = Sampling(n2o,d2o,f2)
  s3o = Sampling(n3o,d3o,f3)
  x1 = rampfloat(f1,d1,n1)
  x2 = rampfloat(f2,d2,n2)
  x3 = rampfloat(f3,d3,n3)
  f = si.getImage()
  if inty=='linear':
    li = TrilinearInterpolator3(x1,x2,x3,f)
  elif inty=='spline':
    cim = TricubicInterpolator3.Method.SPLINE
    li = TricubicInterpolator3(cim,cim,cim,x1,x2,x3,f)
  elif inty=='monotonic':
    cim = TricubicInterpolator3.Method.MONOTONIC
    li = TricubicInterpolator3(cim,cim,cim,x1,x2,x3,f)
  g = li.interpolate(s1o,s2o,s3o)
  return g

def goWrite(si,d2o,d3o,g,file):
  print 'Writing SEG-Y: ',sgyOutFile,'...'
  n3w = len(g); n2w = len(g[0]); n1w = len(g[0][0])
  nhead=3200 # number of bytes in EBCDIC header
  nbhed=400 # number of bytes in binary header
  nthed=240 # number of bytes in trace header
  File(sgyOutFile).delete()
  ais = ArrayFile(sgyInFile,"rw")
  aos = ArrayFile(sgyOutFile,"rw")
  ntth = nthed/4
  yi = zeroint(n1w)
  hi = zeroint(ntth)
  c = si.getCoordInfo()
  cr = si.getCorners()
  # Copy Binary header and EBCIDIC header
  for b in range(nhead):
    aos.writeByte(ais.readByte())
  for b in range(12):
    aos.writeShort(ais.readShort())
  # Change data sample format code to 32-bit IEEE floating-point
  aos.writeShort(5); ais.readShort()
  for b in range(13,200):
    aos.writeShort(ais.readShort())
  # Get data from first trace header
  ais.readInts(hi)
  for i3 in range(n3w):
    for i2 in range(n2w):
      # Write trace header
      inl = i3min+i3
      xln = i2min+i2
      cdp = hi[0]+1
      x = si.floatToIeee(float(getNewX(cr,c,i2,i3,d2o,d3o)))
      y = si.floatToIeee(float(getNewY(cr,c,i2,i3,d2o,d3o)))
      aos.writeInt(cdp)  
      aos.writeInt(xln) 
      aos.writeInt(hi[2])
      aos.writeInt(xln)
      aos.writeInt(inl)
      aos.writeInt(cdp)
      aos.writeInt(xln)
      aos.writeInt(hi[7])
      aos.writeBytes(zerobyte(36))
      aos.writeInt(hi[17])
      aos.writeInt(x)
      aos.writeInt(y)
      aos.writeInt(x)
      aos.writeInt(y)
      aos.writeInt(hi[22])
      aos.writeBytes(zerobyte(20))
      aos.writeInt(hi[28])
      aos.writeInt(hi[29])
      aos.writeBytes(zerobyte(120))
      # Write trace values
      si.floatToIeee(yi,g[i3][i2])
      aos.writeInts(yi)
  ais.close()
  aos.close()
  test = True
  if test:
    print '\n testing...'
    # Set Byte Locations for SEGY
    inlineByte = 17
    xlineByte  = 25
    xByte      = 73
    yByte      = 77
    sio = SegyImage(sgyOutFile)
    sio.setInlineXlineBytes(inlineByte,xlineByte)
    sio.setXYBytes(xByte,yByte)
    sio.printAllInfo()
    #sio.printSummaryInfo()
    #sio.printBinaryHeader()
    sio.printTraceHeader(0)
    sio.printTraceHeader(1)

def getNewX(cr,c,i2,i3,d2,d3):
  return cr[0]+i3*d3*c[9]-i2*d2*c[8]
def getNewY(cr,c,i2,i3,d2,d3):
  return cr[1]+i3*d2*c[6]+i2*d3*c[7]
  
gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def goDisplay(f,clip=None,cm=jet):
  frame = SimpleFrame()
  ipg = frame.addImagePanels(f)
  if clip:
    ipg.setClips(-clip,clip)
  ipg.setColorModel(cm)
  frame.orbitView.setScale(2.0)
  frame.setSize(1000,1000)


#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())