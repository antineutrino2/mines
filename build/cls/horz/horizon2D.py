from imports import *
from wt import WTutils

dataDir = "/Users/amunoz/Home/data/tp/wellties/"

def main(args):
  goHorizon2d()


def goHorizon2d():
  sx,sy,st,cube = getCube()
  nx = sx.getCount()
  ny = sy.getCount()
  dt = st.getDelta()
  nym = ny-1
  xc = int(nx/2)
  ght,ghy = getHorizons(xc,sx,st)
  slc = cube[xc]
  thy = zerofloat(ny)
  thti = zeroint(ny)
  thy = rampfloat(sy.getFirst(),sy.getDelta(),ny)
  
  ml = 20 
  mlp = ml+1
  nl = ml*2+1
  stretchMax = 0.5
  dw = DynamicWarpingF(-ml,ml)
  dw.setStrainMax(stretchMax)

  k=0
  lb = 3 # number of lookback traces
  p = zeroint(lb)
  em = zerofloat(lb)
  thti[0] = int(0.877/dt)#int(ght[0][0]/dt) #best is .877
  for i in range(lb):
    pc = thti[i]
    e = dw.computeErrors(slc[i],slc[i+1])
    d = dw.accumulateForward(e)
    u = dw.findShifts2(e,d)
    thti[i+1] = pc + int(u[pc])
  for i in range(lb,nym):
    for j in range(lb):
      imj = i-j
      pc = thti[imj]
      e = dw.computeErrors(slc[imj],slc[i+1])
      d = dw.accumulateForward(e)
      u = dw.findShifts2(e,d)
      p[j] = pc+int(u[pc])
      em[j] = e[i-j][p[j]-pc]
    emn = em[0]
    for j in range(lb):
      if (em[j]<emn):
        emn = e[j]
	pm = p[j]
    thti[i+1] = pm
    #if (i==k):      
    #  plot2DHorz(slc,st,sy,ght,ghy,mult(thti,dt),thy)
    #  k+=30

  tht = mult(thti,dt)
  plot2DHorz(slc,st,sy,ght,ghy,tht,thy)


###############################################################################################
## utils

"""
Extract slice
"""
def getCube():
  tsdata,st,sy,sx = readTpstData()
  return sx,sy,st,tsdata

"""
Gets tp depth and time horizons in 2D form to plot on a seismic 2D slice
"""
def getHorizons(xc,sx,st):
  tdu = WTutils()
  d = "t"
  DKOT = readHorizonMod("FallRiverDKOT",d)
  CRMT = readHorizonMod("CrowMountainCRMT",d)
  TSAS = readHorizonMod("TensleepASand",d)
  TSBD = readHorizonMod("TensleepBbaseC1Dolo",d)
  tdk,ydk = tdu.horizon2d(DKOT.x1,DKOT.x2,DKOT.x3,sx,st,xc)
  tcr,ycr = tdu.horizon2d(CRMT.x1,CRMT.x2,CRMT.x3,sx,st,xc)
  tta,yta = tdu.horizon2d(TSAS.x1,TSAS.x2,TSAS.x3,sx,st,xc)
  ttb,ytb = tdu.horizon2d(TSBD.x1,TSBD.x2,TSBD.x3,sx,st,xc)
  return [tdk,tcr,tta,ttb],[ydk,ycr,yta,ytb]

"""
Reads in seismic traces
"""
def readTpstData():
  seismicDir = getSeismictDir()
  s1,s2,s3 = getTpstSamplings()
  fileName = seismicDir+"subt_3002_1_0/tpst.dat"
  n1,n2,n3 = s1.count,s2.count,s3.count
  tsdata = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(fileName)
  ais.readFloats(tsdata)
  ais.close()
  return tsdata,s1,s2,s3

def normalize(e):
  emin = min(e)
  emax = max(e)
  return mul(sub(e,emin),1.0/(emax-emin))

def etran(e):
  e = normalize(e)
  return transpose(pow(e,0.25))
  
def mult(x,a):
  y = zerofloat(len(x))
  for i in range(len(x)):
    y[i] = x[i]*a
  return y


###############################################################################################
## plots


"""
Plots a 2-D seismic slice, a well log through that slice, and seismic horizons on the slice
"""
def plot2DHorz(img,st,sy,ghz,ghy,tht,thy,paper=None,slides=None):
  ny = len(img)
  sp = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT) 
  im1 = sp.addPixels(st,sy,img)
  im1.setColorModel(ColorMap.GRAY)#RED_WHITE_BLUE)#GRAY)
  #im1.setInterpolation(PixelsView.Interpolation.NEAREST)
  # Add Horizons
  #pv1 = sp.addPoints(ghz[0],ghy[0])
  #pv1.setStyle("g-")
  #pv1.setLineWidth(2)
  #pv2 = sp.addPoints(ghz[1],ghy[1])
  #pv2.setStyle("c-")
  #pv2.setLineWidth(2)
  #pv3 = sp.addPoints(ghz[2],ghy[2])
  #pv3.setStyle("y-")
  #pv3.setLineWidth(2)
  #pv4 = sp.addPoints(ghz[3],ghy[3])
  #pv4.setStyle("m-")
  #pv4.setLineWidth(2)
  pv5 = sp.addPoints(tht,thy)
  pv5.setStyle("b-")
  pv5.setLineWidth(2)
  sp.setVLimits(0,1.4)
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








#---------------------------------------------------------#

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())












