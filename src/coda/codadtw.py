# Dynamic Time Warping plotter for tests
# This runs all Dtw methods and plots
# @author Andrew Munoz, Colorado School of Mines
# @version 10.20.2011

from imports import *

from wt import WTutils

_pngDir = './pics/'
wtu = WTutils()

def main(args):
	#dname = "synth"
	dname = "real"
	#dname = "filtreal"
	if dname=="synth":
		sf,f = readSig('syn1_131072_2_5',fd=0.001,dd=0.0001)
		sg,g = readSig('syn2_131072_2_5',fd=0.001,dd=0.0001)
		sf,f = decimate(sf,f,1)
		sg,g = decimate(sg,g,1)
	if dname=="real":
		#sf,f = readSig('dat1_17500_4_30',window=3.0) # window the max value in sec
		#sg,g = readSig('dat2_16500_4_53',window=3.8) # window the max value in sec
		sf,f = readSig('dat1_20001_1_-2',cut=6.0) 
		sg,g = readSig('dat2_20001_1_-2',cut=6.0)
	if dname=="filtreal":
		sf,f = readSig('dat3_20001_1_-2',cut=6.0) 
		sg,g = readSig('dat4_20001_1_-2',cut=6.0)
	f = norm(f)
	g = norm(g)
	n1 = len(f) 
	n2 = len(g)
	#goHistogram(f,'f',png='histf')
	#goHistogram(g,'g',png='histg')
	print 'nf=',n1
	print 'ng=',n2
	#snr = computeSNR(0,100,f)
	#plot1D(sg,snr,"snr")
	goWarp(sf,sg,f,g,dname)

def goWarp(sf,sg,f,g,dname,snr=None):
  ng = sg.count
  nf = sf.count
  dt = sf.delta
  print sf.first
  print sg.first
  print sf.last
  print sg.last
  smin = -1.0
  smax =  1.0
  if dname=="synth":
    smin = -0.10
    smax =  0.10
    # Synthetic Params
    dr = 0.0002
    rmin = -0.000
    rmax =  0.006
    ngg = 30
  if dname=="real":
    # Real (unfilt) rough params
    #dr = 0.05
    #rmin = -0.20
    #rmax =  0.20
    #ngg = 31
    # Real (unfilt) smooth params
    smin = -0.04
    smax =  0.04
    dr = 0.004
    rmin = -0.10
    rmax =  0.10
    ngg = 15
    #sg = Sampling(sg.count,sg.delta,sf.first)
  if dname=="filtreal":
    # Real (filt) rough params
    #dr = 0.05
    #rmin = -0.05
    #rmax =  0.05
    #ngg = 125
    # Real (filt) smooth params
    smin = -0.04
    smax =  0.04
    dr = 0.004
    rmin = -0.10
    rmax =  0.10
    ngg = 15
  print 'kmin=',int(floor(rmin/dr))
  print 'kmax=',int(ceil(rmax/dr))
  print '1.0/dr=',1.0/dr
  lc = int(-sf.first/dt)
  print 'lc=',lc
  dw = DynamicWarpingCO(inro(smin/dt),inro(smax/dt),rmin,rmax,dr)
  dw.setInterpolation(DynamicWarpingCO.Interpolation.SPLINE)
  e = dw.computeErrors(f,g)
  #u = dw.findShiftsR(f,g,e,ngg)
  u = dw.findShiftsR(f,g,e)
  #u = dw.findShiftsR(e)
  h = dw.applyShifts(g,u)
  u = mul(u,dt)
  #u = add(u,-smin)
  gri = dw.getG1()
  gr = add(mul(floats(gri),dt),sf.first)
  du = wtu.backwardsDiff(u,dt)
  
  sh = Sampling(len(h),dt,sf.first)
  sl = Sampling(len(e[0]),dt,smin)
  
  #plot(sf,sg,sl,f,g,etran(e))
  #plot(sf,sg,sl,f,g,etran(e),u=u,gr=gr,png='e1')
  #plot(sf,sg,sl,f,h,etran(e),sh=sh,u=u,gr=gr,png='e2')
  plot1D(sf,f,"f")
  plot1D(sg,g,"g")
  plot1D(sh,h,"h")
  plot1D(sf,u,"u")
  plot1D(sf,u,"u",gr=gri)
  plot1D(sf,du,"vr",gr=gri)#,lims=[-0.02,0.02])
  plot1D(sf,du,"vr")#,lims=[-0.02,0.02])
  plot1D(sf,f,"f & g",sg,g)
  plot1D(sf,f,"f & h",sh,h)
  name = str(100*rmax)+"%_max-vr_"
  writeToText(dname+"_out_warped_s1_"+name,h)
  writeToText(dname+"_out_shifts_"+name,u)
  writeToText(dname+"_out_velocity_ratio_"+name,du)
  writeToText(dname+"_out_grid_pts"+name,gr)
	

"""
Computes the signal to noise ratio of a sequence x 
for all times using an sampling interval from [a,b]
by dividing the signal power by the average noise power.
"""
def computeSNR(a,b,x):
	nn = b-a+1
	nx = len(x)
	noise = 0
	for i in range(a,b+1):
		xx = x[i]
		noise += xx*xx
	noise /= nn
	#noise = sqrt(noise)
	sn = zerofloat(nx)
	for i in range(nx):
		xx = x[i]
		sn[i] = xx*xx/noise
	sn = sqrt(sn)
	re = RecursiveExponentialFilter(5)
	re.apply(sn,sn)
	mul(sn,1/max(sn),sn)
	return sn


###########################################################################

def readSig(name,fd=None,dd=None,window=None,cut=None):
	fname = "./data/"+name+".txt"
	parts = name.split("_")
	n = int(parts[1])
	if dd: d = float(parts[2])*dd
	else: d = float(parts[2])*0.001
	if fd: f = float(parts[3])*fd
	else: f = float(parts[3])
	s = Sampling(n,d,f)
	x = zerofloat(n)
	inpu = Scanner(File(fname));
	for i in range(n):
	  x[i] = inpu.nextFloat()
	inpu.close()
	if window:
		ma=0
		for i in range(n):
			if x[i]>ma: 
				ma = x[i]
				mi = i
		wva = window/d
		wst = int(mi-wva)
		if wst<0: wst=0
		wen = int(mi+wva)
		if wen>n: wen=n
		nn = wen-wst
		xw = copy(nn,wst,x)
		sn = Sampling(nn,d,f+wst*d)
		return sn,xw
	if cut:
		nc = int(round(f+cut/d))
		sc = Sampling(nc,d,f)
		xc = copy(nc,x)
		return sc,xc
	return s,x

def counter(x,lim):
  c = 0
  for i in range(len(x)):
    if x[i]>lim: c+=1 
  return c

def clip(x,lim):
  n = len(x)
  xc = zerofloat(n)
  for i in range(n):
    if (x[i]<lim): xc[i] = x[i]
    if (x[i]>-lim): xc[i] = x[i]
    if (x[i]>lim): xc[i] = lim
    if (x[i]<-lim): xc[i] = -lim
  return xc


def norm(x):
  return div(x,max(x))

def etran(e):
  e = norm(e)
  return transpose(pow(e,0.25))

# for n1>=n2
def equivlen(s1,s2,t1,t2):
  n1 = len(t1)
  n2 = len(t2)
  d1 = s1.getDelta()
  d2 = s2.getDelta()
  f1 = s1.getFirst()
  f2 = s2.getFirst()
  t3 = zerofloat(n1)
  s3 = Sampling(n1,d2,f2)
  rf = RandomFloat(23984)
  for i in range(n2):
    t3[i] = t2[i]
  for i in range(n2,n1):
    t3[i] = 0.05*rf.uniform()-0.025
  return s3,t3

def inro(x):
	return int(round(x))

def floats(x):
	n = len(x)
	xd = zerofloat(n)
	for i in range(n):
		xd[i] = float(x[i])
	return xd

def decimate(sx,x,m):
	sy = sx.decimate(m)
	print 'new count=',sy.count
	print 'new delta=',sy.delta
	y = copy(sy.count,0,m,x)
	print 'new len = ',len(y)
	return sy,y

def writeToText(name,data):
	n = len(data)
	out = BufferedWriter(FileWriter("data/out/"+name+".txt"))
	for i in range(n):
		out.write(str(data[i]))
		out.newLine()
	out.close()

###########################################################################


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
  if st.getLast()>st2.getLast(): stg = st.getLast()
  else: stg = st2.getLast()
  p1.setHLimits(0,stg)
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

def plot(sf,sg,sl,f,g,c,sh=None,h=None,s=None,u=None,gr=None,\
		clip=None,paper=None,slide=None,png=None):
	n,nlag = len(c[0]),len(c)
	#slag = Sampling(nlag,1.0,-(nlag-1)/2)
	panel = PlotPanel(2,1,PlotPanel.Orientation.X1RIGHT_X2UP)
	panel.mosaic.setHeightElastic(0,25)
	panel.mosaic.setHeightElastic(1,75)
	panel.setVLimits(0,-1,1)
	panel.setVLimits(1,sl.first,sl.last)
	panel.setHLimits(sf.first,sf.last)
	fv = panel.addPoints(0,0,sf,f)
	if sh: gv = panel.addPoints(0,0,sh,g)
	else: gv = panel.addPoints(0,0,sg,g)
	gv.setLineColor(Color.RED)
	if h:
		hv = panel.addPoints(0,0,sh,h)
		hv.setLineColor(Color.BLUE)
	cv = panel.addPixels(1,0,sf,sl,c)
	cv.setInterpolation(PixelsView.Interpolation.NEAREST)
	if clip:
		cv.setClips(0.0,clip)
	cv.setColorModel(ColorMap.JET)
	if s:
		sv = panel.addPoints(1,0,s)
		sv.setLineColor(Color.WHITE)
		sv.setLineStyle(PointsView.Line.DOT)
		sv.setLineWidth(3)
	if u:
		uv = panel.addPoints(1,0,sf,u)
		uv.setLineColor(Color.WHITE)
		uv.setLineWidth(3)
		minu = min(u); maxu = max(u)
		#panel.setVLimits(1,max(sl.first,min(minu*.75,minu*1.25)),
		#									min(sl.last,max(maxu*.75,maxu*1.25)))
		if gr:
			j=0
			ug = zerofloat(len(gr))
			for i in range(len(u)):
				t = inro((gr[j]-sf.first)/sf.delta)
				if i==t:
					ug[j] = u[i]
					j = j+1
			dg = panel.addPoints(1,0,gr,ug)
			dg.setStyle('rO')
			dg.setLineWidth(5)
	panel.setHLabel("time (s)")
	panel.setVLabel(0,"f & g")
	panel.setVLabel(1,"lag")
	frame = PlotFrame(panel)
	frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
	frame.setFontSize(18)
	frame.setSize(1000,800)
	frame.setVisible(True)
	if paper:
		frame.setSize(600,400)
		frame.setFontSizeForPrint(8.0,240.0)
		frame.paintToPng(720,3.33,_pngDir+paper+'.png')
	if slide:
		frame.setSize(700,500)
		frame.setFontSizeForSlide(1.0,1.0)
		frame.paintToPng(720,7.00,_pngDir+slides+'.png')
	if png:
		frame.paintToPng(720,3.33,_pngDir+png+'.png')
	frame.setVisible(True)


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

def plot1D(s,x,title,s2=None,x2=None,lims=None,gr=None):
  sp = SimplePlot()
  sp.addPoints(s,x)
  if s2 and x2:
    pv = sp.addPoints(s2,x2)
    pv.setLineColor(Color.RED)
  if lims:
    sp.setVLimits(lims[0],lims[1])
  if gr:
    cgs = gr
    cgsv = add(s.first,mul(floats(gr),s.delta))
    cgsx = zerofloat(len(cgs))
    for ic in range(len(cgs)):
      cgsx[ic] = x[cgs[ic]] 
    pv = sp.addPoints(cgsv,cgsx)
    pv.setStyle("rO")
  sp.setTitle(title)
  sp.setSize(1000,300)

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
