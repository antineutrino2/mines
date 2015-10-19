from imports import *
"""
Demo script for structure tensor-guided methods: 
	1. oriented smoothing 
	2. blended neightbor interpolation
Author: Andrew Munoz, Colorado School of Mines
Version: 15/05/2013
"""

# Amplitude
p0 = 0.0
# Linearity
p1 = 1.5
# Planarity (3D only)
p2 = 1.0

# tpst image samplings
#s1 = Sampling(301,0.004,0.250)
#s2 = Sampling(101,0.025,2.500)
#s3 = Sampling(101,0.025,1.000) # (3D only)

# Eshan test
s1 = Sampling(101,0.004,1.5)
s2 = Sampling(701,1.0,0.0)
dataDir = "/Users/amunoz/Downloads/"

def main(args):
	#go2D()
	#go3D()
	ehsanTest()

def ehsanTest():
	f = zerofloat(s1.count,s2.count)
	l = zerofloat(s1.count,s2.count)
	f = readData(f,'seismic')
	l = readData(l,'logs')
	t = makeTensors2(f)
	goInterpLogs2(f,l,t)

def go2D():
	f = getData(51)
	t = makeTensors2(f)
	plot2(f,"f with tensors",tensors=t)
	# increase sigma for smoother 
	sigma = 8.0		
	goSmooth2(f,t,sigma)
	goInterp2(f,t)


def go3D():
	f = getData()
	t = makeTensors3(f)
	plot3(f,"f with tensors",tensors=t)
	# increase sigma for smoother 
	sigma = 4.0		
	#goSmooth3(f,t,sigma)
	#goInterp3(f,t)


##############################################################################
## 2D Methods

"""
Returns a EigenTensors2 class for use in smoothing and interpolation
f: 2D image used to make structure tensors
"""
def makeTensors2(f):
	lof = LocalOrientFilter(4.0)
	d = lof.applyForTensors(f)
	d.invertStructure(p0,p1)
	return d

"""
Applies a tensor guided smoothing using sigma as smoothing control
"""
def goSmooth2(f,t,sigma):
	g = copy(f)
	c = sigma*sigma*0.5
	ls = LocalSmoothingFilter()
	ls.apply(t,c,f,g)
	plot2(f,"f")
	plot2(g,"smoothed f")

"""
Decimates traces in the image f and then attempts to reconstruct it using 
tensor-guided interpolation with blended neighbor interpolation.
"""

def goInterp2(f,t):
	# Parameter for blending
	blnd = 1.5
	# Parameter to control maximum time
	gtmax = 100
	# Decimation factor 
	dec = 5
	# Parameters
	pnull = -999.9
	n1 = s1.count
	n2 = s2.count
	g = zerofloat(n1,n2)
	# decimate f
	gn = decimate2(f,pnull,dec)
	# Get scattered samples and coordinates
	v,x1,x2 = SimpleGridder2.getGriddedSamples(pnull,s1,s2,gn)
	# Interpolate
	bg = BlendedGridder2(t,v,x1,x2)
	bg.setTimeMax(gtmax)
	bg.setSmoothness(blnd)
	gt = bg.gridNearest(pnull,gn)
	bg.gridBlended(gt,gn,g)
	plot2(f,"f")
	plot2(g,"f reconstructed")

def goInterpLogs2(f,l,t):
	# Parameter for blending
	blnd = 0.5
	# Parameter to control maximum time
	gtmax = 1000
	# Parameters
	pnull = 9999
	n1 = s1.count
	n2 = s2.count
	g = zerofloat(n1,n2)
	# Get scattered samples and coordinates
	v,x1,x2 = SimpleGridder2.getGriddedSamples(pnull,s1,s2,l)
	# Interpolate
	bg = BlendedGridder2(t,v,x1,x2)
	bg.setTimeMax(gtmax)
	bg.setSmoothness(blnd)
	gt = bg.gridNearest(pnull,l)
	bg.gridBlended(gt,l,g)
	plot2(l,"logs")
	plot2(g,"logs interp",clips=[4000,10000])
	plot2(f,"seismic with tensors",tensors=t,cmap=gray)

def decimate2(f,pnull,dec):
	n1 = s1.count
	n2 = s2.count
	gn = fillfloat(pnull,n1,n2)
	n2n = n2/dec
	copy(n1,n2n,0,0,1,dec,f,0,0,1,dec,gn)
	return gn

##############################################################################
## 3D Methods

"""
Returns a EigenTensors3 class for use in smoothing and interpolation
f: 3D image used to make structure tensors
"""
def makeTensors3(f):
	lof = LocalOrientFilter(4.0)
	d = lof.applyForTensors(f)
	d.invertStructure(p0,p1,p2)
	return d

"""
Applies a tensor guided smoothing using sigma as smoothing control
"""
def goSmooth3(f,t,sigma):
	g = copy(f)
	c = sigma*sigma*0.5
	ls = LocalSmoothingFilter()
	ls.apply(t,c,f,g)
	plot3(f,"f")
	plot3(g,"smoothed f")

"""
Decimates traces in the image f and then attempts to reconstruct it using 
tensor-guided interpolation with blended neighbor interpolation.
"""
def goInterp3(f,t):
	# Parameter for blending
	blnd = 0.5
	# Parameter to control maximum time
	gtmax = 20
	# Decimation factor 
	dec = 5
	# Parameters
	pnull = -999.9
	n1 = s1.count
	n2 = s2.count
	n3 = s3.count
	g = zerofloat(n1,n2,n3)
	# decimate f
	gn = decimate3(f,pnull,dec)
	# Get scattered samples and coordinates
	v,x1,x2,x3 = SimpleGridder3.getGriddedSamples(pnull,s1,s2,s3,gn)
	# Interpolate
	bg = BlendedGridder3(t,v,x1,x2,x3)
	bg.setTimeMax(gtmax)
	bg.setSmoothness(blnd)
	gt = bg.gridNearest(pnull,gn)
	bg.gridBlended(gt,gn,g)
	plot3(f,"f")
	plot3(g,"f reconstructed")

def decimate3(f,pnull,dec):
	n1 = s1.count
	n2 = s2.count
	n3 = s3.count
	gn = fillfloat(pnull,n1,n2,n3)
	dec = 10
	n2n = n2/dec
	n3n = n3/dec
	copy(n1,n2n,n3n,0,0,0,1,dec,dec,f,0,0,0,1,dec,dec,gn)
	return gn

##############################################################################
## Plots

# Options
gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
xrxu = PlotPanel.Orientation.X1RIGHT_X2UP
xdxr = PlotPanel.Orientation.X1DOWN_X2RIGHT
aplt = PlotPanel.AxesPlacement.LEFT_TOP
inear = PixelsView.Interpolation.NEAREST
eoc = PlotFrame.EXIT_ON_CLOSE

def plot2(f,title,tensors=None,cmap=jet,clips=None):
	pp = PlotPanel(xdxr)
	px = pp.addPixels(s1,s2,f)
	pp.setHLabel("Distance (km)")
	pp.setVLabel("Time (s)")
	px.setColorModel(cmap)
	pp.setTitle(title)
	pp.addColorBar()
	if clips:
		px.setClips(clips[0],clips[1])
	if tensors:
		tv = TensorsView(s1,s2,tensors)
		tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
		tv.setLineColor(Color.YELLOW)
		tv.setLineWidth(2)
		#tv.setEllipsesDisplayed(40)
		fac1 = 30; fac2 =30 
		e1 = Sampling(s1.count/fac1+1,s1.delta*fac1,s1.first)
		e2 = Sampling(s2.count/fac2+1,s2.delta*fac2,s2.first)
		tv.setEllipsesDisplayed(e1,e2)
		tv.setScale(1.00)
		tile = pp.getTile(0,0)
		tile.addTiledView(tv)
	frame = PlotFrame(pp)
	frame.setSize(1000,1000)
	frame.setFontSize(24)
	frame.setDefaultCloseOperation(eoc)
	frame.setVisible(True)

def plot3(f,title,tensors=None):
	cmap = gray
	world = World()
	ipg = ImagePanelGroup(s1,s2,s3,f)
	ipg.setColorModel(cmap)
	world.addChild(ipg)
	if tensors:
		a = 20
		addTensorsInImage(ipg.getImagePanel(Axis.X),tensors,a)
		addTensorsInImage(ipg.getImagePanel(Axis.Y),tensors,a)
		addTensorsInImage(ipg.getImagePanel(Axis.Z),tensors,a)
	sf = SimpleFrame(world)
	ov = sf.getOrbitView()	
	ov.setAxesScale(1.0,1.0,2.0)
	ov.setScale(2.0)
	sf.setTitle(title)
	sf.setVisible(True)

def addTensorsInImage(ipg,et,a):
	tp = TensorsPanel(s1,s2,s3,et)
	tp.setEllipsoidSize(a)
	ipg.getFrame().addChild(tp)
	return tp

##############################################################################
## Data

def getData(a=None):
	xo = zerofloat(s1.count,s2.count,s3.count)
	ais = ArrayInputStream("./data/tpst.dat")
	ais.readFloats(xo)
	ais.close()
	if a:
		return xo[a]
	return xo

def readData(x,fname):
	ais = ArrayInputStream(dataDir+fname+".dat");
	ais.readFloats(x)
	ais.close()
	return x

#------------------------------------------------------------------------------#

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

