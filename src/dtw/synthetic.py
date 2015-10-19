from imports import *
from tp import SeismicModel1D
from wt import WTutils

wtu = WTutils()

# Depth sampling
nz = 1000
dz = 0.1524
fz = 0.0 
sz = Sampling(nz,dz,fz)

# Time sampling 
nt = 700
dt = 0.002
ft = 0.0
st = Sampling(nt,dt,ft)

# Synthetic properties
decay = 0.01
r1 = 0.0 
os = 2.0
nr = 1
rm = 20.0
q = 1000000
fp = 3.0/35

# Layers
q = 10000000
rfv = RandomFloat(432143)
#rfd = RandomFloat(124598)
d = fillfloat(1.0,nz)
v = zerofloat(nz)
for i in range(nz):
	v[i] = rfv.uniform()*5.0

def main(args):
	#goTest1()
	print 'prop-a'
	sm1 = goTest2()
	#print 'prop-d'
	#sm2 = goTest2d()
	#print 'prop-c'
	#sm3 = goTest2c()
	#goTest3()
	#plotCompare(st,sm1,sm2,sm3,"Pa(black), Pd(red), C(blue)",lim1=min(sm2),lim2=max(sm2))

def goTest2():
	ssm = SyntheticSeismogram()
	#ssm.setZones(sz,v,d,q)
	ssm.addZone(0.0,1000,1.0,q)
	ssm.addZone(dz*nz/2.0+fz,2000,1.0,q)
	ssm.addSource(0.0,1.0)
	ssm.addRicker(35)
	ssm.setSRC(r1)
	ssm.setDecay(decay)
	ssm.setOversampling(os)
	#ssm.turnOffMultiples(True)
	#Set recievers
	#for ir in range(nr):
	#	ssm.addReciever(rm*ir+rm)
	ssm.addReciever(0.0)
	seis = ssm.getSeismogram(st)
	sy = seis[0]
	#sy = wtu.addRickerWavelet(fp,seis[0])
	# Plot
	#plotPanels(st,seis,"time (s)","recievers",rm)
	plot(st,sy,"time (s)","propa")
	return sy

def goTest2d():
	ssm = SeismicModel1D()
	ssm.setSourceType(SeismicModel1D.SourceType.ISOTROPIC);
	ssm.setSensorType(SeismicModel1D.SensorType.PRESSURE);
	ssm.addLayer(0.00,v[0],d[0],q)
	for i in range(nz):
		ssm.addLayer(fz+i*dz,v[i],d[i],q)
	ssm.addSource(0.0,1.0)
	ssm.setSurfaceReflectionCoefficient(r1)
	ssm.setDecay(decay)
	ssm.setOversample(os)
	#Set recievers
	#for ir in range(nr):
	#	ssm.addSensor(rm*ir+rm)
	ssm.addSensor(0.0)
	w0 = 0.5/dt
	seis = ssm.makeSeismograms(nt,dt,w0,False)
	sy = wtu.addRickerWavelet(fp,seis[0])
	# Plot
	#plotPanels(st,seis,"time (s)","recievers-d",rm)
	plot(st,sy,"time (s)","propd")
	return sy

def goTest2c():
	r = zerofloat(nz)
	#for i in range(1,nz):
	#	im = i-1
	#	r[i] = (v[i]*d[i]-v1[im]*d[im])/(v[i]*d[i]+v[im]*d[im])
	#r[0] = r[1] 
	r[nt/2] = 1.0/3.0;
	sy = wtu.addRickerWavelet(fp,r)
	plot(st,sy,"time (s)","propc")
	return sy

def goTest1():
	ssm = SyntheticSeismogram()
	ssm.addZone(0.00,1.0,1.0,1000000)
	ssm.addSource(0.0,1.0)
	ssm.setSRC(r1)
	ssm.setDecay(decay)
	ssm.setOversampling(os)
	#Set recievers
	for ir in range(nr):
		ssm.addReciever(rm*ir+rm)
	seis = ssm.getSeismogram(st)
	# Plot
	plotPanels(st,seis,"time (s)","recievers",rm)

def goTest3():
	ssm = SyntheticSeismogram()
	vl = fillfloat(3,nz)
	dl = fillfloat(1.0,nz)
	#ssm.setZones(sz,vl,dl,100000)
	#for i in range(nz):
	#	ssm.addZone(fz+i*dz,vl[0],dl[0],100000)
	ssm.addZone(0.0,3000.0,1.0,100000)
	ssm.addSource(0.0,1.0)
	ssm.setSRC(r1)
	ssm.setDecay(decay)
	ssm.setOversampling(os)
	#Set recievers
	for ir in range(nr):
		ssm.addReciever(rm*ir+rm)
	seis = ssm.getSeismogram(st)
	# Plot
	plotPanels(st,seis,"time (s)","recievers",rm)


def plot(s,x,vaxis,haxis):
	si = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
	si.addPoints(s,x)
	si.setSize(500,1000)
	si.setHLabel(haxis)
	si.setVLabel(vaxis)
	#si.setHLimits(-1,1)

def plotPanels(s,x,vaxis,title,rm):
	nr = len(x)
	pp = PlotPanel(1,nr,PlotPanel.Orientation.X1DOWN_X2RIGHT)
	for ir in range(nr):
		pp.addPoints(0,ir,s,x[ir])
		pp.setHLabel(ir,str(rm*ir+rm)+"m")
	pp.setTitle(title)
	pp.setVLabel(vaxis)
	frame = PlotFrame(pp)
	frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
	width = 500*nr
	if width>2000: width = 2000
	frame.setSize(width,1000)
	frame.setVisible(True)

def plotCompare(s,x1,x2,x3,title,lim1=None,lim2=None):
	pp = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
	pp.setTitle(title)
	pp.setVLabel('time (s)')
	ln1 = pp.addPoints(s,x1)
	ln2 = pp.addPoints(s,x2)
	ln3 = pp.addPoints(s,x3)
	ln2.setLineColor(Color.RED)
	ln3.setLineColor(Color.BLUE)
	if lim1 and lim2:
		pp.setHLimits(lim1,lim2)
	frame = PlotFrame(pp)
	frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
	frame.setSize(500,1000)
	frame.setVisible(True)


#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
