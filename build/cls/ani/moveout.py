from imports import *

sub = 'moveout/'
titles = None
prints = None

def main(args):
	goParams1()
	xe,te = readExactTimes()
	ts1,ts2 = computeTraveltimes(xe,te)
	to = getOptimalVelocities(xe,te)
	plot(xe,te,png='exact')
	plot(xe,ts1,png='unityC')
	plot(xe,ts2,png='optC')
	plot(xe,te,ts1,png='exact_unityC')
	plot(xe,te,ts2,png='exact_optC')
	plot(xe,ts1,ts2,png='unityC_optC')
	plot(xe,te,ts1,ts2,png='exact_unityC_optC')
	plot(xe,percDiff(ts1,ts2),\
				png='pd_unityC_optC',vlabel='% difference')
	plot(xe,percDiff(te,ts1),\
				png='pd_exact_unityC',vlabel='% difference')
	plot(xe,percDiff(te,ts2),\
				png='pd_exact_optC',vlabel='% difference')
	#plot(xe,percDiff(te,to),\
	#			png='pd_exact_optV',vlabel='% difference')
	#plot(xe,te,to,png='exact_optV')

def computeTraveltimes(xe,te):
	nl = len(layers)
	nx = sx.count
	veff = computeVnmoEff()
	neff = computeEtaEff(veff)
	vhor = computeVhor(veff,neff)
	tsl1 = getTimesForLayers(veff,neff,vhor,xe,2.0)
	tsl2 = getTimesForLayers(veff,neff,vhor,xe,1.1655973)
	tsc,minC,ca = getTimesForBestC(2,veff,neff,vhor,xe,te)
	plot(ca[0],ca[1],hlabel='C',vlabel='fit',png='cfit')
	return tsl1[2],tsl2[2]

def getTimesForLayers(veff,neff,vhor,xe,C):
	ts = zerofloat(nx,nl)
	for il in range(nl):
		ton = layers[il][2]
		toto = ton*ton
		vnmo = veff[il]
		vnmo2 = vnmo*vnmo
		for ix in range(nx):
			xx = xe[ix]*xe[ix]
			ts[il][ix] = sqrt(toto+xx/vnmo-\
				((vhor[il]-vnmo)*xx*xx)/\
				(vnmo*(toto*vnmo2+C*vhor[il]*xx)))
	return ts

def getTimesForBestC(il,veff,neff,vhor,xe,te):
	ts = zerofloat(nx,nl)
	fc = 0.9; lc = 1.5;
	dc = 0.001
	nc = int((lc-fc)/dc)+1
	mc = sum(te); minC = 0.0
	ca = zerofloat(nc)
	for ic in range(nc):
		C = ic*dc+fc	
		ton = layers[il][2]
		toto = ton*ton
		vnmo = veff[il]
		vnmo2 = vnmo*vnmo
		for ix in range(nx):
			xx = xe[ix]*xe[ix]
			ts[il][ix] = sqrt(toto+xx/vnmo-\
				((vhor[il]-vnmo)*xx*xx)/\
				(vnmo*(toto*vnmo2+C*vhor[il]*xx)))
		tc = sum(abs(ArrayMath.sub(ts[il],te)))
		ca[ic] = tc
		if tc<mc:
			mc = tc
			minC = C
	print 'min C =',minC
	return ts,minC,[Sampling(nc,dc,fc),ca]

def getOptimalVelocities(xe,te):
	veff = computeVnmoEff()
	neff = computeEtaEff(veff)
	vhor = computeVhor(veff,neff)
	print 'analytic veff =',sqrt(veff[-1])
	print 'analytic vhor =',sqrt(vhor[-1])
	print 'analytic neff =',neff[-1]
	ts = zerofloat(nx)
	to = te[0]
	toto = to*to
	vnf = 1.9; vnl = 3.0;	
	vhb = 0.0; vhl = 3.1;	
	dv = 0.01
	nvn = int((vnl-vnf)/dv)
	nvh = int((vhl-vnf)/dv)
	mv=sum(te);mvnmo=0;mvhor=0;
	for ivn in range(nvn):
		vnmo = vnf+ivn*dv
		vnmo2 = vnmo*vnmo
		vnmo4 = vnmo2*vnmo2
		nvh = int((vhl-vnmo)/dv)
		for ivh in range(nvh):
			vhor = vnmo+ivh*dv
			vhor2 = vhor*vhor
			for ix in range(nx):
				xx = xe[ix]*xe[ix]
				ts[ix] = sqrt(toto+xx/vnmo2-\
					((vhor2-vnmo2)*xx*xx)/\
					(vnmo2*(toto*vnmo4+vhor2*xx)))
			tv = sum(abs(ArrayMath.sub(ts,te)))
			if tv<mv:
				mv = tv
				mvnmo = vnmo
				mvhor = vhor
	print 'optimal vnmo=',mvnmo
	print 'optimal vhor=',mvhor
	print 'optimal eta =',0.5*(pow(mvhor/mvnmo,2)-1)
	to = zerofloat(nx)
	mvnmo2 = mvnmo*mvnmo
	mvhor2 = mvhor*mvhor
	for ix in range(nx):
		xx = xe[ix]*xe[ix]
		to[ix] = sqrt(toto+xx/mvnmo2-\
			((mvhor2-mvnmo2)*xx*xx)/\
			(mvnmo2*(toto*mvnmo2*mvnmo2+mvhor2*xx)))
	return to

def computeVnmoEff():
	vnmo = FloatList()
	vsum=0
	for l in layers:
		toi = l[1]
		ton = l[2]
		vp = l[3]
		vsum += vp*vp*toi
		vnmo.add(vsum/ton)
	return vnmo.trim()

def computeEtaEff(veff):
	eta = FloatList()
	vsum=0; il=0;
	for l in layers:
		toi = l[1]
		ton = l[2]
		vp = l[3]
		et = l[5]
		vnmo = veff[il]
		vsum += vp*vp*vp*vp*(1+8*et)*toi
		eta.add(0.125*(vsum/(vnmo*vnmo*ton)-1))
		il+=1
	return eta.trim()

def computeVhor(veff,neff):
	n = len(veff)
	vhor = zerofloat(n)
	for i in range(n):
		vhor[i] = veff[i]*(1+2*neff[i])
	return vhor

def readExactTimes():
	fname = 'data/PwaveNonHyp.txt'
	f = open(fname)
	x = FloatList()
	t = FloatList()
	for line in f:
		column = line.split(" ")
		x.add(float(column[1]))
		t.add(float(column[3]))
	xa = x.trim()
	ta = t.trim()
	return xa,ta

def goParams1():
	global layers
	global sx
	global nx,nl
	dx = 0.1
	fx = 0.0
	lx = 5.0
	sx = Sampling(int((lx-fx)/dx)+1,dx,fx)
	layers = []
	addLayer(1.480,1.941,0.8,0.230,0.011)
	addLayer(1.980,2.168,0.9,0.247,0.010)
	addLayer(2.244,2.685,1.1,0.264,0.016)
	nl = len(layers)
	nx = sx.count

def addLayer(z,vp,vs,ep,de):
	if len(layers)>0: 
		to = layers[-1][2]
		dz = z-layers[-1][0]
	else: 
		to = 0
		dz = z
	vnmo = vp*sqrt(1+2*de)
	layers.append((z,2*dz/vp,to+2*dz/vp,vnmo,vs,(ep-de)/(1+2*de)))

def findSampling(x):
	n = len(x)
	dx = x[0]
	for i in range(1,n):
		dt = abs(x[i]-x[i-1])
		if dt<dm:
			dm = dt
	return dm

def percDiff(x,y):
	return mul(100,div(ArrayMath.sub(x,y),y))

################################################################################
# plots

def plot(x1,x2,x3=None,x4=None,title=None,vlabel=None,hlabel=None,png=None):
	pp = PlotPanel(PlotPanel.Orientation.X1RIGHT_X2UP)
	l1 = pp.addPoints(x1,x2)
	pp.setVLabel("time (s)")
	pp.setHLabel("offset (km)")
	if x3:
		l2 = pp.addPoints(x1,x3)
		l2.setLineColor(Color.RED)
	if x4:
		l3 = pp.addPoints(x1,x4)
		l3.setLineColor(Color.BLUE)
	if title:
		pp.setTitle(title)
	if titles:
		pp.setTitle(png)
	if vlabel:
		pp.setVLabel(vlabel)
	if hlabel:
		pp.setHLabel(hlabel)
	frame = PlotFrame(pp)
	frame.setFontSize(24)
	frame.setSize(600,600)
	frame.setVisible(True)
	if png:
		if prints: frame.paintToPng(720,3.33,'png/'+sub+png+'.png')

################################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
if __name__=="__main__":
  SwingUtilities.invokeLater(RunMain()) 

