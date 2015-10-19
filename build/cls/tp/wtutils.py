"""
Jython utilities for Teapot Dome.
Author: Andrew Munoz, Colorado School of Mines
Version: 2013.04.09
"""
from dtw import *
from imports import *
from tputils import *
from wtplots import *
from viewer import *

# My Mac
_tpDir = "/Users/amunoz/Home/data/tp/"
# Backus
#_tpDir = "/data/amunoz/tp/"

workDir = "paper2013/"
#workDir = "slides2014/"

seismicDatum = 1.6764-0.003048 # 5500-10 ft
dataDir = _tpDir+"wellties/"
topsDir = _tpDir+"csm/welllogs/welltops/"
flatDir = _tpDir+"csm/flats/"
_pngDir = "/Users/amunoz/Home/pics/tppics/2013/"+workDir
wtu = WTutils()
_gtmax = 50
dz = 0.0001524 # Log depth sampling interval
dz2 = 0.002 # Depth image sampling interval
sigmanorm = 100
fp = 35

#fms = ["F2WC","DKOT","LKOT","CRMT","RDPK","A Sand","C1 Dolo"]
fms = ["DKOT","LKOT","CRMT","RDPK","A Sand","C1 Dolo"]
#fms = ["F2WC","DKOT","A Sand","C1 Dolo"]

#s1 = Sampling(3002,0.001,0.000)
#s1 = Sampling(876,0.002,0.250)
s1 = Sampling(1001,0.002,0.000)
s2 = Sampling(357,0.025,0.000)
s3 = Sampling(161,0.025,0.000)

# tp depth image
s1z = Sampling(2762,0.002,0.000)

# Global complete print file control
_prints = True
geo = True
cwp = False

###############################################################################
###############################################################################
## Warping utilities

"""
Uses Dynamic Time Warping to compute shifts and apply them to f for 
multiple wells
"""
def warpIndividually(tpimg,ssy,sy,sz,x2w,x3w,tzl,vl,smin,smax,rmin,rmax,dr,\
										 wln=None,ph1=None,ph2=None,nr=1,wp=None):
  nw = len(ssy)
  pred = 0.0
  u,h,sh,tz,vint,zt,szt,pv=[],[],[],[],[],[],[],[]
  if ph2:
    traceMulti = []
    for iw in range(nw):
      traceMulti.append(getTraces(tpimg,x2w[iw],x3w[iw],tonly=True,nr=nr))
    phase = rotatePhaseMulti(sy,traceMulti,rmin,rmax,dr,vb=vb,ss=ssy)
    for iw in range(nw):
      sy[iw] = applyPhaseRot(phase,sy[iw])
  for iw in range(nw):
    traces = getTraces(tpimg,x2w[iw],x3w[iw],tonly=True,nr=nr)
    if ph1:
      sy[iw] = rotatePhase(sy[iw],traces,rmin,rmax,dr,ss=ssy[iw],vb=vb)
      #sy[iw] = rotateAllPhase(sy[iw],traces,vrmin,vrmax,dr)
    u0,h0,sh0 = getWarpingShiftsOld(ssy[iw],s1,sy[iw],traces,rmin,rmax,dr)
    #u0,h0,sh0 = getWarpingShifts(ssy[iw],s1,sy[iw],traces,smin,smax,rmin,rmax,dr)
    tzl0 = tzl[iw]
    tz0,vint0,pv0 = updateTDandVelocity(sz[iw],u0,tzl0,ssy[iw])
    #zt0,szt0 = getzt(sz[iw],tz0)
    wn = str(wln[iw])
    #plotCurve(traces[0],s1,sy[iw],ssy[iw],d="t",title=wn+" before")
    #plotCurve(traces[0],s1,h0,sh0,d="t",title=wn+" after")
    #print 'swt R = ',correlation(h0,traces[0],sh0,s1)
    #zt.append(zt0); szt.append(szt0); 
    pv.append(pv0); u.append(u0); h.append(h0); 
    sh.append(sh0); tz.append(tz0); vint.append(vint0)
  #print 'Average Predictability =',pred/nw
  return u,h,sh,tz,vint,zt,szt,pv

"""
Uses Dynamic Time Warping to compute shifts and apply them to f
"""
def getWarpingShifts(sf,sg,f,tr,smin,smax,rmin,rmax,dr,prints=None):
  g = tr[0]
  dt = sf.delta
  dw = DynamicWarpingWT(smin,smax,sg)
  dw.setStrainLimits(rmin,rmax)
  dw.setSmoothing(dr)
  dw.setSyntheticSampling(sf)
  u = dw.findShifts(sf,f,sg,g[0])
  h = dw.applyShifts(f,u)
  sh = Sampling(len(h),dt,u[0]+s1.first)
  #print 'min shift=',sf.first-sh.first
  #print 'max shift=',sh.last-sf.last
  return u,h,sh

"""
Uses Dynamic Time Warping to compute shifts and apply them to f
"""
def getWarpingShiftsOld(sf,sg,f,tr,rmin,rmax,dr,prints=None):
  print 'using old swt class...'
  g = tr[0]
  dt = sf.delta
  nl = sg.count-sf.count+1
  if prints:
    print 'kmax=',inro(rmax/dr)
    print 'kmin=',inro(rmin/dr)
    print 'dr= ',dr
  dw = DynamicWarpingSWT(nl,rmin,rmax,dr)
  e = dw.computeErrorsMulti(f,tr)
  u = dw.findShiftsR(e)
  h = dw.applyShifts(u,f)
  u = mul(u,dt)
  #h,fh = dw.getWarpedSyntheticSeismogram(sf,u,f)
  #sh = Sampling(len(h),dt,fh[0])#u[0]+s1.first)
  sh = Sampling(len(h),dt,u[0]+s1.first)
  #print 'min shift=',sf.first-sh.first
  #print 'max shift=',sh.last-sf.last
  return u,h,sh



###############################################################################
###############################################################################
## 2D/3D Multiple-tie utilities

"""
Sorts the well list by x2 values
"""
def sortWells(wells):
	nw = len(wells)
	wells2,ysort = [],[]
	for i in wells:
		_set,_ID = getWellFromSet(i) 
		yv0 = getLogs(_set,_ID,getyv=True)
		ysort.append((i,yv0))
	ysort = sorted(ysort,key=itemgetter(1))
	for i in range(nw):
		wells2.append(ysort[i][0])
	print 'sorted wells:',wells2
	return wells2

"""
Reduces arrays for well list (w1) where elements of logs are length nw1
w1 and logs must be sorted to same order
"""
def reduceWellList(w1,w2,logs):
	nw1 = len(w1)
	nw2 = len(w2)
	nl  = len(logs)
	rl = [[] for i in range(nl)]
	for il in range(nl):
		log = logs[il]	
		for iw in range(nw1):
			if w1[iw] in w2:
				rl[il].append(log[iw])
	return rl

"""
Computes multiple synthetic seismograms and returns them in a list along
with the samplings, time-depths, and depth-times 
"""
def getSeismograms(sz,vl,dl,q,phase=None,wav=None,simple=None,wells=None):
	nw = len(sz)
	ssy,sy,tzl,ztl,sztl = [],[],[],[],[]
	normalize=True
	wavt,phaset=None,None
	for iw in range(nw):
		if wav: wavt=wav[iw]
		if phase: phaset = phase[iw]
		if simple:
			ssy0,sy0,tzl0 = getSimpleSeismogram(\
										sz[iw],vl[iw],mul(vl[iw],dl[iw]),phaset,normalize,wavt)
		else:	
			ssy0,sy0,tzl0 = getPropagatorSeismogram(\
										sz[iw],vl[iw],dl[iw],q,phaset,normalize,wavt)#,taper=True)
		#tzlt = wtu.tzVlog(vl[iw],sz[iw].first,sz[iw].delta)
		#zt0,szt0 = getzt(sz[iw],tzl0)
		ssy.append(ssy0); sy.append(sy0); tzl.append(tzl0)
		#ztl.append(zt0); sztl.append(szt0);
	return ssy,sy,tzl,ztl,sztl

"""
Gets all logs from specified well set and applys despike or smooth for certian
samples in the top of the logs
"""
def getMultipleLogs(wells,datums,conly=None,id=None):
	nw = len(wells); j=0
	ail,vl,dl,gl,sz,x2,x3,ids=[],[],[],[],[],[],[],[]
	for i in wells:
		_set,_ID = getWellFromSet(i) 
		if id: print 'well',i,'has UWI',_ID
		ids.append(_ID)
		# smooth samples, sigma, cut
		# TODO fix smooth methods
		smooth=[2,1,100]
		ail0,yv0,xv0,zv0,vl0,dl0,sz0,gll = getLogs(\
														_set,_ID,g=True,smooth=smooth,cut=[0,10],dspk=[40,2.0])
		#if i==0 or i==2:
		#	ail0,yv0,xv0,zv0,vl0,dl0,sz0,gll = getLogs(_set,_ID,g=True,dspk=[40,2.0],smooth=smooth)
		#if i==3:
			#ail0,yv0,xv0,zv0,vl0,dl0,sz0,gll = getLogs(_set,_ID,g=True,cut=[0,10],smooth=smooth)
		#if i==9:
		#	ail0,yv0,xv0,zv0,vl0,dl0,sz0,gll = getLogs(_set,_ID,g=True,dspk=[40,2.0],smooth=smooth)
		#if i==10 or i==12:
			#ail0,yv0,xv0,zv0,vl0,dl0,sz0,gll = getLogs(_set,_ID,g=True,smooth=smooth,cut=[0,40])
		x2.append(yv0); x3.append(xv0); 
		if i==5: zv0-=3.0
		sz.append(sz0)
		ail.append(ail0); vl.append(vl0); dl.append(dl0); 
		gl.append([gll[0],Sampling(len(gll[0]),dz,gll[1]+datums[j])])
		#plotLogPanel(sz[j],vl[j],dl[j],reflectivity(ail[j]))
		j+=1
	if conly:
		return x2,x3
	showMaxWellDist(x2,x3)
	#for iw in range(len(wells)):
	#	print str(wells[iw])+" has coords: "+str(x2w[iw])+","+str(x3w[iw])
	return ail,vl,dl,gl,sz,x2,x3,ids


###############################################################################
###############################################################################
## 2D Multi-tie utilities

"""
Gets finely sampled curve through points (well locations)
"""
def getCurveThruWells(x2s,x3s,dy,edg=0.3):
	x2o = x2s
	x3o = x3s
	ns = len(x2s)
	#put first well at 0,0
	for i in range(2): # 2 iterations should be sufficient
		ds = zerofloat(ns)
		ds[0] = 0.0
		for js in range(1,ns):
			ds[js] = ds[js-1]+hypot(x2s[js]-x2s[js-1],x3s[js]-x3s[js-1])
		if i==0: xsp = copy(ds)
		ci2 = CubicInterpolator(ds,x2s)
		ci3 = CubicInterpolator(ds,x3s)
		dsy = ds[-1]
		smin,smax = ds[0],ds[-1]
		smin -= edg
		smax += edg
		xsp = add(xsp,edg)
		ns = 1+int((smax-smin)/dy)
		ds = (smax-smin)/(ns-1)
		sj = rampfloat(smin,ds,ns)
		x2s = zerofloat(ns)
		x3s = zerofloat(ns)
		ci2.interpolate(sj,x2s)
		ci3.interpolate(sj,x3s)
	sn2 = Sampling(ns,ds,0.0)#x2s[0])
	sn3 = Sampling(ns,ds,0.0)#x3s[0])
	return x2s,x3s,sn2,xsp

"""
Gets 2D seismic image along specified curve.
"""
def getImageAlongCurve(f,x2s,x3s,sd=None):
	if sd: sz = sd
	else: sz = s1
	n1 = sz.count
	ns = len(x2s)
	g = zerofloat(n1,ns)
	si = SincInterp()
	for js in range(ns):
		x2j = x2s[js]
		x3j = x3s[js]
		for j1 in range(n1):
			x1j = sz.getValue(j1)
			g[js][j1] = si.interpolate(sz,s2,s3,f,x1j,x2j,x3j)
	return g

"""
Interpolates synthetic seismograms in 2D using image-guided tensors
"""
def interpolateSynthetics2(sy,ssy,s2x,x2c,bgs=0.5,ts=None,scat=None,sz=None):
	name = kind+"_"+str(bgs)+"_"+str(n1)+"_"+str(n2)+"_"+str(n3)+".dat"
	n1 = s1.count
	if sz:
		n1 = sz.count
	n2 = s2x.count
	nw = len(sy)
	dt = s1.delta
	ft = s1.first
	if sz:
		dt = sz.delta
		ft = sz.first
	pnull = -999.9
	g = zerofloat(n1,n2)
	gv = fillfloat(pnull,n1,n2)
	for i in range(nw):
		x2i = inro((x2c[i]-s2x.first)/s2x.delta)
		ssyi = infl((ssy[i].first-ft)/s1.delta)
		ns = len(sy[i])
		gv = wtu.syntheticImageFill(x2i,ns,ssyi,sy[i],gv);
	#plotSlice(s2x,s1,gp,title="pre-interpolated wells")
	wa,xa1,xa2 = SimpleGridder2.getGriddedSamples(pnull,s1,s2x,gv)
	if ts:
		bg = BlendedGridder2(ts,wa,xa1,xa2)
	else:
		bg = BlendedGridder2(wa,xa1,xa2)
	bg.setTimeMax(_gtmax)
	bg.setSmoothness(bgs)
	gt = bg.gridNearest(pnull,gv)
	bg.gridBlended(gt,gv,g)
	return g

"""
Interpolates depth logs in 2D using image-guided tensors
"""
def interpolateLogs2(log,slog,s2x,x2c,tz,bgs=0.5,ts=None,scat=None):
	n1 = s1.count
	n2 = s2x.count
	nw = len(log)
	dt = s1.delta
	ft = s1.first
	nt = s1.count
	pnull = -999.9
	ta = floats(s1.getValues())
	g = zerofloat(n1,n2)
	gv = fillfloat(pnull,n1,n2)
	for i in range(nw):
		x2i = inro((x2c[i]-s2x.first)/s2x.delta)
		nz = slog[i].count
		gv = wtu.logImageFill(x2i,nt,nz,pnull,log[i],tz[i],ta,gv)
	#plotSlice(s2x,s1,gp,title="pre-interpolated wells",cmap=jet)
	wa,xa1,xa2 = SimpleGridder2.getGriddedSamples(pnull,s1,s2x,gv)
	# Get distances
	if ts:
		bg = BlendedGridder2(ts,wa,xa1,xa2)
	else:
		bg = BlendedGridder2(wa,xa1,xa2)
	bg.setTimeMax(_gtmax)
	bg.setSmoothness(bgs)
	gt = bg.gridNearest(pnull,gv)
	bg.gridBlended(gt,gv,g)
	return g

def makeVot(sz,v,tz,sh):
  nw = len(sz)
  f1 = s1.first
  dt = s1.delta
  n1 = s1.count
  ta = floats(s1.getValues())
  vot = []
  for iw in range(nw):
    nt = sh[iw].count
    nz = sz[iw].count
    ft = sh[iw].first
    if1 = inro((ft-f1)/dt)
    if (nt+if1>n1): nt -= (nt+if1)-n1
    vn = zerofloat(n1)
    vn = MultipleWellTies.logImageFill(n1,nz,0,v[iw],tz[iw],ta,vn)
    vnc = copy(nt,if1,vn)
    vot.append(vnc)
  return vot


###############################################################################
###############################################################################
## 3D Multi-tie utilities

"""
Interpolates synthetic seismograms in 3D using image-guided tensors
"""
def interpolateSynthetics3(sy,ssy,x2c,x3c,bgs=0.5,ts=None,kind=None,cname=None):
	print 'Interpolate synthetic image: '+kind
	n1 = s1.count; n2 = s2.count; n3 = s3.count
	d1 = s1.delta; d2 = s2.delta; d3 = s3.delta
	f1 = s1.first; f2 = s2.first; f3 = s3.first
	nw = len(sy)  
	print 'n1=',n1,' n2=',n2,' n3=',n3
	name = kind+"_"+str(bgs)+"_"+str(n1)+"_"+str(n2)+"_"+str(n3)+".dat"
	fpath = _tpDir+"csm/interps/"
	fname = fpath+name
	if cname:
		fname = fpath+cname
	if path.exists(fname) and path.isfile(fname):
		print kind+" interp exists..."
		g = zerofloat(n1,n2,n3)
		ais = ArrayInputStream(fname)
		ais.readFloats(g)
		ais.close()
	else:
		print "make "+kind+" interps..."
		pnull = -999.9
		ta = floats(s1.getValues())
		g = zerofloat(n1,n2,n3)
		gx = zerofloat(n1,n2,n3)
		gv = fillfloat(pnull,n1,n2,n3)
		for i in range(nw):
			x2i = inro((x2c[i]-f2)/d2)
			x3i = inro((x3c[i]-f3)/d3)
			ssyi = inro((ssy[i].first-f1)/d1)
			ns = len(sy[i])
			gv = wtu.syntheticImageFill(x2i,x3i,ns,ssyi,sy[i],gv)
		wa,xa1,xa2,xa3 = SimpleGridder3.getGriddedSamples(pnull,s1,s2,s3,gv)
		print len(wa),'samples to interpolate'
		if ts:
			bg = BlendedGridder3(ts,wa,xa1,xa2,xa3)
		else:
			bg = BlendedGridder3(wa,xa1,xa2,xa3)
		bg.setSmoothness(bgs)
		bg.setTimeMax(_gtmax)
		print '  gridding nearest....'
		gt = bg.gridNearest(pnull,gv)
		print '  grid blended....'
		bg.gridBlended(gt,gv,g)
		aos = ArrayOutputStream(fname)
		aos.writeFloats(g)
		aos.close()
	return g

"""
Interpolates depth logs in 3D using image-guided tensors
"""
def interpolateLogs3(log,slog,x2c,x3c,tz,bgs=0.5,ts=None,kind=None,cname=None):
	n1 = s1.count; n2 = s2.count; n3 = s3.count
	d1 = s1.delta; d2 = s2.delta; d3 = s3.delta
	f1 = s1.first; f2 = s2.first; f3 = s3.first
	nw = len(log)  
	name = kind+"_"+str(bgs)+".dat"
	fpath = _tpDir+"csm/interps/"
	fname = fpath+name
	if cname:
		fname = fpath+cname
	if path.exists(fname) and path.isfile(fname):
		print kind+" interp exists..."
		g = zerofloat(n1,n2,n3)
		ais = ArrayInputStream(fname)
		ais.readFloats(g)
		ais.close()
	else:
		print "make "+kind+" interps..."
		pnull = -999.9
		g = zerofloat(n1,n2,n3)
		gx = zerofloat(n1,n2,n3)
		gv = fillfloat(pnull,n1,n2,n3)
		ta = floats(s1.getValues())
		for i in range(nw):
			x2i = inro((x2c[i]-f2)/d2)
			x3i = inro((x3c[i]-f3)/d3)
			nz = slog[i].count
			gv = wtu.logImageFill(x2i,x3i,n1,nz,pnull,log[i],tz[i],ta,gv)
		wa,xa1,xa2,xa3 = SimpleGridder3.getGriddedSamples(pnull,s1,s2,s3,gv)
		if ts:
			bg = BlendedGridder3(ts,wa,xa1,xa2,xa3)
		else:
			bg = BlendedGridder3(wa,xa1,xa2,xa3)
		bg.setSmoothness(bgs)
		bg.setTimeMax(_gtmax)
		print '  gridding nearest....'
		gt = bg.gridNearest(pnull,gv)
		print '   grid blended....'
		bg.gridBlended(gt,gv,g)
		aos = ArrayOutputStream(fname)
		aos.writeFloats(g)
		aos.close()
	return g


###############################################################################
###############################################################################
## Synthetic Seismograms

"""
Uses propagator matrix method to generate a synthetic seismogram
Peak frequency for ricker wavelet
Can apply a phase rotation
"""
def getPropagatorSeismogram(sz,v,d,q,\
		phase=None,normalize=True,wav=None,taper=None,\
		nomult=None,anisotropic=None,pressure=None,nocut=None):
	dt = s1.delta
	f1 = s1.first
	fz = sz.first
	dz = sz.delta
	nz = sz.count
	tzl = wtu.tzVlog(v,fz,dz)
	nt = int(ceil(tzl[nz-1]/dt))
	ft = tzl[0];
	fref = 10000*dt
	r1 = 1.0
	samp = 1.0
	st = Sampling(nt,dt,ft)
	#tzl0 = wtu.tzVlog(v,0.0,dz)
	wlth = 0
	if wav is not None: wlth = len(wav)
	fname = dataDir+'PropSeis/Prop_'+str(q)+'_'+str(fp)+'_'+str(wlth)\
								 +'_'+xstr(nomult)+'_'+xstr(anisotropic)+'_'+xstr(pressure)\
								 +'_'+str(nt)+'_'+str(ft)+'_'+str(fref)+'_'+str(r1)+'_'+str(samp)+'.dat'
	if path.exists(fname) and path.isfile(fname) and access(fname, R_OK):
		print 'Seismogram exists...'
		sy = zerofloat(nt)
		ais = ArrayInputStream(fname)
		ais.readFloats(sy)
		ais.close()
	else:
		print 'Making seismogram...'
		ssm = SeismicModel1D()
		if nomult:
			ssm.turnOffMultiples(True)
			print 'remove multiples'
		if anisotropic:
			ssm.setSourceType(SeismicModel1D.SourceType.ANISOTROPIC) # Land dyn
		else:
			ssm.setSourceType(SeismicModel1D.SourceType.ISOTROPIC) # Airgun, Vib, Marine dyn
		if pressure:
			ssm.setSensorType(SeismicModel1D.SensorType.PRESSURE) # Hydrophone
		else:
			ssm.setSensorType(SeismicModel1D.SensorType.VELOCITY) # Geophone
		ssm.setSurfaceReflectionCoefficient(r1)
		ssm.setLayers(sz,v,d,q)
		ssm.addSensor(0.0)
		ssm.addSource(0.0,samp)
		if not wav:
			ssm.setRickerWavelet(fp)
			ssm.removeSourceSignature(True)
		sy = ssm.makeSeismograms(nt,dt,fref)[0]
		if wav:
			sy = wtu.addWavelet(wav,sy)
		aos = ArrayOutputStream(fname)
		aos.writeFloats(sy)
		aos.close()
	if phase:
		sy = applyPhaseRot(phase,sy) 
	syo = sy; sto = st
	if not nocut:
		# Remove the data above ft
		#rmft = 2.0*fz/v[0]-1.0/fp
		rmft = ft
		fts = infl((rmft-f1)/dt)
		nt -= fts; 
		if nt>0:
			syt = copy(nt,fts,sy)
			ic=0
			for i in range(nt):
				if syt[i]<0.005: ic=i
				else: break
			nt -= ic
			ft += ic*dt
			sy = copy(nt,ic,syt); 
			st = Sampling(nt,dt,sto.first)
		else: print 'nt<=0'
	if taper:
		#plotCurve(sy,st,d="t",title="b4 taper")
		lt  = ft+nt*dt;  lp = lt*0.1;
		tb1 = ft+lp; 
		te1 = lt-lp;
		sy = wtu.applyDampner(st,tb1,te1,sy)
		#plotCurve(sy,st,d="t",title="tapered")
	if normalize:
		#plotCurve(sy,st,d="t",title="prop b4 norm")
		sy = normalizeRMSlocal(sy,sigmanorm)
		#plotCurve(sy,st,d="t",title="prop norm")
	if nomult: title='Prop- nm'
	else: title='Prop- m'
	#plotCurve(mul(syo,20.0),sto,sy,st,d="t",title='init prop')
	##print 'q=',q
	#print 'nprop=',nt
	return st,sy,tzl

def getPropDir():
  return dataDir+'PropSeis/'

def getCsmDir():
  return _tpDir+"csm/"

"""
Uses convolution or summation delayed by t(z) and scaled by r(z)
to make a simple synthetic seismogram
"""
def getSimpleSeismogram(sz,v,ail,phase=None,normalize=True,wav=None):
	nz = sz.count
	dz = sz.delta
	fz = sz.first
	dt = s1.delta
	rf = reflectivity(ail)
	tzl = wtu.tzVlog(v,fz,dz)
	#tzl0 = wtu.tzVlog(v,0.0,dz)
	if wav:
		sy = wtu.convolutionSeismogram(rf,tzl,wav,s1,sz)
	else:
		sy = wtu.simpleSeismogram(sz,s1,fp,rf,tzl)
	sf = Sampling(len(sy),dt,tzl[0])
	if phase:
		sy = applyPhaseRot(phase,sy) 
	if normalize:
		#plotCurve(sy,sf,d="t",title="smpl b4 norm")
		sy = normalizeRMSlocal(sy,sigmanorm)
		#plotCurve(sy,sf,d="t",title="smpl norm")
	#plotCurve(sy,sf,d="t",title='Simple')
	print 'nsimple=',sf.count
	return sf,sy,tzl


###############################################################################
###############################################################################
## Time-depth utilities

"""
Converts a seismic time image to depth using t(x,y,z)
"""
def imageTimeToDepth3(timg,tzs):
	return wtu.imageTimeToDepth3(s1,s2,s3,timg,tzs)

"""
Computes a time-depth image in 3D by:
1. resample the time-depth curves at each well to 2 m
2. extrapolate the time-depth curves to t=z=0
3. compute scattered tz image with resampled logs
4. do a isotropic tensor interpolation of the scattered samples
"""
def computetz3(kind,tzl,sz,x2w,x3w,bgs):
	nw = len(sz)
	n1 = s1.count; n2 = s2.count; n3 = s3.count
	d1 = s1.delta; d2 = s2.delta; d3 = s3.delta
	f1 = s1.first; f2 = s2.first; f3 = s3.first
	print 'n1=',n1,' n2=',n2,' n3=',n3
	name = kind+"_"+str(bgs)+"_"+str(n1)+"_"+str(n2)+"_"+str(n3)
	fpath = _tpDir+"csm/interps/"
	fname = fpath+name+".dat"
	sfname = fpath+name+"_samplings.dat"
	if path.exists(fname) and path.isfile(fname):
		print kind+" interp exists..."
		si = zerodouble(3)
		ais = ArrayInputStream(sfname)
		ais.readDoubles(si)
		ais.close()
		sd = Sampling(int(si[0]),si[1],si[2])
		tzf = zerofloat(sd.count,n2,n3)
		ais = ArrayInputStream(fname)
		ais.readFloats(tzf)
		ais.close()
	else:
		nr = 100 # regression off the end of the time depth curves
		pnull = -999.9
		mt = 0
		#lz = 2.0 #km
		lz = sz[0].last
		for i in range(nw):
			if (sz[i].last>lz):
				lz = sz[i].last
		print 'max z =',lz
		nz = ince(lz/dz2)
		sd = Sampling(nz,dz2,0.0)
		tzf = zerofloat(nz,n2,n3)
		tzs = fillfloat(pnull,nz,n2,n3)
		########### resample the time-depth curves to 2 m ############
		si = LinearInterpolator()
		atzi,aszi = [],[]
		for i in range(nw):
			# resample
			szl = sz[i]
			nzi = ince((szl.last-szl.first)/dz2)
			szi = Sampling(nzi,dz2,szl.first)
			tzi = zerofloat(nzi)
			si.setUniform(szl.count,szl.delta,szl.first,tzl[i])
			si.interpolate(szi.count,szi.delta,szi.first,tzi)
			atzi.append(tzi)
			aszi.append(szi)
			mtzi = max(tzi)
			if mtzi>mt:	
				mt = mtzi
			#SimplePlot.asPoints(szl,tzl[i])
		print 'max t =',mt
		############## convert to vavg and extrapolate #################
		#b1,b0 = getRegTerms(atzi,aszi,nr)
		for i in range(nw):
			vavg = div(mul(2.0,floats(aszi[i].getValues())),atzi[i])
			tzi = extrapolate(vavg,aszi[i],sd,mt,lz)#,b1,b0)
			# fill depth image
			x2i = inro((x2w[i]-f2)/d2)
			x3i = inro((x3w[i]-f3)/d3)
			tzs = wtu.syntheticImageFill(x2i,x3i,sd.count,0,tzi,tzs)
		# interpolate depth image
		wa,xa1,xa2,xa3 = SimpleGridder3.getGriddedSamples(pnull,s1,s2,s3,tzs)
		bg = BlendedGridder3(wa,xa1,xa2,xa3)
		bg.setTimeMax(_gtmax)
		bg.setSmoothness(bgs)
		gt = bg.gridNearest(pnull,tzs)
		bg.gridBlended(gt,tzs,tzf)
		#xa1 = floats(sd.getValues())
		#xa2 = x2w
		#bic = BicubicInterpolator2(BicubicInterpolator2.Method.SPLINE,\
		#													 BicubicInterpolator2.Method.SPLINE,xa1,xa2,tzs)
		#bic = BilinearInterpolator2(xa1,xa2,tzs)
		#bic.interpolate(sd,sn,tzf)
		aos = ArrayOutputStream(fname)
		aos.writeFloats(tzf)
		aos.close()
		aos = ArrayOutputStream(sfname)
		smp = zerodouble(3)
		smp[0]=sd.count; smp[1]=sd.delta; smp[2]=sd.first;
		aos.writeDoubles(smp)
		aos.close()
	return sd,tzf 


"""
Converts a seismic time image to depth using t(x,z)
"""
def imageTimeToDepth2(timg,tzs,sn):
	nz = len(tzs[0])
	n2 = sn.count
	dimg = zerofloat(nz,n2)
	si = SincInterp()
	snv = sn.getValues()
	for i2 in range(n2):
		for i1 in range(nz):
			dimg[i2][i1] = si.interpolate(s1,sn,timg,tzs[i2][i1],snv[i2])
	return dimg

"""
Computes a time-depth image in 2D by:
1. resample the time-depth curves at each well to 2 m
2. extrapolate the time-depth curves to t=z=0
3. compute scattered tz image with resampled logs
4. do a isotropic tensor interpolation of the scattered samples
"""
def computetz2(tzl,sz,sn,x2w,bgs,scat=None):
	nw = len(sz)
	nr = 100 # regression off the end of the time depth curves
	pnull = -999.9
	mt = 0
	#lz = 2.0 #km
	lz = sz[0].last
	for i in range(nw):
		if (sz[i].last>lz):
			lz = sz[i].last
	print 'max z =',lz
	nz = ince(lz/dz2)
	sd = Sampling(nz,dz2,0.0)
	tzf = zerofloat(nz,sn.count)
	tzs = fillfloat(pnull,nz,sn.count)
	#tzs = zerofloat(nz,nw)
	# resample the time-depth curves to 2 m
	si = LinearInterpolator()
	atzi,aszi = [],[]
	for i in range(nw):
		# resample
		szl = sz[i]
		nzi = ince((szl.last-szl.first)/dz2)
		szi = Sampling(nzi,dz2,szl.first)
		tzi = zerofloat(nzi)
		#si.interpolate(szl,tzl[i],szi,tzi)
		si.setUniform(szl.count,szl.delta,szl.first,tzl[i])
		si.interpolate(szi.count,szi.delta,szi.first,tzi)
		atzi.append(tzi)
		aszi.append(szi)
		mtzi = max(tzi)
		if mtzi>mt:	
			mt = mtzi
		#SimplePlot.asPoints(szl,tzl[i])
	print 'max t =',mt
	# extrapolate
	#b1,b0 = getRegTerms(atzi,aszi,nr)
	if scat: tsc = tzf
	for i in range(nw):
		tzi = extrapolate(atzi[i],aszi[i],sd,mt,lz)#,b1,b0)
		#SimplePlot.asPoints(sd,tzi)
		# fill depth image
		x2i = inro((x2w[i]-sn.first)/sn.delta)
		tzs = wtu.syntheticImageFill(x2i,sd.count,0,tzi,tzs)
		#tzs = wtu.syntheticImageFill(i,sd.count,0,tzi,tzs)
		if scat:
			tsc = wtu.syntheticImageFill(x2i,sd.count,0,tzi,tsc)
	if scat: return sd,tsc
	#SimplePlot.asPixels(sd,sn,tzs)
	# interpolate depth image
	wa,xa1,xa2 = SimpleGridder2.getGriddedSamples(pnull,sd,sn,tzs)
	tens = makeTensors(tzf)
	bg = BlendedGridder2(tens,wa,xa1,xa2)
	bg.setTimeMax(_gtmax)
	bg.setSmoothness(bgs)
	gt = bg.gridNearest(pnull,tzs)
	bg.gridBlended(gt,tzs,tzf)
	#xa1 = floats(sd.getValues())
	#xa2 = x2w
	#bic = BicubicInterpolator2(BicubicInterpolator2.Method.SPLINE,\
	#													 BicubicInterpolator2.Method.SPLINE,xa1,xa2,tzs)
	#bic = BilinearInterpolator2(xa1,xa2,tzs)
	#bic.interpolate(sd,sn,tzf)
	return sd,tzf 

"""
Extrapolates a function to zero and uses simple regression of the last 
specified (nr) number of samples to extrapolate to the end of the depth image
"""
def extrapolate(tz,stz,sd,mt,mz,b1,b0):
  dze = stz.delta
  fze = stz.first
  nze = inro(sd.last/dze)
  ntz = len(tz)
  nb = infl(fze/dze)
  tze = zerofloat(nze)
  # extrapolate beggining
  br0 = (tz[0]-0.3)
  br1 = 0.3/nb
  if nze<(nb+ntz): ntz -= (ntz+nb)-nze
  for i in range(nb):
    tze[i] = br0+i*br1
  # fill in known log values
  for i in range(ntz):
    tze[i+nb] = tz[i]
  # extrapolate ends using last depth ###simple regression terms
  ne = nze-(ntz+nb)
  if ne>0:
    b0 = tz[-1]
    b2 = (mt-b0)/ne 
    for i in range(ne):
      tze[i+nb+ntz] = i*b2+b0
  sze = Sampling(nze,dze,0)
  return tze,sze

"""
Solves for simple regression term for all wells by averaging the terms
in a least squares sense.
"""
def getRegTerms(tz,stz,nr):
	nw = len(stz)
	b1=0; b0=0;
	for i in range(nw):
		ntz = stz[i].count
		zvl = stz[i].getValues()
		nre = ntz-nr
		ym = meanW(tz[i],nre,ntz)	
		xm = meanW(zvl,nre,ntz)	
		num=0; den=0;
		for j in range(nre,ntz):
			xi = zvl[j]-xm
			num += (tz[i][j]-ym)*xi
			den += xi*xi
		c1 = num/den
		c0 = ym - c1*xm
		b1 += c1
		b0 += c0
	b1 /= nw
	b0 /= nw
	print 'b1 = ',b1
	print 'b0 = ',b0
	return b1,b0

"""
Computes mean in a specified sample window [w1,w2]
"""
def meanW(x,w1,w2):
	n = w2-w1
	xw = copy(n,w1,x)
	m = sum(xw)/n
	return m
		
"""
Compute TD curves and update velocity
"""
def updateTDandVelocity(sz,u,tzl,sf):
  dt = s1.delta 
  dz = sz.delta 
  nu = len(u) 
  nz = len(tzl) 
  x = rampfloat(s1.first,dt,nu)
  x2 = rampfloat(sf.first,dt,nu)
  ttau = add(x,u)
  tz = zerofloat(nz)
  # compute t(z) = t(tau=tau(z))
  ci = CubicInterpolator(\
  	CubicInterpolator.Method.LINEAR,nu,x2,ttau)
  ci.interpolate(nz,tzl,tz)
  print 'u0=',u[0]
  print 'ttau0=',ttau[0]
  print 'tz00=',tzl[0]
  #v1 = div(2,wtu.forwardDiff(tz,dz))
  #v0 = div(2,wtu.forwardDiff(tzl,dz))
  v1 = div(2,wtu.backwardsDiff(tz,dz))
  v0 = div(2,wtu.backwardsDiff(tzl,dz))
  #v1 = div(2,wtu.centerDiff(tz,dz))
  #v0 = div(2,wtu.centerDiff(tzt,dz))
  pvz = mul(div(sub(v1,v0),v0),100)
  return tz,v1,pvz 

"""
Gets z(x,t) from the t(x,z)
"""
def getzt2(sz,tz):
	#plotCurve(tz,sz,title='tz',d="z",hlabel='time (s)')
	dt = s1.delta
	fz = sz.first
	n2 = len(tz)
	tm = max(tz)
	nt = ince(tm/dt)
	szt = Sampling(nt,dt,0.0)
	zt = zerofloat(nt,n2)
	for i in range(n2):
		zt[i] = wtu.getzt(tz[i],szt,sz)
	#plotCurve(zt,sth,title='zt',d="t",hlabel='depth (km)')
	return zt,szt

"""
Gets z(t) from the log velocity and warping shifts
"""
def getzt(sz,tz):
	#plotCurve(tz,sz,title='tz',d="z",hlabel='time (s)')
	dt = s1.delta
	fz = sz.first
	ftau = tz[0]
	ntau = ince((tz[-1]-tz[0])/dt)
	stau = Sampling(ntau,dt,ftau)
	zt = wtu.getzt(tz,stau,sz)
	nth = len(zt)
	szt = Sampling(nth,dt,ftau)
	#plotCurve(zt,sth,title='zt',d="t",hlabel='depth (km)')
	return zt,szt

"""
Extrapolates t(z) to the surface where t = z = 0 
"""
def extraptz(sz,tzl,lz,lt):
	nz = sz.count
	dz = sz.delta
	nf = inro(sz.first/dz)
	nl = inro(lz/dz)-nf
	nzn = nz + nf + nl
	m1 = tzl[0]/nf
	m2 = (lt-tzl[-1])/nl
	tzn = zerofloat(nzn)
	for i in range(nf):
		tzn[i] = m1*i
	for i in range(nz):
		tzn[i+nf] = tzl[i]
	for i in range(nl):
		tzn[i+nf+nz] = tzl[-1]+m2*i
	szn = Sampling(nzn,dz,0.0)
	return tzn,szn

def getFirstZ(t,tz):
	n1 = len(tz[0][0])
	n2 = len(tz[0])
	n3 = len(tz)
	for i1 in range(n1):
		t1 = tz[n3/2][n2/2][i1]
		if t1>=t:
			z = i1
			break
	return z

def getDz():
  return dz
def getDz2():
  return dz2

###############################################################################
###############################################################################
## Seismic utilities

"""
Extract traces around inline/crossline of log
Traces (looking down z-axis) (nr=1):
  
	 3 5 8 
   2 0 7 
   1 4 6

(Well located at trace 0)
"""
def getTraces(tsdata,y,x,nr=1,tonly=None):
	yc = inro((y-s2.first)/s2.delta)
	xc = inro((x-s3.first)/s3.delta)
	#print "ntrace=",st.count
	#dev = maxdev(x,y)
	#tdev = max(int(dev*s2.delta),int(dev*s3.delta))
	#print '\nmaxdev is '+str(tdev)+' traces'
	ntr = 9+12*(nr-1)+(nr-1)*(nr-1)*4
	traces = zerofloat(s1.count,ntr)
	it = 1
	for ix in range(-nr,nr+1):
		for iy in range(-nr,nr+1):
			if ix!=0 or iy!=0:
				traces[it] = tsdata[xc+ix][yc+iy]
				it += 1
	traces[0] = tsdata[xc][yc]
	if tonly: return traces
	return traces,yc,xc

"""
Get the 3D seismic image and normalize each trace if desired
"""
def getImage(normalize=None,smooth=None,cut1=None,cut2=None,cut3=None,rnrm=None,z=None):
  tpimg,st,sy,sx = readTpstData()
  if z:
    fz,s1z,s2z,s3z = readTpszData()
  if smooth:
    p1 = 1.0
    p2 = 1.0
    sig = 4.0
    c = sig*sig*0.5
    seismicDir = getSeismictDir()
    fname = seismicDir+\
      "subt_"+str(st.count)+\
      "_"+str(int(st.delta*1000))+\
      "_"+str(int(st.first*1000))+\
      "/tpsts_"+str(p1)+"_"+str(p2)+"_"+str(sig)+".dat"
    if path.exists(fname) and path.isfile(fname):
      print 'smooth image exists'
      ais = ArrayInputStream(fname)
      ais.readFloats(tpimg)
      ais.close()
    else:
      print 'compute smooth image'
      plot3P(st,sy,sx,tpimg,"unsmooth image")
      tens = makeTensors(tpimg,p1=p1,p2=p2)
      ls = LocalSmoothingFilter()
      ls.apply(tens,c,tpimg,tpimg)
      aos = ArrayOutputStream(fname)
      aos.writeFloats(tpimg)
      aos.close()
      plot3P(st,sy,sx,tpimg,"smooth image")
  if not cut1:
    cut1=[st.first,st.last]; nn1 = st.count
  if not cut2:
    cut2=[sy.first,sy.last]; nn2 = sy.count
  if not cut3:
    cut3=[sx.first,sx.last]; nn3 = sx.count
  if cut1:
    x1fi = cut1[0]
    x1li = cut1[1]
    fn1 = inro(x1fi/st.delta)
    nn1 = inro(x1li/st.delta)-fn1+1
  if cut2:
    x2fi = cut2[0]
    x2li = cut2[1]
    fn2 = inro(x2fi/sy.delta)
    nn2 = inro(x2li/sy.delta)-fn2+1
  if cut3:
    x3fi = cut3[0]
    x3li = cut3[1]
    fn3 = inro(x3fi/sx.delta)
    nn3 = inro(x3li/sx.delta)-fn3+1
  if normalize:
    tpimgu = tpimg
    tpimg = wtu.localRMSnorm(tpimg,sigmanorm)
  if cut1 or cut2 or cut3:
    tpimg = copy(nn1,nn2,nn3,fn1,fn2,fn3,tpimg)
    if z: fz = copy(s1z.count,nn2,nn3,inro(s1z.first/s1z.delta),fn2,fn3,fz)
    if rnrm: tpimgu = copy(nn1,nn2,nn3,fn1,fn2,fn3,tpimgu)
    st = Sampling(nn1,st.delta,fn1*st.delta)
    sy = Sampling(nn2,sy.delta,fn2*sy.delta)
    sx = Sampling(nn3,sx.delta,fn3*sx.delta)
    # reset global params
    global s1,s2,s3
    s1 = st
    s2 = sy
    s3 = sx
    #plot3P(s1,s2,s3,tpimg,"cut image")
  if rnrm: 
    if z: 
      return tpimg,sx,sy,st,fz,s1z,tpimgu
    else: 
      return tpimg,sx,sy,st,tpimgu
  return tpimg,sx,sy,st

"""
Reads in seismic image 
"""
def readTpstData():
	seismicDir = getSeismictDir()
	name = "subt_"+str(s1.count)\
						+"_"+str(int(s1.delta*1000))\
						+"_"+str(int(s1.first*1000))
	fileName = seismicDir+name+"/tpst.dat"
	n1,n2,n3 = s1.count,s2.count,s3.count
	tsdata = zerofloat(n1,n2,n3)
	ais = ArrayInputStream(fileName)
	ais.readFloats(tsdata)
	ais.close()
	return tsdata,s1,s2,s3

"""
Reads in seismic image 
"""
def readTpszData():
	seismicDir = getSeismiczDir()
	name = "subz_"+str(s1z.count)\
						+"_"+str(int(s1z.delta*1000))\
						+"_"+str(int(s1z.first*1000))
	fileName = seismicDir+name+"/tpsz.dat"
	n1,n2,n3 = s1z.count,s2.count,s3.count
	fz = zerofloat(n1,n2,n3)
	ais = ArrayInputStream(fileName)
	ais.readFloats(fz)
	ais.close()
	return fz,s1z,s2,s3


"""
Computes 2D or 3D structure-oriented tensors from the seismic image
specify p2 to compute 3D tensors, which can be written to file
"""
def makeTensors(g,p0=0.0,p1=1.0,p2=None):
	lof = LocalOrientFilter(4.0)
	if p2:
		fname = _tpDir+"csm/tensors/"+"tpts_"+str(p0)+"_"+str(p1)+"_"+str(p2)+\
			      "_"+str(s3.count)+"_"+str(s2.count)+"_"+str(s1.count)
		if path.exists(fname+".dat") and path.isfile(fname+".dat"):
			print 'Tensors exist...'
			d = readTensors(fname)
		else:
			print 'Making tensors...'
			d = lof.applyForTensors(g)
			d.invertStructure(p0,p1,p2)
			writeTensors(fname,d)
		return d
	d = lof.applyForTensors(g)
	d.invertStructure(p0,p1)
	return d

"""
Computes predictability by White (1980)
"""
def predictability(h,tr,nh,sh):
	fg = zerofloat(nh)
	ff = zerofloat(nh)
	gg = zerofloat(nh)
	Conv.xcor(nh,0,tr,nh,0,h,nh,0,fg)
	Conv.xcor(nh,0,h,nh,0,h,nh,0,ff)
	Conv.xcor(nh,0,tr,nh,0,tr,nh,0,gg)
	fgfg = mul(fg,fg)
	ffgg = mul(ff,gg)
	num = sum(fgfg)
	den = sum(ffgg)
	pred = num/den
	pred/=1.214878
	return pred

def goSlopes2(f,sn):
	n2 = len(f)
	n1 = len(f[0])
	zm = ZeroMask(0.1,1.0,1.0,f)
	sigma1 = 8.0
	sigma2 = 8.0
	print "LSF: sigma1=%g, sigma2=%g"%(sigma1,sigma2)
	pmax = 2.0
	lsf = LocalSlopeFinder(sigma1,sigma2,pmax)
	p2 = zerofloat(n1,n2)
	ep = zerofloat(n1,n2)
	lsf.findSlopes(f,p2,ep)
	zero = 0.00;
	tiny = 0.01;
	zm.apply(zero,p2);
	zm.apply(tiny,ep);
	writeImage(flatDir+"p2.dat",p2)
	writeImage(flatDir+"ep.dat",ep)
	s1f = Sampling(n1,s1.delta,s1.first)
	plotSlice(sn,s1f,p2,title="P2",cbar="Slope",cmap=jet)
	plotSlice(sn,s1f,ep,title="Planarity",cbar="Planarity",cmap=jet)

#def goFlat2(f,dw):
#  ne1 = dw.getPPErrorLength()
#  s1f = Sampling(ne1)
#  s2f = Sampling(n2)
#  s3f = Sampling(n3)
#  p2 = readImage(threeDDir,"p2",ne1,n2,n3)
#  p3 = readImage(threeDDir,"p3",ne1,n2,n3)
#  ep = readImage(threeDDir,"ep",ne1,n2,n3)
#  ep = pow(ep,2)
#  fl = Flattener3()
#  fl.setWeight1(1.0)
#  fl.setIterations(0.1,1000)
#  fm = fl.getMappingsFromSlopes(s1f,s2f,s3f,p2,p3,ep)
#  ff = fm.flatten(f)
#  h = fm.unflatten(ff)
#  s = fm.getShiftsS()
#  x1 = fm.x1
#  y1 = fm.u1
#  writeImage(threeDDir,"x1",x1)
#  writeImage(threeDDir,"y1",y1)
#  writeImage(threeDDir,"ff",ff)
#  writeImage(threeDDir,"fs",s)
#  sfs = Sampling(ne1,d1,f1)
#  plotPP3(ff,title="PP flat",s2=s2,s3=s3,label1="Tau",slices=slices)
#  plotPP3(h,title="PP unflattened",s1=sfs,s2=s2,s3=s3,label1="PP time (s)",
#          slices=slices,)
#  plotPP3(s,title="PP flattening shifts",s1=sfs,s2=s2,s3=s3,slices=slices,
#          label1="PP time (s)",cbar="Shift (samples)",cmap1=jet,clips1=None)
#  print "average shift =",sum(s)/(n1*n2),"samples"

def printError(st,sh,h,tr):
	nh = len(h)
	fh = int(round(sh.first/st.delta))
	tre = copy(nh,fh,tr)
	tavg = sum(tre)/nh
	havg = sum(h)/nh
	sser=0.0;stot=0.0;htot=0.0;mtot=0.0;hrms=0.0;trms=0.0;enre=0.0;entr=0.0;
	for it in range(nh):
		trit = tre[it]
		hit = h[it]
		tmh = trit-hit
		tma = trit-tavg
		hma = hit -havg
		#sser += tmh*tmh
		#enre += tmh*tmh
		#entr += trit*trit
		stot += tma*tma
		htot += hma*hma
		mtot += tma*hma
		#hrms += hit*hit
		#trms += trit*trit
	#rsq = 1.0 - sser/stot
	prs = mtot/sqrt(stot*htot)
	#nrms = 2.0*sqrt(sser/nh)/(sqrt(hrms/nh)+sqrt(trms/nh))
	pred = predictability(h,tre,nh,sh)
	#pred = 1.0 - (enre/entr)
	#print 'R^2=',rsq
	print 'R  =',prs
	print 'PRED =',pred
	#plotSequences(tr,h,st,sh,title='f & h')
	return pred

"""
Statisctically extracts a wavelet for a number wells
"""
def extractWavelets(sz,x2w,x3w,vl):
	dt = s1.delta
	wavelets = []
	nw = len(vl)
	wvl  = 0.500 # wavelet length
	for iw in range(nw):
		# Extract a statistical wavelet
		traces = getTraces(x2w[iw],x3w[iw],tonly=True)
		wav = extractWavelet(vl[iw],wvl,sz[iw],dt,traces)
		wavelets.append(wav)
	return wavelets

"""
Extract a statistical wavelet from the traces within a window 
"""
def extractWavelet(vl,wtvl,sz,traces,noprints=None):
	dt = s1.delta
	wmin = int(2.0*(sz.first)/medianFinder(vl,30)/dt) # window min
	wmax = int(2.0*(sz.last)/medianFinder(vl,100)/dt) # window max
	wvl  = int(wtvl/dt) # wavelet length
	wav = we.getWavelet(wvl,wmin,wmax,traces)
	if not noprints:
		print 'wavelet window min= ',wmin*dt
		print 'wavelet window max= ',wmax*dt
		print 'wavelet length    = ',wvl*dt
		plotCurve(wav,Sampling(len(wav),dt,0.0),d="t",hlabel="wavelet")
	return wav

def getTpstSamplings():
	return s1,s2,s3

###############################################################################
###############################################################################
## Log utilities

"""
List of all deeps well ID's with velocity and density curves
"""
def deepWellSet():
	dws = [490251095000, 490251091800, 490251105400, 490252305400, 490251061000, 
	 			 490251090200, 490251087700, 490251104600, 490251113400, 490251116100, 
	 			 490251094400, 490251096100, 490251094600, 490251104700, 490251099200, 
	 			 490251091600, 490251097300, 490251106400]
	return dws

def getWellFromSet(x):
	_set = "d"
	dws = deepWellSet()
	_ID = dws[x]
	return _set,_ID
def getIDFromDeepSet(x):
	_set = "d"
	dws = deepWellSet()
	_ID = dws[x]
	return _ID
def getWell1():
  _set = "d"
  _ID = 490252305400
  return _set,_ID

def getUwis(uwis,wells):
  uwia = []
  for i in range(len(wells)):
    uwia.append(uwis[wells[i]])
  return uwia

"""
Gets well logs (velocity, density, gamma ray) 
for a specified tp well. 
"""
def getLogs(_set,_ID,g=None,getyv=None,smooth=None,dspk=None,cut=None):
  """
   smooth: smooth samples, sigma, cut samples
   dspk: samples for median, stddevs for spike
  """
  wdata = readLogDataset(_set) # in tputils.py
  wlog = wdata.get(_ID) # in tputils.py
  #wlog = setCoords(_ID,wlog) # deviated wells
  s1,s2,s3 = getTpstSamplings() 
  ls = LogSamples(); ls.fixLogsForSynthetics(wlog)
  fv = ls.vl; xv = ls.xl;
  fd = ls.dl; yv = ls.yl;
  fg = ls.gl; zv = ls.zl;
  if getyv:
    return yv
  #fd = gardnersRelation(fv)
  #fv,fd,zv = fixLogs2(fv,fd,zv,dz,smooth,cut)
  #plotCurve(fv,Sampling(len(fv),dz,zv[0]))
  if dspk:
    fv = despike(fv,m=dspk[0],l=dspk[1])
    fd = despike(fd,m=dspk[0],l=dspk[1])
  ail = mul(fv,fd)
  sz = Sampling(len(ail),dz,zv)
  if g:
    return ail,yv,xv,zv,fv,fd,sz,[fg,zv]
  return ail,yv,xv,zv,fv,fd,sz

"""
Prints impedance log lengths and IDs
"""
def printWellLengths():
	dws = deepWellSet()
	for i in range(len(dws)):
		_set,_ID = getWellFromSet(i) 
		ail0,yv0,xv0,zv0,vl0,dl0,sz0 = getLogs(_set,_ID)
		print i
		print _ID
		print len(ail0)

"""
Makes log lengths and zmins equal
"""
def fixLogs2(f1,f2,zf,dz,smooth=None,cuts=None):
	n = len(f1)
	# init variables for copy
	n1n = n
	n2n = n
	o1n = 0
	o2n = 0
	zfn = zf
	# fix the starts:
	if smooth:
		cutb,cute=0,0
		if cuts:
			cutb = cuts[0]
			cute = cuts[1]
		pad = smooth[0]
		sig = smooth[1]
		cut = smooth[2]
		n1n = n1n-cut-cute
		n2n = n2n-cut-cute
		o1n = o1n+cut
		o2n = o2n+cut
		ft1 = copy(pad,o1n,f1)
		ft2 = copy(pad,o2n,f2)
		ExponentialSmoother(sig).apply(ft1,ft1)
		ExponentialSmoother(sig).apply(ft2,ft2)
		copy(pad,0,ft1,o1n,f1)
		copy(pad,0,ft2,o2n,f2)
		zfn = zf+cut*dz
	f1n = zerofloat(n1n)
	f2n = zerofloat(n2n)
	copy(n1n,o1n,f1,0,f1n); copy(n2n,o2n,f2,0,f2n);
	return f1n,f2n,zfn

"""
prints the maximum distance between wells in sampling for the maximum time
option on the blended gridder.
"""
def showMaxWellDist(x2w,x3w):
	def ir(x):
		return int(round(x))
	nw = len(x2w)
	ds2,ds3 = s2.delta,s3.delta
	dmax = 0
	for i in range(nw):
		for j in range(nw):
			d = hypot(x2w[i]-x2w[j],x3w[i]-x3w[j])
			if d>dmax:
				dmax = d
				ds = hypot(ir(x2w[i]/ds2)-ir(x2w[j]/ds2),ir(x3w[i]/ds3)-ir(x3w[j]/ds3))
	print 'Maximum distance between wells is '+str(ds)+' samples'
	return ds

"""
prints the maximum distance between wells in sampling for the maximum time
option on the blended gridder.
"""
def showMinWellDist(x2w,x3w):
	def ir(x):
		return int(round(x))
	nw = len(x2w)
	ds2,ds3 = s2.delta,s3.delta
	dmin = s2.last+s3.last
	for i in range(nw):
		for j in range(nw):
			d = hypot(x2w[i]-x2w[j],x3w[i]-x3w[j])
			if d<dmin and d>0.0:
				dmin = d
				ds = hypot(ir(x2w[i]/ds2)-ir(x2w[j]/ds2),ir(x3w[i]/ds3)-ir(x3w[j]/ds3))
	print 'Minimum distance between wells is '+str(ds)+' samples'
	return ds

"""
makes points grid from g1 and g2
"""
def getGS(d1,d2,f2,g1,g2):
	ng1 = len(g1)
	ng2 = len(g2)
	n = ng1*ng2
	g1n = zerofloat(n)
	g2n = zerofloat(n)
	j=0
	for i2 in range(ng2):
		for i1 in range(ng1):
			g1n[j] = float(g1[i1])*d1
			g2n[j] = float(g2[i2])*d2+f2
			j+=1
	return [g1n,g2n]

"""
Used to over sample a signal
@param x1 the signal
@param xs1 the uniform Sampling
@param xs2 the new uniform Sampling
"""
def overSample(x1,x2,xs1,xs2):
  n1 = len(x1)
  n2 = xs2.getCount() 
  si = SincInterp()
  si.interpolate(xs1,x1,xs2,x2)
  return x2

"""
Converts x,y,z coords from kilometers to meters
"""
def km2m(x,y,z):
  x = mul(x[i],1000.0)
  y = mul(y[i],1000.0)
  z = mul(z[i],1000.0)
  return x,y,z

"""
Clips a signal based on the average to remove the outliers
"""
def despike(a,m=3,l=0.5):
  return WTutils().despike(m,l,a)

"""
Computes acoustic relectivity from an acoustic impedance log.
"""
def reflectivity(ai):
  n = len(ai)
  r = zerofloat(n)
  for i in range(1,n):
    r[i] = (ai[i-1]-ai[i])/(ai[i-1]+ai[i])
  r[0]=r[1]
  #refl = Clip(refl)
  return r

"""
Gardners relation to convert a velocity log to a density log
"""
def gardnersRelation(v):
  n = len(v)
  d = zerofloat(n)
  vft = mul(v,3.2808399) # m to ft
  d = mul(0.23,pow(v,0.25))  # constants proposed by Gardner
  return d

"""
Prints wells that contain a combination of specified logs
"""
def printWellIDList():
  #(set,type1,type2,...)
  printLogIDCombo("d","v","d","g")

"""
Returns the maximum x,y deviation distance from the surface
"""
def maxdev(x,y):
  return WTutils().maxDeviation(x,y)

"""
Fixes an input log to a reference log
"""
def fixGR(fv,fg):
  return WTutils().fixLog(fv,fg)


###############################################################################
###############################################################################
## Tops and horizon utilities

"""
Gets tp depth and time horizons in 2D form to plot on an arbitray seismic 2D slice
"""
def getXSliceHorizons(xc,sz):
	ght,ghy,csty = [],[],[]
	nf = len(fms)
	d = "t"
	for i in range(nf):
		hrz,clr = getHorzName(fms[i])
		NAME = readHorizonMod(hrz,d)
		tf,yf = wtu.horizon2d(NAME.x1,NAME.x2,NAME.x3,s3,s1,xc)
		ght.append(tf); ghy.append(yf); csty.append(clr)
	return [ght,ghy,csty]

"""
Gets tp depth and time horizons in 2D form to plot on a constant seismic 2D slice
"""
def get2DHorizons(y2s,y3s):
	ght,ghy,csty = [],[],[]
	nf = len(fms)
	d = "t"
	for i in range(nf):
		hrz,clr = getHorzName(fms[i])
		NAME = readHorizonMod(hrz,d)
		tf,yf = wtu.horizon2d(NAME.x1,NAME.x2,NAME.x3,s2,s1,y2s,y3s)
		ght.append(tf); ghy.append(yf); csty.append(clr)
	return [ght,ghy,csty]

"""
Gets tp depth and time horizons in 3D 
"""
def get3DHorizons():
	x1s,x2s,x3s,clrs = [],[],[],[]
	nf = len(fms)
	d = "t"
	for i in range(nf):
		hrz,clr = getHorzName(fms[i])
		NAME = readHorizonMod(hrz,d)
		x1s.append(NAME.x1); x2s.append(NAME.x2); 
		x3s.append(NAME.x3); clrs.append(clr);
	return [x1s,x2s,x3s,clrs]

"""
Horizon names colors and styles
"""
def getHorzName(fmn):
	sty = "-"
	fmx,clr = "",""
	if fmn=="F2WC":
		fmx="KF2F2WC"
		clr="r"+sty
	if fmn=="DKOT":
		fmx="FallRiverDKOT"
		clr="g"+sty
	if fmn=="LKOT":
		fmx="LakotaMorrison"
		clr="b"+sty
	if fmn=="CRMT":
		fmx="CrowMountainCRMT"
		clr="c"+sty
	if fmn=="RDPK":
		fmx="RedPeak"
		clr="k"+sty
	if fmn=="A Sand":
		fmx="TensleepASand"
		clr="y"+sty
	if fmn=="C1 Dolo":
		fmx="TensleepBbaseC1Dolo"
		clr="m"+sty
	return fmx,clr

"""
Get well tops from file
"""
def getWellTops(wells):
	ids,fmsp,datums,clr,i1 = [],[],[],"",0
	for i in range(len(wells)):
		s,id1 = getWellFromSet(wells[i]) 
		ids.append(id1)
	wtops = WellTops(ids,fms,_tpDir+"doe/WellLogs/")
	nw = len(ids)
	nf = len(fms)
	tops = [[] for i in range(nw)]
	for idn in ids:	
		# Get the datum
		datum = wtops.getTop(idn,"datum")
		datumdiff = seismicDatum-datum
		datums.append(datumdiff)
		# Then get the tops and apply the datum shift for the seismic
		#print id
		for fm in fms:
			fmv = wtops.getTop(idn,fm)
			#print fm
			#print fmv
			if fmv>0.0: fmv = fmv+datumdiff#-0.55
			tops[i1].append(fmv)
		i1+=1
	# Fm tops for plotting with color format
	for i in range(nf):
		fmn = fms[i]
		if fmn=="F2WC":
			clr="r-"
		if fmn=="DKOT":
		  clr="g-"
		if fmn=="LKOT":
			clr="b-"
		if fmn=="CRMT":
			clr="c-"
		if fmn=="RDPK":
			clr="k-"
		if fmn=="A Sand":
			clr="y-"
		if fmn=="C1 Dolo":
			clr="m-"
		fmsp.append((fmn,clr))
	return tops,fmsp,datums

"""
returns the indexes for time horizons and well tops for dtw constraints
"""
def constraintsForTops(tops,gh,tzl,x2w,sz,ncl):
	nc = 1
	# Choose the top/horizon to use:
	wc = 0
	nw = len(tzl)
	nf = len(tops[0])
	dt = s1.delta
	topc = zerofloat(nw)
	const = [[0,0,0] for r in range(nw)]
	for iw in range(nw):
		wci = wc
		for ic in range(nf):
			tpd = tops[iw][ic]
			if tpd>=sz[iw].first:
				wci = ic
				break
			else: tpd=0.0
		ght = gh[0][wci]
		tzl0 = sub(tzl[iw],tzl[iw][0])
		# Find horizon time sample
		x2c = x2w[iw]
		ghtm = gh[1][wci][0]
		for j in range(1,len(gh[1][wci])):
			ghtp = gh[1][wci][j]
			if (ghtm<=x2c and x2c<ghtp):
				hc = j
			ghtm = ghtp
		# Get constraints in samples for the well and the trace
		iz = sz[iw].indexOfNearest(tpd)
		const[iw][0] = int(round(ght[hc]/dt))
		const[iw][1] = int(round(tzl0[iz]/dt))
		const[iw][2] = ncl 
	return const 

"""
returns the indexes for time horizons and well tops for dtw constraints
"""
def constraintsForTops3(tops,gh3,tzl,x2w,x3w,sz,ncl):
	nc = 1
	# Choose the top/horizon to use:
	wc = 0
	nw = len(tzl)
	nf = len(tops[0])
	dt = s1.delta
	topc = zerofloat(nw)
	const = [[0,0,0] for r in range(nw)]
	for iw in range(nw):
		wci = wc
		for ic in range(nf):
			tpd = tops[iw][ic]
			if tpd>=sz[iw].first:
				wci = ic
				break
			else: tpd=0.0
		ght = gh3[0][wci]
		tzl0 = sub(tzl[iw],tzl[iw][0])
		# Find horizon time sample
		x2c = x2w[iw]
		x3c = x3w[iw]
		gh2m = gh3[1][wci][0]
		gh3m = gh3[2][wci][0]
		for j in range(1,len(gh3[1][wci])):
			gh2p = gh3[1][wci][j]
			gh3p = gh3[2][wci][j]
			if (gh2m<=x2c and x2c<gh2p and gh3m<=x3c and x3c<gh3p):
				hc = j
			gh2m = gh2p
			gh3m = gh3p
		# Get constraints in samples for the well and the trace
		iz = sz[iw].indexOfNearest(tpd)
		const[iw][0] = int(round(ght[hc]/dt))
		const[iw][1] = int(round(tzl0[iz]/dt))
		const[iw][2] = ncl
	return const

"""
prints the difference between a formation top and it's seismic time horizon
"""
def printTopDiffs(wells,gh3,tops,tz,sz,x2w,x3w,wc=0):
	nw = len(wells)	
	hc = 0
	totaldiff = 0
	for iw in range(nw):
		x2c = x2w[iw]
		x3c = x3w[iw]
		gh2m = gh3[1][wc][0]
		gh3m = gh3[2][wc][0]
		for j in range(1,len(gh3[1][0])):
			gh2p = gh3[1][wc][j]
			gh3p = gh3[2][wc][j]
			if (gh2m<=x2c and x2c<gh2p and gh3m<=x3c and x3c<gh3p):
				hc = j
			gh2m = gh2p
			gh3m = gh3p
		tpd = tops[iw][wc]
		if tpd>=sz[iw].first:
			iz = sz[iw].indexOfNearest(tpd)
			ttpd = tz[iw][iz]
			diff = gh3[0][wc][hc]-ttpd
			print "Well "+str(wells[iw])+" with a diff of "+str(diff)+" sec"
			totaldiff += abs(diff)
	print "The total differences are",totaldiff

"""
Prints the difference between certian well tops and their 
correspoding time horizons at the well. 
"""
def printperc(hv,wv):
  print 'F2WC diff= '+str(hv[0]-wv[0])
  print 'DKOT diff= '+str(hv[1]-wv[1])
  print 'LKOT diff= '+str(hv[2]-wv[2])
  print 'CRMT diff= '+str(hv[3]-wv[3])
  print 'RDPK diff= '+str(hv[4]-wv[4])
  print 'TSAS diff= '+str(hv[5]-wv[5])
  print 'TSBD diff= '+str(hv[6]-wv[6])


###############################################################################
###############################################################################
## Phase estimation 

"""
Rotates the phase of the input synthetic seismogram by a constant specified 
in degrees.
"""
def applyPhaseRot(ph,xr):
  n = len(xr)
  xim = zerofloat(n)
  y = zerofloat(n)
  hb = HilbertTransformFilter()
  hb.apply(n,xr,xim)
  phr = (FLT_PI/180)*ph 
  y = add(mul(cos(phr),xr),mul(sin(phr),xim))
  return y

"""
Finds constant phase rotation of synthetic more efficiently.
"""
def rotatePhase(sy,tr,vrmin,vrmax,dr,vb=None,ss=None,wp=None):
	n1 = len(sy)
	n2 = len(tr[0])
	nl = 1+(n2-n1)
	dw = DynamicWarpingWT(nl,vrmin-1.0,vrmax-1.0,dr)
	nh = 360
	mset = [30,5,1]
	mphr = zerofloat(nh); 
	pmrms=1; sa=0; en=nh; ihm=0;
	if vb:
		taumin = ss.first
		dt = ss.delta
		smin = -taumin*vb[0] - 1/fp
		smax =  taumin*vb[1] + 1/fp
		tmini = inro((taumin+smin)/dt)
		tmaxi = inro((ss.last+smax)/dt)
	else:
		tmini = 0
		tmaxi = n2
	if tmaxi>n2: tmaxi=n2
	if tmini<0:  tmini=0
	for mu in mset:
		for ih in range(sa,en,mu):
			#print ih
			syr = applyPhaseRot(ih,sy)
			dmin = dw.findDmin(syr,tr,tmini,tmaxi)
			mphr[ih] = dmin  
			#mi = d[n1-1].index(dmin)
			if ih==0: pmrms=mphr[ih]
			if (mphr[ih]<=pmrms): 
				ihm=ih
				pmrms = mphr[ih]
		sa = ihm-mu
		en = sa+2*mu  
		if sa<0: sa=0
		if en>nh: en=nh
	print '  Optimum Phase is  ',ihm
	return applyPhaseRot(ihm,sy)

"""
Finds constant phase rotation of multiple synthetics.
"""
def rotatePhaseMulti(sy,tr,vrmin,vrmax,dr,vb=None,ss=None):
	nw = len(sy)
	nh = 361
	mset = [1]
	mphr = zerofloat(nh); 
	pmrms=1; sa=0; en=nh; ihm=0; 
	dt = s1.delta
	f1 = s1.first
	for mu in mset:
		for ih in range(sa,en,mu):
			for iw in range(nw):
				n1 = len(sy[iw])
				n2 = len(tr[iw][0])
				nl = 1+(n2-n1)
				if vb:
					taumin = ss[iw].first
					smin = -taumin*vb[0] - 1/fp
					smax =  taumin*vb[1] + 1/fp
					tmini = inro((taumin+smin)/dt)
					tmaxi = inro((ss[iw].last+smax)/dt)
				else:
					tmini = 0
					tmaxi = n2
				if tmaxi>n2: tmaxi=n2
				if tmini<0:  tmini=0
				dw = DynamicWarpingWT(nl,vrmin-1.0,vrmax-1.0,dr)
				if ih%10==0 and iw==0: print ih
				syr = applyPhaseRot(ih,sy[iw])
					#dmin = dw.findDmin(syr,tr[iw])
				e = dw.computeErrorsMulti(syr,tr[iw])
				u = dw.findShiftsR(e)
				h = dw.applyShifts(u,syr)
				sh = Sampling(len(h),dt,u[0]*dt+f1)
					#dmin = 1.0-predictability(h,tr[iw][0],len(h),sh)
					#if dmin<0.0: dmin = 0.0
					#mphr[ih] += dmin
				#dmin = dw.findDmin(syr,tr[iw],tmini,tmaxi)
				r = correlation(h,tr[iw][0],sh,s1)
				#mphr[ih] += dmin*n1
				mphr[ih] += r/nw
				#mi = d[n1-1].index(dmin)
			if ih==0: pmrms=mphr[ih]
			if (mphr[ih]>=pmrms): 
				ihm=ih
				pmrms = mphr[ih]
		sa = ihm-mu
		en = sa+2*mu  
		if sa<0: sa=0
		if en>nh: en=nh
	print 'Optimum Phase is',ihm
	plotPhaseRMS(mphr,nh,1,slides='multiphaseplot4',paper='multiphaseplot4')
	return ihm

"""
Finds constant phase rotation of synthetic by rotating through and checking
360 degrees using dtw.
"""
def rotateAllPhase(sy,tr,vrmin,vrmax,dr):
  n1 = len(sy)
  n2 = len(tr[0])
  # Constrain n2 below n1-n2
  nl = 1+(n2-n1)
  dw = DynamicWarpingWT(nl,vrmin-1.0,vrmax-1.0,dr)
  nh = 361
  mphr = zerofloat(nh); 
  pmrms = 1; ihm=0;
  for ih in range(nh):
    #print ih
    syr = applyPhaseRot(ih,sy)
    dmin = dw.findDmin(syr,tr)
    mphr[ih] = dmin 
    #mi = d[n1-1].index(dmin)
    if ih==0: pmrms=mphr[ih]
    if (mphr[ih]<=pmrms): 
      ihm=ih
      pmrms = mphr[ih]
    #print 'min for '+str(ih)+' is '+str(mphr[ih])+' and mi is '+str(mi)+' and nd is '+str(nd)
  print 'Optimum Phase is',ihm
  plotPhaseRMS(mphr,nh,1)#,slides='sphaseplot',paper='phaseplot')
  return applyPhaseRot(ihm,sy)

"""
Finds constant phase rotation and dr.
"""
def reduceAll(sy,tr,vrmin,vrmax,dr,drstep):
	n1 = len(sy)
	n2 = len(tr[0])
	lnl = 1+(n2-n1)
	ndr = int(ceil(dr/drstep))
	dw = DynamicWarpingWT(nl,vrmin-1.0,vrmax-1.0,dr)
	nh = 360
	mset = [30,5,1]
	mphr = zerofloat(nh)
	pmrms= float('inf'); ihm=0;
	mdr=0;
	for idr in range(ndr):
		drt = dr-idr*drstep
		dw.setSmoothing(drt)
		sa=0; en=nh; iht=0;
		for mu in mset:
			for ih in range(sa,en,mu):
				syr = applyPhaseRot(ih,sy)
				dmin = dw.findDmin(syr,tr)
				mphr[ih] = dmin  
				if (mphr[ih]<=pmrms): 
					iht = ih
					pmrms = mphr[ih]
					ihm=ih
					mdr = drt
			sa = iht-mu
			en = sa+2*mu  
			if sa<0: sa=0
			if en>nh: en=nh
	print 'minimum d is', pmrms
	print 'Optimum Phase is',ihm
	print 'Optimum dr    is',mdr
	return applyPhaseRot(ihm,sy),mdr


###############################################################################
#############################################################################
# General utilities

"""
Computes nrms for a 1D, 2D, or 3D sequences/images
"""
def nrms(x1,x2):
	return wtu.computeNrms(x1,x2)

"""
writes out image 
"""
def writeImage(x,fpath):
	aos = ArrayOutputStream(fpath)	
	aos.writeFloats(x)
	aos.close()

"""
reads in image 
"""
def readImage(x,fpath):
	ais = ArrayInputStream(fpath)	
	ais.readFloats(x)
	ais.close()
	return x

"""
Converts 1D float array to 1D double array
"""
def double(x):
  xd = zerodouble(len(x))
  xd = add(x,0.0)
  return xd

"""
Converts 1D array to 1D float array
"""
def floats(x):
	n = len(x)
	xd = zerofloat(n)
	for i in range(n):
		xd[i] = float(x[i])
	return xd

"""
Converts 1D double array to 1D float array
"""
def afloat(x):
  xf = zerofloat(len(x))
  xf = add(x,0.)
  return xf

"""
Find the median of the top nmx values
"""
def medianFinder(x,nmx):
	xm = zerofloat(nmx);
	copy(nmx,x,xm);
	mf = MedianFinder(nmx);
	md = mf.findMedian(xm);
	return md

"""
Returns the mean of a 1D array
"""
def mean(a):
	return sum(a)/float(len(a))

"""
Sorts array a1 and array a2 by a2
"""
def sortarrays(a1,a2):
  n = len(a1)
  b1 = zerofloat(n)
  x = rampint(0,1,n)
  quickIndexSort(a2,x)
  quickSort(a2) 
  wtu.asort(a1,x,b1)
  return b1,a2

"""
Computes and returns pearson's correlation coefficient
"""
def correlation(x,y,cx=None,cy=None):
	nx = len(x)
	ny = len(y)
	if not cx: cx = Sampling(nx,1.0,0.0)
	if not cy: cy = Sampling(ny,1.0,0.0)
	# x must have first sample greater than y and x must be shorter than y
	if nx>ny:
		print 'switched!'
		nt = nx
		nx = ny
		ny = nt
		ct = cx
		cx = cy
		cy = ct
	d = cx.delta;
	f1 = cx.first 
	f2 = cy.first 
	f1i = infl((f1-f2)/d)
	def cc(x,y):
		mx = mean(x)
		my = mean(y)
		sx = sub(x,mx)
		sy = sub(y,my)
		sxsx = mul(sx,sx)
		sysy = mul(sy,sy)
		return sum(mul(sx,sy))/sqrt(sum(sxsx)*sum(sysy))
	return cc(copy(nx,0,x),copy(nx,f1i,y))

"""
Prints correlations for sets of wells and traces
"""
def printCorrelations(f,x,cx,x2,x3,which=''):
	nw = len(x)
	for i in range(nw):
		y = getTraces(f,x2[i],x3[i],tonly=True)
		r = correlation(x[i],y[0],cx[i],s1)
		print which,'R =',r

"""
Trims an array and sampling
"""
def trim(sx,x,v=0.0):
	ni = len(x)
	c1 = 0
	c2 = 0
	st = False
	for i in range(ni):
		if x[i]>v:
			if not st:
				st = True
		else:
			if st:
				c2+=1
			else:
				c1+=1
	nf = ni-c1-c2
	y = zerofloat(nf)
	for i in range(nf):
		y[i] = x[i+c1]
	sy = Sampling(nf,sx.delta,sx.first+c1*sx.delta)
	return sy,y

"""
Checks if the type is None and returns a string of the input
"""
def xstr(s):
	if s is None or False:
		return 'N'
	if s is True:
		return 'T'
	return str(s)

"""
normalizes max to 1
"""
def normalize(e):
  emin = min(e)
  emax = max(e)
  return mul(sub(e,emin),1.0/(emax-emin))

"""
Local RMS normalization.
See java method documentation. 
"""
def normalizeRMSlocal(x,sig):
  return WTutils().localRMSnorm(x,sig)

"""
Global RMS normalization
"""
def normalizeRMS(x):
  n = len(x)
  rms = sqrt(sum(mul(x,x))/n)
  div(x,rms)
  return x

"""
Global normalization
"""
def globNorm(a):
  amax = max(a)
  div(a,amax)
  return a

def inro(x):
	return int(round(x))

def ince(x):
	return int(ceil(x))

def infl(x):
	return int(floor(x))

def etran(e):
  #e = normalize(e)
  return transpose(pow(e,0.25))

def mtran(e):
  return transpose(e)

def isEven(number):
	return number % 2 == 0

def getTimePrint(t):
	# seconds
	if t<=60.0:
		return str(t)+" seconds"
	# minutes
	if t<=(60.0*60.0):
		return str(t/60.0)+" minutes"
	# hours
	if t<=(60.0*60.0*24.0):
		return str(t/(60.0*60.0))+" hours"
	# else
	return str(t/(60.0*60.0*24.0))+" days"

def getTpDir():
	return _tpDir

def like3f(x):
	return zerofloat(len(x[0][0]),len(x[0]),len(x))
def like2f(x):
	return zerofloat(len(x[0]),len(x))
def like1f(x):
	return zerofloat(len(x))

"""
Plots the amplitude spectrum of the input sequence using
the sampling, the values, and the title of the plot (s,x,l)
"""
def plotSpectrum(s,x,l):
  Spectrum().apply(s,x,l)

"""
3D exponential smoother in (x2,x3) lateral dims
niter=16 approximates gaussian
"""
def smooth23(x,s=5,niter=16):
	y = like3f(x)
	RecursiveExponentialFilter(s).apply2(x,y)
	RecursiveExponentialFilter(s).apply3(y,y)
	for i in range(niter-1):
		RecursiveExponentialFilter(s).apply2(y,y)
		RecursiveExponentialFilter(s).apply3(y,y)
	return y

def smooth1(x,s=5,niter=16):
	y = like1f(x)
	RecursiveExponentialFilter(s).apply1(x,y)
	for i in range(niter-1):
		RecursiveExponentialFilter(s).apply1(y,y)
	return y

def smooth2(x,s=5,niter=16):
	y = like2f(x)
	RecursiveExponentialFilter(s).apply2(x,y)
	for i in range(niter-1):
		RecursiveExponentialFilter(s).apply2(y,y)
	return y




###############################################################################
#############################################################################
# Plots

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
xrxu = PlotPanel.Orientation.X1RIGHT_X2UP
xdxr = PlotPanel.Orientation.X1DOWN_X2RIGHT
aplt = PlotPanel.AxesPlacement.LEFT_TOP
inear = PixelsView.Interpolation.NEAREST
eoc = PlotFrame.EXIT_ON_CLOSE

def plotMatrix(c,sf,slag,u=None,cp=None,lim=None,title=None,png=None,paper=None,slides=None):
  n1,nlag = len(c[0]),len(c)
  #slag = Sampling(nlag,1.0,-(nlag-1)/2)
  panel = PlotPanel(xrxu)
  cv = panel.addPixels(sf,slag,c)
  cv.setInterpolation(inear)
  cv.setColorModel(jet)
  panel.setVLimits(slag.first,slag.last)
  if u:
    uv = panel.addPoints(sf,u)
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
    panel.setVLimits(min(u)*0.5,max(u)*1.5)
  if cp:
    cv.setPercentiles(0,cp)
  if lim:
  	panel.setVLimits(lim[0],lim[1])
  panel.setHLimits(sf.first,sf.last)
  panel.setHLabel("time (s)")
  panel.setVLabel("time lag (s)")
  if title:
    panel.setTitle(title)
  panel.addColorBar()
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(eoc)
  frame.setFontSize(18)
  frame.setSize(1000,800)
  frame.setVisible(True)
  if png:
    if _prints: frame.paintToPng(720,3.33,_pngDir+png+'.png')
  if paper:
    frame.setSize(732,641)
    frame.setFontSizeForPrint(8.0,240.0)
    if _prints: frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    frame.setSize(550,500)
    frame.setFontSizeForSlide(1.0,1.0)
    if _prints: frame.paintToPng(720,3.50,_pngDir+slides+'.png')
  
def plotSequences(f,g,sx=None,sy=None,u=None,title=None):
  n1 = len(f)
  n2 = len(g)
  if sx is None:
    sx = Sampling(n1,1.0,0.0)
  if sy is None:
    sy = Sampling(n2,1.0,0.0)
  panel = PlotPanel(xdxr)
  fv = panel.addPoints(sx,f)
  gv = panel.addPoints(sy,g)
  if u:
    gv = panel.addPoints(u,g)
  gv.setLineColor(Color.RED)
  panel.setVLabel("time (s)")
  if title:
    panel.setTitle(title)
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(eoc)
  frame.setFontSize(18)
  frame.setSize(300,800)
  frame.setVisible(True)

"""
Plots a panel of well logs
"""
def plotLogPanel(sz,v,d,rf,paper=None,slides=None):
  p1 = PlotPanel(1,3,xdxr)
  p1.setVLabel("Depth (km)")
  vp = p1.addPoints(0,0,sz,v)
  dp = p1.addPoints(0,1,sz,d)
  rf = p1.addPoints(0,2,sz,rf)
  p1.setHLabel(0,"Velocity (km/s)")
  p1.setHLabel(1,"Density (g/cc)")
  p1.setHLabel(2,"Reflectivity")
  p1.setHLimits(0,min(v),max(v))
  p1.setHLimits(1,1.51,max(d))
  p1.setHInterval(1,0.5)
  p1.setHLimits(2,-0.15,0.15)
  p1.setHInterval(2,0.1)
  #dp.setLineColor(BLUE)
  #gp.setLineColor(GREEN)
  frame = PlotFrame(p1)
  #frame.setSize(800,1500)
  #frame.setFontSize(28)
  #frame.setDefaultCloseOperation(eoc);
  if paper:
    frame.setSize(470,643)
    frame.setFontSizeForPrint(8.0,222.0)
    if _prints: frame.paintToPng(720,3.08,_pngDir+paper+'.png')
  if slides:
    frame.setSize(690,561)
    frame.setFontSizeForSlide(1.0,1.0)
    if _prints: frame.paintToPng(720,6.00,_pngDir+slides+'.png')
  frame.setVisible(True)

"""
Plots a panel of 3-5 seismograms
"""
def plotSeismogramPanel(ssy1,ssy2,ssy3,sy1,sy2,sy3,ssy4=None,ssy5=None,sy4=None,sy5=None,hlim1=None,hlim2=None,hlim3=None,hlim4=None,hlim5=None,vlim=None,rm3=None,paper=None,slides=None,size=None):
  i=0
  if rm3:
    e = 2
  else:
    e = 3
  if ssy4 and sy4: i+=1
  if ssy5 and sy5: i+=1
  p1 = PlotPanel(1,e+i,xdxr)
  p1.setVLabel("Time (s)")
  p1.addPoints(0,0,ssy1,sy1)
  p1.addPoints(0,1,ssy2,sy2)
  if not rm3: p1.addPoints(0,2,ssy3,sy3)
  if ssy4 and sy4: p1.addPoints(0,3,ssy4,sy4)
  if ssy5 and sy5: p1.addPoints(0,4,ssy5,sy5)
  #p1.setHLabel(0,"simple")
  #p1.setHLabel(1,"multiples")
  #p1.setHLabel(2,"multiples and Q")
  if hlim1: p1.setHLimits(0,hlim1[0],hlim1[1])
  if hlim2: p1.setHLimits(1,hlim2[0],hlim2[1])
  if hlim3: p1.setHLimits(2,hlim3[0],hlim3[1])
  if hlim4: p1.setHLimits(3,hlim4[0],hlim4[1])
  if hlim5: p1.setHLimits(4,hlim5[0],hlim5[1])
  if vlim:  p1.setVLimits(vlim[0],vlim[1])
  #dp.setLineColor(BLUE)
  #gp.setLineColor(GREEN)
  frame = PlotFrame(p1)
  frame.setSize(800,1200)
  #frame.setFontSize(28)
  frame.setDefaultCloseOperation(eoc);
  if paper:
    frame.setSize(515,650)
    if size: frame.setSize(size[0],size[1])
    w = 3.33
    if geo: 
      frame.setFontSizeForPrint(8.0,240.0)
      w = 3.33
    elif cwp: 
      frame.setFontSizeForPrint(8.0,222.0)
      w = 3.08
    #override to fit axis numbers
    if not rm3: frame.setFontSizeForPrint(8.0,270.0)
    if _prints: frame.paintToPng(720,w,_pngDir+paper+'.png')
  if slides:
    frame.setSize(690,561+187*i)
    frame.setFontSizeForSlide(1.0,1.0)
    if _prints: frame.paintToPng(720,6.00,_pngDir+slides+'.png')
  frame.setVisible(True)

"""
Plots a panel of 3 velocities
"""
def plotVelPanel(sz,sy1,sy2,sy3,sy12=None,sy22=None,sy32=None,avel=None,pvel=None,hlim=None,title=None,paper=None,slides=None):
	p1 = PlotPanel(1,3,xdxr)
	p1.setVLabel("Depth (km)")
	vp = p1.addPoints(0,0,sz,sy1)
	dp = p1.addPoints(0,1,sz,sy2)
	rf = p1.addPoints(0,2,sz,sy3)
	if sy12:
		l1 = p1.addPoints(0,0,sz,sy12)
		l1.setLineColor(RED)
	if sy22:
		l2 = p1.addPoints(0,1,sz,sy22)
		l2.setLineColor(RED)
	if sy32:
		l3 = p1.addPoints(0,2,sz,sy32)
		l3.setLineColor(RED)
	if avel:
		p1.setHLabel(0,"Velocity (km/s)")
		p1.setHLabel(1,"Velocity (km/s)")
		p1.setHLabel(2,"Velocity (km/s)")
	elif pvel:
		p1.setHLabel(0,"% difference")
		p1.setHLabel(1,"% difference")
		p1.setHLabel(2,"% difference")
	else:
		p1.setHLabel(0,"Time (s)")
		p1.setHLabel(1,"Velocity (km/s)")
		p1.setHLabel(2,"% difference")
	if hlim:
		p1.setHLimits(2,hlim[0],hlim[1])
	#p1.setHLabel(0,"Velocity (km/s)")
	#p1.setHLabel(1,"Velocity (km/s)")
	#p1.setHLabel(2,"Velocity (km/s)")
	#p1.setHLimits(0,-0.51,0.51)
	#p1.setHLimits(1,-0.29,0.20)
	#p1.setHLimits(2,-0.51,0.51)
	#dp.setLineColor(BLUE)
	#gp.setLineColor(GREEN)
	if title: 
		p1.setTitle(title)
	frame = PlotFrame(p1)
	frame.setSize(800,1200)
	frame.setFontSize(28)
	frame.setDefaultCloseOperation(eoc);
	if paper:
		frame.setSize(711,600)
		frame.setFontSizeForPrint(8.0,240.0)
		if _prints: frame.paintToPng(720,3.33,_pngDir+paper+'.png')
	if slides:
		frame.setSize(711,600)
		frame.setFontSizeForSlide(1.0,1.0)
		if _prints: frame.paintToPng(720,6.00,_pngDir+slides+'.png')
	frame.setVisible(True)
	return frame


"""
Plots 3 curves on top of eachother
"""
def plotCurvePanel3(sf,sg,f,g,sh=None,h=None,tlim=None,paper=None,slides=None):
	if h:
		p1 = PlotPanel(3,1,xrxu,aplt)
		pf = p1.addPoints(0,0,sf,f)
		pg = p1.addPoints(1,0,sg,g)
		pg2f = p1.addPoints(0,0,sh,h)
		pg2g = p1.addPoints(2,0,sh,h)
		pg2f.setLineColor(RED)
		pg2g.setLineColor(RED)
		#p1.setVLabel(0,"f & gw")
		#p1.setVLabel(1,"g & gw")
	else: 
		p1 = PlotPanel(2,1,xrxu,aplt)
		pf = p1.addPoints(0,0,sf,f)
		pg = p1.addPoints(1,0,sg,g)
	p1.setHLabel(0,"Time (s)")
	#p1.setVLabel(0,"f")
	#p1.setVLabel(1,"g")
	if tlim:
		p1.setHLimits(0,tlim)
	else: p1.setHLimits(0,sf.last)
	frame = PlotFrame(p1)
	frame.setDefaultCloseOperation(eoc)
	frame.setFontSize(28)
	frame.setSize(1500,350)
	if paper:
		frame.setSize(1100,413)
		frame.setFontSizeForPrint(8.0,240.0)
		if _prints: frame.paintToPng(720,3.33,_pngDir+paper+'.png')
	if slides:
		#frame.setSize(1000,500)
		frame.setSize(1000,350)
		frame.setFontSizeForSlide(1.0,1.0)
		if _prints: frame.paintToPng(720,6.00,_pngDir+slides+'.png')
		#if g2:
		#  frame.setSize(965,455)
		#else: frame.setSize(965,349)
		#frame.setFontSizeForSlide(1.0,1.0)
		#frame.paintToPng(720,6.00,_pngDir+slides+'.png')
	frame.setVisible(True)

def plotCurvePanel2(sf,sg,sh,f,g,h,lims=None,paper=None,slides=None):
	p1 = PlotPanel(2,1,xrxu,aplt)
	p1.addPoints(0,0,sg,g)
	p1.addPoints(1,0,sg,g)
	pf = p1.addPoints(0,0,sf,f)
	ph = p1.addPoints(1,0,sh,h)
	pf.setLineColor(RED)
	ph.setLineColor(RED)
	#p1.setVLabel(0,"f & g")
	#p1.setVLabel(1,"f & h")
	p1.setHLabel(0,"Time (s)")
	if lims:
		p1.setHLimits(lims[0],lims[1])
	frame = PlotFrame(p1)
	frame.setDefaultCloseOperation(eoc)
	frame.setFontSize(28)
	frame.setSize(1500,350)
	if paper:
		frame.setSize(712,418)
		frame.setFontSizeForPrint(8.0,222.0)
		if _prints: frame.paintToPng(720,3.08,_pngDir+paper+'.png')
	if slides:
		#frame.setSize(1000,500)
		frame.setSize(1000,350)
		frame.setFontSizeForSlide(1.0,1.0)
		if _prints: frame.paintToPng(720,6.00,_pngDir+slides+'.png')
		#if g2:
		#  frame.setSize(965,455)
		#else: frame.setSize(965,349)
		#frame.setFontSizeForSlide(1.0,1.0)
		#frame.paintToPng(720,6.00,_pngDir+slides+'.png')
	frame.setVisible(True)

"""
Plots comparison of reflectvity to its resulting synthetic seismogram
"""
def plotCurveW(sz,rf,ss,sy,paper=None):
  p1 = PlotPanel(1,1,xdxr)
  p2 = PlotPanel(1,1,xdxr)
  p1.addPoints(0,0,sz,rf)
  p2.addPoints(0,0,ss,sy)
  p1.setVLabel(0,"Depth (km)")
  p2.setVLabel(0,"t (s)")
  p1.setHLimits(0,-0.13,0.13)
  p2.setHLimits(0,-0.5,0.5)
  frame = PlotFrame(p1,p2,PlotFrame.Split.HORIZONTAL)
  frame.setSize(300,1500)
  frame.setDefaultCloseOperation(eoc);
  frame.setVisible(True)
  if paper:
    frame.setSize(584,757)
    frame.setFontSizeForPrint(8.0,240.0)
    if _prints: frame.paintToPng(720,3.33,_pngDir+paper+'.png')

"""
Plots comparison of synthetic seismogram before and after normalization
"""
def plotCurveN(st,tr,trn,ss,sy,syn,lim=None,paper=None):
  p1 = PlotPanel(1,2,xdxr)
  p2 = PlotPanel(1,2,xdxr)
  p1.addPoints(0,0,st,tr)
  p1.addPoints(0,1,st,trn)
  p2.addPoints(0,0,ss,sy)
  p2.addPoints(0,1,ss,syn)
  p1.setVLabel(0,"Time (s)")
  p2.setVLabel(0,"t (s)")
  p2.setHLimits(0,-1.1,1.1)
  frame = PlotFrame(p1,p2,PlotFrame.Split.HORIZONTAL)
  frame.setSize(300,1500)
  frame.setDefaultCloseOperation(eoc);
  frame.setVisible(True)
  if paper:
    frame.setSize(584,757)
    frame.setFontSizeForPrint(8.0,240.0)
    if _prints: frame.paintToPng(720,3.33,_pngDir+paper+'.png')

"""
Plots time-depth curves before and after warping
"""
def plotTDCurves(sz,tauz,tz,tzn=None,paper=None,slides=None):
  p1 = PlotPanel(xdxr)
  l1 = p1.addPoints(sz,tz)
  l2 = p1.addPoints(sz,tauz)
  l1.setLineColor(RED)
  if tzn:
    l3 = p1.addPoints(sz,tzn)
    l3.setLineColor(RED)
    l1.setLineColor(BLUE)
  p1.setHLabel('Time (s)')
  p1.setVLabel('Depth (km)')
  frame = PlotFrame(p1)
  frame.setSize(700,500)
  frame.setFontSizeForPrint(8.0,240.0)
  frame.setDefaultCloseOperation(eoc);
  if paper:
		#colm = 222.0
    colm = 148.0
    frame.setSize(513,379) 
    frame.setFontSizeForPrint(8.0,colm)
    if _prints: frame.paintToPng(720,colm/72.0,_pngDir+paper+'.png')
  if slides:
    frame.setSize(1000,500)
    frame.setFontSizeForSlide(1.0,1.0)
    if _prints: frame.paintToPng(720,6.00,_pngDir+slides+'.png')
    #frame.setSize(513,379)
    #frame.setFontSizeForSlide(1.0,1.0)
    #frame.paintToPng(720,6.00,_pngDir+slides+'.png')
  frame.setVisible(True)

def plotTDandV(sz,v1,v2,tauz,tz,h2label=None,paper=None,slides=None):
	p1 = PlotPanel(1,2,xdxr)  
	p1.mosaic.setWidthElastic(0,67)
	p1.mosaic.setWidthElastic(1,33)
	l1 = p1.addPoints(0,0,sz,tz)                              
	l2 = p1.addPoints(0,0,sz,tauz)                          	
	l1.setLineColor(RED)                                      
	p1.setHLabel(0,'Time (s)')                                
	p1.setVLabel(0,'Depth (km)')                            	
	if h2label:                                             	
		p1.setHLabel(1,h2label)                               	
	l1 = p1.addPoints(0,1,sz,v1)                            	
	l1.setLineColor(BLACK)                                  	
	l2 = p1.addPoints(0,1,sz,v2)                            	
	l2.setLineColor(RED)                                    	
	frame = PlotFrame(p1)                                   	
	frame.setSize(1000,1000)	                              	
	if slides:                                              	
		frame.setSize(800,500)                                	
		frame.setFontSizeForSlide(1.0,1.0)                    	
		if _prints: frame.paintToPng(720,6.00,_pngDir+slides+'.png')      	
	if paper:                                               	
		colm = 222.0                                          	
		frame.setSize(712,524)                                	
		frame.setFontSizeForPrint(8.0,colm)                   	
		if _prints: frame.paintToPng(720,colm/72.0,_pngDir+paper+'.png')  	
	frame.setVisible(True)                                  	
	frame.setDefaultCloseOperation(eoc)

def plotAllTDCurves(sz,tz,flz=None,paper=None,slides=None):
  p1 = PlotPanel(xdxr)  
  lz = 0
  nw = len(sz)
  for i in range(nw):
    p1.addPoints(sz[i],tz[i])                              
    z = sz[i].last
    if z>lz: lz = z
  p1.setHLabel('Time (s)')                                
  p1.setVLabel('Depth (km)')                            	
  p1.setVLimits(0,lz)
  if flz: p1.setVLimits(0,flz)
  p1.setHLimits(0,s1.last)
  frame = PlotFrame(p1)                                   	
  frame.setSize(1000,1000)	                              	
  if slides:                                              	
    frame.setSize(800,800)                                	
    frame.setFontSizeForSlide(1.0,1.0)                    	
    if _prints: frame.paintToPng(720,6.00,_pngDir+slides+'.png')      	
  if paper:                                               	
    colm = 222.0                                          	
    frame.setSize(712,524)                                	
    frame.setFontSizeForPrint(8.0,colm)                   	
    if _prints: frame.paintToPng(720,colm/72.0,_pngDir+paper+'.png')  	
  frame.setVisible(True)                                  	
  frame.setDefaultCloseOperation(eoc)

def plotAllVavgCurves(sz,v,flz=None,paper=None,slides=None):
  p1 = PlotPanel(xdxr)  
  lz = 0
  nw = len(sz)
  for i in range(nw):
    p1.addPoints(sz[i],v[i])                              
    z = sz[i].last
    if z>lz: lz = z
  p1.setHLabel('Velocity (km/s)')                                
  p1.setVLabel('Depth (km)')                            	
  p1.setVLimits(0,lz)
  if flz: p1.setVLimits(0,flz)
  p1.setHLimits(2.0,4.5)
  frame = PlotFrame(p1)                                   	
  frame.setSize(1000,1000)	                              	
  if slides:                                              	
    frame.setSize(800,800)                                	
    frame.setFontSizeForSlide(1.0,1.0)                    	
    if _prints: frame.paintToPng(720,6.00,_pngDir+slides+'.png')      	
  if paper:                                               	
    colm = 222.0                                          	
    frame.setSize(712,524)                                	
    frame.setFontSizeForPrint(8.0,colm)                   	
    if _prints: frame.paintToPng(720,colm/72.0,_pngDir+paper+'.png')  	
  frame.setVisible(True)                                  	
  frame.setDefaultCloseOperation(eoc)

"""
Plots a curve horizontally and optionally another as a red in the same panel
"""
def plotCurveHorz(c1,s1,c2=None,s2=None,title=None,lim=None,\
		hlabel=None,png=None,d=None,slides=None,paper=None):
	p1 = PlotPanel(xrxu)
	l1 = p1.addPoints(s1,c1)
	l1.setLineColor(BLACK)
	if c2 and s2:
		l2 = p1.addPoints(s2,c2)
		l2.setLineColor(RED)
	if d=="z":
		p1.setHLabel("Depth (km)")
	if d=="t":
		p1.setHLabel("Time (s)")
	if lim:
		p1.setVLimits(lim[0],lim[1])
	if title: 
		p1.setTitle(title)
	if hlabel:
		p1.setVLabel(hlabel)
	frame = PlotFrame(p1)
	frame.setSize(1500,300)
	if png:
		if _prints: frame.paintToPng(720,3.33,_pngDir+png+'.png') 
	if slides:
		frame.setSize(800,300)
		frame.setFontSizeForSlide(1.0,1.0)
		if _prints: frame.paintToPng(720,6.00,_pngDir+slides+'.png')
	if paper:
		#colm = 90.0
		colm=222.0
		frame.setSize(581,266)
		frame.setFontSizeForPrint(8.0,colm)
		if _prints: frame.paintToPng(720,colm/72.0,_pngDir+paper+'.png')
	frame.setVisible(True)
	frame.setDefaultCloseOperation(eoc);

"""
Plots 4 curve horizontally in two panels 
"""
def plotCurveHorz2(c1,s1,c2,s2,c3,s3,c4,s4,title=None,\
		lim=None,vlabel=None,hlabel=None,paper=None,slides=None):
  p1 = PlotPanel(2,1,xrxu,aplt)
  l1 = p1.addPoints(0,0,s1,c1)
  l1.setLineColor(BLACK)
  l2 = p1.addPoints(0,0,s2,c2)
  l2.setLineColor(RED)
  l3 = p1.addPoints(1,0,s3,c3)
  l3.setLineColor(BLACK)
  l4 = p1.addPoints(1,0,s4,c4)
  l4.setLineColor(RED)
  if lim:
    p1.setVLimits(lim[0],lim[1])
  if title: 
    p1.setTitle(title)
  if hlabel:
    p1.setVLabel(hlabel)
  if vlabel:
    p1.setHLabel(vlabel)
  else:
    p1.setHLabel("Time (s)")
  frame = PlotFrame(p1)
  frame.setSize(700,600)
  if paper:
 	  frame.setSize(711,440)
 	  frame.setFontSizeForPrint(8.0,240.0)
 	  if _prints: frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    frame.setSize(711,500)
    frame.setFontSizeForSlide(1.0,1.0)
    if _prints: frame.paintToPng(720,6.00,_pngDir+slides+'.png')
  frame.setVisible(True)
  frame.setDefaultCloseOperation(eoc);
  return frame

"""
Plots 2 PlotPanels in the same frame
"""
def plotMultiPanels(p1,p2,horz=None,paper=None,slides=None):
	pp1 = p1.getPlotPanel()
	pp2 = p2.getPlotPanel()
	frame = PlotFrame(pp1,pp2,PlotFrame.Split.VERTICAL)
	if horz:
		frame = PlotFrame(pp1,pp2,PlotFrame.Split.HORIZONTAL)
	se1 = p1.getSize()
	se2 = p2.getSize()
	width1,width2  = int(se1.getWidth()),int(se2.getWidth())
	height1,height2 =  int(se1.getHeight()),int(se2.getHeight())
	frame.setSize(width1,width2,height1,height2)
	if paper:
		frame.setFontSizeForPrint(8.0,240.0)
		if _prints: frame.paintToPng(720,3.33,_pngDir+paper+'.png')
	if slides:
		frame.setFontSizeForSlide(1.0,1.0)
		if _prints: frame.paintToPng(720,7.00,_pngDir+slides+'.png')
	frame.setVisible(True)
	frame.setDefaultCloseOperation(eoc);

"""
Plots a curve vertically and optionally another as a red in the same panel
"""
def plotCurve(c1,s1,c2=None,s2=None,title=None,lim=None,hlabel=None,png=None,\
		d=None,slides=None,paper=None,psize=None,pfill=None,nfill=None,vlim=None):
	p1 = PlotPanel(xdxr)
	l1 = p1.addPoints(s1,c1)
	l1.setLineColor(BLACK)
	if c2 and s2:
		l2 = p1.addPoints(s2,c2)
		l2.setLineColor(RED)
	if pfill:
		l1.setPositiveFill(pfill)
	if nfill:
		l1.setNegativeFill(nfill)
	p1.setVLabel("Time (s)")
	if d=="z":
		p1.setVLabel("Depth (km)")
	if d=="t":
		p1.setVLabel("Time (s)")
	if lim:
		p1.setHLimits(-lim,lim)
	if vlim:
		p1.setVLimits(vlim[0],vlim[1])
	if title: 
		p1.setTitle(title)
	if hlabel:
		p1.setHLabel(hlabel)
	frame = PlotFrame(p1)
	frame.setSize(300,1500)
	if png:
		if _prints: frame.paintToPng(720,3.33,_pngDir+png+'.png') 
	if slides:
		frame.setSize(300,800)
		frame.setFontSizeForSlide(1.0,1.0)
		if _prints: frame.paintToPng(720,6.00,_pngDir+slides+'.png')
	if paper:
		colm = 90.0
		frame.setSize(200,379)
		frame.setFontSizeForPrint(8.0,colm)
		if _prints: frame.paintToPng(720,colm/72.0,_pngDir+paper+'.png')
	frame.setVisible(True)
	frame.setDefaultCloseOperation(eoc);


"""
Plots a 2-D seismic slice, and well logs through that slice, and seismic horizons on the slice
"""
def plot2DLogs(s1,s2,img,sh,h,x2w,wwidth=3,tops=None,fmsp=None,tz=None,sz=None,gh=None,\
    title=None,nolog=None,paper=None,slides=None,cmap=gray,lims=None,cbar=None,\
		psize=None,lcmap=gray,clips=None,lclips=None,rmb=None,rml=None,velc=None):
  nw = len(sh)
  g = copy(img)
  d2 = s2.delta
  if isEven(wwidth): wwidth -= 1
  sp = PlotPanel(xdxr) 
  im1 = sp.addPixels(s1,s2,g)
  im1.setColorModel(cmap)
  im1.setInterpolation(inear)
  if clips:
    im1.setClips(clips[0],clips[1])
  # remove background
  if rmb:
    sp.remove(im1)
  if gh:
    # Add Horizons
    for i in range(len(gh[0])):
      pv = sp.addPoints(gh[0][i],floats(s2.getValues()))#gh[1][i])
      pv.setStyle(gh[2][i])
      pv.setLineWidth(2)
  # Add wells
  for i in range(nw):
    hi = h[i]
    shi = sh[i]
    if velc:
      shi,hi = trim(sh[i],h[i],v=1.0)
    nj = len(hi)
    la = zerofloat(nj,wwidth)
    sl2 = Sampling(wwidth,d2,x2w[i]-d2*(wwidth-1)/2)
    for w in range(wwidth):
      for j in range(nj):
        la[w][j] = hi[j]
    log = sp.addPixels(shi,sl2,la)
    if clips:
      log.setClips(clips[0],clips[1])
    else:
      log.setClips(min(g),max(g))
    if lclips:
      log.setClips(lclips[0],lclips[1])
    log.setColorModel(lcmap)
    log.setInterpolation(inear)
    if rml:
      sp.remove(log)
  if tops:
    # Add tops
    wd = wwidth
    for iw in range(nw):
      nfms = len(fmsp)
      for ic in range(nfms):
        tpd = tops[iw][ic]
        if tpd>=sz[iw].first:
          iz = sz[iw].indexOfNearest(tpd)
          ttpd = tz[iw][iz]
          ttpd = fillfloat(ttpd,wd)
          xvs = zerofloat(wd)
          for ij in range(wd):
            xvs[ij] = x2w[iw]-d2*(wd-1)/2+ij*d2
          tcolor = fmsp[ic][1]
          pvt = sp.addPoints(ttpd,xvs)
          pvt.setStyle(tcolor); pvt.setLineWidth(3)
  sp.setHLabel("Distance (km)")
  sp.setVLabel("Time (s)")
  if lims:
    sp.setLimits(lims[0],lims[1],lims[2],lims[3])
  else: 
    sp.setLimits(s2.first,s1.first,s2.last,s1.last)
  if title:
    sp.setTitle(title)
  if paper or slides:
    sp.removeTitle()
  if cbar:
    sp.addColorBar(cbar)
    #sp.setColorBarWidthMinimum(100)
  frame = PlotFrame(sp)
  frame.setSize(1000,1000)
  frame.setFontSize(24)
  if paper:
    frame.setSize(711,325)
    if psize:
      frame.setSize(psize[0],psize[1])
    if cwp:
      frame.setFontSizeForPrint(8.0,222.0)
      if _prints: frame.paintToPng(720,3.08,_pngDir+paper+'.png')
    if geo:
      frame.setFontSizeForPrint(8.0,240.0)
      if _prints: frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    frame.setSize(1000,700)
    #frame.setSize(800,700)
    frame.setFontSizeForSlide(1.0,1.0)
    if _prints: frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)
  frame.setDefaultCloseOperation(eoc);

def plotLogPanelsTops(st,sz,a,g,wells,fmsp,tops,zt,tz,gh,x2w):
	nw = len(a); 
	alim = [0,20]
	#alim = [1,8]
	glim = [0,200]
	p1 = PlotPanel(1,nw*2,xdxr)
	p1.setVLabel("Depth (km)")
	nfms = len(fmsp)
	fmtps,fmhrz = [],[]
	for ic in range(nfms):
		for clr in fmsp[ic][1]:
			fmtps.append((fmsp[ic][0],clr+"O"))
			fmhrz.append((fmsp[ic][0],clr+"S"))
			break
	hc = zeroint(nfms,nw)
	for iw in range(nw):
		x2c = x2w[iw]
		for ic in range(nfms):
			ghtm = gh[1][ic][0]
			for j in range(1,len(gh[1][ic])):
				ghtp = gh[1][ic][j]
				if (ghtm<=x2c and x2c<ghtp):
					hc[iw][ic] = j
				ghtm = ghtp
	iwp = 0
	def addTops(x,iwp,iw,mrk):
		dist = 1.0
		for ic in range(nfms):
			fmpt = fillfloat(-99.9,sz[iw].count)
			tpd = tops[iw][ic]
			if tpd>sz[iw].first:
				tcolor = fmtps[ic][1]
				iz = sz[iw].indexOfNearest(tpd)
				fmpt[iz] = x[iz]
				pvt = p1.addPoints(0,iwp,sz[iw],fmpt)
				pvt.setStyle(tcolor); pvt.setMarkSize(mrk)
			fmpt = fillfloat(-99.9,sz[iw].count)
			hrt = int((gh[0][ic][hc[iw][ic]]-tz[iw][0])/st.delta)
			if hrt<len(zt[iw]) and hrt>0:
				hpd = zt[iw][hrt]
				iz = sz[iw].indexOfNearest(hpd)
				fmpt[iz] = x[iz]
				pvh = p1.addPoints(0,iwp,sz[iw],fmpt)
				tcolorh = fmhrz[ic][1]
				pvh.setStyle(tcolorh); pvh.setMarkSize(mrk)
				#if tpd>0.0 and ic==0 and iw==3: 
				#	dist = abs(hpd-tpd)
	for iw in range(nw):
		ap = p1.addPoints(0,iwp,sz[iw],a[iw])
		p1.setHLabel(iwp,"I"+str(wells[iw]))
		p1.setHLimits(iwp,alim[0],alim[1])
		dist = addTops(a[iw],iwp,iw,10)
		iwp+=1
	for iw in range(nw):
		gp = p1.addPoints(0,iwp,g[iw][1],g[iw][0])
		p1.setHLabel(iwp,"G"+str(wells[iw]))
		p1.setHLimits(iwp,glim[0],glim[1])
		addTops(g[iw][0],iwp,iw,10)
		iwp+=1
	frame = PlotFrame(p1)
	frame.setSize(1500,1500)
	frame.setFontSize(28)
	frame.setDefaultCloseOperation(eoc);
	frame.setVisible(True)



"""
Plots a 2-D seismic slice
"""
def plotSlice(sy,st,img,v=None,title=None,lims=None,png=None,slides=None,paper=None,linear=None,\
		cmap=gray,tens=None,cbar=None,psize=None,clips=None,clips2=None,pts=None,vlabel=None,cont=None):
	p1 = PlotPanel(xdxr)
	px = p1.addPixels(st,sy,img)
	if v:
		px2 = p1.addPixels(st,sy,v)
		p1.addColorBar("Velocity (km/s)")
		px2.setColorModel(ColorMap.getJet(0.75))
		if clips2:
			px2.setClips(clips2[0],clips2[1])
	if clips:
		px.setClips(clips[0],clips[1])
	p1.setHLabel("Distance (km)")
	p1.setVLabel("Time (s)")
	if linear:
		px.setInterpolation(inear)
	if vlabel:
		p1.setVLabel(vlabel)
	if title:
		p1.setTitle(title)
	px.setColorModel(cmap)
	if cbar:
		p1.addColorBar(cbar)
		#p1.setColorBarWidthMinimum(75)
	p1.setLimits(sy.first,st.first,sy.last,st.last)
	if lims:
		p1.setLimits(lims[0],lims[1],lims[2],lims[3])
	if cont:
		cv = ContoursView(st,sy,img)
		cv.setOrientation(ContoursView.Orientation.X1DOWN_X2RIGHT)
		cv.setLineColor(Color.WHITE)
		cv.setLineWidth(3)
		tile = p1.getTile(0,0)
		tile.addTiledView(cv)
	if tens:
		tv = TensorsView(st,sy,tens)
		tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
		tv.setLineColor(Color.YELLOW)
		tv.setLineWidth(2)
		#tv.setEllipsesDisplayed(40)
		fac1 = 40; fac2 = 20
		e1 = Sampling(st.count/fac1+1,st.delta*fac1,st.first)
		e2 = Sampling(sy.count/fac2+1,sy.delta*fac2,sy.first)
		tv.setEllipsesDisplayed(e1,e2)
		tv.setScale(0.55)
		#tv.setScale(10.00)
		tile = p1.getTile(0,0)
		tile.addTiledView(tv)
	if pts:
		ptv = p1.addPoints(pts[0],pts[1])
		ptv.setStyle("rO")
		ptv.setMarkSize(8)
	if slides: p1.removeTitle()
	if paper: p1.removeTitle()
	frame = PlotFrame(p1)
	frame.setSize(1000,1000)
	frame.setFontSize(24)
	frame.setDefaultCloseOperation(eoc);
	frame.setVisible(True)
	if png:
		if _prints: frame.paintToPng(1000,5.00,_pngDir+png+'.png')  
	if slides:
		frame.setSize(1000,700)
		#frame.setSize(720,700)
		frame.setFontSizeForSlide(1.0,1.0)
		if _prints: frame.paintToPng(720,7.00,_pngDir+slides+'.png')
	if paper:
		frame.setSize(711,325)
		if psize:
			frame.setSize(psize[0],psize[1])
		if cwp:
			frame.setFontSizeForPrint(8.0,222.0)
			if _prints: frame.paintToPng(720,3.08,_pngDir+paper+'.png')
		if geo:
			frame.setFontSizeForPrint(8.0,240.0)
			if _prints: frame.paintToPng(720,3.33,_pngDir+paper+'.png')
	frame.setVisible(True)

def plotCoordSlice(s2,s3,x2s,x3s,x2w,x3w,tpimg,pt2=None,line=None,dash=None,slides=None,paper=None):
  pp = PlotPanel()
  pp.setHLabel("Crossline (km)")
  pp.setVLabel("Inline (km)")
  pp.addPixels(s2,s3,tpimg)
  pp.setLimits(s2.first,s3.first,s2.last,s3.last)
  if line:
    ln = pp.addPoints(x2s,x3s)
    ln.setStyle("w-")
    ln.setLineWidth(5)
  pt = pp.addPoints(x2w,x3w)
  pt.setStyle("wO")
  pt.setMarkSize(10)
  if pt2:
    x2s,x3s,x2w,x3w = pt2 
    pt = pp.addPoints(x2w,x3w)
    pt.setStyle("wO")
    pt.setMarkSize(10)
    if dash:
      ln = pp.addPoints(x2s,x3s)
      ln.setStyle("w--")
      ln.setLineWidth(5)
  frame = PlotFrame(pp)
  frame.setSize(1000,1000)
  frame.setFontSize(24)
  frame.setDefaultCloseOperation(eoc)
  if slides:
    frame.setSize(1000,700)
    frame.setFontSizeForSlide(1.0,1.0)
    if _prints: frame.paintToPng(720,6.00,_pngDir+slides+'.png')
  if paper:
    frame.setSize(556,314)
    if cwp:
      frame.setFontSizeForPrint(8.0,222.0)
      if _prints: frame.paintToPng(720,3.08,_pngDir+paper+'.png')
    if geo:
      frame.setFontSizeForPrint(8.0,240.0)
      if _prints: frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  frame.setVisible(True)

def plotHSlice(s2,s3,tpimg,x2w=None,x3w=None,slides=None,paper=None):
	pp = PlotPanel()
	pp.setHLabel("Distance (km)")
	pp.setVLabel("Distance (km)")
	pp.addPixels(s2,s3,tpimg)
	if x2w and x3w:
		pt = pp.addPoints(x2w,x3w)
		pt.setStyle("yO")
		pt.setMarkSize(5)
	frame = PlotFrame(pp)
	frame.setSize(1000,1000)
	frame.setFontSize(24)
	frame.setDefaultCloseOperation(eoc)
	if slides:
		frame.setSize(1000,700)
		frame.setFontSizeForSlide(1.0,1.0)
		if _prints: frame.paintToPng(720,6.00,_pngDir+slides+'.png')
	if paper:
		frame.setSize(556,314)
		frame.setFontSizeForPrint(8.0,222.0)
		if _prints: frame.paintToPng(720,3.08,_pngDir+paper+'.png')
	frame.setVisible(True)

def plotPanel3(s1,s2,s3,x,clrpix=None,maps=None,mpoints=None,\
		label1=None,label2=None,label3=None,cbar=gray,coord=None,paper=None,slides=None):
	svp = Viewer3P(s1,s2,s3,x)
	if maps:	
		svp.getPP().pixelsView23.tile.addTiledView(maps)
	if mpoints:
		svp.getPP().pixelsView23.tile.addTiledView(mpoints)
	svp.setColorModel1(cbar)
	if label1:
		svp.setLabel1(label1)
	if label2:
		svp.setLabel2(label2)
	if label3:
		svp.setLabel3(label3)

	svp.show()

"""
Plots phase rotation errors
"""
def plotPhaseRMS(mrms,nph,phf,png=None,paper=None,slides=None):
  sp = SimplePlot()
  sph = Sampling(nph/phf,phf,0.0)
  sp.addPoints(mrms)
  sp.setHLabel("phase (degrees)")
  sp.setVLabel("error")
  sp.setSize(1500,1000)
  if png:
    if _prints: sp.paintToPng(1200,6.00,_pngDir+png+'.png')
  if paper:
    sp.removeTitle()
    sp.setSize(1000,800)
    sp.setFontSizeForPrint(8.0,240.0)
    if _prints: sp.paintToPng(720,3.33,_pngDir+paper+'.png')  
  if slides:
    sp.removeTitle()
    sp.setSize(950,800)
    sp.setFontSizeForSlide(1.0,1.0)
    if _prints: sp.paintToPng(720,6.00,_pngDir+slides+'.png')

"""
Plots u0 vs dmin
"""
def plotuo(us,sl,png=None,paper=None,slides=None):
  sp = SimplePlot()
  sp.addPoints(sl,us)
  sp.setHLabel("u0 time shift (s)")
  sp.setVLabel("minimum D")
  sp.setSize(1200,1000)
  if png:
    if _prints: sp.paintToPng(1200,6.00,_pngDir+png+'.png')
  if paper:
    sp.removeTitle()
    sp.setSize(1000,800)
    sp.setFontSizeForPrint(8.0,240.0)
    if _prints: sp.paintToPng(720,3.33,_pngDir+paper+'.png')  
  if slides:
    sp.removeTitle()
    sp.setSize(950,800)
    sp.setFontSizeForSlide(1.0,1.0)
    if _prints: sp.paintToPng(720,6.00,_pngDir+slides+'.png')

def plotCompare(s1,s2,x1,x2,title):
	pp = PlotPanel(xdxr)
	pp.setTitle(title)
	pp.setVLabel('time (s)')
	ln1 = pp.addPoints(s1,x1)
	ln2 = pp.addPoints(s2,x2)
	ln2.setLineColor(Color.RED)
	#pp.setVLimits(max(s1.first,s2.first),max(s1.last,s2.last))
	frame = PlotFrame(pp)
	frame.setDefaultCloseOperation(eoc)
	frame.setSize(500,1000)
	frame.setVisible(True)

"""
Plots just a gaussian
"""
def gaussianPlot(sig):
  nf = 801
  f1 = zerofloat(nf)
  f2 = zerofloat(nf)
  f1[int(nf/2)] = 100000.0
  RecursiveGaussianFilter(sig).apply2(f1,f2)
  f2 = mul(-1,f2)
  si = SimplePlot()
  si.addPoints(f2)
  si.setSize(700,500)
  if _prints: si.paintToPng(900,3.33,_pngDir+"RICKER"+".png")

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
    if _prints: sp.paintToPng(360,3.33,_pngDir+png+'.png') 


"""
Plots a wavelet impulse response
"""
def waveletTest():
  fp = 35 
  nt = 200; nz = 700;
  st = Sampling(nt,0.001,0.0)
  sz = Sampling(nz,0.0001524,0.0)
  rf = zerofloat(nz)
  #rf[0] = 1.0
  #rf[(nz-1)/3] = 1.0
  rf[(nz-1)/2] = 1.0
  #rf[2*(nz-1)/3] = 1.0
  #rf[(nz-2)] = 1.0
  vl = fillfloat(2.0,nz)
  tzl = WTutils().tzVlog(vl,sz)
  #sy = WTutils().syntheticSeismogram(fp,rf,tzl,st,sz,"ricker")
  sy = WTutils().syntheticSeismogram(fp,rf,tzl,st,sz,"morlet")
  #sy = applyPhaseRot(92,sy)
  ss = Sampling(len(sy),st.delta,tzl[0])
  pp = PlotPanel(xdxr)
  pp.addPoints(ss,sy)
  pp.setVLabel("Time (s)")
  #pp.setVLimits(0,0.1);
  frame = PlotFrame(pp)
  frame.setDefaultCloseOperation(eoc)
  frame.setSize(350,600)
  frame.setFontSizeForSlide(1.0,1.0)
  #frame.paintToPng(720,4.00,_pngDir+'swavelet.png')
  #frame.setFontSizeForPrint(8.0,240.0)
  #frame.paintToPng(720,6.00,_pngDir+'pwavelet.png')
  frame.setVisible(True)


"""
Plots teaser figure for CWP report 2012 of synthetic and trace before and after alignment
"""
def plotTeaser(tr,st,y,ys,sst,ss,tm=None,lims=None,paper=None,slides=None,paper2=None,paper3=None):
  if tm:
    p1 = PlotPanel(1,3,xdxr)
    p1.setHLabel(2,"tau")
    stau = Sampling(ss.count,ss.delta,tm[0])
    p1.addPoints(0,2,stau,tm)
  else: 
    p1 = PlotPanel(2,1)#,xdxr)
  l1 = p1.addPoints(0,0,st,tr)
  l1.setLineWidth(2)
  p1.setHLabel("Time (s)")
  #ss2 = Sampling(ss.count,ss.delta,sst.first)
  l2 = p1.addPoints(0,0,ss,y)
  l2.setStyle('r-')
  l2.setLineWidth(2)
  l3 = p1.addPoints(1,0,st,tr)
  l3.setLineWidth(2)
  l4 = p1.addPoints(1,0,sst,ys)
  l4.setStyle('r-')
  l4.setLineWidth(2)
  if lims:
    p1.setVLimits(0,-lims,lims)
    p1.setVLimits(1,-lims,lims)
    if tm: p1.setVLimits(2,-lims,lims)
  p1.setHLimits(0,sst.last+0.5)
  frame = PlotFrame(p1)
  frame.setSize(900,1500)
  frame.setDefaultCloseOperation(eoc);
  frame.setVisible(True)
  if paper:
    frame.setSize(1360,354)
    frame.setFontSizeForPrint(8.0,504.0)
    if _prints: frame.paintToPng(720,6.00,_pngDir+paper+'.png')
  if paper2:
    frame.setSize(904,379)
    frame.setFontSizeForPrint(8.0,240.0)
    if _prints: frame.paintToPng(720,3.33,_pngDir+paper2+'.png')
  if paper3:
    frame.setSize(904,379)
    frame.setFontSizeForPrint(8.0,312.0)
    if _prints: frame.paintToPng(720,4.33,_pngDir+paper3+'.png')
  if slides:
    frame.setSize(500,900)
    frame.setFontSizeForSlide(1.0,1.0)
    if _prints: frame.paintToPng(720,6.00,_pngDir+slides+'.png')
  frame.setVisible(True)


def plot3P(s1,s2,s3,x,title,cmap=gray,p=None,m=None,d=None,v=None,z=None,\
						clips1=None,clips2=None,ss=None,lims=None,slc=None,con=None,\
            cbar=None,paper=None,slides=None):
  svp = Viewer3P(s1,s2,s3,x)
  if v:
    svp.addPixels(v)
    svp.addColorBar("Velocity (km/s)")
    svp.setColorModel2(ColorMap.getJet(0.75))
    if clips2:
      svp.setClips2(clips2[0],clips2[1])
  if d:
    svp.addPixels(d)
    svp.addColorBar("Depth (km)")
    svp.setColorModel2(ColorMap.getJet(0.75))
    if clips2: 
      svp.setClips2(clips2[0],clips2[1])
    if con:
      svp.addContours(s1,s2,s3,d)
  if m:
    svp.getPP().pixelsView23.tile.addTiledView(m)
  if p:
    #svp.addPnts23(p[0],p[1])
    svp.addPts(p)
  #svp.setTitle(title)
  #svp.setSize(1200,1000)
  if cmap:
    svp.setColorModel1(cmap)
  if lims:
    svp.setLimits1(lims[0],lims[1])
  if clips1:
    svp.setClips1(clips1[0],clips1[1])
  if cbar:
    svp.addColorBar(cbar)
  svp.setLabel1("Time (s)")
  svp.setLabel2("Crossline (km)")
  svp.setLabel3("Inline (km)")
  if z:
    svp.setLabel1("Depth (km)")
  if slc:
    svp.setSlices(inro((slc[0]-s1.first)/s1.delta),\
  	              inro((slc[1]-s2.first)/s2.delta),\
  								inro((slc[2]-s3.first)/s3.delta))
  #svp.setSlices(80,160,343)
  #svp.setSlices(343,160,80)
  #svp.setSlices(563,174,84)
  #svp.setSlices(575,228,84)
  svp.setSize(600,600)
  if paper:
    svp.setFontSizeForPrint(8.0,240.0)
    if _prints: svp.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    svp.setFontSizeForSlide(1.0,1.0)
    if _prints: svp.paintToPng(720,7.00,_pngDir+slides+'.png')
  svp.show()

def plot3(f,title,g=None,sd=None,tensors=None,logs=None,\
          cmap=gray,cmap2=gray,plots=None,lmap=None,slides=None,paper=None):
  world = World()
  ipg = ImagePanelGroup(s1,s2,s3,f);
  if sd: 
    ipg = ImagePanelGroup(sd,s2,s3,f); 
  ipg.setColorModel(cmap)
  world.addChild(ipg)
  if g: 
    ipg2 = ImagePanelGroup(s1,s2,s3,g); 
    if sd: 
      ipg2 = ImagePanelGroup(sd,s2,s3,g); 
    ipg2.setColorModel(cmap2)
    world.addChild(ipg2)
  if tensors:
    a = 15
    addTensorsInImage(ipg.getImagePanel(Axis.X),tensors,a)
    addTensorsInImage(ipg.getImagePanel(Axis.Y),tensors,a)
    addTensorsInImage(ipg.getImagePanel(Axis.Z),tensors,a)
  if logs:
    samples = makeLogSamples(logs)
    lg = makeLogPoints(samples,min(f),max(f),size=8)
    world.addChild(lg)
  if lmap:
    lm = makeLmap(lmap)
    world.addChild(lm)
  sf = SimpleFrame(world)
  ov = sf.getOrbitView()
  ov.setAxesScale(1.0,1.0,2.0)
  ov.setScale(2.0)
  if plots:
    # For Geophysics Paper:
    azimuth=219.35944
    elevation=32.0
    scale=1.00
    scx,scy,scz = 1.0,1.0,3.5
    ix=inro((2.425-s3.first)/s3.delta)
    iy=inro((6.375-s2.first)/s2.delta)
    iz=inro((0.954-s1.first)/s1.delta)
    tx,ty,tz = -0.025472915230442886,-0.16815276025007686,-0.10889203275845535
    # For Geophysics Paper Tensors
    #azimuth=-34.94662
    #elevation=23.705036
    #scale=1.9671422
    #scx,scy,scz = 1.0,1.0,4.0
    #tx,ty,tz = -0.06336199061374731,-0.06187212077664625,-0.0131456674616428
    #ix=inro((2.475-s3.first)/s3.delta)
    #iy=inro((5.000-s2.first)/s2.delta)
    #iz=inro((0.894-s1.first)/s1.delta)
    ov.setScale(scale)
    ov.setAxesScale(scx,scy,scz)
    ov.setElevation(elevation)
    ov.setAzimuth(azimuth)
    ov.setTranslate(Vector3(tx,ty,tz))
    ipg.setSlices(iz,iy,ix)
  sf.setTitle(title)
  sf.setVisible(True)
  if slides:
  	if _prints: sf.paintToFile(_pngDir+slides+'.png')
  if paper:
  	if _prints: sf.paintToFile(_pngDir+paper+'.png')
  
def makeLmap(lmap):
	gr = Group()
	x2p,x3p,x2s,x3s = lmap
	pg = makeTopPointGroup(x2p,x3p,Color.RED)
	lg = makeTopPointGroup(x2s,x3s,Color.RED,line=3)
	makeTopPointGroupBox(gr,Color.BLUE)
	wl = makeWellLine(x2p,x3p,Color.RED)
	gr.addChild(pg)
	gr.addChild(lg)
	gr.addChild(wl)
	return gr

def makeWellLine(x2,x3,color):
	nw = len(x2)
	s1d = s1.decimate(4)
	f = zerofloat(s1d.count,nw)
	s1g = [s1d for i in range(nw)]
	logs = f,s1g,add(x2,-s2.delta),add(x3,-s3.delta)
	samples = makeLogSamples(logs)
	lg = makeLogPoints(samples,0,0,size=1)
	states = StateSet()
	cs = ColorState()
	cs.setColor(color)
	states.add(cs)
	lg.setStates(states)
	return lg

def makeTopPointGroup(x2,x3,color,line=None):
	n = len(x2)
	x1 = fillfloat(s1.first-s1.delta,n)	
	xyz = zerofloat(3*n)
	copy(n,0,1,x3,0,3,xyz)
	copy(n,0,1,x2,1,3,xyz)
	copy(n,0,1,x1,2,3,xyz)
	rgb = None
	if line:
		pg = LineGroup(xyz,rgb)
		ls = LineState()
		ls.setWidth(line)
		ss = StateSet()
		ss.add(ls)
	else:
		pg = PointGroup(xyz,rgb)
		ps = PointState()
		ps.setSize(8)
		ps.setSmooth(False)
		ss = StateSet()
		ss.add(ps)
	cs=ColorState()
	cs.setColor(color)
	ss.add(cs)
	pg.setStates(ss)
	return pg

def makeTopPointGroupBox(gr,color):
	x2 = floats(s2.getValues())
	x3 = floats(s3.getValues())
	x2f = fillfloat(s3.first,len(x2))
	x3f = fillfloat(s2.first,len(x3))
	x2l = fillfloat(s3.last, len(x2))
	x3l = fillfloat(s2.last, len(x3))
	bg1 = makeTopPointGroup(x2,x2f,color,line=1)
	bg2 = makeTopPointGroup(x2,x2l,color,line=1)
	bg3 = makeTopPointGroup(x3f,x3,color,line=1)
	bg4 = makeTopPointGroup(x3l,x3,color,line=1)
	gr.addChild(bg1)
	gr.addChild(bg2)
	gr.addChild(bg3)
	gr.addChild(bg4)
	return gr

def makeLogSamples(logs):
	f,sl,x2,x3 = logs
	fl,x1l,x2l,x3l = [],[],[],[]
	nw = len(f)
	for i in range(nw):
		ns = len(f[i])
		fl.append(f[i])
		x1l.append(floats(sl[i].getValues()))
		x2l.append(fillfloat(x2[i],ns))
		x3l.append(fillfloat(x3[i],ns))
	return fl,x1l,x2l,x3l

def makeLogPoints(samples,cmin,cmax,cbar=None,cmap=gray,size=None):
  lg = Group()
  fl,x1l,x2l,x3l = samples
  for i,f in enumerate(fl):
    f = fl[i]
    x1 = x1l[i]
    x2 = x2l[i]
    x3 = x3l[i]
    pg = makePointGroup(f,x1,x2,x3,cmin,cmax,cbar,cmap,size)
    lg.addChild(pg)
  return lg

def makePointGroup(f,x1,x2,x3,cmin,cmax,cbar,cmap,size):
  n = len(x1)
  xyz = zerofloat(3*n)
  copy(n,0,1,x3,0,3,xyz)
  copy(n,0,1,x2,1,3,xyz)
  copy(n,0,1,x1,2,3,xyz)
  rgb = None
  if cmin<cmax:
    cmap = ColorMap(cmin,cmax,cmap)
    if cbar:
      cmap.addListener(cbar)
    rgb = cmap.getRgbFloats(f)
  pg = PointGroup(xyz,rgb)
  ps = PointState()
  ps.setSize(size)
  ps.setSmooth(False)
  ss = StateSet()
  ss.add(ps)
  pg.setStates(ss)
  return pg

def addTensorsInImage(ipg,et,a):
	tp = TensorsPanel(s1,s2,s3,et)
	tp.setEllipsoidSize(a)
	ipg.getFrame().addChild(tp)
	return tp

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

"""
	if dmap:
		dist = zerofloat(n1,n2)
		nrst = zerofloat(n1,n2)
		dmap1 = zerofloat(n2)
		wds = zerofloat(nw+1)
		x2cp = -s2x.last/10.0
		for i in range(nw):
			wds[i] = x2c[i]-x2cp
			x2cp = x2c[i]
		wds[nw] = wds[0]
		wds = smooth(wds,n2)
		wds = sub(1.0,div(wds,max(wds)))
		ng = NearestGridder2(wa,xa1,xa2)
		ng.computeDistancesAndValues(s1,s2x,dist,nrst)
		for i in range(n2):
			dmap1[i] = pow(min(dist[i]),wds[i])
		#SimplePlot().addPoints(dmap1)
		return g,dmap1,dist
def smooth(x,n):
	y = zerofloat(n)
	li = LinearInterpolator()
	li.setUniform(len(x),1.0,0.0,x)
	li.interpolate(n,float(len(x)-1)/float(n),0.0,y)
	return y

"""
