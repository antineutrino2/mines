"""
Acoustic wavefield modeling
"""
from java.lang import *
from javax.swing import *
from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util.ArrayMath import *

from fd import *

#############################################################################
# Parameters

ffile = 'marmousi'
sz = Sampling(176,16.0,0.0)
sx = Sampling(421,16.0,1400.0)
#st = Sampling(5101,0.0012,0.0)
st = Sampling(3001,0.0012,0.0)
nz,nx,nt = sz.count,sx.count,st.count
dz,dx,dt = sz.delta,sx.delta,st.delta
fpeak = 10.0 # Ricker wavelet peak frequency
nit = 51 

#ffile = 'layered'
#sz = Sampling(151,16.0,0.0)
#sx = Sampling(301,16.0,1400.0)
##st = Sampling(2001,0.0012,0.0)
#st = Sampling(1001,0.002,0.0)
#nz,nx,nt = sz.count,sx.count,st.count
#dz,dx,dt = sz.delta,sx.delta,st.delta
#fpeak = 50.0 # Ricker wavelet peak frequency

# source
kzs = 0
kxs = (nx-1)/2
kxsa = zerofloat(nx)

# receivers
kzr = 0
kxr = rampint(0,1,nx)
kzra = fillint(kzr,len(kxr))

saveDir = '../../../../../data/migclass/'
dataDir = saveDir
pngDir = '../dtw/dtw1d2d/'

#############################################################################
# Functions

def main(args):
  #forwardBackModel()
  #forwardModelMM()
  #forwardModelDAR()
	plot(readModel())

def forwardBackModel():
  v1 = readModelLat()
  #v2 = readModel2()
  u = AcousticWavefield(sz,sx,st)
  u.forwardPropagate(fpeak,kzs,kxs,v1)
  data = u.wavefield(kzr,kxr) # data
  w1 = u.wavefield()

  #u.forwardPropagate(fpeak,kzs,kxs,v2)
  #data2 = u.wavefield(kzr,kxr)
  #data = sub(data1,data2)

  u2 = AcousticWavefield(sz,sx,st)
  u2.backPropagate(data,kzra,kxr,v1)
  w2 = u2.wavefield()
  img = u.correlateSR(w1,w2)
  #plot(data,cmap=gray,perc=100)#,png='marmousi-sg_'+str(kxs))
  plot(img,cmap=gray,perc=100)#,png='marmousi-sg_'+str(kxs)+'_IC')
#  dsum = zerofloat(len(data[0]),len(data))
#  dsum = add(dsum,data)
#  for i in range(1,dl):
#    kkk = i*4
#    u.forwardPropagate(fpeak,kzs,kkk,v1)
#    data1 = u.wavefield(kzr,kxr)
#    u.forwardPropagate(fpeak,kzs,kkk,v2)
#    data2 = u.wavefield(kzr,kxr)
#    data = sub(data1,data2)
#    dsum = add(dsum,data)
#    plot(data,cmap=gray,perc=98)
  plot(v1)#,png='marmousi-model')
#  plot(dsum,cmap=gray,perc=98)
#  writeImage('mmtest',data)

def forwardModelMM():
  u1 = AcousticWavefield(sz,sx,st)
  for it in range(1):
    v1 = readModel()
    vw = zerofloat(len(v1[0]),len(v1))
    for i2 in range(len(v1)):
      for i1 in range(len(v1[0])):
        if (v1[i2][i1]==1500.0): vw[i2][i1]=1500.0
        else: v1[i2][i1] *= (1.25+it*0.01)
    ExponentialSmoother(100,16).apply(v1,v1)
    for i2 in range(len(v1)):
      for i1 in range(len(v1[0])):
        if (vw[i2][i1] == 1500.0):
          v1[i2][i1] = vw[i2][i1]
    u1.forwardPropagate(fpeak,kzs,kxs,v1)
    data1 = u1.wavefield(kzr,kxr)

    v2 = readModel2()
    u1.forwardPropagate(fpeak,kzs,kxs,v2)
    data2 = u1.wavefield(kzr,kxr)

    data = sub(data1,data2)
    plot(v1)#,png='marm-model-smoothed-125perc')
    nm = str(it)
    plot(data,cmap=gray,perc=99)#,png='marm-125perc-smoothed')
    plot(data1,cmap=gray,perc=99)#,png='marm-125pc-dr')
    plot(data2,cmap=gray,perc=99)#,png='marm-dr_only')
    #writeImage('marmousiTest-v'+nm,data)

"""
Direct arrival removal forward model layered
"""
def forwardModelDAR():
  u = AcousticWavefield(sz,sx,st)
  #f = u.getSource(fpeak)
  v1 = readModelLat()
  for it in range(51):
    u.forwardPropagate(fpeak,kzs,kxs,v1)
    data1 = u.wavefield(kzr,kxr)
    plot(data,cmap=gray,perc=100)#,png='LayeredModel-sine_src-v'+str(it))
    writeImage('layered-sineSourceTapered2-'+str(it),data)
  plot(v1)
  #plot(data,cmap=gray,perc=98)
  #writeImage('layered-sineSource-'+str(it),data)

"""
Modeling
"""
_v1 = 2500.0
_v2 = 4000.0
_vg = 1500.0
_nz2 = int(nz/6)
_nx2 = int(nx/10)
_sig = 1

def readModel():
  if ffile=='layered':
    v = fillfloat(_v1,nz,nx)
    for ix in range(nx):
      for iz in range(nz/2,nz):
        v[ix][iz] = _v2
  else:
    v2 = zerofloat(nz,nx)
    nz2 = nz#+50
    v = zerofloat(nz2,nx)
    readImage('marmousi',v2)
    for iz in range(nz):
      for ix in range(nx):
        v[ix][iz+(nz2-nz)] = v2[ix][iz]
    for iz in range(nz2-nz):
      for ix in range(nx):
        v[ix][iz] = 1500.0
  return v

def readModel2():
  if ffile=='layered':
    v = fillfloat(_v1,nz,nx)
  else:
    v = fillfloat(1500,nz,nx)
  return v

def readModelI(v1,v2):
  if ffile=='layered':
    v = fillfloat(v1,nz,nx)
    for ix in range(nx):
      for iz in range(nz/2,nz):
        v[ix][iz] = v2
  return v

def readModel2I(v1):
  if ffile=='layered':
    v = fillfloat(v1,nz,nx)
  return v

"""
Multiple vertical layering
"""
def readModelVert():
  if ffile=='layered':
    nl = 3 # number of interfaces
    v = zerofloat(nz,nx)
    for ix in range(nx):
      for iz in range(nz/100):
        v[ix][iz] = _v1
    block = int(nz/nl)
    ble = nz%nl
    for bl in range(1,nl+1):
      end = bl*block
      if (bl == nl): end += ble - 1 
      for ix in range(nx):
        for iz in range((nz/100)+(bl-1)*block,end):
            v[ix][iz] = _v1+bl*80
	    v[ix][end] = _v1+bl*80
  return v

def readModelVert2():
  if ffile=='layered':
    v = fillfloat(_v1,nz,nx)
  return v

"""
Lateral v change
"""
def readModelLat():
  if ffile == 'layered':
    v = zerofloat(nz,nx)
    for ix in range(nx):
      for iz in range(nz):
        v[ix][iz] = _v1+5*ix
    #for ix in reversed(range(nx)):
    #  for iz in range(nz/2,nz):
    #    v[ix][iz] = _v2+5*ix
  return v

def readModelLat2():
  if ffile == 'layered':
    v = zerofloat(nz,nx)
    for ix in range(nx):
      for iz in range(nz/2):
        v[ix][iz] = _v1+5*ix 
  return v

"""
Dipping v layers
"""
#def readModelDipLayers():
#	nr = 2
#	nl = nz/nr
#	v = zerofloat(nz,nx)
#	fv = 2500.0
#	dv = 1500.0/nr
#	nxf1 = 20
#	nxf2 = nx-nfx1
#	# Do straight part first
#	il = 0
#	tl = nl
#	for ix in range(nxf1):
#		for iz in range(nz):
#			v[ix][iz] = fv+il*dv
#			if iz==tl: 
#				tl+=nl	
#			  il+=1
#	il = 0
#	tl = nl
#	for ix in range(nxf1):
#	for ix in range(nxf2,nx):
#		for iz in range(nz):
#			v[ix][iz] = fv+il*dv
#			if iz==tl: 
#				tl+=nl	
#			  il+=1



"""
Gaussian
"""
def readModelG():
  if ffile=='layered':
    v = fillfloat(_v1,nz,nx)
    v2 = fillfloat(_v1,_nz2,_nx2)
    v2[_nx2/2][_nz2/2] = _vg
    RecursiveGaussianFilter(_sig).apply00(v2,v2)
    for ix in range(nx):
      for iz in range(2*nz/3,nz):
        v[ix][iz] = _v2
    for ix in range(_nx2):
      for iz in range(_nz2):
        v[ix+nx/2-_nx2/2][iz+nz/3-_nz2/2] = v2[ix][iz]
  return v

def readModelG2():
  if ffile=='layered':
    v = fillfloat(_v1,nz,nx)
    v2 = fillfloat(_v1,_nz2,_nx2)
    v2[_nx2/2][_nz2/2] = _vg
    #RecursiveGaussianFilter(_sig).apply00(v2,v2)
    for ix in range(_nx2):
      for iz in range(_nz2):
        v[ix+nx/2-_nx2/2][iz+nz/3-_nz2/2] = v2[ix][iz]
  return v


#############################################################################
# Plots/IO

def depthAxis():
  x = zerofloat(nz)
  for i in range(nz):
    x[i] = (i*dz)/1000
  return x

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
def plot(x,cmap=jet,perc=100,axis=None,png=None):
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.addColorBar()
  sp.setSize(1200,500)
  sz1 = Sampling(len(x[0]),sz.delta,sz.first)
  pv = sp.addPixels(sz1,sx,x)
  pv.setColorModel(cmap)
  sp.setHLabel("Distance (m)")
  sp.setVLabel("Depth (m)")
  #pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if perc<100:
    pv.setPercentiles(100-perc,perc)
  if png:
    sp.paintToPng(360,3.33,pngDir+png+'.png')

def plot1(x,png=None):
  sp = SimplePlot()
  sp.addPoints(x)
  sp.setSize(1200,500)
  sp.setHLabel("Time (ms)")
  sp.setVLabel("Amplitude")
  if png:
    sp.paintToPng(360,3.33,pngDir+png+'.png')
    

def readImage(name,image):
  fileName = dataDir+name+".dat"
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()

def writeImage(name,image):
  fileName = saveDir+name+'.dat'
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()

#############################################################################
# Do everything on Swing thread.

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
