
from imports import *

dataDir = '/Users/amunoz/Home/data/coda/'
_pngDir = './pics/'
fpath = './data/Inversion_data/'
dpath = './data/invdat/'
outDir = dataDir+'invout/'
s1 = Sampling(6001,0.001,0.0)

def main(args):
  data = readDats()  
  pdata = pairData(data)
  go2Dtest1(pdata)  

def go2Dtest1(pdata):
  dt = s1.delta
  smin = -0.2
  smax =  0.2
  rmin = -0.02
  rmax =  0.02
  dr = 0.005
  ngg = 15
  flv = FloatList()
  print 'number of pairs is',len(pdata)
  pi = 0
  for pair in pdata:
    fd,gd = pair
    f = fd[0]
    g = gd[0]
    dw = DynamicWarpingCO(inro(smin/dt),inro(smax/dt),rmin,rmax,dr)
    dw.setInterpolation(DynamicWarpingCO.Interpolation.SPLINE)
    e = dw.computeErrors(f,g)
    #SimplePlot.asPixels(e)
    #u = dw.findShiftsR(f,e)#,ngg)
    #ae = dw.accumulate(f,e)
    #SimplePlot.asPixels(ae)
    #u = dw.findShiftsR(e)
    #u = dw.findShiftsR(f,g,e,ngg)
    u = dw.findShiftsR(f,g,e)
    h = dw.applyShifts(g,u)
    u = mul(u,dt)
    gri = dw.getG1()
    gr = add(mul(floats(gri),dt),s1.first)
    #du = wtu.backwardsDiff(u,dt)
    sh = Sampling(len(h),dt,s1.first+u[0])
    sl = Sampling(len(e[0]),dt,smin)
    vr = bderivative(u,dt)
    flv.add(vr)
    #plot1D(s1,f,"f")
    #lot1D(s1,g,"g")
    #lot1D(sh,h,"h")
    #plot1D(s1,u,"u")
    #plot1D(s1,vr,"vr")
    #plot1D(s1,f,"f",gr=gri)
    #plot1D(s1,u,"u",gr=gri)
    #plot1D(s1,f,"f & g",s1,g)
    #plot1D(s1,f,"f & h",sh,h)
    stb = str(fd[2]); stc = str(fd[3])
    writeText('spline_shifts_'+stb+'_'+stc,u)
    writeText('spline_velocity_'+stb+'_'+stc,vr)
    writeText('spline_h_'+stb+'_'+stc,vr)
  vra = flv.trim(0.005)
  goHistogram(vra,"velocity ratio")

################################################################################
## Utils

def pairData(data):
  pdata = []
  for sig1 in data:
    a1 = sig1[1] 
    if a1>0: continue
    for sig2 in data:
      a2 = sig2[1]
      if sig1[2]==sig2[2] and sig1[3]==sig2[3] and a1!=a2:
        pdata.append((sig1,sig2))
  return pdata

def inro(x):
	return int(round(x))

def floats(x):
	n = len(x)
	xd = zerofloat(n)
	for i in range(n):
		xd[i] = float(x[i])
	return xd

def bderivative(x,dx):
  n = len(x)
  y = zerofloat(n)
  for i in range(1,n):
    y[i] = (x[i]-x[i-1])/(dx);
  y[0] = y[1];
  return y;

################################################################################
## Plots

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

################################################################################
## I/O

def readDats():
  files = File(dpath)
  data = []
  for fil in files.listFiles():
    name = fil.getName()
    parts = (name.replace(".dat","")).split("_")
    a = int(parts[0])
    b = int(parts[1])
    c = int(parts[2])
    n = int(parts[3])
    fa = zerofloat(n)
    ais = ArrayInputStream(dpath+name)
    ais.readFloats(fa)
    ais.close()
    data.append((fa,a,b,c))
  return data 

def readTxtWriteDat():
  pfx = 'tmp_data_bh_sr_f15hz_rec_'
  files = File(fpath)
  for fil in files.listFiles():
    fl = FloatList()  
    br = BufferedReader(FileReader(fil))
    name = fil.getName()
    if pfx not in name: continue
    lns = "0"
    while lns:
      lns=br.readLine()
      if lns==None: break
      lnf = Float.parseFloat(lns)
      fl.add(lnf)
    br.close()
    fa = fl.trim()
    n = str(len(fa))
    name = name.replace(".asc","")
    parts1 = name.split(pfx)
    parts2 = parts1[1].split("_")
    a = str(int(float(parts2[0]))) # source number
    b = str(int(float(parts2[1]))) # baseline (0) or time-lapse (1)
    c = str(int(float(parts2[2]))) # receiver number
    dname = b+"_"+a+"_"+c+"_"+n+".dat"
    aos = ArrayOutputStream(dpath+dname)
    aos.writeFloats(fa)
    aos.close()
    print 'read',name,'write',dname
    #SimplePlot.asPoints(fa) 

def writeText(name,data):
	n = len(data)
	out = BufferedWriter(FileWriter(outDir+name+'.txt'))
	for i in range(n):
		out.write(str(data[i]))
		out.newLine()
	out.close()

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
