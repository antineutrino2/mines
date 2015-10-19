"""
Automated semblance velocity picking with dynamic programming
@author Andrew Munoz, CSM	
"""

from bputils import *
from imports import *

papers=False
slides=False
setvis=True

if papers:
  pngDir = "/Users/amunoz/Home/pics/sem/p_"
elif slides:
  pngDir = "/Users/amunoz/Home/pics/sem/s_"
# cmp number -803 (SP1 and x=0)
#gix2 = 0 
gix2 = 251
fpeak = 30
#xmax = -7.0

def main(args):
  #data = getData(gix2)#xmax=xmax)#,stack=1)
 #hypTest(data)
  #nonHypTest(data)
  nonHypTestLoop()
## others
  #goHypTest()
  #goNonHypTest(data)
  #goHyp(data)
  #goNonHyp(data)
  #plotDq(data)
  #goNishantTest()

def nonHypTestLoop():
  fcmp = 0#803
  n = 803#4799-fcmp
  eer = zerofloat(n)
  ver = zerofloat(n)
  eir = zerofloat(n)
  vir = zerofloat(n)
  print 'n=',n
  sw = Stopwatch()
  sw.start()
  for i in range(n):
    data = getData(i+fcmp)
    pdata = nonHypLoop(data)
    rv,re,ev,ee,rvi,rei,evi,eei = pdata
    eer[i] = sum(abs(sub(re,ee)))
    ver[i] = sum(abs(mul(div(sub(ev,rv),rv),100)))
    eir[i] = sum(abs(sub(rei,eei)))
    vir[i] = sum(abs(mul(div(sub(evi,rvi),rvi),100)))
    print 'i=',i
  sw.stop()
  print 'took',sw.time()/60,'minutes to complete'
  plot1special(eer,ver)
  plot1special(eir,vir)
  eei = zeroint(1)
  vei = zeroint(1)
  eii = zeroint(1)
  vii = zeroint(1)
  print 'min eta error',min(eer,eei),'is at cmp',eei[0]
  print 'min velocity error',min(ver,vei),'is at cmp',vei[0]
  print 'min eta interval error',min(eir,eii),'is at cmp',eii[0]
  print 'min velocity interval error',min(vir,vii),'is at cmp',vii[0]

def nonHypLoop(data):
  s1,so,g,vnmo,eta,wb,vpt,etat = data
  n1 = s1.count
  d1 = s1.delta
  vpl = [1.492,3.00,0.004]
  epl = [0.000,0.15,0.004]
  vpli = [1.400,4.50]
  epli = [-0.00,0.40]
  wbi = int(wb/d1)
  dr = 0.03
  r1min,r1max = -0.05,0.55
  r2min,r2max = -0.05,0.15
  DvDt = d1/vpl[2]
  DeDt = d1/epl[2]
  pl1,pu1 = 0.5,2.0
  pl2,pu2 = 1.0,1.0
  nv = int((vpl[1]-vpl[0])/vpl[2])
  ne = int((epl[1]-epl[0])/epl[2])
  sv = Sampling(nv,vpl[2],vpl[0])
  se = Sampling(ne,epl[2],epl[0])
  vl1,vl2 = makeIntconsts(vpli[0],vpli[1],s1,wbi,pl1,pu1)
  el1,el2 = fillfloat(epli[0],n1),fillfloat(epli[1],n1)
  p1min,p1max,p2min,p2max,vel,veu,eel,eeu = getTimeVaryingConstraints(s1,wbi,vl1,vl2,el1,el2,dr)
  p1max,p1min = mul(p1min,DvDt),mul(p1max,-DvDt)
  p2max,p2min = mul(p2min,DeDt),mul(p2max,DeDt)
  s = Semblance(s1,so,fpeak)
  s.setStretchMute(0.75)
  ps = s.applyNonHypV(g,sv,se,vel,veu,eel,eeu)
  #ps = s.applyNonHypV(g,sv,se)
  #stk = stackGats2(g,s1)
  #ps = mul(ps,ps)
  #ps = mul(ps,ps)
  vln1,vln2 = makeIntconsts(vpl[0],vpl[1],s1,wbi,0.8,1.2)
  stk = stackSem2(pow(ps,4),vln1,vln2,sv)
  ds = DynamicSolver(sv.count,se.count,r1min,r1max,r2min,r2max,dr,wbi)
  ds.setTimeVaryingStrainLimits(p1min,p1max,p2min,p2max)
  #ds.setWbParam(sv.indexOfNearest(1.492),se.indexOfNearest(0.0))
  #up1,up2 = ds.findSolution(stk,ps,ng)
  up1,up2 = ds.findSolution(stk,ps)
  #up1,up2 = ds.findSolution(ps)
  #cgs = ds.getCgs()
  #plot1(stk,s1,name="stk-pts",cgs=cgs)
  up1 = add(sv.first,mul(up1,sv.delta))
  up2 = add(se.first,mul(up2,se.delta))
  uv = up1
  ue = up2
  #up1m,up2m = s.pickMaxSemblance(ps,se,sv)
  #q  = s.flattenGatherNonHypE(g,up1,up2)
  #pduv = mul(div(sub(uv,vnmo),vnmo),100)
  vnmoi = dixInversionVnmo(s1,uv)
  etai = dixInversionEta(s1,ue,uv,vnmoi)
  pdata = vnmo,eta,uv,ue,vnmoi,etai,vpt,etat
  return pdata

def plotDq(data):
  s1,so,g,vnmo,eta,wb,vp = data
  n1 = s1.count
  vrmax = mul(vnmo,1.2)
  vrmin = mul(vnmo,0.8)
  vpmax = mul(vp,1.2)
  vpmin = mul(vp,0.8)
  def plot(s1,so,g,vnmo,wb,vp,eta):
    vrms2 = mul(vnmo,vnmo)
    vp2 = mul(vp,vp)
    dvdt = zerofloat(n1)
    dt = s1.delta
    dts = wb
    wbi = int(wb/dt)
    for i in range(wbi,n1):
      dvdt[i] = (vnmo[i]-vnmo[i-1])/(dt)
      #dvdt[i] = dt*(vp2[i]/dts-vrms2[i]/(2*sqrt(dts)))/(2*vnmo[i])/1000
      #dts += dt
    dq1dt = div(mul(-2,dvdt),mul(vrms2,vnmo))
    plot1(dvdt,s1,name="dvdt",hlabel="dv/dt (km/s^2)")
    plot1(dq1dt,s1,name="dq1dt",hlabel="1/km")
    dedt = zerofloat(n1)
    for i in range(wbi,n1):
      dedt[i] = (eta[i]-eta[i-1])/(dt)
    dq2dt = add(mul(2,mul(eta,dq1dt)),mul(div(2,vrms2),dedt))
    to = floats(s1.getValues())
    xmax = so.getFirst()
    xmax2 = xmax*xmax
    toxm = mul(2,sqrt(mul(to,to)+div(xmax2,vrms2)))
    ddtvdt = sub(1,div(sub(mul(2,to),mul(dvdt,div(xmax2,mul(2,mul(vrms2,vnmo))))),toxm))
    plot1(dedt,s1,name="dedt",hlabel="1/s")
    plot1(dq2dt,s1,name="dq2dt",hlabel="1/km")
    plot1(ddtvdt,s1,name="ddtvdt")
  #plot(s1,so,g,vrmax,wb,vpmax,eta)
  #plot(s1,so,g,vrmin,wb,vpmin,eta)
  plot(s1,so,g,vnmo,wb,vp,eta)
  vmax = max(vnmo); vmin = min(vnmo)
  dtmin = -sqrt(100/(vmin*vmin))
  dtmax = 1-sqrt(1+100/(vmax*vmax))
  print vmin
  print vmax
  print dtmin
  print dtmax
  ndt = int((dtmax-dtmin)/s1.delta)
  dta = rampfloat(dtmin,s1.delta,ndt)
  plot1(dta,name="dta 1 sec")


def goHyp(data):
  s1,so,g,vnmo,eta,wb,vpt = data
  vpl = [1.40,3.00,0.005]
  r1min,r1max = -0.50,0.10
  dr = 0.02
  nv = int((vpl[1]-vpl[0])/vpl[2])
  sv = Sampling(nv,vpl[2],vpl[0])
  # Test Data
  #g = Gather(s1,so,sv)
  #vnmo = g.makeLinearParameter(vpl[0],vpl[1])
  #g = g.modelGatherHyp(30,1e6,vnmo)
  s = Semblance(s1,so,20)
  s.setStretchMute(0.50)
  #ps = s.applyHypDt(g)
  #ps = s.applyHypQ(g)
  ps = s.applyHypV(g,vpl[0],vpl[1],vpl[2])
  wbi = int(wb/s1.delta)
  stk = stackSem(ps)
  dw = DynamicSolver(sv.count,r1min,r1max,dr,wbi)
  ac = dw.accumulate(transpose(ps))
  uv = dw.findSolution(transpose(ps))
  #uv = dw.findSolution(stk,transpose(ps))
  uv = add(sv.first,mul(uv,sv.delta))
  q = s.flattenGatherHyp(g,uv)
  #cv = s.pickMaxSemblance(ps)
  pduv = mul(div(sub(uv,vnmo),vnmo),100)
  #pdcv = mul(div(sub(cv,vnmo),vnmo),100)
  print 'plotting'
  cmplims=[-6.5,s1.first,so.last,s1.last]
  s1c = Sampling(len(ac),s1.delta,s1.first)
  plotModelsCmp(g,s1,so,vnmo,eta,'modelscmp',clims=[cmplims[0],cmplims[2]])
  plot2(transpose(ac),s1c,sv,cmap=jet,name='accumulated semblance',sem=True)
  #plot1(vnmo,s1,uv,cv,name="vnmo")
  #plot2(ps,s1,sv,u=cv,name='semblance max pick',sem=True)
  plot1(pduv,s1,name='vnmo % difference dp')
  #plot1(pdcv,s1,name='vnmo % difference max')
  plot1(vnmo,s1,uv,name="vnmo and vnmo estimated")
  plot2(ps,s1,sv,u=vnmo,u2=uv,name='semblance comp',sem=True,size=[728,753])
  plot2(ps,s1,sv,u=uv,name='semblance est',sem=True,size=[728,753])
  plot2(ps,s1,sv,name='semblance',sem=True,size=[728,753])
  plot2(q,s1,so,name='flat cmp',cmp=True,perc=99)
  plot2(g,s1,so,name='cmp',cmp=True,perc=99,lims=cmplims)

def hypTest(data):
  s1,so,g,vnmo,eta,wb,vpt,etat = data
  pdata = hypTestV(data)
  hypPlots(pdata)
  #pdata = hypTestDt(data)
  #hypPlots(pdata)
  #pdata = hypTestQ(data)
  #hypPlots(pdata)

def hypPlots(pdata):
  s1,so,g,vnmo,eta,wb,pduv,uv,up,vnmop,sv,q,ps,cv = pdata
  print 'plotting'
  cmplims=[so.first,s1.first,so.last,s1.last]
  plotModelsCmp(g,s1,so,vnmo,eta,'modelscmp',clims=[cmplims[0],cmplims[2]])
  plot1(pduv,s1,name='vnmo % difference dp')
  plot1(vnmo,s1,uv,name="vnmo and vnmo estimated")
  plot1comp(s1,vnmo,uv,pduv,name="vnmo comparison",vel=True)
  plot2(ps,s1,sv,u=vnmop,u2=up,name='semblance comp',sem=True,size=[728,753])
  plot2(ps,s1,sv,u2=up,name='semblance est',sem=True,size=[728,753])
  plot2(ps,s1,sv,u=vnmop,u2=cv,name='semblance max comp',sem=True,size=[728,753])
  plot2(ps,s1,sv,u2=cv,name='semblance max',sem=True,size=[728,753])
  plot2(ps,s1,sv,name='semblance',sem=True,size=[728,753])
  plot2(q,s1,so,name='flat cmp zoom',cmp=True,perc=99,lims=[so.first,5.5,so.last,7.0])
  plot2(q,s1,so,name='flat cmp',cmp=True,perc=99)
  plot2(g,s1,so,name='cmp',cmp=True,perc=99,lims=cmplims)

def hypTestV(data):
  s1,so,g,vnmo,eta,wb,vpt,etat = data
  r1min,r1max,dr1 = -0.10,0.60,0.02;
  vpl = [1.402,3.00,0.005]
  nv = int((vpl[1]-vpl[0])/vpl[2])
  sv = Sampling(nv,vpl[2],vpl[0])
  s = Semblance(s1,so,20)
  s.setStretchMute(0.75)
  ps = s.applyHypV(g,vpl[0],vpl[1],vpl[2])
  wbi = int(wb/s1.delta)
  stk = stackSem(ps)
  dw = DynamicSolver(sv.count,r1min,r1max,dr1,wbi)
  dw.setWbParam(sv.indexOfNearest(1.492))
  #uv = dw.findSolution(transpose(ps))
  uv = dw.findSolution(stk,transpose(ps))
  uv = add(sv.first,mul(uv,sv.delta))
  up = uv
  vnmop = vnmo
  q = s.flattenGatherHyp(g,uv)
  #q = s.flattenGatherHyp(g,vnmo)
  cv = s.pickMaxSemblance(ps,sv)
  pduv = mul(div(sub(uv,vnmo),vnmo),100)
  pdata = s1,so,g,vnmo,eta,wb,pduv,uv,up,vnmop,sv,q,ps,cv
  return pdata

def hypTestDt(data):
  s1,so,g,vnmo,eta,wb,vpt = data
  r1min,r1max,dr1 = -0.50,0.50,0.02; 
  vpl = [1.40,2.80,0.005]
  s = Semblance(s1,so,20)
  s.setStretchMute(0.50)
  vmin = vpl[0]
  vmax = vpl[1]
  xmax = so.first
  tmax = s1.last
  dtmax =  -sqrt(xmax*xmax/(vmin*vmin));
  #dtmin = -sqrt(xmax*xmax/(vmax*vmax));
  dtmin = 0
  dt = s1.delta;
  nv = inro((dtmin-dtmax)/dt);
  sv = Sampling(nv,dt,dtmax);
  ps = s.applyHypDt(g,sv)
  plot2(ps,s1,sv,name='semblance',sem=True,size=[728,753])
  wbi = int(wb/s1.delta)
  stk = stackSem(ps)
  dw = DynamicSolver(sv.count,r1min,r1max,dr1,wbi)
  ac = dw.accumulate(transpose(ps))
  #up = dw.findSolution(transpose(ps))
  up = dw.findSolution(stk,transpose(ps))
  up = add(sv.first,mul(up,sv.delta))
  to = floats(s1.getValues())
  uv = sqrt(div(xmax*xmax,sub(pow(sub(to,up),2),to)))
  vnmop = mul(-1,sqrt(div(xmax*xmax,mul(vnmo,vnmo))))
  q = s.flattenGatherHyp(g,uv)
  pduv = mul(div(sub(uv,vnmo),vnmo),100)
  pdata = s1,so,g,vnmo,eta,wb,pduv,uv,up,vnmop,sv,q,ps
  return pdata

def hypTestQ(data):
  s1,so,g,vnmo,eta,wb,vpt = data
  r1min,r1max,dr1 = -0.05,0.10,0.02;
  vpl = [1.40,3.00,0.005]
  s = Semblance(s1,so,20)
  s.setStretchMute(0.50)
  vmin = vpl[0]; vmax = vpl[1]; dq = vpl[2];
  qmin = 1.0/(vmax*vmax);
  qmax = 1.0/(vmin*vmin);
  nq = inro((qmax-qmin)/dq)
  sq = Sampling(nq,dq,qmin);
  ps = s.applyHypQ(g,vmin,vmax,dq)
  wbi = int(wb/s1.delta)
  stk = stackSem(ps)
  dw = DynamicSolver(sq.count,r1min,r1max,dr1,wbi)
  ac = dw.accumulate(transpose(ps))
  #uv = dw.findSolution(transpose(ps))
  up = dw.findSolution(stk,transpose(ps))
  up = add(sq.first,mul(up,sq.delta))
  uv = div(1,sqrt(up))
  vnmop = div(1,mul(vnmo,vnmo))
  q = s.flattenGatherHyp(g,uv)
  pduv = mul(div(sub(uv,vnmo),vnmo),100)
  pdata = s1,so,g,vnmo,eta,wb,pduv,uv,up,vnmop,sq,q,ps
  return pdata

def normalizeRMSlocal(x,sig):
  return WTutils().localRMSnorm(x,sig)

################################################################################
# Nonhyperbolic methods

def nonHypTest(data):
  s1,so,g,vnmo,eta,wb,vpt,etat = data
  pdata = nonHypV(data)
  nonHypPlots(pdata)
  #pdata = nonHypDt(data)
  #nonHypPlots(pdata)
  #pdata = nonHypQ(data)
  #nonHypPlots(pdata)
  return pdata

def nonHypQ(data):
  s1,so,g,vnmo,eta,wb,vpt = data
  vpl = [1.492,3.00,0.005]
  epl = [0.00,0.15,0.005]
  r1min,r1max = -0.10,0.50
  r2min,r2max = -0.10,0.50
  dr = 1.0
  dq1 = 0.002
  dq2 = dq1
  vmin = vpl[0]; vmax = vpl[1];
  emin = epl[0]; emax = epl[1];
  q1min = 1.0/(vmax*vmax)
  q1max = 1.0/(vmin*vmin)
  q2min = 2.0*q1min*emin
  q2max = 2.0*q1max*emax
  nq1 = inro((q1max-q1min)/dq1)
  nq2 = inro((q2max-q2min)/dq2)
  sq1 = Sampling(nq1,dq1,q1min)
  sq2 = Sampling(nq2,dq2,q2min)
  # Test Data
  nv = int((vpl[1]-vpl[0])/vpl[2])
  ne = int((epl[1]-epl[0])/epl[2])
  sv = Sampling(nv,vpl[2],vpl[0])
  se = Sampling(ne,epl[2],epl[0])
  gt = Gather(s1,so,sv,se)
  vnmo = gt.makeLinearParameter(vpl[0],vpl[1])
  eta = gt.makeLinearParameter(epl[0],epl[1])
  g = gt.modelGatherNonHyp(30,1e6,vnmo,eta)
  s = Semblance(s1,so,30)
  s.setStretchMute(0.50)
  ps = s.applyNonHypQ(g,sq1,sq2)
  ps = mul(ps,ps)
  wbi = int(wb/s1.delta)
  wbi=0;
  #stk = stackSem2(ps,vl1,vl2,sv)
  ds = DynamicSolver(sq1.count,sq2.count,r1min,r1max,r2min,r2max,dr,wbi)
  #ds.setWbParam(sq1.indexOfNearest(1/(1.492*1.492)),sq2.indexOfNearest(0.0))
  #up1,up2 = ds.findSolution(stk,ps)
  up1,up2 = ds.findSolution(ps)
  cgs = ds.getCgs()
  up1 = add(sq1.first,mul(up1,sq1.delta))
  up2 = add(sq2.first,mul(up2,sq2.delta))
  uv = div(1,sqrt(up1))
  ue = mul(0.5,div(up2,up1))
  q  = s.flattenGatherNonHypE(g,uv,ue)
  qe = s.flattenGatherNonHypE(g,vnmo,eta)
  pduv = mul(div(sub(uv,vnmo),vnmo),100)
  vnmop = div(1,mul(vnmo,vnmo))
  etap = mul(2,mul(eta,vnmop))
  pdata = s1,so,g,vnmo,eta,wb,pduv,uv,ue,up1,up2,vnmop,etap,sq1,sq2,q,qe,ps,cgs
  return pdata

#GOTO
def nonHypV(data):
  s1,so,g,vnmo,eta,wb,vpt,etat = data
  n1 = s1.count
  d1 = s1.delta
  vpl = [1.492,3.00,0.004]
  epl = [0.000,0.15,0.004]
  vpli = [1.400,4.50]
  epli = [-0.00,0.40]
  wbi = int(wb/d1)
  dr = 0.03
  #r1min,r1max = -0.50,0.50
  #r2min,r2max = -0.50,0.50
  r1min,r1max = -0.05,0.55
  r2min,r2max = -0.05,0.15
  DvDt = d1/vpl[2]
  DeDt = d1/epl[2]
  pl1,pu1 = 0.5,2.0
  pl2,pu2 = 1.0,1.0
  nv = int((vpl[1]-vpl[0])/vpl[2])
  ne = int((epl[1]-epl[0])/epl[2])
  sv = Sampling(nv,vpl[2],vpl[0])
  se = Sampling(ne,epl[2],epl[0])
  dedt = zerofloat(n1)
  dvdt = zerofloat(n1)
  for i in range(n1):
    dedt[i] = ((eta[i]-eta[i-1])/d1)*DeDt
    dvdt[i] = ((vnmo[i]-vnmo[i-1])/d1)*DvDt
  dedt[0] = dedt[1]
  dvdt[0] = dvdt[1]
  plot1(dedt,s1,name="exact dedt slopes")
  plot1(dvdt,s1,name="exact dvdt slopes")
  vl1,vl2 = makeIntconsts(vpli[0],vpli[1],s1,wbi,pl1,pu1)
  #el1,el2 = makeIntconsts(epli[0],epli[1],s1,wbi,pl2,pu2)
  el1,el2 = fillfloat(epli[0],n1),fillfloat(epli[1],n1)
  p1min,p1max,p2min,p2max,vel,veu,eel,eeu = getTimeVaryingConstraints(s1,wbi,vl1,vl2,el1,el2,dr)
  p1max,p1min = mul(p1min,DvDt),mul(p1max,-DvDt)
  p2max,p2min = mul(p2min,DeDt),mul(p2max,DeDt)
  plot1(vl1,s1,vl2,name="vnmo interval bounds")
  plot1(el1,s1,el2,name="eta interval bounds",hlim=[epli[0]-0.1,epli[1]+0.1])
  plot1(floats(p1min),s1,floats(p1max),name="vnmo derivative bounds")
  plot1(floats(p2min),s1,floats(p2max),name="eeta derivative bounds")
  plot1(floats(p1min),s1,floats(p1max),dvdt,name="vnmo derivative bounds comp")
  plot1(floats(p2min),s1,floats(p2max),dedt,name="eeta derivative bounds comp")
  # Test Data
  #gt = Gather(s1,so,sv,se)
  #vnmo = gt.makeLinearParameter(vpl[0],vpl[1]-0.1)
  #eta = gt.makeLinearParameter(epl[0],epl[1]-.03)
  #g = gt.modelGatherNonHyp(30,1e6,vnmo,eta)
  s = Semblance(s1,so,fpeak)
  s.setStretchMute(0.75)
  #ps = s.applyNonHypV(g,sv,se,vel,veu,eel,eeu)
  ps = s.applyNonHypV(g,sv,se)
  #stk = stackGats2(g,s1)
  #ps = mul(ps,ps)
  #ps = mul(ps,ps)
  vln1,vln2 = makeIntconsts(vpl[0],vpl[1],s1,wbi,0.8,1.2)
  stk = stackSem2(pow(ps,4),vln1,vln2,sv)
  #ng = 40
  #print 'ng=',ng
  #wbi=0;
  ds = DynamicSolver(sv.count,se.count,r1min,r1max,r2min,r2max,dr,wbi)
  ds.setTimeVaryingStrainLimits(p1min,p1max,p2min,p2max)
  #ds.setWbParam(sv.indexOfNearest(1.492),se.indexOfNearest(0.0))
  #up1,up2 = ds.findSolution(stk,ps,ng)
  up1,up2 = ds.findSolution(stk,ps)
  #up1,up2 = ds.findSolution(ps)
  cgs = ds.getCgs()
  plot1(stk,s1,name="stk-pts",cgs=cgs)
  up1 = add(sv.first,mul(up1,sv.delta))
  up2 = add(se.first,mul(up2,se.delta))
  uv = up1 
  ue = up2 
  up1m,up2m = s.pickMaxSemblance(ps,se,sv)
  #q  = s.flattenGatherNonHypE(g,up1m,up2m)
  q  = s.flattenGatherNonHypE(g,up1,up2)
  #qe = s.flattenGatherNonHypE(g,vnmo,eta)
  qe = s.flattenGatherNonHypE(g,fillfloat(1.492,len(vnmo)),fillfloat(0.0,len(eta)))
  pduv = mul(div(sub(uv,vnmo),vnmo),100)
  vnmop = vnmo
  etap = eta
  vnmoi = dixInversionVnmo(s1,uv)
  etai = dixInversionEta(s1,ue,uv,vnmoi)
  pdata = s1,so,g,vnmo,eta,wb,pduv,uv,ue,up1,up2,vnmop,etap,sv,se,q,qe,ps,cgs,vel,veu,eel,eeu,up1m,up2m,vnmoi,etai,vpt,etat,vl1,vl2,el1,el2
  up2 = div(sub(up2,se.first),se.delta)
  dedt = zerofloat(n1)
  eeta = div(sub(eta,se.first),se.delta)
  deedt = zerofloat(n1)
  for i in range(n1):
    dedt[i] = (up2[i]-up2[i-1])
    deedt[i] =(eeta[i]-eeta[i-1])
  dedt[0] = dedt[1]
  deedt[0] = deedt[1]
  plot1(dedt,s1,deedt,name="computed dedt slopes")
  return pdata

def dixInversionVnmo(s1,vnmo):
  n = s1.count 
  vnmoi = zerofloat(n)
  vnmo2 = mul(vnmo,vnmo)
  for i in range(1,n):
    t0 = s1.getValue(i-1)
    t1 = s1.getValue(i)
    vnmoi[i] = (vnmo2[i]*t1 - vnmo2[i-1]*t0)/(t1-t0)
  vnmoi[0] = vnmoi[1]
  return sqrt(vnmoi)

def dixInversionEta(s1,eta,vnmo,vnmoi):
  n = s1.count
  etai = zerofloat(n)
  vnmoi4 = mul(mul(mul(vnmoi,vnmoi),vnmoi),vnmoi)
  vnmo4 = mul(mul(mul(vnmo,vnmo),vnmo),vnmo)
  for i in range(1,n):
    t0 = s1.getValue(i-1)
    t1 = s1.getValue(i)
    etai[i] = ((vnmo4[i]*(1.0+8.0*eta[i])*t1-vnmo4[i-1]*(1.0+8.0*eta[i-1])*t0)/(t1-t0)-vnmoi4[i])/(8.0*vnmoi4[i])
  etai[0] = etai[1]
  return etai

def getTimeVaryingConstraints(s1,wbi,vl1,vl2,el1,el2,dr):
  nt = s1.count
  dt = s1.delta
  ne = nt-wbi
  p1min = zerodouble(nt)
  p1max = zerodouble(nt)
  p2min = zerodouble(nt)
  p2max = zerodouble(nt)
  leeta = zerofloat(nt)
  ueeta = zerofloat(nt)
  to = s1.getValues()
  def vnmo2(vp,dt,n,to):
    vnmo2 = zerodouble(n)
    vps = mul(vp,vp)
    vnmo2[0] = vps[0]*dt
    for i in range(1,n):
     vnmo2[i] = vnmo2[i-1]+vps[i]*dt
    dts = dt
    for i in range(n):
      vnmo2[i] /= dts
      dts += dt
    return vnmo2
  # compute dvdt and dedt min and max
  lvnmo2 = vnmo2(vl1,dt,nt,to)
  uvnmo2 = vnmo2(vl2,dt,nt,to)
  lvnmo1 = sqrt(lvnmo2) 
  uvnmo1 = sqrt(uvnmo2)
  vl1s = mul(vl1,vl1)
  vl2s = mul(vl2,vl2)
  #plot1(vl1,s1,floats(lvnmo1),name="vnmo vs Vnmo")
  for i in range(wbi,nt):
    p1min[i] = (vl1s[i] - lvnmo2[i])/(2*to[i]*lvnmo1[i])
    p1max[i] = (vl2s[i] - uvnmo2[i])/(2*to[i]*uvnmo1[i])
    p2min[i] = (p1min[i]/lvnmo1[i] - 1/to[i] + vl1s[i]*vl1s[i]*(1+8*el1[i]))/(8*to[i]*vl1s[i]*vl1s[i])
    p2max[i] = (p1max[i]/uvnmo1[i] - 1/to[i] + vl2s[i]*vl2s[i]*(1+8*el2[i]))/(8*to[i]*vl2s[i]*vl2s[i])
  p1minm,p2minm = mul(p1max,1),mul(p2max,1)
  p1maxm,p2maxm = mul(p1min,1),mul(p2min,-1)#filldouble(dr,nt),filldouble(dr,nt)
  vpq1 = mul(vl1s,vl1s)
  vpq2 = mul(vl2s,vl2s)
  leeta[0] = vpq1[0]*(1+8*el1[0])*dt
  ueeta[0] = vpq2[0]*(1+8*el2[0])*dt
  for i in range(1,nt):
    leeta[i] = leeta[i-1]+vpq1[i]*(1+8*el1[i])*dt
    ueeta[i] = ueeta[i-1]+vpq2[i]*(1+8*el2[i])*dt
  v1rms4 = zerofloat(nt);
  v2rms4 = zerofloat(nt); dts = dt
  for i in range(nt):
    v1rms4[i] = lvnmo2[i]*lvnmo2[i]*dts
    v2rms4[i] = uvnmo2[i]*uvnmo2[i]*dts
    dts += dt
  leeta = mul(sub(div(leeta,v1rms4),1),0.125)
  ueeta = mul(sub(div(ueeta,v2rms4),1),0.125)
  return p1minm,p1maxm,p2minm,p2maxm,floats(lvnmo1),floats(uvnmo1),leeta,ueeta
 
def makeIntconsts(vpl,vpu,s1,wbi,lp,up):
  nt = s1.count
  ne = nt-wbi
  vpl1 = lp*vpl
  vpl2 = up*vpl
  vpu1 = lp*vpu
  vpu2 = up*vpu
  dv1 = (vpu1-vpl1)/ne
  dv2 = (vpu2-vpl2)/ne
  vp1 = fillfloat(vpl1,nt)
  vp2 = fillfloat(vpl2,nt)
  for i in range(ne):
    vp1[i+wbi] = vpl1+i*dv1
    vp2[i+wbi] = vpl2+i*dv2
  return vp1,vp2

def nonHypPlots(pdata):
  s1,so,g,vnmo,eta,wb,pduv,uv,ue,up1,up2,vnmop,etap,\
    sv,se,q,qe,ps,cgs,vel,veu,eel,eeu,up1m,up2m,vnmoi,etai,vpt,etat,vl1,vl2,el1,el2 = pdata
  vve=True
  print 'plotting'
  ut = floats(s1.getValues())
  cmplims=[so.first,s1.first,so.last,s1.last]
  pduvm = mul(div(sub(up1m,vnmo),vnmo),100)
  #plotModelsCmp(g,s1,so,vnmo,eta,'modelscmp',clims=[cmplims[0],cmplims[2]])
  plot1(sub(eta,up2),s1,name='eta dp difference')
  #plot1(vnmo,s1,uv,cv,name="vnmo")
  plot1(eel,s1,eeu,name="eta eff bounds")
  plot1(vel,s1,veu,name="vnmo eff bounds")
  plot1(pduv,s1,name='vnmo % difference dp')
  plot1(eta,s1,ue,name="eta and eta estimated",colm=3)
  plot1(vnmo,s1,uv,name="vnmo and vnmo estimated",colm=3)
  plot1(vnmo,s1,uv,name="vnmo and vnmo estimated pts",cgs=cgs)
  plot1(eta,s1,ue,name="eta and eta estimated pts",cgs=cgs)
  plot1comp(s1,vnmo,uv,pduv,name="vnmo nh comp",vel=True)
  plot1comp(s1,eta,ue,sub(ue,eta),name="eta nh comp",eta=True,clims2=[-0.025,0.025])
  plot1comp(s1,vnmo,up1m,pduvm,name="vnmo and vnmo max comp",vel=True)
  plot1comp(s1,eta,up2m,sub(up2m,eta),name="eta and eta max comp",eta=True)
  plot1(mul(vpt,0.001),s1,vnmoi,name="vnmo interval vs dix",bounds=[vl1,vl2])
  plot1(etat,s1,etai,name="eta interval vs dix",bounds=[el1,el2])
  plot1(sub(mul(vpt,0.001),vnmoi),s1,name="vnmoi errors")
  plot1(sub(etat,etai),s1,name="etai errors")
  ll = [3.0,4.0,5.0,6.0,7.0,8.0,9.0] # 3-9 sec.
  #ll = [2.5,3.5,4.5,5.5,6.5,7.5,8.5] # 3-9 sec.
  #ll = cgs
  llims = [0.0,0.60]
  for l in ll:
    ls = str(l)
    l  = round(l/s1.delta)
    plot2(ps[l],s1=sv,s2=se,name='vnmo vs eta t='+ls,vve=vve,\
          mx=[up1[l],up2[l],vnmop[l],etap[l]],\
          sem=True,cmin=llims[0],cmax=llims[1])
  plot2(q,s1,so,name='flat nh cmp',cmp=True,perc=99)
  plot2(q,s1,so,name='flat nh cmp zoom',cmp=True,perc=99,lims=[so.first,5.5,so.last,7.0])
  plot2(qe,s1,so,name='flat real nh cmp',cmp=True,perc=99)
  plot2(g,s1,so,name='cmp',cmp=True,perc=99,lims=cmplims)
  ps2 = ps
  ps = zerofloat(len(ps2),len(ps2[0]),len(ps2[0][0]))
  Transpose.transposeP(ps2,ps)
  plot3(ps,cmap=jet,s1=s1,s2=se,s3=sv,u=[up1,up2,ut],name='eta semblance est')
  plot3(ps,cmap=jet,s1=s1,s2=se,s3=sv,name='eta semblance clip',perc=99)
  plot3(ps,cmap=jet,s1=s1,s2=se,s3=sv,name='eta semblance')
  plot3(ps,cmap=jet,s1=s1,s2=se,s3=sv,u=[vnmop,etap,ut],name='eta semblance exact')


def goNonHyp(data):
  s1,so,g,vnmo,eta,wb = data
  vpl = [1.40,3.00,0.005]
  epl = [0,0.15,0.005]
  r1min,r1max = -0.50,0.10
  r2min,r2max = -0.40,0.08
  dr = 0.02
  nv = int((vpl[1]-vpl[0])/vpl[2])
  ne = int((epl[1]-epl[0])/epl[2])
  sv = Sampling(nv,vpl[2],vpl[0])
  se = Sampling(ne,epl[2],epl[0])
  # Test Data
  #g = Gather(s1,so,sv)
  #vnmo = g.makeLinearParameter(vpl[0],vpl[1])
  #g = g.modelGatherHyp(30,1e6,vnmo)
  s = SemblanceOld(s1,so,sv,se,20)
  s.setStretchMute(0.50)
  ps = s.applyNonHyp(g)
  wbi = int(wb/s1.delta)
  stk = stackSem2(ps)
  ds = DynamicSolver(sv.count,se.count,r1min,r1max,r2min,r2max,dr,wbi)
  uv,ue = ds.findSolution(stk,ps)
  ue = add(se.first,mul(ue,se.delta))
  uv = add(sv.first,mul(uv,sv.delta))
  ut = floats(s1.getValues())
  q  = s.flattenGatherNonHypE(g,uv,ue)
  qe = s.flattenGatherNonHypE(g,vnmo,eta)
  cv,ce = s.pickMaxSemblance(ps,se)
  pduv = mul(div(sub(uv,vnmo),vnmo),100)
  pdcv = mul(div(sub(cv,vnmo),vnmo),100)
  vve=True
  print 'plotting'
  cmplims=[-6.5,s1.first,so.last,s1.last]
  plotModelsCmp(g,s1,so,vnmo,eta,'modelscmp',clims=[cmplims[0],cmplims[2]])
  plot1(stk,s1,name="sem stack")
  plot1(sub(eta,ue),s1,name='eta dp difference')
  plot1(vnmo,s1,uv,cv,name="vnmo")
  plot1(pduv,s1,name='vnmo % difference dp')
  plot1(pdcv,s1,name='vnmo % difference max')
  plot1(eta,s1,ue,name="eta and eta estimated")
  plot1(vnmo,s1,uv,name="vnmo and vnmo estimated")
  plot2(ps[500],s1=sv,s2=se,name='vnmo vs eta t=500',vve=vve)
  plot2(q,s1,so,name='flat cmp',cmp=True,perc=98)
  plot2(qe,s1,so,name='flat real cmp',cmp=True,perc=98)
  plot2(g,s1,so,name='cmp',cmp=True,perc=99,lims=cmplims)
  plot3(ps,cmap=jet,s1=sv,s2=se,s3=s1,u=[uv,ue,ut],perc=99,\
                                      name='eta semblance',paper='bpnhsu')
  plot3(ps,cmap=jet,s1=sv,s2=se,s3=s1,name='eta semblance')

def stackSem(s):
  n2 = len(s)
  n1 = len(s[0])
  stk = zerofloat(n1)
  for i in range(n2):
    stk = add(stk,s[i])
  return stk

def stackSem2(s,vl1,vl2,sv):
  n3 = len(s)
  n2 = len(s[0])
  n1 = len(s[0][0])
  stk = zerofloat(n3)
  for i3 in range(n3):
    il = sv.indexOfNearest(vl1[i3])
    iu = sv.indexOfNearest(vl2[i3])
    for i2 in range(n2):
      for i1 in range(il,iu):
        stk[i3] += s[i3][i2][i1]
  return stk

def stackGats2(g,s1):
  f = g[-1]
  es = ExponentialSmoother(1.0/(fpeak*s1.delta))
  es.apply(f,f)
  nt = len(f)
  ht = HilbertTransformFilter()
  ft = zerofloat(nt)
  ht.apply(nt,f,ft)
  h = sqrt(add(pow(f,2.0),pow(ft,2.0)));
  return h

def getData(ix2,stack=None,xmax=None):
  makeSubset(ix2)
  s1,so,g,vp,de,ep,th,wb = readSubset(ix2,xmax)
  if stack:
    makeSubset(ix2-1)
    makeSubset(ix2+1)
    dat1 = readSubset(ix2-1)
    dat2 = readSubset(ix2+1)
    g1 = dat1[2]; g2 = dat2[2];
    g = add(g1,add(g2,g))
  eta = div(sub(ep,de),add(1,mul(2,de)))
  #plot1(mul(vp,sqrt(add(1,mul(2,de)))),s1,name="vnmo exact interval")
  #plot1(eta,s1,name="eta exact interval")
  vrms,eeta = makeEffectiveParams(s1,vp,eta,de)
  #plot1(vp,s1,hlabel="vp")
  #plot1(ep,s1,hlabel="epsilon")
  #plot1(de,s1,hlabel="delta")
  #plot1(eta,s1,hlabel="eta")
  #plot1(vrms,s1,hlabel="vrms")
  #plot1(eeta,s1,hlabel="eeta")
  #plot2(g,s1,so,cmap=gray,perc=99,cmp=True)
  #vpz = getModel("vp")
  #epz = getModel("epsilon")
  #dez = getModel("delta")
  #thz = getModel("theta")
  #s1z,s2 = getSamplings()
  #lims = [s2.first,s1z.first,20,s1z.last]
  #vpz = mul(vpz,0.001)
  #plot2(vpz,s1z,s2,name='vpz',vpz=True,paper="vpz",lims=lims,cbar="Velocity (km/s)",cmppts=ix2)
  #plot2(epz,s1z,s2,name='epz',vpz=True,paper="epz",lims=lims,cbar="Epsilon",cmppts=ix2)
  #plot2(dez,s1z,s2,name='dez',vpz=True,paper="dez",lims=lims,cbar="Delta",cmppts=ix2)
  #plot2(thz,s1z,s2,name='thz',vpz=True,paper="thz",lims=lims,cbar="Theta",cmppts=ix2,cmap=rwb,perc=98)
  return s1,so,g,vrms,eeta,wb,vp,eta

def makeEffectiveParams(s1,vp,eta,de):
  n = len(vp)
  dt = s1.delta
  vrms = zerofloat(n)
  eeta = zerofloat(n)
  vpa = mul(vp,sqrt(add(1,mul(2,de))))
  vps = mul(vpa,vpa)
  vpq = mul(vps,vps)
  to = floats(s1.getValues())
  vrms[0] = vps[0]*dt
  for i in range(1,n):
    vrms[i] = vrms[i-1]+vps[i]*dt
  dts = dt
  for i in range(n):
    vrms[i] /= dts
    dts += dt
  eeta[0] = vpq[0]*(1+8*eta[0])*dt
  for i in range(1,n):
    eeta[i] = eeta[i-1]+vpq[i]*(1+8*eta[i])*dt
  vrms4 = zerofloat(n); dts = dt
  for i in range(n):
    vrms4[i] = vrms[i]*vrms[i]*dts
    dts += dt
  eeta = mul(sub(div(eeta,vrms4),1),0.125)
  return mul(sqrt(vrms),0.001),eeta

################################################################################
# Tests


def goHypTest():
  goTestData()
  fpeak = 30
  snr = 1.0e6
  #snr = 5.0
  #snr = 0.1
  rmin,rmax = -0.10,0.10
  dr = 0.01
  print 'vcmin=',rmin*sv.delta,',vcmax=',rmax*sv.delta
  print 'sparse grid=',st.delta/dr,'seconds'
  print 'sigma=',1.0/(fpeak*st.delta)
  print 'model a gather'
  g = Gather(st,sx,sv)
  vnmo = g.makeLinearParameter(2.1,3.5)
  p = g.modelGatherHyp(fpeak,snr,vnmo)
  print 'compute semblance spectrum'
  s = Semblance(st,sx,sv)
  s.setStretchMute(1.0)
  ps = s.applyHyp(p)
  print 'pick velocities with dynamic programming'
  ds = DynamicSolver(sv.count,rmin,rmax,dr)
  ac = ds.accumulate(transpose(ps))
  uv = ds.findSolution(transpose(ps))
  uv = add(sv.first,mul(uv,sv.delta))
  q = s.flattenGatherHyp(p,uv)
  cv = s.pickMaxSemblance(ps)
  pduv = mul(div(sub(uv,vnmo),vnmo),100)
  pdcv = mul(div(sub(cv,vnmo),vnmo),100)
  print 'plotting'
  #plot2(transpose(ac),cmap=jet,u=vs,s1=st,s2=sv,name='accumulated semblance',sem=True)
  plot1(vnmo,st,uv,cv,name="vnmo")
  plot2(ps,u=cv,s1=st,s2=sv,name='semblance max pick',sem=True)
  plot1(pduv,st,name='vnmo % difference dp')
  plot1(pdcv,st,name='vnmo % difference max')
  plot1(vnmo,st,uv,name="vnmo and vnmo estimated")
  plot2(ps,u=uv,s1=st,s2=sv,name='semblance',sem=True)
  plot2(ps,s1=st,s2=sv,name='semblance',sem=True)
  plot2(q,s1=st,s2=sx,name='flat cmp',cmp=True)
  plot2(p,s1=st,s2=sx,name='cmp',cmp=True)

def goNonHypTest(data):
  s1,so,g,vnmo,eta,wb,vpt,etat = data
  st = s1
  sx = so
  #goTestData()
  fpeak = 30
  #snr = 1.0e6
  snr = 5.0
  #snr = 0.1
  rmin1,rmax1 = -0.05,0.55
  rmin2,rmax2 = -0.05,0.15
  dr = 0.05
  vpl = [1.492,3.00,0.004]
  epl = [0.000,0.15,0.004]
  nv = int((vpl[1]-vpl[0])/vpl[2])
  ne = int((epl[1]-epl[0])/epl[2])
  sv = Sampling(nv,vpl[2],vpl[0])
  se = Sampling(ne,epl[2],epl[0])
  print 'vcmin=',rmin1*sv.delta,',vcmax=',rmax1*sv.delta
  print 'ecmin=',rmin2*se.delta,',ecmax=',rmax2*se.delta
  print 'sparse grid=',st.delta/dr,'seconds'
  print 'sigma=',1.0/(fpeak*st.delta)
  print 'model a gather'
  g = Gather(st,sx,sv,se)
  #vnmo = g.makeLinearParameter(2.1,3.5)
  #eta = g.makeLinearParameter(0.05,0.4)
  p = g.modelGatherNonHyp(fpeak,snr,vnmo,eta)
  rs = [g,vnmo,eta,p,fpeak,snr,rmin1,rmax1,rmin2,rmax2,dr,sv,se,st,sx]
  goNonHypTest1(rs)
  #goNonHypTest2(rs)


def goNonHypTest1(rs):
  g,vnmo,eta,p,fpeak,snr,rmin1,rmax1,rmin2,rmax2,dr,sv,se,st,sx = rs
  print 'compute semblance spectrum'
  sw = Stopwatch(); sw.restart();
  s = SemblanceOld(st,sx,sv,se,fpeak)
  s.setStretchMute(1.0)
  ps = s.applyNonHyp(p)
  sw.stop(); print 'total time =',sw.time(),'seconds'
  print 'pick velocities with dynamic programming'
  ds = DynamicSolver(sv.count,se.count,rmin1,rmax1,rmin2,rmax2,dr)
  uv,ue = ds.findSolution(ps)
  ue = add(se.first,mul(ue,se.delta))
  uv = add(sv.first,mul(uv,sv.delta))
  ut = floats(st.getValues())
  q  = s.flattenGatherNonHypE(p,uv,ue)
  qe = s.flattenGatherNonHypE(p,vnmo,eta)
  cv,ce = s.pickMaxSemblance(ps,se)
  plots(vnmo,eta,uv,ue,ut,p,q,ps,cv,ce,se,qe,st,sx,sv,se,"eta")

def goNonHypTest2(rs):
  g,vnmo,eta,p,fpeak,snr,rmin1,rmax1,rmin2,rmax2,dr = rs
  vhor = cvhor(vnmo,eta)
  print 'compute semblance spectrum'
  s = Semblance(st,sx,sv,se,fpeak)
  s.setStretchMute(1.0)
  ps = s.applyNonHyp(p,sh)
  print 'pick velocities with dynamic programming'
  ds = DynamicSolver(sv.count,sh.count,rmin1,rmax1,rmin2,rmax2,dr)
  uv,uh = ds.findSolution(ps)
  uh = add(sh.first,mul(uh,sh.delta))
  uv = add(sv.first,mul(uv,sv.delta))
  ut = floats(st.getValues())
  q  = s.flattenGatherNonHypH(p,uv,uh)
  cv,ch = s.pickMaxSemblance(ps,sh)
  ce = ceta(uv,uh)
  plot1(eta,st,ce,name="eta from vhor")
  plots(vnmo,vhor,uv,uh,ut,p,q,ps,cv,ch,sh,v2="vhor")

def plots(vnmo,eta,uv,ue,ut,p,q,ps,cv,ca,sa,qe,st,sx,sv,se,v2):
  vve = None; vvh = None;
  if v2=="eta": vve=True
  if v2=="vhor": vvh=True
  print 'plot '+v2
  pduv = mul(div(sub(uv,vnmo),vnmo),100)
  pdcv = mul(div(sub(cv,vnmo),vnmo),100)
  plot1(eta,st,ue,ca,name=v2)
  plot1(sub(eta,ue),st,name=v2+' dp difference')
  plot1(sub(eta,ca),st,name=v2+' max difference')
  plot1(eta,st,ue,name=v2,hlabel="eta")
  plot1(vnmo,st,uv,cv,name='vnmo',hlabel="velocity (km/s)")
  plot1(pduv,st,name='vnmo % difference dp')
  plot1(pdcv,st,name='vnmo % difference max')
  plot1(vnmo,st,uv,name='vnmo',hlabel="velocity (km/s)",paper='vnmoestm1')
  plot1(eta,st,ue,name='eta',hlabel="eta",paper='estaestm1')
  plot2(ps[900],s1=sv,s2=sa,name='vnmo vs '+v2+'t=900',vve=vve,vvh=vvh)
  plot2(ps[700],s1=sv,s2=sa,name='vnmo vs '+v2+'t=700',vve=vve,vvh=vvh)
  plot2(ps[500],s1=sv,s2=sa,name='vnmo vs '+v2+'t=500',vve=vve,vvh=vvh)
  plot2(ps[300],s1=sv,s2=sa,name='vnmo vs '+v2+'t=300',vve=vve,vvh=vvh)
  plot2(ps[100],s1=sv,s2=sa,name='vnmo vs '+v2+'t=100',vve=vve,vvh=vvh)
  plot2(ps[10 ],s1=sv,s2=sa,name='vnmo vs '+v2+'t=10',vve=vve,vvh=vvh)
  plot3(ps,cmap=jet,s1=sv,s2=sa,s3=st,u=[uv,ue,ut],perc=99,\
                                      name=v2+' semblance',paper='nhsu')
  plot3(ps,cmap=jet,s1=sv,s2=sa,s3=st,name=v2+' semblance')
  plot2(q,s1=st,s2=sx,name=v2+' flat cmp',cmp=True)
  plot2(qe,s1=st,s2=sx,name=v2+' flat exact cmp',cmp=True)
  plot2(p,s1=st,s2=sx,name=v2+' cmp',cmp=True)

def cvhor(vnmo,eta):
	return mul(vnmo,sqrt(add(mul(2,eta),1)))
def ceta(vnmo,vhor):
	return mul(0.5,sub(pow(div(vhor,vnmo),2),1))
  

def goNishantTest():
  omin = 0.025
  omax = 2.5
  do = 0.025
  no = inro((omax-omin)/do)
  so = Sampling(no,do,omin)
  lt = 3
  dt = 0.004
  nt = inro(lt/dt)
  st = Sampling(nt,dt,0.0)
  fname = "/Users/amunoz/Home/data/bp/dat/test/pp_test.txt" 
  g = zerofloat(nt,no)
  fil = File(fname)
  br = BufferedReader(FileReader(fil))
  for io in range(no):
    for it in range(nt):
      lns = br.readLine()
      if lns: 
        lnf = Float.parseFloat(lns)
        g[io][it] = lnf
  br.close()
  r1min,r1max = -0.50,0.50
  r2min,r2max = -0.50,0.50
  dr = 0.5
  sv = Sampling(300,0.01,1.0)
  se = Sampling(80,0.005,0)
  #s = Semblance(st,so,15)
  #ps = s.applyNonHypV(g,sv,se)
  #ds = DynamicSolver(sv.count,se.count,r1min,r1max,r2min,r2max,dr)
  #up1,up2 = ds.findSolution(ps)
  #up1 = add(sv.first,mul(up1,sv.delta))
  #up2 = add(se.first,mul(up2,se.delta))
  #q  = s.flattenGatherNonHypE(g,up1,up2)
  #plot3(ps,cmap=jet,s1=sv,s2=se,s3=st,name='semblance')
  #plot2(q,s1=st,s2=so,name='flat cmp',cmp=True)
  #plot2(q,s1=st,s2=so,name='cmp',cmp=True)
  s = SemblanceOld(st,so,sv,15)
  ps = s.applyHyp(g)
  plot2(ps,s1=st,s2=sv,name='semblance',sem=True)
  plot2(g,s1=st,s2=so,name='cmp',cmp=True)
  return
  ds = DynamicSolver(sv.count,r1min,r1max,dr)
  uv = ds.findSolution(transpose(ps))
  uv = add(sv.first,mul(uv,sv.delta))
  q = s.flattenGatherHyp(g,uv)
  plot2(ps,u=uv,s1=st,s2=sv,name='semblance',sem=True)
  plot2(ps,s1=st,s2=sv,name='semblance',sem=True)
  plot2(q,s1=st,s2=so,name='flat cmp',cmp=True)
  plot2(g,s1=st,s2=so,name='cmp',cmp=True)


###############################################################################
# Data Samplings

def goTestData():
  global st,sx,sv,se,sh
  global vp,ep,hp
  vp = [1.8,4.0,0.02]
  ep = [0.0,0.6,0.01]
  hp = csvhor(vp,ep)
  st = Sampling(1001,0.004,0.0)
  sx = Sampling(251, 0.050,0.0)
  nv = int((vp[1]-vp[0])/vp[2])
  ne = int((ep[1]-ep[0])/ep[2])
  nh = int((hp[1]-hp[0])/hp[2])
  sv = Sampling(nv,vp[2],vp[0])
  se = Sampling(ne,ep[2],ep[0])
  sh = Sampling(nh,hp[2],hp[0])
  print 'nv=',nv
  print 'ne=',ne
  print 'nh=',nh
  print 'fv,lv=',vp[0],',',vp[1]
  print 'fe,le=',ep[0],',',ep[1]
  print 'fh,lh=',hp[0],',',hp[1]

def csvhor(vp,ep):
  v1 = vp[0]*sqrt(2*ep[0]+1)
  v2 = vp[1]*sqrt(2*ep[1]+1)
  return [v1,v2,vp[2]]

###############################################################################
# Plots

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE

def plot1(x,s1=None,u=None,v=None,name=None,\
          hlabel=None,cgs=None,paper=None,colm=1,hlim=None,bounds=None):
  pp = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  pv = pp.addPoints(x)
  if s1: pv.set(s1,x)
  if u: addLine(pp,u,s1,RED)
  if v: addLine(pp,v,s1,BLUE)
  if cgs and s1:
    l = x
    if u:
      l = u
    cgsv = add(s1.first,mul(floats(cgs),s1.delta))
    cgsx = zerofloat(len(cgs))
    for ic in range(len(cgs)):
      cgsx[ic] = l[cgs[ic]] 
    pv = pp.addPoints(cgsv,cgsx)
    pv.setStyle("rO")
  if bounds:
    ul,ll = bounds
    addLine(pp,ul,s1,BLUE)
    addLine(pp,ll,s1,BLUE)
  pp.setVLabel("Time (s)")
  if hlim: pp.setHLimits(hlim[0],hlim[1])
  if hlabel: pp.setHLabel(hlabel)
  frame(pp,name,[500,500],paper,colm)

def addLine(p,u,s=None,color=BLACK):
  if s: pv = p.addPoints(s,u)
  else:	pv = p.addPoints(u)
  pv.setLineColor(color)

def plot2(x,s1=None,s2=None,u=None,u2=None,cmap=None,\
         cmin=0,cmax=0,perc=100,cbar=None,name=None,\
				 nrt=None,cmp=None,sem=None,vve=None,vvh=None,\
         vpz=None,lims=None,cmppts=None,paper=None,\
         size=[673,640],colm=1,mx=None):
  if sem: cbar = "Semblance"
  pan = panel(cbar)
  pix = pan.addPixels(x)
  if s1 and s2:
    pix.set(s1,s2,x)
  if cmp:
    cmapo = gray
    pan.setHLabel("Offset (km)")
    pan.setVLabel("Time (s)")
    size=[800,1000]
    colm=2
  if sem:
    cmapo = jet 
    pan.setHLabel("Vnmo (km/s)")
    pan.setVLabel("Time (s)")
    colm=2
  if vve:
    cmapo = jet 
    pan.setHLabel("Eta")
    pan.setVLabel("Vnmo (km/s)")
    colm=2
  if vvh:
    cmapo = jet 
    pan.setHLabel("Vhor (km/s)")
    pan.setVLabel("Vnmo (km/s)")
  if vpz:
    cmapo = jet 
    pan.setVLabel("Depth (km)")
    pan.setHLabel("Distance (km)")
    size = [569,537]
    colm=3
  if cmap: cmapo = cmap
  pix.setColorModel(cmapo)
  if nrt:
    pix.setInterpolation(PixelsView.Interpolation.NEAREST)
  if u:
    if s1:
      pt = pan.addPoints(s1,u)
    else:
      pt = pan.addPoints(u)
    pt.setLineColor(WHITE)
    pt.setLineWidth(3)
  if u2:
    if s1:
      pt2 = pan.addPoints(s1,u2)
    else:
      pt2 = pan.addPoints(u2)
    pt2.setLineColor(WHITE)
    pt2.setLineWidth(3)
    pt2.setLineStyle(PointsView.Line.DASH)
  if mx:
    mx1,mx2,mxe1,mxe2 = mx
    # est
    amx1 = zerofloat(1)
    amx2 = zerofloat(1)
    amx1[0],amx2[0] = mx1,mx2
    # max
    index = zeroint(2)
    mv = max(x,index)
    ipx1,ipx2 = index
    px1,px2 = ipx1*s1.delta+s1.first,ipx2*s2.delta+s2.first
    apx1 = zerofloat(1)
    apx2 = zerofloat(1)
    apx1[0],apx2[0] = px1,px2
    # exact 
    aex1 = zerofloat(1)
    aex2 = zerofloat(1)
    aex1[0],aex2[0] = mxe1,mxe2
    pv1 = pan.addPoints(amx1,amx2)
    pv2 = pan.addPoints(apx1,apx2)
    pv3 = pan.addPoints(aex1,aex2)
    pv1.setStyle("wo")
    pv2.setStyle("wx")
    pv3.setStyle("wO")
    pv1.setMarkSize(20)
    pv2.setMarkSize(20)
    pv3.setMarkSize(15)
  if cmppts:
    ix = cmppts
    cmpl = fillfloat(s2.getValue(ix),1)
    x1 = fillfloat(s1.delta*15,1)
    pts = pan.addPoints(x1,cmpl)
    pts.setMarkSize(10)
    pts.setStyle("rO")
  if cmin<cmax:
    pix.setClips(cmin,cmax)
  if perc<100:
    pix.setPercentiles(100-perc,perc)
  if lims:
    pan.setLimits(lims[0],lims[1],lims[2],lims[3])
  elif s1 and s2:
    pan.setLimits(s2.first,s1.first,s2.last,s1.last)
  frame(pan,name,paper=paper,size=size,colm=colm)

def plot3(x,cmap=jet,s1=None,s2=None,s3=None,u=None,\
         cmin=0,cmax=0,perc=100,cbar=None,name=None,paper=None):
  world = World()
  if s1 and s2 and s3:
    ipg = ImagePanelGroup(s1,s2,s3,x)
  else:
    ipg = ImagePanelGroup(x)
  ipg.setColorModel(cmap)
  if cmin<cmax:
    ipg.setClips(cmin,cmax)
  if perc<100:
    ipg.setPercentiles(100-perc,perc)
  world.addChild(ipg)
  if u:
    lg = makePoints(u,YELLOW,8)
    world.addChild(lg)
    if s1:
      nt = s1.count
      us = zerofloat(nt)
      for i in range(nt):
        us[i] = x[s3.indexOfNearest(u[0][i])][s2.indexOfNearest(u[1][i])][i]
      plot1(us,s1,name="semblance along path "+name)
  sf = SimpleFrame(world)
  ov = sf.getOrbitView()	
  ov.setAxesScale(0.5,2.0,0.1)
  ov.setScale(14.843407)
  ov.setAzimuth(123.91459)
  ov.setElevation(26.942446)
  ipg.setSlices(\
    round((7.008-s1.first)/s1.delta),\
    round((0.112-s2.first)/s2.delta),\
    round((1.659-s3.first)/s3.delta))
  if name:
    sf.setTitle(name)
  sf.setVisible(setvis)
  if paper and papers:
    sf.paintToFile(pngDir+paper+".png")
  if slides:
    sf.paintToFile(pngDir+name+".png")

def makePoints(u,color,size):
  x1,x2,x3 = u
  n = len(x1)
  xyz = zerofloat(3*n)
  copy(n,0,1,x1,0,3,xyz)
  copy(n,0,1,x2,1,3,xyz)
  copy(n,0,1,x3,2,3,xyz)
  rgb = None
  pg = LineGroup(xyz,rgb)
  ls = LineState()
  ls.setWidth(size)
  ls.setSmooth(False)
  ss = StateSet()
  ss.add(ls)
  #pg = PointGroup(xyz,rgb)
  #ps = PointState()
  #ps.setSize(8)
  #ps.setSmooth(False)
  #ss = StateSet()
  #ss.add(ps)
  cs=ColorState()
  cs.setColor(color)
  ss.add(cs)
  pg.setStates(ss)
  return pg

def panel(cbar=None):
  p = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  cb = p.addColorBar()
  #cb.setWidthMinimum(100)
  if cbar:
    cb.setLabel(cbar)
  return p

def frame(panel,name,size=[600,800],paper=None,colm=None):
  #panel.setVLabel('time (s)')
  frame = PlotFrame(panel)
  #frame.setBackground(Color(204,204,204,255))
  #frame.setFontSizeForSlide(1.0,1.0)
  frame.setSize(size[0],size[1])
  if name:
    frame.setTitle(name)
  frame.setVisible(setvis)
  if paper and papers:
    if paper: name=paper
    if colm==2:
      frame.setFontSizeForPrint(8.0,469.0) # 2 column
      frame.paintToPng(720,6.51,pngDir+name+'.png')
    elif colm==3:
      frame.setFontSizeForPrint(8.0,165.3) # 1/3 column
      frame.paintToPng(720,2.17,pngDir+name+'.png')
    else:
      frame.setFontSizeForPrint(8.0,222.0) # 1 column
      frame.paintToPng(720,3.08,pngDir+name+'.png')
  if slides:
    ar = 16.0/9.0; 
    fw = 0.9; fh = 0.90; 
    if colm==2:
      fw = 0.45; fh = 0.90; 
    if colm==3:
      fw = 0.25; fh = 0.90; 
    printSlide(frame,pngDir+name,fw,fh,ar)

"""
Method to print slides.
@param frame a PlotFrame class
@param fname the png file name
@param fw    the fraction of the slide width that the figure occupies
@param fh    the fraction of the slide height that the figure occupies
@param ar    the aspect ratio of the slide (use keynote default pixel widths)
@param sc    scalar for frame size that is dependent on screen resolution 
              (my MBP 15" needs sc=2)
"""
def printSlide(frame,fname,fw,fh,ar,sc=2):
  swp = 1920 # 16/9 keynote slide default 
  if ar==4.0/3.0: swp = 1024 # 4/3 keynote slide default
  fwi = int(swp*fw/sc)+1
  fhi = int(swp/ar*fh/sc)+1
  frame.setSize(fwi,fhi)
  frame.setFontSizeForSlide(fw,fh,ar) 
  frame.paintToPng(swp*fw,1.0,fname+".png")

def plotModelsCmp(g,s1,so,v,e,name,clims=None,perc=99):
  pp = PlotPanel(1,3,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  pg = pp.addPixels(0,0,s1,so,g)
  pv = pp.addPoints(0,1,s1,v)
  pe = pp.addPoints(0,2,s1,e)
  pp.setVLabel(0,"Time (s)")
  pp.setHLabel(0,"Offset (km)")
  pp.setHLabel(1,"Velocity (km/s)")
  pp.setHLabel(2,"Eta")
  pg.setPercentiles(100-perc,perc)
  if clims:
    pp.setHLimits(0,clims[0],clims[1])
  frame = PlotFrame(pp)
  frame.setSize(800,400)
  frame.setTitle(name)
  frame.setVisible(setvis)
  if papers:
    frame.setFontSizeForPrint(8.0,469.0) # 2 column
    frame.paintToPng(720,6.51,pngDir+name+'.png')
  if slides:
    frame.setFontSizeForSlide(1.0,0.9,16.0/9.0) # 2 column
    frame.paintToPng(720,2.4,pngDir+name+'.png')

def plot1comp(s1,v1,v2,dv,name,clims=None,clims2=None,vel=None,eta=None):
  pp = PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  pv1 = pp.addPoints(0,0,s1,v1)
  pv2 = pp.addPoints(0,0,s1,v2)
  pv2.setStyle("r-")
  pd = pp.addPoints(0,1,s1,dv)
  pp.setVLabel(0,"Time (s)")
  if vel:
    pp.setHLabel(0,"Velocity (km/s)")
    pp.setHLabel(1,"% Difference")
  if eta:
    pp.setHLabel(0,"Eta")
    pp.setHLabel(1,"Difference")
  if clims:
    pp.setHLimits(0,clims[0],clims[1])
  frame = PlotFrame(pp)
  if clims2:
    pp.setHLimits(1,clims2[0],clims2[1])
  #frame.setSize(646,600)
  frame.setSize(658,555)
  frame.setTitle(name)
  frame.setVisible(setvis)
  if papers:
    frame.setFontSizeForPrint(8.0,469.0) # 2 column
    frame.paintToPng(720,6.51,pngDir+name+'.png')
  if slides:
    frame.setFontSizeForSlide(0.65,0.9,16.0/9.0) # 2 column
    frame.paintToPng(720,2.4,pngDir+name+'.png')

def plot1special(ver,eer):
  pp = PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  pv1 = pp.addPoints(0,0,ver)
  pv2 = pp.addPoints(0,1,eer)
  pp.setHLabel(0,"Total eta error")
  pp.setHLabel(1,"Total vnmo error")
  pp.setVLabel(0,"CMP")
  frame = PlotFrame(pp)
  frame.setSize(658,555)
  frame.setVisible(True)
  frame.setFontSizeForSlide(0.65,0.9,16.0/9.0) # 2 column
  frame.paintToPng(720,2.4,'./rms.png')


"""
Converts 1D array to 1D float array
"""
def floats(x):
  n = len(x)
  xd = zerofloat(n)
  for i in range(n):
    xd[i] = float(x[i])
  return xd

def inro(x):
	return int(round(x))

def ince(x):
	return int(ceil(x))

def infl(x):
	return int(floor(x))


#------------------------------------------------------------------------------#

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
