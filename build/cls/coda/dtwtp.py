# Dynamic Time Warping plotter for teapot dome data
# This runs all Dtw methods and plots
# @author Andrew Munoz, Colorado School of Mines
# @version 10.20.2011

from imports import *
from tputils import *
from tp import *
from dtw import Dtw1 

_pngDir = "./tppics/"

def main(args):
  goDtwData()
  #printWellIDList()

  """
  For testing
  """
  #goDtwTest()
  #goDtwTestData()
  #getAILog()
  #getSeismicTraces(35)
  #getReflectivity(35)
  #waveletTest()
  
def goDtwData():
 #Params
  ail,zl,yl,xl,vl = getAILog()
  datum = 0.4041648;
  sz = Sampling(len(zl),0.0001524,zl[0]+datum)
  fp = 35 #(Hz) peak frequency
  trace,st,sy,sx,tsdataxc,yc,xc,yv,xv = getSeismicTraces(fp)
  tr = 0
  binterval = 0 # in samples
  bst = 2 #int(len(trace[tr])*0.1) # sample at which paths can start

  # Warp well to trace directly
  #warpData1(sx,sy,st,trace[tr],tsdataxc,yc,fp,binterval,bst)

  # Stretch-Convolve-Unstretch well then warp
  #warpData2(sx,sy,sz,st,trace[tr],tsdataxc,yc,fp,binterval,bst,ail,vl)

  # Stretch-Convolve well then warp
  warpData3(sx,sy,sz,st,trace[tr],tsdataxc,yc,xc,yv,xv,fp,ail,vl)

"""
Smooth the reflectivity in depth and then dtw with the trace and warp 
"""
def warpData1(sx,sy,st,trace,tsdataxc,yc,fp,binterval,bst):
  well,sz,wello,tli,tui,v1,v2 = getReflectivity(fp)
  #wellr = applyPhaseRot(iph,well)
  cm  = Dtw1().cm(trace,well)
  acm = Dtw1().acm(cm) 
  pe,tpi,zpi = Dtw1().paths(binterval,bst,cm,tui,tli)
  ri,mrms = minimums(rms2(pe))
  ip=ri
  
  tp2,zp2 = reverseInputThrow(tpi[ip],zpi[ip],st,sz)
  #tp2,zp2 = smoothTDcurve(tp2,zp2,2)
  toz,zot = interpolations(tp2,zp2,st,sz)
  toz2,zot2 = smoothTDcurve(toz,zot,6)
  vz,vt = velocities(toz2,zot2,sz,st)
  wellt = stretch(zot,well,st)
  sts = Sampling(len(wellt),st.getDelta(),zot[0])
 
  #Plots 
  #plotImgCurves(st,sz,cm,trace,well)
  #plotImgCurves(st,sz,acm,trace,well)
  plotImgCurves(st,sz,acm,trace,well,cn=len(tpi),pt=tpi,pd=zpi)#paper="dtwddc")
  #plotImgCurves(st,sz,acm,trace,well,cn=len(tpi),pt=tpi,pd=zpi,v1=v1,v2=v2)
  #plotImgCurves(st,sz,acm,trace,well,cn=len(tpi),pt=tpi,pd=zpi,tlim1=tli,tlim2=tui)
  #plotImgCurves(st,sz,acm,trace,well,c1=ri,pt=tpi,pd=zpi)
  #plotCurvePanel2(trace,st,well,zot,vt,"t",Hlims=3)
  #plot2DLog(tsdataxc,wellt,yc,st,sy,sts,paper="teaser")  

"""
Convert the well to time, convolve with a wavelet then apply dtw in time
"""
def warpData3(sx,smy,sz,st,tro,tsdataxc,yc,xc,yv,xv,fp,ail,vl):
  dz = sz.getDelta(); dt = st.getDelta();
  sigma = 1/(2*FLT_PI*fp*dt);
  sigf = sigma*8
  r = 5 # number of lowest minumum paths to keep
  # Make synthetic
  rf = reflcalc(ail)
  #rf = Dtw1().makeEvents(len(ail), 30859) 
  nr = len(rf)
  tauz = TDutils().maketauz(vl,sz)
  syo = TDutils().applyTWavelet(rf,st,sz,fp,tauz)
  ss = Sampling(len(syo),dt,tauz[0])

  tui,tli = 0,int(st.getLast()/dt) #int((2*sz.getFirst()/2.5)/dt)
  sy = normalizeRMSlocal(syo,sigf)
  tr = normalizeRMSlocal(tro,sigf)


  phase=False
  ephase=False
  # DTW
  if phase:
    nh = 361
    mphr = zerofloat(nh); 
    pmrms = 1; ihm=0;
    #pset = [0,84,200];bn = zerofloat(len(tr),3);ct=0;
    for ih in range(nh):
      syr = applyPhaseRot(ih,sy)
      cm  = Dtw1().cm(tr,syr)
      acm  = Dtw1().acm(cm)
      pe,tpi,zpi = Dtw1().paths(cm)
      pe1,tpi1,zpi1 = Dtw1().pathRemoval(pe,tpi,zpi,tui,tli)
      pe3,tpi3,zpi3 = Dtw1().pathRemoval(pe1,tpi1,zpi1,acm[len(acm)-1],1)
      mphr[ih] = pe3[0][0]/len(pe3[0])
      # For B plot for multiple phases
      #plotImgCurves2(st,ss,acm,tr,sy,cn=len(tpi3),pt=tpi3,pd=zpi3)
      #for i in range(len(tr)):
      #  bn[ct][i] = acm[len(acm)-1][i]/len(pe[i])
      #ct+=1
      if (mphr[ih]<=pmrms): 
        ihm=ih
        pmrms = mphr[ih]
      print ih
    print 'Optimum Phase is',ihm
    #plotDphase(bn,st,paper='pdnphases')#,slides='Dnphases')
    plotPhaseRMS(mphr,nh,1,slides='sphaseplot',paper='phaseplot')
  elif ephase:
    nh = 360
    mset = [30,5,1]
    mphr = zerofloat(nh); 
    pmrms = 1
    sa = 0; en = nh; ihm=0;
    for mu in mset:
      for ih in range(sa,en,mu):
        print ih
        syr = applyPhaseRot(ih,sy)
        cm  = Dtw1().cm(tr,syr)
        acm  = Dtw1().acm(cm)
        pe,tpi,zpi = Dtw1().paths(cm)
        pe1,tpi1,zpi1 = Dtw1().pathRemoval(pe,tpi,zpi,tui,tli)
        #pe2,tpi2,zpi2 = Dtw1().pathRemoval(pe1,tpi1,zpi1,acm[len(acm)-1],r) # best r paths
        pe3,tpi3,zpi3 = Dtw1().pathRemoval(pe1,tpi1,zpi1,acm[len(acm)-1],1)
        mphr[ih] = pe3[0][0]/len(pe3[0])
        #mphr[ih] = pe2[5][0]/len(pe2[5])
        # For B plot for multiple phases
        #plotImgCurves2(st,ss,acm,tr,sy,cn=len(tpi3),pt=tpi3,pd=zpi3)
        #for i in range(len(tr)):
        #  bn[ct][i] = acm[len(acm)-1][i]/len(pe[i])
        #ct+=1
        if (mphr[ih]<=pmrms): 
          ihm=ih
          pmrms = mphr[ih]
      sa = ihm-mu
      en = sa+2*mu  
      if sa<0: sa=0
      if en>nh: en=nh;
    print 'Optimum Phase is',ihm
  else: ihm= 92

  # With minimum D phase
  syr = applyPhaseRot(ihm,sy)
  cm  = Dtw1().cm(tr,syr)
  acm = Dtw1().acm(cm) 
  pe,tpi,zpi = Dtw1().paths(cm) # all paths from every point
  pe1,tpi1,zpi1 = Dtw1().pathRemoval(pe,tpi,zpi,tui,tli) # paths ending w/in trace
  pe2,tpi2,zpi2 = Dtw1().pathRemoval(pe1,tpi1,zpi1,acm[len(acm)-1],r) # best r paths
  pe3,tpi3,zpi3 = Dtw1().pathRemoval(pe1,tpi1,zpi1,acm[len(acm)-1],1) # best path
   
  ttau,taut,tz = TDutils().tdcurve(sz,st,ss,tpi3[0],zpi3[0],tauz)
  #ttau,taut,tz = TDutils().tdcurve(sz,st,ss,tpi2[1],zpi2[1],tauz)
  syrt = TDutils().stretch(ttau,syr,ss)
  sst = Sampling(len(syrt),ss.getDelta(),ttau[0])
  # Bulk shifted tauz
  tauzbs = add(tauz,tz[0]-tauz[0])
  sstbs = Sampling(len(sy),ss.getDelta(),tauzbs[0])
 
  # Selected horizons and well tops for well
  ght,ghz,ghy = getHorizons(xc,sx,sz,st,tz)
  print 'Horizon Times'
  printds(xc,yc,sx,smy,"t")
  print '\n'
  print 'Horizon Depths'
  printds(xc,yc,sx,smy,"z")
  print '\n'
  print 'Well Top Depths'
  gw = getWellTops()
  wtops = zerofloat(sz.getCount())
  for i in range(len(gw)):
    wtops[sz.indexOfNearest(gw[i])] = 1

#  nz=sz.getCount(); dz=sz.getDelta(); fz=sz.getFirst();
#  snrgf = zerofloat(nz)
#  snrs = zerofloat(nz)
#  RecursiveGaussianFilter(128).apply1(tz,snrgf)
#  RgSmoother(128).apply1(tz,snrs)
#  
#  for i in range(10):
#    print i
#    print tz[i]
#  for i in range(10):
#    print i
#    print tz[nz-1-i]
#  vnrgf = div(2*dz,snrgf)
#  vnrs = div(2*dz,snrs)
#  plotCurve(tz,"tz",s=sz)
#  plotCurve(tauz,"tauz",s=sz)
#  #plotCurve(tz2,"smoothed tz",s=sz,d='d')
#  plotCurve(vl,"velocity log",s=sz,d='d')
#  #plotCurve(sn,"new slowness",s=sz,d='d')
#  plotCurve(vnrgf,"new velocity rgf",s=sz,d='d')
#  plotCurve(vnrs,"new velocity rs",s=sz,d='d')
#  plotCurve2(tauz,tz,sz,sz,"tz vs. tauz",d='d')
#  plotCurve2(vl,vnrgf,sz,sz,"new velocity rgf vs. log",d='d')
#  plotCurve2(vl,vnrs,sz,sz,"new velocity rs vs. log",d='d')
 
  #plotCurve(sub(abs(syrt2),abs(syrt)),"sysub")
  #plotCurve(sub(tm2,tm),"tmsub")
  #plotCurve(syrt2,"sy2",s=sst2)
  #plotCurve(syrt,"sy",s=sst)
  #plotCurve(rf,"",s=sz,d="d")#,slides='srf')
  #plotCurve(syo,"",s=ss,intv=0.2,lim=0.25,slides='ssyo')#,paper='psyo2')#
  #plotCurve(tro,"",s=st,ss=ss,slides='stro') #paper='ptro2')#
  #plotCurve(syrt,"",s=sst,intv=2.0,slides='ssyrt')
  #plotCurve(sy,"",s=ss,intv=2.0,slides='ssy')#paper='psy')#,fill="p")
  #plotCurve(tr,"",s=st,ss=ss,slides='strn')#paper='ptr')#,fill="p")
  #plotCurveN(st,tro,tr,ss,syo,sy,paper='pnorms')
  #plotCurveW(sz,rf,ss,syo,paper='pwavsy')
  #plotTDCurves(sz,tauz,tz,slides='stdcurve')#paper='tdplot23')

  #plot2D(st,ss,cm,cbwmin=85,slides='sdtwde')#paper='pdtwde')
  #plot2D(st,ss,acm,cbintv=300,slides='sdtwdd')#paper='pdtwdd')
  #plot2D(st,ss,acm,cn=len(pe),pt=tpi,pd=zpi,cbintv=300,paper='pdtwdca')
  #plot2D(st,ss,acm,cn=len(pe1),pt=tpi1,pd=zpi1,cbintv=300,paper='pdtwdcs')
  #plot2D(st,ss,acm,cn=len(pe2),pt=tpi2,pd=zpi2,cbintv=300,paper='pdtwdcr')
  #plot2D(st,ss,acm,cn=len(pe3),pt=tpi2,pd=zpi2,cbintv=300,slides='sdtwdcm')#,paper='pdtwdcm')
  #plotCurvePanelT(tr,st,syr,syrt,sst,ss,paper3='pteaser43')#,slides='sdtwdcpo')
  #plotCurvePanel3(tr,st,syr,syrt,sst,ss,ttau,slides='sdtwdcpobs')#paper='pwarps1')#

  plotImgCurves2(st,ss, cm,tr,sy)#,slides='sdtwde')#,paper='dtwdte')
  plotImgCurves2(st,ss,acm,tr,sy)#,slides='sdtwdd')#,paper='dtwdtd')
  #plotImgCurves2(st,ss,acm,tr,sy,cn=len(tpi),pt=tpi,pd=zpi)#,slides='sdtwddca')#,paper='dtwdtdc')
  #plotImgCurves2(st,ss,acm,tr,sy,cn=len(tpi1),pt=tpi1,pd=zpi1)#,slides='sdtwddct')#,paper='dtwdtdc')
  #plotImgCurves2(st,ss,acm,tr,sy,cn=len(tpi2),pt=tpi2,pd=zpi2)#,slides='sdtwddcr')#,paper='dtwdtdc')
  plotImgCurves2(st,ss,acm,tr,sy,cn=len(tpi3),pt=tpi3,pd=zpi3)#,slides='sdtwddcm')#,paper='dtwdtedm')

  plotImgB(st,ss,acm,pe=pe)#,slides='sdtwddb')
  #plotImgB(st,ss,acm,cn=len(tpi),pt=tpi,pd=zpi,slides='sdtwddba')
  #plotImgB(st,ss,acm,pe=pe,cn=len(tpi1),pt=tpi1,pd=zpi1,pt1=tpi1)#,slides='sdtwddbt')
  #plotImgB(st,ss,acm,pe=pe,cn=len(tpi2),pt=tpi2,pd=zpi2,pt1=tpi1)#,slides='sdtwddbr')
  plotImgB(st,ss,acm,pe=pe,cn=len(tpi3),pt=tpi3,pd=zpi3)#,slides='sdtwdcm2')
  
  #plotb(st,ss,acm,pe,paper='pdtwdbn')
  #plot2DLog(tsdataxc,syrt,yc,st,smy,sst,nolog="n",slides='sseismicwithoutlog')#,paper="teaser")  
  plot2DLog(tsdataxc,syrt,yc,st,smy,sst,tops=wtops,tz=tz,ghz=ght,ghy=ghy)#,slides='sseismicwithlogth')#,paper="teaser") 
  plot2DLog(tsdataxc,syrt,yc,st,smy,sst,tops=wtops,tz=tz,ghz=ghz,ghy=ghy)#,slides='sseismicwithlogzh')#,paper="teaser")  
  #plot2DLog(tsdataxc,sy,yc,st,smy,sstbs,tops=wtops,tz=tauzbs,ghz=ghz,ghy=ghy,slides='sseismicwithlogzhbs')#,paper="teaser")  

  #plotSubsq(st,ss,tr,syr,st3=sst,g2=syrt,tlim=2,slides='sdtwdcs')
  #plotSubsq(st,ss,tr,syr,tlim=2,slides='sdtwdc')




def rmsdiff(st,sst,sy,tr):
  fi = int(sst.getFirst()/sst.getDelta())
  n = sst.getCount()
  ms = 0
  for i in range(n):
    m = (tr[i+fi]-sy[i])
    ms += m*m
  rms = sqrt(ms/n)
  return rms


############################################################################################################
############################################################################################################


def getAILog():
  _set = "d"
  _ID = 490252305400
  #_ID = 490251094600
  #_ID = 490251106400
  wdata = readLogDataset(_set)
  wlog = wdata.get(_ID)
  wlog = setCoords(_ID,wlog)
  s1,s2,s3 = getTpstSamplings()
  fv,zv,yv,xv = wlog.getSamples("v",None,s2,s3) # (f,x1,x2,x3) in CSM coords
  fd,zd,yd,xd = wlog.getSamples("d",None,s2,s3)
  #fg,zg,yg,xg = wlog.getSamples("g",None,s2,s3)
  #fv = mul(fv,1/3.2808399) # ft to km
  #fd = gardnersRelation(fv)
  if (len(fv)>len(fd)): 
    fd,zd,yd,xd = addPad(fv,fd,zd,yd,xd)
  if (len(fv)<len(fd)): 
    fv,zv,yv,xv = addPad(fd,fv,zv,yv,xv)
  #if (len(fv)!=len(fg)):
  #  fg = fixGR(fv,fg)
  ailog = mul(fv,fd)
  #ailog,zv,yv,xv = endCut(ailog),endCut(zv),endCut(yv),endCut(xv)
  #zv,yv,xv = km2m(zv,yv,xv)
  #plotCurve(ailog,"Acoustic Impedance (kg/(s*m^2))",s=zv,d=1,png="AILOG")
  #plotCurve(ailog,"AI in samples")
  #plotCurve(fv,"Velocity (km/s)",s=zv,d=1,png="VLOG")
  #plotCurve(fd,"Density (g/cc)",s=zv,d=1,png="DLOG")
  #rf = reflcalc(ailog) # make reflectivity log
  #plotLogPanel(zv,fv,fd,rf,paper='plogpanel')#slides="saipanel")
  return ailog,zv,yv,xv,fv

"""
Extract traces around inline/crossline of log
Traces (looking down z-axis):

   2 3 4
   1 0 5
   8 7 6

(Well at Trace 0)
"""
def getSeismicTraces(fp):
  tsdata,st,sy,sx = readTpstData()
  ai,z,y,x,fv = getAILog()
  yc = int(y[0]/sy.getDelta())
  xc = int(x[0]/sx.getDelta())
  dev = maxdev(x,y)
  tdev = int(dev*sx.getDelta())
  print '\nmaxdev is '+str(tdev)+' traces'
  trace = zerofloat(st.getCount(),9)
  trace[0] = tsdata[xc  ][yc  ]
  trace[1] = tsdata[xc-1][yc  ]
  trace[2] = tsdata[xc-1][yc+1]
  trace[3] = tsdata[xc  ][yc+1]
  trace[4] = tsdata[xc+1][yc+1]
  trace[5] = tsdata[xc+1][yc  ]
  trace[6] = tsdata[xc+1][yc-1]
  trace[7] = tsdata[xc  ][yc-1]
  trace[8] = tsdata[xc-1][yc-1]
  sigft = 4/(FLT_PI*fp*st.getDelta())
  for i in range(1):
    trace[i] = normalizeRMSlocal(trace[i],sigft)
  #plotCurve2(trace[0],otrace[0],ts,ots,"Orig (black) Oversampled (red)",lim=2) # Tests oversampling
  #plotCurve(trace[0],"Trace 0 RMS local normalized",s=Sampling(ts))
  #plotSlice(sy,st,tsdata[xc],slides="tpslice")#png="Seismic Slice at well")
  return trace,st,sy,sx,tsdata[xc],yc,xc,y[0],x[0]

def getHorizons(xc,sx,sz,st,tz):
  d = "z"
  DKOT = readHorizonMod("FallRiverDKOT",d)
  CRMT = readHorizonMod("CrowMountainCRMT",d)
  TSAS = readHorizonMod("TensleepASand",d)
  TSBD = readHorizonMod("TensleepBbaseC1Dolo",d)
  zdk,ydk = TDutils().horizon2d(DKOT.x1,DKOT.x2,DKOT.x3,sx,sz,xc,tz,d)
  zcr,ycr = TDutils().horizon2d(CRMT.x1,CRMT.x2,CRMT.x3,sx,sz,xc,tz,d)
  zta,yta = TDutils().horizon2d(TSAS.x1,TSAS.x2,TSAS.x3,sx,sz,xc,tz,d)
  ztb,ytb = TDutils().horizon2d(TSBD.x1,TSBD.x2,TSBD.x3,sx,sz,xc,tz,d)
  d = "t"
  DKOT = readHorizonMod("FallRiverDKOT",d)
  CRMT = readHorizonMod("CrowMountainCRMT",d)
  TSAS = readHorizonMod("TensleepASand",d)
  TSBD = readHorizonMod("TensleepBbaseC1Dolo",d)
  tdk,ydk = TDutils().horizon2d(DKOT.x1,DKOT.x2,DKOT.x3,sx,st,xc,tz,d)
  tcr,ycr = TDutils().horizon2d(CRMT.x1,CRMT.x2,CRMT.x3,sx,st,xc,tz,d)
  tta,yta = TDutils().horizon2d(TSAS.x1,TSAS.x2,TSAS.x3,sx,st,xc,tz,d)
  ttb,ytb = TDutils().horizon2d(TSBD.x1,TSBD.x2,TSBD.x3,sx,st,xc,tz,d)
  return [zdk,zcr,zta,ztb],[tdk,tcr,tta,ttb],[ydk,ycr,yta,ytb]

def printds(xc,yc,sx,sy,d,convd=None):
  DKOT = readHorizonMod("FallRiverDKOT",d)
  CRMT = readHorizonMod("CrowMountainCRMT",d)
  TSAS = readHorizonMod("TensleepASand",d)
  TSBD = readHorizonMod("TensleepBbaseC1Dolo",d)
  TDutils().printds(DKOT.x1,DKOT.x2,DKOT.x3,sx,sy,xc,yc,"DKOT")
  TDutils().printds(CRMT.x1,CRMT.x2,CRMT.x3,sx,sy,xc,yc,"CRMT")
  TDutils().printds(TSAS.x1,TSAS.x2,TSAS.x3,sx,sy,xc,yc,"TSAS")
  TDutils().printds(TSBD.x1,TSBD.x2,TSBD.x3,sx,sy,xc,yc,"TSBD")
  if convd:
    TDutils().printds(DKOT.x1,DKOT.x2,convd[0],sx,sy,xc,yc,"DKOT")
    TDutils().printds(CRMT.x1,CRMT.x2,convd[1],sx,sy,xc,yc,"CRMT")
    TDutils().printds(TSAS.x1,TSAS.x2,convd[2],sx,sy,xc,yc,"TSAS")
    TDutils().printds(TSBD.x1,TSBD.x2,convd[3],sx,sy,xc,yc,"TSBD")


def getWellTops():
  # For 490252305400 
    # to match to seismic horizons, add 102m in depth,
    #  which is the difference between the seismic datum and 
    #  the KB for the well.
  zdkot = 1.173629352; 
  zcrmt = 1.3596366;   
  ztsas = 1.668877536; 
  ztsbd = 1.701323496; 
  print 'wDKOT=',zdkot
  print 'wCRMT=',zcrmt
  print 'wTSAS=',ztsas
  print 'wTSBD=',ztsbd
  return [zdkot,zcrmt,ztsas,ztsbd]

def getReflectivity(fp):
  ailog,z,y,x,fv = getAILog()
  n = len(ailog)
  sigma = getSigma(fp)
  sigf = getSigmaF(fp)
  print 'sigma is '+str(sigma)+"\n"
  synt = zerofloat(n) 
  refl = reflcalc(ailog) # make reflectivity log
  RecursiveGaussianFilter(sigma).apply2(refl,synt)
  synto = mul(-1,synt)
  synt = normalizeRMSlocal(synto,sigf)
  dz = 0.0001524 # 6 in log sampling
  dt = 0.001
  sz = Sampling(len(z),dz,z[0])
  m = 11
  dz = m*sz.getDelta()
  v1 = 4*dz/dt
  v2 = 1.5#dz/dt 
  tl = 2*sz.getFirst()/v2 #lower v bound set by water.
  tu = 2*sz.getFirst()/v1
  tli = int(tl/dt)
  tui = int(tu/dt)
  print 'New dz is '+str(dz)+' km'
  print 'Velocity Bounds are '+str(v2)+' to '+str(v1)+' km/s'
  print 'Deepest time is '+str(tl)+' s'+' and it\'s index is '+str(tli)
  print 'Shallowest time is '+str(tu)+' s'+' and it\'s index is '+str(tui)
  synt = decimate(synt,m)
  synto = decimate(synto,m)
  sz = sz.decimate(m)
  #gaussianPlot(sigma)
  #plotCurve2(y2,yd,ysu,ysd,"Reflectivity (black) vs. Decimated (red)",lim=2 # Tests oversampling)
  #plotCurve(refl,"Reflectivity",s=z,d=1,png="RLOG")
  #plotCurve(synt,"Synthetic",s=sz,d=1,png="Synthetic")
  #plotCurve2(y2,"reflectivity RMS local normalized (AGC)")#,s=z,d=1)
  #plotCurve(yd,"reflectivity decimated",s=ys,d=1)
  return synt,sz,synto,tli,tui,v1,v2

def reflcalc(ai):
  n = len(ai)
  refl = zerofloat(n)
  for i in range(1,n):
    if ((ai[i]+ai[i-1])==0): refl[i]=0
    else: refl[i] = (ai[i]-ai[i-1])/(ai[i]+ai[i-1])
  refl[0]=refl[1]
  refl = Clip(refl)
  return refl


"""
Get sigma by using seismic input  
velocity scalar (to scale v to depth)
peak frequency of seismic
log sampling interval
"""
def getSigma(fp):
  vs = 2000 #(m/s) velocity scalar
  dl = 0.1524 #(m) = 6 in log sampling (original sampling)
  sig = vs/(2*FLT_PI*fp*dl)
  return sig
"""
Sigma for local RMS normalization using
a two-sided exponential filter. Sigma is 
at least 2 wavelengths
"""
def getSigmaF(fp):
  sig = getSigma(fp)
  sigf = sig*8 
  return sigf


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
  ftoz = TDutils().maketauz(vl,sz)
  sy = TDutils().applyTWavelet(rf,st,sz,fp,ftoz)
  sy = applyPhaseRot(92,sy)
  ss = Sampling(len(sy),st.getDelta(),ftoz[0])
  pp = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  pp.addPoints(ss,sy)
  pp.setVLabel("Time (s)")
  pp.setVLimits(0,0.1);
  frame = PlotFrame(pp)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setSize(350,600)
  frame.setFontSizeForSlide(1.0,1.0)
  frame.paintToPng(720,4.00,_pngDir+'swavelet.png')
  #frame.setFontSizeForPrint(8.0,240.0)
  #frame.paintToPng(720,7.00,_pngDir+'pwavelet.png')
  frame.setVisible(True)


  #si2 = SimplePlot()
  #si2.addPoints(rf)

###########################################################################
# Utils

def rms2(pe):
  np = len(pe)
  rms = zerofloat(np)
  for ip in range(np):
    rms[ip] = pe[ip][0]/len(pe[ip])
  return rms

def rmstest(tr,t):
  dif = sub(tr,t)
  sq = mul(dif,dif)
  sums = sum(sq)
  rms = sums/len(t)
  return rms

def applyPhaseRot(ph,xr):
  n = len(xr)
  xim = zerofloat(n)
  y = zerofloat(n)
  hb = HilbertTransformFilter()
  hb.apply(n,xr,xim)
  phr = (FLT_PI/180)*ph 
  y = add(mul(cos(phr),xr),mul(sin(phr),xim))
  return y

def fixGR(fv,fg):
  fg2 = zerofloat(len(fv))
  j=0
  for i in range(len(fg)-4):
    if j==len(fv):
      break
    c1 = fg[i+1]
    c2 = fg[i+2]
    c3 = fg[i+3]
    if (fg[i]!=c1 and fg[i]!=c2 and fg[i]!=c3): 
      fg2[j] = fg[i]
      j+=1
  return fg2

"""
Smooths the time-depth curves using an
two-sided exponential smoother with a specified
sigma.
@param t times
@param d depths
@param s sigma
"""
def smoothTDcurve(t,z,s):
  nt,nz = len(t),len(z)
  t2,z2 = zerofloat(nt),zerofloat(nz)
  ExponentialSmoother(s).apply(t,t2)
  ExponentialSmoother(s).apply(z,z2)
  return t2,z2

"""
Finds the T2D velocity and D2T velocity
"""
def velocities(zot,toz,sz,st):
  return TDutils().velocities(zot,toz,sz,st)

"""
Finds the RMS error of each path
"""
def RMS(tpe):
  return TDutils().RMS(tpe)

"""
Check java documentation
"""
def reverseInputThrow(tpi,zpi,st,sz):
  return TDutils().reverseInputThrow(tpi,zpi,st,sz)

"""
Check java documentation
"""
def interpolations(tp2,zp2,st,sz):
  return TDutils().interpolation(tp2,zp2,st,sz)

"""
Check java documentation
"""
def stretch(gof,f,sg):
  return TDutils().stretch(gof,f,sg)

def makeSynthetic():
  wavelet = waveletExtractor()
  refl,rs,rsu = getReflectivity()
  synthetic = zerofloat(len(refl))
  Conv.conv(len(refl),0,refl,len(wavelet),0,wavelet,len(refl),0,synthetic)
  #plotCurve(synthetic,"synthetic")
  #plotCurve(synthetic,"synthetic vs depth",s=z)
  return synthetic

def waveletExtractor():
  trace = getSeismicTraces()
  wmin = 0
  wmax = len(trace[0])
  wvl = 100  # length of wavelet
  wavelet = zerofloat(wvl)
  wavelet = WaveletExtractor.getWavelet(wmin,wmax,wvl,trace)
  plotCurve(wavelet,"wavelet")
  return wavelet

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

"""
Prints wells that contain a combination of specified logs
"""
def printWellIDList():
  #(set,type1,type2,...)
  printLogIDCombo("d","v","d","g")

"""
Test for strictly increasing or decreasing
"""
def monotonicityTest(x):
  n = len(x)
  for i in range(n-1):
    #print x[i]
    if (x[i]>=x[i+1]): print 'not montonic between '+str(i)+' and '+str(i+1)+' val of '+str(x[i])+' and '+str(x[i+1])

"""
Test for v=0,infinty in path slopes
"""
def slopeTest(pt,pd):
  n,nb = len(pt[0]),len(pt)
  print '\nSLOPE TEST:'
  for b in range(nb):
    for i in range(1,n):
      ti = pt[b][i]
      tm1 = pt[b][i-1]
      di = pd[b][i]
      dm1 = pd[b][i-1]
      if (ti==tm1 and ti!=0.0 and ti!=1.0): 
        print 'Time  match- i='+str(i)+' b='+str(b)+' ti='+str(ti)+' tm1='+str(tm1)
      if (di==dm1 and di!=0.0 and di!=1.0): 
        print 'Depth match- i='+str(i)+' b='+str(b)+' di='+str(di)+' dm1='+str(dm1)
  print '\nEND SLOPE TEST'


"""
Find the minimum mean and rms and their indicies 
for an array of values and returns the index
"""
def minimums(rms):
  mrms,mr = findMin(rms)
  return mr,mrms

def findMin(x):
  n = len(x)
  m = x[0]
  d = 0
  for i in range(1,n):
    if (x[i]<m):
      m = x[i]
      d = i
  return m,d

"""
Decimates the specified curve x 
by a specified factor of m. 
"""
def decimate(x,m):
  n = len(x)
  #m = int(n/(4*fp)) # decimation factor
  n2 = 1+int((n-1)/m)
  y = zerofloat(n2)
  j=0;
  for i in range(len(y)):
    y[i] = x[j]
    j+=m
  return y

"""
Creates a double array with uniform sampling
"""
def uniformSampling(n,sr,fs=None):
  ts = zerodouble(n)
  if fs: sc = fs
  else: sc = 0.0
  for i in range(n):
    ts[i] = sc + i*sr
  return ts

"""
Gardners relation to convert a velocity log to a density log
"""
def gardnersRelation(v):
  n = len(v)
  d = zerofloat(n)
  for i in range(n):
    v[i] *= 3.2808399  # m to ft
    d[i] = 0.23*pow(v[i],0.25) # constants proposed by Gardner
  return d

"""
Local RMS normalization.
Square the input array, and then use a local exponential smoother
which weights the values based on the surrounding local values.
Then divide the original values by the sqrt of the smoothed values 
at each point to get the local normalization. 
"""
def normalizeRMSlocal(x,sig):
  n = len(x)
  y,xx,yy = zerofloat(n),zerofloat(n),zerofloat(n)
  xx = mul(x,x)
  ExponentialSmoother(sig,16).apply(xx,yy)
  for i in range(n):
    if (yy[i]==0): y[i]=0
    else: y[i] = x[i]/sqrt(yy[i])
  return y

"""
Global RMS normalization
"""
def normalizeRMS(x):
  n = len(x)
  rms = sqrt(sum(mul(x,x))/n)
  for i in range(n):
    x[i] /= rms
  return x

"""
Global normalization
"""
def globNorm(a):
  amax = max(a)
  for i in range(len(a)):
    a[i] /= amax
  return a

"""
Clips a signal based on the average to remove the outliers
"""
def Clip(a):
  avg = sum(abs(a))/len(a)
  for i in range(len(a)):
    if (abs(a[i])>10*avg): 
      a[i] = avg
  return a

"""
Finds the end of the nonzero values in an array and
cuts the zeroes off the end of the array 
"""
def endCut(a):
  n = len(a)
  l = 0
  t = 0
  ma = min(a)
  for i in range(n):
    if (a[i]<=1.001*ma):
      for j in range(i,n-i):
        if (a[j]<=1.001*ma):
	  t+=1
      if (t==n-i-1):
        l=t
      else: t=0
  if (l>0):
    b = zerofloat(l)
    for i in range(l):
      b[i] = a[i]
    return b
  else: return a

"""
Pad 2nd array if
len(f1)>len(f2) 
"""
def addPad(f1,f2,z2,y2,x2):
  n1,n2 = len(f1), len(f2)
  ft,zt,yt,xt = zerofloat(n1),zerofloat(n1),zerofloat(n1),zerofloat(n1)
  for i in range(n2):
    ft[i] = f2[i]
    xt[i] = x2[i]
    yt[i] = y2[i]
    zt[i] = z2[i]
  return ft,xt,yt,zt

def addPad2(x1,x2):
  n1,n2,n=len(x1),len(x2),0
  if (n1>n2):
    x3,x4 = zerofloat(n1), zerofloat(n1)
    n = n1
  if (n2>n1):
    x3,x4 = zerofloat(n2), zerofloat(n2)
    n = n2
  if (n2==n1):
    return
  for i in range(n):
    x3[i]=x1[i]
    x4[i]=x2[i]
  return x3,x4

"""
Used to over sample a signal
@param x1 the signal
@param xs1 the uniform Sampling
@param xs2 the new uniform Sampling
"""
def overSample(x1,x2,xs1,xs2):
  n1 = len(x1)
  n2 = xs2.getCount() 
  si = SincInterpolator()
  si.setUniform(n1,xs1.getDelta(),xs1.getFirst(),x1)
  si.interpolate(n2,xs2.getDelta(),xs2.getFirst(),x2)
  return x2

"""
Returns the maximum x,y deviation distance from the surface
"""
def maxdev(x,y):
  xmax,ymax = 0,0
  for i in range(1,len(x)):
    if (x[i]-x[0] > x[i-1]-x[0]):
      xmax = x[i]-x[0]
  for i in range(1,len(y)):
    if (y[i]-y[0] > y[i-1]-y[0]):
      ymax = y[i]-y[0]
  return sqrt(xmax*xmax+ymax*ymax)

"""
Converts x,y,z coords from kilometers to meters
"""
def km2m(x,y,z):
#FIX, USE mul( )
  for i in range(len(z)):
    x[i] = x[i]*1000.0
    y[i] = y[i]*1000.0
    z[i] = z[i]*1000.0
  return x,y,z

"""
Converts 1D float array to 1D double array
"""
def double(x):
  xd = zerodouble(len(x))
  for i in range(len(x)):
    xd[i] = x[i]
  return xd

"""
Converts 1D double array to 1D float array
"""
def afloat(x):
  xf = zerofloat(len(x))
  for i in range(len(x)):
    xf[i] = x[i]
  return xf


"""
To view the sampling rates of logs in a histogram
"""
def getSz(x,y,z):
  n = len(z)
  sx = zerofloat(n)
  sy = zerofloat(n)
  sz = zerofloat(n)
  sa = zerofloat(n)
  for i in range(1,n):
    sa[i] = sqrt(pow(x[i]-x[i-1],2)+pow(y[i]-y[i-1],2)+pow(z[i]-z[i-1],2))
    #sz[i] = sqrt(pow(x[i],2)+pow(y[i],2)+pow(z[i],2))
    sx[i] = sqrt(pow(x[i]-x[i-1],2))
    sy[i] = sqrt(pow(y[i]-y[i-1],2))
    sz[i] = sqrt(pow(z[i]-z[i-1],2))
    #sx[i] = sqrt(pow(x[i],2))
    #sy[i] = sqrt(pow(y[i],2))
    #sz[i] = sqrt(pow(z[i],2))
  sx[0]=sx[1]
  sy[0]=sy[1]
  sz[0]=sz[1]
  sa[0]=sa[1]
  goHistogram(sx,"sx")#,png="sx-hist")
  goHistogram(sy,"sy")#,png="sy-hist")
  goHistogram(sz,"sz")#,png="sz-hist")
  goHistogram(sa,"sa")#,png="sa-hist")
  return sz

###########################################################################


def plot2D(st,sz,img,title=None,pt=None,pd=None,v1=None,v2=None,tlim1=None,tlim2=None,cn=None,c1=None,cbintv=None,cbwmin=None,paper=None,slides=None):
  # For v1 and v2
  dt = st.getDelta()
  dz = sz.getDelta()
  ft = st.getFirst()
  fz = sz.getFirst()
  lt = st.getLast()
  lz = sz.getLast()
  t1,t2,z0 = zerofloat(2),zerofloat(2),zerofloat(2)
  z0[0] = fz; z0[1] = lz;
  pp = PlotPanel()
  pv = pp.addPixels(st,sz,img)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.JET)
  if pd and pt:
    cp = [0]
    if c1: 
      cp = [c1]
    if cn: 
      cp = range(cn)
    for i in cp:
        pt2 = mul(pt[i],dt)
        pd2 = mul(pd[i],dz)
        pd2 = add(pd2,fz)
        t1[0] = pt2[len(pt[i])-1] 
        ln = pp.addPoints(pt2,pd2)
        ln.setLineColor(WHITE)
        ln.setLineWidth(2)
    if v1 and v2:
      t2[0] = t1[0]
      t1[1] = t1[0]+2*lz/v1; 
      t2[1] = t2[0]+2*lz/v2; 
      vlim1 = pp.addPoints(t1,z0)
      vlim2 = pp.addPoints(t2,z0)
      vlim1.setLineColor(WHITE)
	    #vlim1.setLineStyle(PointsView.Line.DASH)
      vlim1.setLineWidth(3)
      vlim2.setLineColor(WHITE)
      vlim2.setLineWidth(3)
	    #vlim2.setLineStyle(PointsView.Line.DASH)
    if tlim1:
      tl1,tl2 = zerofloat(2),zerofloat(2)
      tl1[0] = tlim1*dt+ft; tl1[1] = tlim1*dt+ft; 
      tl2[0] = tlim2*dt+ft; tl2[1] = tlim2*dt+ft; 
      tlims1 = pp.addPoints(tl1,z0)
      tlims2 = pp.addPoints(tl2,z0)
      tlims1.setLineColor(BLACK)
      tlims1.setLineStyle(PointsView.Line.DASH)
      tlims1.setLineWidth(3)
      tlims2.setLineColor(BLACK)
      tlims2.setLineStyle(PointsView.Line.DASH)
      tlims2.setLineWidth(3)
  #pp.setHLabel("Time (s)")
  #pp.setVLabel("Tau (s)")
  if title: pp.setTitle(title)
  pp.setLimits(ft,fz,lt,lz)
  cb = pp.addColorBar()#"Error")
  if cbwmin: pp.setColorBarWidthMinimum(cbwmin)
  if cbintv:
    cb.setInterval(cbintv)
  pp.setVLimits(sz.getFirst(),sz.getLast())
  pp.setVInterval(0.2)
  frame = PlotFrame(pp)
  frame.setSize(1169,344) # Perfect value to frame ratio
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  if paper:
    frame.setSize(1002,329)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    frame.setSize(700,350)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)

def plotImgCurves2(st,sz,img,t,z,c1=None,cn=None,pt=None,pd=None,paper=None,slides=None):
  p1 = PlotPanel(2,2,PlotPanel.Orientation.X1RIGHT_X2UP,PlotPanel.AxesPlacement.LEFT_TOP)
  p1.mosaic.setWidthElastic(0,0)
  p1.mosaic.setHeightElastic(0,0)
  p1.mosaic.setWidthElastic(1,100)
  p1.mosaic.setHeightElastic(1,100)
  pp1 = p1.addPoints(0,1,st,t)
  pp1.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
  p1.setVLabel(0,"trace")
  p1.setVLabel(1,"Tau (s)")
  p1.setHLabel(0,"synth")
  p1.setHLabel(1,"Time (s)")
  pp2 = p1.addPoints(1,0,sz,z)
  pp2.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
  im = p1.addPixels(1,1,st,sz,img)
  im.setInterpolation(PixelsView.Interpolation.NEAREST)
  im.setColorModel(ColorMap.JET)
  im.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
  #p1.addColorBar()
  #p1.setColorBarWidthMinimum(100)
  # ADDED PATHS AND CURVES
  dz=sz.getDelta(); dt=st.getDelta();
  fz=sz.getFirst(); ft=st.getDelta();
  if pd and pt:
    cp = [0]
    if c1: 
      cp = [c1]
    if cn: 
      cp = range(cn)
    for i in cp:
      pt2 = mul(pt[i],dt)
      pd2 = mul(pd[i],dz)
      pd2 = add(pd2,fz)
      ln = p1.addPoints(1,1,pt2,pd2)
      ln.setLineColor(WHITE)
      ln.setLineWidth(3)
      ln.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
  p1.setHLimits(1,st.getFirst(),st.getLast())
  p1.setVLimits(1,sz.getFirst(),sz.getLast())
  p1.addColorBar()
  frame = PlotFrame(p1)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setSize(1500,1000)
  if paper:
    frame.setSize(1000,700)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    frame.setSize(1115,525)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)

def plotImgB(st,sz,img,pe=None,c1=None,cn=None,pt=None,pd=None,pt1=None,paper=None,slides=None):
  p1 = PlotPanel(2,1,PlotPanel.Orientation.X1RIGHT_X2UP,PlotPanel.AxesPlacement.LEFT_TOP)
  p1.mosaic.setHeightElastic(0,25)
  p1.mosaic.setHeightElastic(1,100)
  p1.setVLabel(1,"Tau (s)")
  p1.setHLabel(0,"Time (s)")
  p1.setVLabel(0,"D")
  b = add(img[len(img)-1],0)
  if pe:
    b2 = zerofloat(len(b))
    if pt1:
      for i in range(len(b)):
        if (i>=pt1[0][0] and i<=pt1[len(pt1)-1][0]):
          b2[i] = b[i]/len(pe[i])
        else:
       	  b2[i] = b[i]
    else:
      for i in range(len(b)):
        b2[i] = b[i]/len(pe[i])
    p1.setVLabel(0,"Dn")
    p1.addPoints(0,0,st,b2)
    p1.setVLimits(0,min(b2)*0.5,min(b2)*3.0)
    p1.setVFormat(0,"%.1f")
    p1.setVInterval(0,0.25)
  else: 
    p1.addPoints(0,0,st,b)
    p1.setVInterval(0,300)
    p1.setVLimits(0,0,max(b))
    p1.setVFormat(0,"%1.0f")
  im = p1.addPixels(1,0,st,sz,img)
  im.setInterpolation(PixelsView.Interpolation.NEAREST)
  im.setColorModel(ColorMap.JET)
  im.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
  #p1.addColorBar()
  #p1.setColorBarWidthMinimum(100)
  # ADDED PATHS AND CURVES
  dz=sz.getDelta(); dt=st.getDelta();
  fz=sz.getFirst(); ft=st.getDelta();
  if pd and pt:
    cp = [0]
    if c1: 
      cp = [c1]
    if cn: 
      cp = range(cn)
    for i in cp:
      pt2 = mul(pt[i],dt)
      pd2 = mul(pd[i],dz)
      pd2 = add(pd2,fz)
      ln = p1.addPoints(1,0,pt2,pd2)
      ln.setLineColor(WHITE)
      ln.setLineWidth(3)
      ln.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
  p1.setHLimits(0,st.getFirst(),st.getLast())
  p1.setVLimits(1,sz.getFirst(),sz.getLast())
  p1.addColorBar()
  frame = PlotFrame(p1)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setSize(1500,1000)
  if paper:
    frame.setSize(1000,700)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    frame.setSize(900,560)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)


def plotCurvePanel2(x,xs,y,ys,v,unit,title=None,Hlims=None,png=None):
  cs = Sampling(double(ys))
  p1 = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  p2 = PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  p1.addPoints(0,0,xs,x)
  lred = p1.addPoints(0,0,cs,y)
  lred.setLineColor(RED)
  p2.addPoints(0,1,cs,ys)
  if (unit=="d"):
    p1.setVLabel("Depth (km)")
    p1.setVLimits(0,xs.getFirst(),xs.getLast())
    p1.setHLabel(0,"Trace (Red) and Synthetic (Black)")
    p2.addPoints(0,0,cs,v)
    p2.setVLabel("Depth (km)")
    p2.setHLabel(0,"Velocity (km/s)")
    p2.setHLabel(1,"Time (s)")
    p2.setVLimits(0,xs.getFirst(),xs.getLast())
    p2.setHLimits(0,0,max(v))
  if (unit=="t"):
    p1.setVLabel("Time (s)")
    p1.setVLimits(0,cs.getFirst(),cs.getLast())
    p1.setHLabel(0,"Trace (Black) and Synthetic (Red)")
    p2.addPoints(0,0,cs,v)
    p2.setVLabel("Time (s)")
    p2.setHLabel(0,"Velocity (km/s)")
    p2.setHLabel(1,"Depth (km)")
    p2.setHLimits(0,0,max(v))
  if title:
    p1.setTitle(title)
    p2.setTitle(title)
  if Hlims:
    p1.setHLimits(0,-Hlims,Hlims)
  frame = PlotFrame(p1,p2,PlotFrame.Split.HORIZONTAL)
  frame.setSize(900,1500)
  frame.setFontSize(24)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True)
  if png:
    frame.paintToPng(1000,5.00,_pngDir+png+'.png')

def plotImgCurves(st,sz,img,t,z,c1=None,cn=None,pt=None,pd=None,v1=None,v2=None,tlim1=None,tlim2=None,tax=None,paper=None,slides=None):
  p1 = PlotPanel(2,2,PlotPanel.Orientation.X1RIGHT_X2UP,PlotPanel.AxesPlacement.LEFT_TOP)
  p1.mosaic.setWidthElastic(1,25)
  p1.mosaic.setHeightElastic(1,25)
  p1.mosaic.setWidthElastic(0,100)
  p1.mosaic.setHeightElastic(0,100)
  pp1 = p1.addPoints(1,0,st,t)
  pp1.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
  p1.setHLabel(1,"Synthetic")
  p1.setVLabel(0,"Depth (km)")
  p1.setVLabel(1,"Trace")
  p1.setHLabel(0,"Time (s)")
  if tax:
    p1.setVLabel(0,"Time (s)")
  pp2 = p1.addPoints(0,1,sz,z)
  pp2.setOrientation(PointsView.Orientation.X1DOWN_X2RIGHT)
  im = p1.addPixels(0,0,st,sz,img)
  im.setInterpolation(PixelsView.Interpolation.NEAREST)
  im.setColorModel(ColorMap.JET)
  im.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP)
  p1.addColorBar()
  #p1.setColorBarWidthMinimum(100)
  # ADDED PATHS AND CURVES
  dz=sz.getDelta(); fz=sz.getFirst(); lz=sz.getLast();
  dt=st.getDelta(); ft=st.getFirst(); lt=st.getLast();
  z0 = zerofloat(2); z0[0] = fz; z0[1] = lz;
  if pd and pt:
    cp = [0]
    if c1: 
      cp = [c1]
    if cn: 
      cp = range(cn)
    for i in cp:
      pt2 = mul(pt[i],dt)
      pd2 = mul(pd[i],dz)
      pd2 = add(pd2,fz)
      ln = p1.addPoints(0,0,pt2,pd2)
      ln.setLineColor(WHITE)
      ln.setLineWidth(3)
      ln.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
      if v1 and v2:
        t1,t2 = zerofloat(2),zerofloat(2)
        t1[0] = pt2[len(pt[i])-1] 
        t2[0] = t1[0]
        t1[1] = t1[0]+2*lz/v1; 
        t2[1] = t2[0]+2*lz/v2; 
        vlim1 = p1.addPoints(0,0,z0,t1); 
        vlim2 = p1.addPoints(0,0,z0,t2);
        vlim1.setLineColor(WHITE); vlim1.setLineStyle(PointsView.Line.DASH); vlim1.setLineWidth(3); 
        vlim2.setLineColor(WHITE); vlim2.setLineStyle(PointsView.Line.DASH); vlim2.setLineWidth(3);
        vlim1.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
        vlim2.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
  if tlim1:
    tl1,tl2 = zerofloat(2),zerofloat(2)
    tl1[0] = tlim1*dt+ft; tl1[1] = tlim1*dt+ft; 
    tl2[0] = tlim2*dt+ft; tl2[1] = tlim2*dt+ft; 
    tlims1 = p1.addPoints(0,0,z0,tl1); 
    tlims2 = p1.addPoints(0,0,z0,tl2);
    tlims1.setLineColor(BLACK); tlims1.setLineStyle(PointsView.Line.DASH); tlims1.setLineWidth(3);
    tlims2.setLineColor(BLACK); tlims2.setLineStyle(PointsView.Line.DASH); tlims2.setLineWidth(3);
    tlims1.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
    tlims2.setOrientation(PointsView.Orientation.X1RIGHT_X2UP)
  p1.setLimits(ft,fz,lt,lz)
  frame = PlotFrame(p1)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(28)
  frame.setSize(1500,1000)
  if paper:
    frame.setSize(1000,700)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    frame.setSize(1500,1000)
    frame.setFontSizeForSlide(8.0,504.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)

def plotCurvePanel(lx,lsx,tx,tsx,vt,vd,cb,title1=None,title2=None,lims=None,png=None):
  p1 = PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  p1.addPoints(0,0,lsx,lx)
  p1.setVLabel("Depth (km)")
  if title1:
    p1.setHLabel(0,title1)
  #p1.setTitle("Log and TD Curve")
  p2 = PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  p2.addPoints(0,1,tsx,tx)
  p2.setVLabel("Time (s)")
  if title2:
    p2.setHLabel(1,title2)
  else:
    p2.setHLabel("Trace")
  p1.addPoints(0,1,lsx,vd[cb])
  p2.addPoints(0,0,tsx,vt[cb])
  p1.setHLabel(1,"Velocity km/s")
  p2.setHLabel(0,"Slowness s/km")
  if lims:
    p1.setHLimits(0,-lims,lims)
    p2.setHLimits(1,-lims,lims)
  frame = PlotFrame(p2,p1,PlotFrame.Split.HORIZONTAL)
  frame.setSize(900,1500)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True)
  if png:
    frame.paintToPng(500,3.33,_pngDir+png+'.png')

def plotCurvePanelT(tr,st,y,ys,sst,ss,tm=None,lims=None,paper=None,slides=None,paper2=None,paper3=None):
  if tm:
    p1 = PlotPanel(1,3,PlotPanel.Orientation.X1DOWN_X2RIGHT)
    p1.setHLabel(2,"tau")
    stau = Sampling(ss.getCount(),ss.getDelta(),tm[0])
    p1.addPoints(0,2,stau,tm)
  else: 
    p1 = PlotPanel(2,1)#,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  l1 = p1.addPoints(0,0,st,tr)
  l1.setLineWidth(2)
  p1.setHLabel("Time (s)")
  #ss2 = Sampling(ss.getCount(),ss.getDelta(),sst.getFirst())
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
  p1.setHLimits(0,sst.getLast()+0.5)
  frame = PlotFrame(p1)
  frame.setSize(900,1500)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True)
  if paper:
    frame.setSize(1360,354)
    frame.setFontSizeForPrint(8.0,504.0)
    frame.paintToPng(720,7.00,_pngDir+paper+'.png')
  if paper2:
    frame.setSize(904,379)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper2+'.png')
  if paper3:
    frame.setSize(904,379)
    frame.setFontSizeForPrint(8.0,312.0)
    frame.paintToPng(720,4.33,_pngDir+paper3+'.png')
  if slides:
    frame.setSize(500,900)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)


def plotCurvePanel3(tr,st,y,ys,sst,ss,tm=None,lims=None,paper=None,slides=None):
  if tm:
    p1 = PlotPanel(1,3,PlotPanel.Orientation.X1DOWN_X2RIGHT)
    p1.setHLabel(2,"tau")
    stau = Sampling(ss.getCount(),ss.getDelta(),tm[0])
    p1.addPoints(0,2,stau,tm)
  else: 
    p1 = PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  l1 = p1.addPoints(0,0,st,tr)
  l1.setLineWidth(2)
  p1.setVLabel("Time t (s)")
  #ss2 = Sampling(ss.getCount(),ss.getDelta(),sst.getFirst())
  l2 = p1.addPoints(0,0,ss,y)
  l2.setStyle('r-')
  l2.setLineWidth(2)
  l3 = p1.addPoints(0,1,st,tr)
  l3.setLineWidth(2)
  l4 = p1.addPoints(0,1,sst,ys)
  l4.setStyle('r-')
  l4.setLineWidth(2)
  if lims:
    p1.setHLimits(0,-lims,lims)
    p1.setHLimits(1,-lims,lims)
    if tm: p1.setHLimits(2,-lims,lims)
  p1.setVLimits(0,sst.getLast()+0.5)
  frame = PlotFrame(p1)
  frame.setSize(900,1500)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True)
  if paper:
    frame.setSize(500,900)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    frame.setSize(668,811)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)

def plotLogPanel(sz,v,d,rf,paper=None,slides=None):
  p1 = PlotPanel(1,3,PlotPanel.Orientation.X1DOWN_X2RIGHT)
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
  p1.setHLimits(2,-0.13,0.13)
  p1.setHInterval(2,0.1)
  #dp.setLineColor(BLUE)
  #gp.setLineColor(GREEN)
  frame = PlotFrame(p1)
  frame.setSize(800,1500)
  frame.setFontSize(28)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  if paper:
    frame.setSize(548,738)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    frame.setSize(690,561)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)


def plotCurve(c,title,fill=None,s=None,d=None,lim=None,intv=None,ss=None,png=None,paper=None,slides=None):
  p1 = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT,PlotPanel.AxesPlacement.LEFT_TOP)
  #p1 = PlotPanel()
  if s:
    pv = p1.addPoints(s,c)
    if d:
      p1.setVLabel("Depth (km)")
    else:
      p1.setVLabel("Time (s)")
  else:
    pv = p1.addPoints(c)
    p1.setVLabel("Samples")
  if lim:
    p1.setHLimits(-lim,lim)
  if intv:
    p1.setHInterval(intv)
  if fill:
    f = zerofloat(len(c))
    for i in range(0,len(c),2):
      if c[i]>0:
        f[i] = c[i]
    for i in range(0,len(c)):
      if c[i]<0:
        f[i] = c[i]
    fl = p1.addPoints(s,f)
    fl.setStyle("-")
    fl.setLineWidth(3)
  if ss:
    p1.setVLimits(ss.getFirst(),ss.getLast())
  p1.setHLabel(title)
  frame = PlotFrame(p1)
  #frame.setSize(1500,300)
  frame.setSize(500,1500)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  if png:
    frame.paintToPng(900,3.33,_pngDir+png+'.png') 
  if paper:
    frame.setSize(548,738)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    frame.setSize(260,700)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)

def plotCurveN(st,tr,trn,ss,sy,syn,lim=None,paper=None):
  p1 = PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  p2 = PlotPanel(1,2,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  p1.addPoints(0,0,st,tr)
  p1.addPoints(0,1,st,trn)
  p2.addPoints(0,0,ss,sy)
  p2.addPoints(0,1,ss,syn)
  p1.setVLabel(0,"Time (s)")
  p2.setVLabel(0,"t (s)")
  p2.setHLimits(0,-1.1,1.1)
  frame = PlotFrame(p1,p2,PlotFrame.Split.HORIZONTAL)
  frame.setSize(300,1500)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True)
  if paper:
    frame.setSize(584,757)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')

def plotCurveW(sz,rf,ss,sy,paper=None):
  p1 = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  p2 = PlotPanel(1,1,PlotPanel.Orientation.X1DOWN_X2RIGHT)
  p1.addPoints(0,0,sz,rf)
  p2.addPoints(0,0,ss,sy)
  p1.setVLabel(0,"Depth (km)")
  p2.setVLabel(0,"t (s)")
  p1.setHLimits(0,-0.13,0.13)
  p2.setHLimits(0,-0.5,0.5)
  frame = PlotFrame(p1,p2,PlotFrame.Split.HORIZONTAL)
  frame.setSize(300,1500)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True)
  if paper:
    frame.setSize(584,757)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')

def plotTDCurves(sz,tauz,tz,paper=None,slides=None):
  p1 = PlotPanel()
  l1 = p1.addPoints(sz,tz)
  l2 = p1.addPoints(sz,tauz)
  l1.setLineColor(RED)
  p1.setVLabel('Time (s)')
  p1.setHLabel('Depth (km)')
  frame = PlotFrame(p1)
  frame.setSize(700,500)
  frame.setFontSizeForPrint(8.0,240.0)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  if paper:
    frame.setSize(513,379) #700,500)
    frame.setFontSizeForPrint(8.0,192.0)
    frame.paintToPng(720,2.165,_pngDir+paper+'.png')
  if slides:
    frame.setSize(513,379)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)

def plotCurve2(c1,c2,s1,s2,title,lim=None,png=None,d=None):
  p1 = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  l1 = p1.addPoints(s1,c1)
  l2 = p1.addPoints(s2,c2)
  l1.setLineColor(BLACK)
  l2.setLineColor(RED)
  p1.setVLabel("Samples")
  if d:
    p1.setVLabel("Depth (km)")
  if lim:
    p1.setHLimits(-lim,lim)
  p1.setTitle(title)
  frame = PlotFrame(p1)
  frame.setSize(300,1500)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True)
  if png:
    frame.paintToPng(1080,3.33,_pngDir+png+'.png') 

def plot2DLog(img,x,yc,st,sy,zot,tops=None,tz=None,ghz=None,ghy=None,nolog=None,paper=None,slides=None):
  nx = len(x)
  ny = len(img)
  x2 = zerofloat(nx,3)
  copy(x,x2[0])         
  copy(x,x2[1])
  copy(x,x2[2])         
  sy2 = Sampling(3,sy.getDelta(),(yc-1)*sy.getDelta())
  sp = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT) 
  im1 = sp.addPixels(st,sy,img)
  im1.setColorModel(ColorMap.GRAY)#RED_WHITE_BLUE)#GRAY)
  if ghz:
    # Add Horizons
    pv1 = sp.addPoints(ghz[0],ghy[0])
    pv1.setStyle("g-")
    pv1.setLineWidth(2)
    pv2 = sp.addPoints(ghz[1],ghy[1])
    pv2.setStyle("c-")
    pv2.setLineWidth(2)
    pv3 = sp.addPoints(ghz[2],ghy[2])
    pv3.setStyle("y-")
    pv3.setLineWidth(2)
    pv4 = sp.addPoints(ghz[3],ghy[3])
    pv4.setStyle("m-")
    pv4.setLineWidth(2)
  if tops:
    for i in range(len(tops)):
      if (tops[i]==1):
        it = int((tz[i]-tz[0])/st.getDelta())
	x2[0][it] = 100; x2[1][it] = 100; x2[2][it] = 100;
	x2[0][it-1] = 100; x2[1][it-1] = 100; x2[2][it-1] = 100;
	x2[0][it+1] = 100; x2[1][it+1] = 100; x2[2][it+1] = 100;
  log = sp.addPixels(zot,sy2,x2)
  log.setColorModel(ColorMap.GRAY)#RED_WHITE_BLUE)#GRAY)
  log.setInterpolation(PixelsView.Interpolation.NEAREST)
  log.setClips(min(img),max(img))
  if nolog: sp.remove(log)
  if ((zot.getFirst()-0.2)<0): zotgf = 0
  else: zotgf = zot.getFirst()-0.2
  sp.setLimits(yc*sy.getDelta()-1.2,zotgf,yc*sy.getDelta()+1.2,zot.getLast()+0.2)
  #sp.setLimits(yc*sy.getDelta()-1.2,0.45,yc*sy.getDelta()+1.2,1.4)
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

  
def plotPhaseRMS(mrms,nph,phf,png=None,paper=None,slides=None):
  sp = SimplePlot()
  sph = Sampling(nph/phf,phf,0.0)
  sp.addPoints(mrms)
  sp.setHLabel("Phase (degrees)")
  sp.setVLabel("Minimum Dn")
  sp.setTitle("Minumum RMS vs. Phase")
  sp.setSize(1500,1000)
  if png:
    sp.paintToPng(1200,6.00,_pngDir+png+'.png')
  if paper:
    sp.removeTitle()
    sp.setSize(1000,800)
    sp.setFontSizeForPrint(8.0,240.0)
    sp.paintToPng(720,3.33,_pngDir+paper+'.png')  
  if slides:
    sp.removeTitle()
    sp.setSize(950,800)
    sp.setFontSizeForSlide(1.0,1.0)
    sp.paintToPng(720,7.00,_pngDir+slides+'.png')

def plotDphase(bn,st,slides=None,paper=None):
  p1 = PlotPanel()
  p1.setVLabel("Dn")  
  p1.setHLabel("Time (s)")
  ln1 = p1.addPoints(st,bn[0])
  ln2 = p1.addPoints(st,bn[1])
  ln3 = p1.addPoints(st,bn[2])
  ln2.setLineColor(RED)
  ln3.setLineColor(BLUE)
  mbn = min(bn[0])
  if mbn>0:
    p1.setVLimits(0,mbn*0.5,mbn*2.0)
  p1.setVFormat(0,"%.1f")
  p1.setVInterval(0,0.1)
  frame = PlotFrame(p1)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  if paper:
    frame.setSize(1248,572)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    frame.setSize(700,400)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)

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
  si.paintToPng(900,3.33,_pngDir+"RICKER"+".png")

def plotb(sx,sy,img,pe,paper=None,slides=None):
  n1 = len(img[0]) #long
  n2 = len(img)    #short
  
  p2 = PlotPanel(2,1,PlotPanel.Orientation.X1RIGHT_X2UP,PlotPanel.AxesPlacement.LEFT_TOP)
  #img2 = transpose(img)
  ba = add(img[n2-1],0.)
  b2 = zerofloat(len(ba))
  for i in range(len(ba)):
    b2[i] = ba[i]/len(pe[i])
  p2.addPoints(0,0,sx,b2)
  p2.setVLimits(0,min(b2)*0.5,min(b2)*3.0)
  p2.setVFormat(0,"%.1f")
  p2.setVInterval(0,0.2)
  p2.mosaic.setHeightElastic(0,25)
  p2.mosaic.setWidthElastic(0,100)
  p2.mosaic.setHeightElastic(1,100)
  #p2.setVLabel(1,"f(x)")
  #p2.setHLabel(1,"g(x)")
  p2.addPixels(1,0,sx,sy,img)
  p2.addColorBar()
  frame = PlotFrame(p2)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setSize(1500,300)
  if paper:
    frame.setSize(1002,700)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    frame.setSize(1500,1000)
    frame.setFontSizeForSlide(8.0,504.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)

def plotSubsq(st,st2,f,g,st3=None,g2=None,tlim=None,paper=None,slides=None):
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
  if tlim:
    p1.setHLimits(0,tlim)
  else: p1.setHLimits(0,st.getLast())
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


"""
Histogram plot 
@author Farhad
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
"""
def plot2P(st,sz,x,y,img,title,png=None):
  #p1 = PlotPanel() 
  #ln1 = p1.addPoints(x)
  #ln2 = p1.addPoints(y)
  #ln1.setLineColor(BLACK)
  #ln2.setLineColor(BLUE)
  #p1.setHLabel("Samples")
  #p1.setVLabel("Amplitude")
  #p1.setTitle(title)
  
  p2 = PlotPanel(PlotPanel.Orientation.X1RIGHT_X2UP)
  pv = p2.addPixels(st,sz,img)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  pv.setColorModel(ColorMap.JET)
  #pv.setColorModel(ColorMap.GRAY)
  #p2.addColorBar()
  #p2.setColorBarWidthMinimum(100)
  p2.setHLabel("Time (s)")
  p2.setVLabel("Depth (km)")
  p2.setTitle(title)
  #frame = PlotFrame(p1,p2,PlotFrame.Split.VERTICAL)
  frame = PlotFrame(p2)
  frame.setSize(1500,1000)
  frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
  frame.setVisible(True)
  if png:
    frame.paintToPng(360,3.33,_pngDir+png+'.png') 
"""

##########################################################################################
"""
Test signals
@author Luming Liang
"""

def sig1L(n1):
  seed   = 23232
  seedN1 = 31415 
  fpeak  = 0.125
  rsmNoise = 0.1
  X = Dtw1().makeEvents(n1, seed) 
  X = Dtw1().addRickerWavelet(fpeak, X)
  X = mul(1.0/max(abs(X)),X)
  X = Dtw1().addNoise(rsmNoise, seedN1, X)
  return X

def sig1S(n1,n2):
  _si = SincInterpolator()
  _si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT)
  seedN2 = 10000
  rsmNoise = 0.5
  amp = 30
  Y = zerofloat(n2)
  t = zerofloat(n2)
  X = sig1L(n1);
  w = float((2*PI/n2))
  start = 300
  for i in range (0,n2):
    t[i]=start+i-amp*sin(i*w)    #shift between the signals. 'start' corresponds to where Y starts on X. 
  _si.setUniform(n1,1.0,0.0,X) #(length of X, uniform sampling interval, the value of the first unform sample of X, array)
  _si.interpolate(n2,t,Y)     #(num of output samples, array of t values at which to interpolate Y, array of interpolated output Y
  Y = Dtw1().addNoise(rsmNoise, seedN2, Y)
  return Y

#---------------------------------------------------------#

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

 

"""
# Phase Stuff... 
  
  #plotCurve2(well,wellps,sz,sz,"Orig (black) Phase Rotated (red) by "+str(ps))#,lim=2)   
  nph = 361 # number of phases
  phf = 1 # phase skip factor
  minrms = zerofloat(nph/phf)
  im,imm = 0,0
  #phaseset = [0,136,189,282]
#  for iph in range(0,nph,phf):
#    
#    wellr = applyPhaseRot(iph,well)
#    cm  = Dtw1().cm(trace[tr],wellr)
#    acm = Dtw1().acm(cm) 
#    pe,tpi,zpi = Dtw1().paths(binterval,bst,cm,tui,tli)
#    rms = rms2(pe)
#    #rms = RMS(pe)
#    ri,mrms = minimums(rms)
#    minrms[im] = mrms
#    if (minrms[im]<minrms[imm]):
#      miniph = iph
#      imm = im
#    im+=1
#  plotPhaseRMS(minrms,nph,phf,png="Phase plot")

    #Spectrum().apply(sz,well,"Well")
    #Spectrum().apply(st,trace[tr],"Trace")


  #Plots
#  plot2D(st,sz,acm," Smoothed Path-"+str(ip)+" with a RMS of "+str(minrms[miniph])+" Phase Rot- "+str(iph),pt=tp2,pd=zp2,val="v",curve=ip+1,v1=v1,v2=v2) 
#  plotCurvePanel2(trace[tr],st,wellr,zot,vt,"t",Hlims=3,title='D2T for Phase- '+str(iph))
#  #plotCurvePanel2(well,sz,trace[tr],toz,vz,"d",Hlims=3,title='T2D for Phase- '+str(iph))
#  #plot2DLog(tsdataxc,wellr,yc,st,sy,Sampling(double(zot)),"Log in Seismic Slice- "+str(iph))  
#  #for i in range(np):
#  plot2D(st,sz,acm," Path-"+str(ip)+" with a RMS of "+str(minrms[miniph])+" Phase Rot- "+str(iph),pt=tpi,pd=zpi,val="i",curve=ip+1,v1=v1,v2=v2) 
#  plot2D(st,sz,acm,"Minimum Distance Paths"+" Phase Rot- "+str(iph),pt=tpi,pd=zpi,val="i",v1=v1,v2=v2,tlim1=tli,tlim2=tui)
#  plot2D(st,sz,acm,"Accumulated Distance"+" Phase Rot- "+str(iph))
#  plot2D(st,sz,cm,"Errors "+" Phase Rot- "+str(iph),pt=tpi,pd=zpi,v1=v1,v2=v2,val="i")
#  plot2D(st,sz,cm,"Errors"+" Phase Rot- "+str(iph))



def goDtwTestData():
  rf,tr = getTestData()
  sz = Sampling(len(rf),0.0001524,0.0)
  st = Sampling(len(tr),0.001,0.0)
  # First test the convolve then stretch method
  #testCS(rf,tr,sz,st)
  testSC(rf,tr,sz,st)

def testLims(m,st,sz):
  dt = st.getDelta()
  dz = m*sz.getDelta()
  v1 = 4*dz/dt
  v2 = dz/dt 
  tl = 2*sz.getFirst()/v2 #lower v bound set by water.
  tu = 2*sz.getFirst()/v1
  tli = int(tl/dt)
  tui = int(tu/dt)
  print 'New dz is '+str(dz)+' km'
  print 'Velocity Bounds are '+str(v2)+' to '+str(v1)+' km/s'
  print 'Deepest time is '+str(tl)+' s'+' and it\'s index is '+str(tli)
  print 'Shallowest time is '+str(tu)+' s'+' and it\'s index is '+str(tui)
  return v1,v2,tli,tui

def getTestData():
  ais1 = ArrayInputStream("./data/test-reflectivity.dat")
  ais2 = ArrayInputStream("./data/test-trace.dat")
  rf = zerofloat(3280)
  tr = zerofloat(731)
  ais1.readFloats(rf)
  ais2.readFloats(tr)
  ais1.close(); ais2.close();
  return rf,tr



def goDtwTest():
  n1,n2 = 1000, 300     # n1 > n2 always... 
  binterval = 40 #int(tn2*0.25)  # make local minumum search window half n2
  bst = int(0.1*n1)
  x = sig1L(n1)
  y = sig1S(n1,n2)
  xs = Sampling(n1,1.0,0.0)
  ys = Sampling(n2,1.0,0.0)
 #DTW
  costm  = Dtw1().cm(x,y)
  acostm = Dtw1().acm(costm) 
  xvals,yvals,pathx,pathy,mean,rms = Dtw1().wpath1(binterval,bst,costm)
  #h = align(minb,yvals[minb],x,y)
 #Plots
  #plotCurvePanel(h,ys,x,xs,xvals,yvals,minb,"H","X",2,png="X-H-TD-CurvePanel")
  plotCurvePanel(y,ys,x,xs,xvals,yvals,minb,"Y","X",2,png="X-Y-TD-CurvePanel")
  plotCurve(h,"H")
  plotCurve(y,"Y")
  plotCurve(x,"X")
  for i in range(len(rms)):
    plot2D(xs,ys,acostm,"Path-"+str(i)+" with a RMS of "+str(rms[i]),pt=xvals,pd=yvals,curve=i+1) 
  plot2D(xs,ys,acostm,"Accumulated Distance and Paths",pt=xvals,pd=yvals)
  plot2D(xs,ys,acostm,"Accumulated Distance")
  plot2D(xs,ys,costm,"|f-g|^2 'The Error Matrix'")


Convert the well to time, convolve with a wavelet, convert back to depth, then apply dtw

def warpData2(sx,sy,sz,st,trace,tsdataxc,yc,fp,binterval,bst,ail,vl):
  dz = sz.getDelta(); dt = st.getDelta();

  # Find sigmas
  vs = 2 #(km/s) velocity scalar
  sigma = vs/(2*FLT_PI*fp*dz)
  print 'sigma= ',sigma
  sigf = 8*sigma
  sigft = 4/(FLT_PI*fp*dt)
 
  # Constant stationary wavelet applied in time
  rf = reflcalc(ail)
  nr = len(rf)
  vzot = zerofloat(nr)
  RecursiveGaussianFilter(nr*0.1).apply(vl,vzot) #smooth the log velocity
  fzot = TDutils().makezot(vzot,sz)
  rt,ftoz = stretch(fzot,rf,sz) # Stretch to time with the same sampling in seconds 
  rtw = TDutils().applyTWavelet(rt,sz,fp)
  rzw,dm = stretch(ftoz,rtw,sz) # Unstretch back to depth after applying the wavelet
  nz2 = len(rzw)
  sz = Sampling(nz2,dz,0.0)
  synt = normalizeRMSlocal(rzw,sigf)

  tr = normalizeRMSlocal(tr,sigft)

  # Get lims and decimate
  m = 10 
  v1,v2,tli,tui = testLims(m,st,sz)
  synth = decimate(synt,m)
  sz = sz.decimate(m)
  slope = float(len(synth))/float(len(tr))
  print 'slope= ',slope

  # DTW
  cm  = Dtw1().cm(tr,synth)
  acm = Dtw1().acm(cm) 
  pe,tpi,zpi = Dtw1().paths(binterval,bst,cm,tui,tli)
  # DT conversion
  minc = zerofloat(len(pe)); ri=0;
  for ip in range(len(pe)):
    tp2,zp2 = reverseInputThrow(tpi[ip],zpi[ip],st,sz)
    #tp2,zp2 = smoothTDcurve(tp2,zp2,2)
    toz,zot = interpolations(tp2,zp2,st,sz)
    #toz2,zot2 = smoothTDcurve(toz,zot,6)
    #vz,vt = velocities(toz2,zot2,sz,st)
    wellt,dm = stretch(zot,synth,st)
    sts = Sampling(len(wellt),st.getDelta(),zot[0])
    crms = TDutils().correlation(tr,wellt,st,sts)
    print "%d with %f rms"%(ip,crms)
    if ip==0: minc[0]=crms
    if (crms>minc[ip]):
      ri=ip 

  tp2,zp2 = reverseInputThrow(tpi[ri],zpi[ri],st,sz)
  #tp2,zp2 = smoothTDcurve(tp2,zp2,2)
  toz,zot = interpolations(tp2,zp2,st,sz)
  #toz2,zot2 = smoothTDcurve(toz,zot,6)
  #vz,vt = velocities(toz2,zot2,sz,st)
  wellt = stretch(zot,synth,st)
  sts = Sampling(len(wellt),st.getDelta(),zot[0])

  # Plots
  plotImgCurves(st,sz, cm,tr,synth)
  plotImgCurves(st,sz,acm,tr,synth)
  plotImgCurves(st,sz,acm,tr,synth,cn=len(tpi),pt=tpi,pd=zpi)
  #plotImgCurves(st,sz,acm,tr,synth,cn=len(tpi),pt=tpi,pd=zpi,v1=v1,v2=v2)
  #plotImgCurves(st,sz,acm,tr,synth,cn=len(tpi),pt=tpi,pd=zpi,tlim1=tli,tlim2=tui)
  plotImgCurves(st,sz,acm,tr,synth,c1=ri,pt=tpi,pd=zpi)
  #plotCurvePanel2(tr,st,synth,zot,vt,"t",Hlims=3)




"""
