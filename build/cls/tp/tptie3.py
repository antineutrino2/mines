# Multiple 3D well ties
# @author Andrew Munoz, Colorado School of Mines
# @version 04.10.2013

from dtw import *
from tputils import *
from wtutils import *
from imports import *

#_pres = '_rsch'
_pres = 'paper_'

#printWellLengths()
# Params:
fp = 35
q  = 100
# Image warping
dr1   = 0.02 # 2D smoothness vertically 
dr2   = 0.50 # 2D smoothness horizontally
dr3   = 0.50 # 3D smoothness horizontally
rmin1,rmax1 = -0.06,0.12 # vertical constraint # use -0.06,0.12,0.02,0.5,-1.0,1.0
rmin2,rmax2 = -1.00,1.00 # horizontal constraint (s/km)
rmin3,rmax3 = -1.00,1.00 # horizontal constraint (s/km)
smin,smax = -0.10,0.10 # min and max time shifts for 2D/3D warping
print 'dr1  = ',dr1,  'dr2  = ',dr2,  'dr3  = ',dr3
print 'rmin1=', rmin1,'rmin2=', rmin2,'rmin3=', rmin3
print 'rmax1= ',rmax1,'rmax2= ',rmax2,'rmax3= ',rmax3

#3D interp	
vinterps3=False
dinterps3=False
sinterps3=False
#2D interp
vinterps2=False
dinterps2=False
sinterps2=False
#plots
singleplots=True

def main(args):
  #printWellLengths()
  goMultiWell3D()

def goMultiWell3D():
  # 1. Get the 3D seismic image
  global s1,s2,s3
  global smin,smax
  global wells,ids
  cut1=[0.10,1.25];cut2=[3.40,6.65];cut3=[0.6,2.7];
  #f,s3,s2,s1 = getImage(normalize=True,smooth=True,cut1=cut1,cut2=cut2,cut3=cut3)
  f,s3,s2,s1,fpun = getImage(normalize=True,cut1=cut1,cut2=cut2,cut3=cut3,rnrm=True)
  #f,s3,s2,s1 = getImage(normalize=True)
  dt = s1.delta
  
  # 2. Get the wells, coords, and tops
  #wells = [0,1,2,3,4,6,7,8,9,10,11,12,13,15,16,17]
  #wells = [0,1,  3,4,      9,10,11,12,   15,16   ] #>=100 #15
  #wells = [0,3,4,9,15,16,12,17,2,11] #seis6
  #wells = [3,15,16,11,4]; nm="1"; #seis1
  #phases = [71]*len(wells) #nm="1"
  #wells = [3,15,16,11,10]; nm="2"; #seis2
  #phases = [100]*len(wells) #15
  #wells = [3,12,16,15,4,11]; nm="3" # 3D
  #phases = [30]*len(wells) #nm=3
  wells = [3,12,16,15,7,4,2]; nm="7" # 3D take out 0,1,11
  phases = [44]*len(wells) #nm=6,7
  wells = sortWells(wells)
  tops,fmsp,datums = getWellTops(wells)
  ail,vl,dl,gl,sz,x2w,x3w,ids = getMultipleLogs(wells,datums,id=True)
  
  # 3. Get the synthetics
  ssy,sy,tzl,ztl,sztl = getSeismograms(sz,vl,dl,q,phase=phases)
  printCorrelations(f,sy,ssy,x2w,x3w,'initial')

  plotNewCoords(f,[3,16,15,7,2],[3,12,4,11])
  return
  
  # 3. Apply single warp to all wells
  smin1 = -1.0; smax1 = 1.0
  u,h,sh,tz,vint,zt,szt,pv = warpIndividually(\
  			f,ssy,sy,sz,x2w,x3w,tzl,vl,smin1,smax1,rmin1,rmax1,dr1,wells)#,ph2=True)
  
  print '---3D warp refinement---'
  # 7. Interpolate and warp the synthetics
  bgs = 2.0
  tensors = makeTensors(f,p2=1.0)
  #plot3(f,"3D Tensors",tensors=tensors,paper='3d tensors',plots=True)
  
  print 'Interpolate 1D synthetics...'
  kind = "swt_synimg"+nm+"_"+str(len(sy))
  sw = Stopwatch()
  sw.start()
  #gi = interpolateSynthetics3(sy,ssy,x2w,x3w,bgs,kind=kind,ts=tensors)
  #gt = interpolateSynthetics3(h,sh,x2w,x3w,bgs,kind=kind,ts=tensors)
  sw.stop()
  pts = [x3w,x2w]
  #plot3P(s1,s2,s3,gi,"synthetic image",paper='isynimg')
  print 'interpolation took '+(getTimePrint(sw.time()))
  
  #print 'f3 = ',len(f) ,' f2 = ',len(f[0]) ,' f1 = ',len(f[0][0])
  #print 'g3 = ',len(gi),' g2 = ',len(gi[0]),' g1 = ',len(gi[0][0])
  
  #h3d,u3d,hn,shn,tzn,vintn,ztn,sztn,pvn,name = warpAllin3D(gi,f,x2w,x3w,sy,ssy,vl,tzl,sz,kind)#,new=True)
  #h3d,u3d,hn,shn,tzn,vintn,ztn,sztn,pvn,name = warpAllin3D(gt,f,x2w,x3w,h,sh,vl,tz,sz,kind+"_swtf_")#,new=True)
  #readAndMakePlots(nm,f,sy,ssy,tzl,sz,x2w,x3w,bgs,gi,h3d,u3d,hn,shn,tz,vintn,pvn,name)

  # Plotting
  #lconmw = [ail,vl,dl,sz,tzl,tzn,vintn,sy,ssy,hn,shn,pvn,phases]
  lconsw = [ail,vl,dl,sz,tzl,tz,vint,sy,ssy,h,sh,pv,phases]
  #line1
  #wellline1 = [0,1,7,2];2
  wellline2 = [3,16,15,7,2]; 
  #getProfile(wells,wellline2,lconsw,f,fg,h3d,u3d,edg=0.1,wln="4",plot1Donly=True)
  #line2
  wellline1 = [3,16,15,7,4,11]; 
  #getProfile(wells,wellline2,lconmw,f,gt,h3d,u3d,edg=0.1,wln="_swtf",plot1Dtie=[h,sh,tz])
  #line3
  wellline3 = [3,12,4,11];
  #getProfile(wells,wellline3,lconmw,f,gi,h3d,u3d,edg=0.1,wln="3")
 
  wellsl = [15]#[3,4,15,16]#,11,7,12,2]
  if singleplots:
    for i in wellsl:
      goSingleWellPlots(f,x2w,x3w,lconsw,well=i,swt=True,vpl=True,syn=fpun)
      #goSingleWellPlots(f,x2w,x3w,lconmw,well=i)

  if vinterps3:
    print 'Interpolate updated velocities...'
    vtensors = makeTensors(f,p0=0.0,p1=1.50,p2=1.50)
    bgs = 1.0
    sw.restart()
    vname = "nvel"+nm
    veln = interpolateLogs3(vintn,shn,x2w,x3w,tzn,bgs,\
    												 kind=vname+"_"+name,ts=vtensors) 
    vname = "svel"+nm
    vels = interpolateLogs3(vl,ssy,x2w,x3w,tzn,bgs,\
    												 kind=vname+"_"+name,ts=vtensors)
    vname = "ivel"+nm
    veli = interpolateLogs3(vl,ssy,x2w,x3w,tzl,bgs,\
    												 kind=vname+"_"+name,ts=vtensors)
    sw.stop()
    v21 = getImageAlongCurve(veln,x2s,x3s)
    v22 = getImageAlongCurve(vels,x2s,x3s)
    v23 = getImageAlongCurve(veli,x2s,x3s)
    plotSlice(sn,s1,v21,title='3D new velocity',cmap=jet,cbar="Velocity (km/s)")
    plotSlice(sn,s1,v22,title='3D shifted velocity',cmap=jet,cbar="Velocity (km/s)")
    plotSlice(sn,s1,v23,title='3D inital velocity',cmap=jet,cbar="Velocity (km/s)")
    print 'interpolation took '+(getTimePrint(sw.time()))
   
  if sinterps3:
    print 'Interpolate updated seismograms...'
    bgs = 2.0
    sw.restart()
    sname = "nsynimg"+nm
    gn = interpolateSynthetics3(hn,shn,x2w,x3w,bgs,\
    												 kind=sname+"_"+name,ts=tensors)
    sw.stop()
    plot3(gn,"updated synthetic image",logs=[hn,shn,x2w,x3w],plots=True,lmap=[x2w,x3w,x2s,x3s],paper ='ipnsynimg')
    gi2 = getImageAlongCurve(gi,x2s,x3s)
    h2n = getImageAlongCurve(h3d,x2s,x3s)
    gn2 = getImageAlongCurve(gn,x2s,x3s)
    plotSlice(sn,s1,gn2,title='3D tied synthetic image')#,paper ='2dpnsynimg')
    plotSlice(sn,s1,h2n,title='3D warped synthetic image',paper ='2dpwarpsynimg')
    plotSlice(sn,s1,gi2,title='3D init synthetic image',paper ='2dpinitsynimg')
    print 'interpolation took '+(getTimePrint(sw.time()))
  
  if dinterps3:
    print 'Interpolate updated time-depths...'
    bgs = 0.5
    sw.restart()
    dname = "ntz"+nm
    sd,tz3n = computetz3(dname+"_"+name,tzn,sz,x2w,x3w,bgs)
    plot3(tz3n,"time-depths",sd=sd,cmap=jet,paper ='iptz',plots=True)
    sw.stop()
    tz3n = smooth23(tz3n)
    plot3(tz3n,"time-depths smoothed",sd=sd,cmap=jet,paper ='iptz',plots=True)
    print 'interpolation took '+(getTimePrint(sw.time()))
    print 'Depth-convert 3D seismic image...'
    gz = imageTimeToDepth3(f,tz3n)
    plot3(gz,"depth-converted seismic image",sd=sd,paper ='ipzimg',plots=True)
    gz2 = getImageAlongCurve(gz,x2s,x3s,sd=sd)
    tlim = [sn.first,s1.first,sn.last,max(tz3n)]
    dlim = [sn.first,sd.delta*getFirstZ(s1.first,tz3n),sn.last,sd.last]
    plotSlice(sn,sd,gz2,lims=dlim,title='depth-converted image',vlabel="Depth (km)",paper ='2dpnzimg')
    plotSlice(sn,s1,cf,lims=tlim,title='seismic image',paper ='2dpimg')
    plot3(f,"seismic image",paper ='iptimg',plots=True)
    plot3(gz,"depth-converted seismic image",sd=sd,paper ='ipzimg',plots=True)
  return	
  if vinterps2 or dinterps2 or sinterps2:
    lcon = [ail,vl,dl,sz,tzl,tzn,vintn,sy,ssy,hn,shn,vint,tz]
    ail2,vl2,dl2,sz2,tzl2,tzn2,vintn2,sy2,ssy2,hn2,shn2,vint2,tz2=reduceWellList(wells,wellline1,lcon)
  if vinterps2:
    print 'Interpolate velocities in 2D...'
    bgs = 1.5
    cv = [2.74,5.5]
    vtensors = makeTensors(cf,p1=1.0)
    vlmg = interpolateLogs2(vl2,sz2,sn,x2wn,tzl2,bgs,ts=vtensors)
    vumg = interpolateLogs2(vintn2,sz2,sn,x2wn,tzn2,bgs,ts=vtensors)
    vsmg = interpolateLogs2(vl2,sz2,sn,x2wn,tzn2,bgs,ts=vtensors)
    vimg = interpolateLogs2(vint2,sz2,sn,x2wn,tz2,bgs,ts=vtensors)
    vt#mg = interpolateLogs2(vl2,sz2,sn,x2wn,tz2,bgs,ts=vtensors)
    plotSlice(sn,s1,vlmg,cmap=jet,cbar='Velocity (km/s)',title='init velocities',clips=cv)
    plotSlice(sn,s1,vsmg,cmap=jet,cbar='Velocity (km/s)',title='shifted velocities',clips=cv)
    plotSlice(sn,s1,vumg,cmap=jet,cbar='Velocity (km/s)',title='tied velocities',clips=cv)
    plotSlice(sn,s1,vimg,cmap=jet,cbar='Velocity (km/s)',title='indiv warp velocities',clips=cv)
    plotSlice(sn,s1,vtmg,cmap=jet,cbar='Velocity (km/s)',title='indiv shift velocities',clips=cv)
    plotSlice(sn,s1,cf,v=vlmg,title='init velocities',clips2=cv,paper='2dpvelinit')
    plotSlice(sn,s1,cf,v=vsmg,title='shifted velocities',clips2=cv,paper='2dpvelshift')
    plotSlice(sn,s1,cf,v=vumg,title='tied velocities',clips2=cv,paper='2dpvelwarp')
    plotSlice(sn,s1,cf,v=vimg,title='indiv warp velocities',clips2=cv,paper='2dpvel1Dwarp')
    plotSlice(sn,s1,cf,v=vtmg,title='indiv shift velocities',clips2=cv,paper='2dpvel1Dshift')
  if sinterps2:
    print 'Interpolate synthetics in 2D...'
    bgs = 2.0
    stensors = makeTensors(cf,p1=1.0)
    #gi2 = getImageAlongCurve(gi,x2s,x3s)
    #h2n = getImageAlongCurve(h3d,x2s,x3s)
    #gn2 = getImageAlongCurve(gn,x2s,x3s)
    gi2 = interpolateSynthetics2(sy2,ssy2,sn,x2wn,bgs,ts=stensors)
    gn2 = interpolateSynthetics2(hn2,shn2,sn,x2wn,bgs,ts=stensors)
    plotSlice(sn,s1,gn2,title='tied synthetic image')#,paper ='2dpnsynimg')
    #plotSlice(sn,s1,h2n,title='warped synthetic image',paper ='2dpwarpsynimg')
    plotSlice(sn,s1,gi2,title='init synthetic image',paper ='2dpinitsynimg')
    #plot2DLogs(s1,sn,cf,shn2,hn2,x2wn,title='tied synthetics',paper ='2dpwarpsyn')
    #plot2DLogs(s1,sn,cf,ssy2,sy2,x2wn,title='init synthetics',paper ='2dpinitsyn')


def warpAllin3D(f,g,x2w,x3w,sy,ssy,vl,tzl,sz,kind,new=None):
  nw = len(vl)
  dt = s1.delta
  f1 = s1.first
  n1 = s1.count; n2 = s2.count; n3 = s3.count
  name = str(nw)+"_"+str(rmin1)+"_"+str(rmin2)+"_"+str(rmin3)+"_"+str(rmax1)+"_"+\
  			 str(rmax2)+"_"+str(rmax3)+"_"+str(dr1)+"_"+str(dr2)+"_"+str(dr3)+"_"+\
  			 str(smin)+"_"+str(smax)+"_"+str(n3)+"_"+str(n2)+"_"+str(n1)
  tpDir = getTpDir()
  fname = tpDir+"csm/warps/"+kind+"_3dwarpu_"+name+".dat"
  dwc = DynamicWarpingWT(smin,smax,s1,s2,s3)
  #dwo = DynamicWarpingWTM(g,f,infl(smin/dt),ince(smax/dt))
  if path.exists(fname) and path.isfile(fname) and not new:
    print "u exists..."
    u = zerofloat(n1,n2,n3)
    readImage(u,fname)
  else:
    print "find u..."
    dwc.setStrainLimits(rmin1,rmax1,rmin2,rmax2,rmin3,rmax3)
    dwc.setSmoothing(dr1,dr2,dr3)
    u = dwc.findShifts(s1,f,s1,g)
    #u = dwo.findShifts(rmin1,rmax1,dr1,rmin2,rmax2,dr2,rmin3,rmax3,dr3)
    writeImage(u,fname)
  h = dwc.applyShifts(f,u)
  #h,s = dwo.applyShifts(u)
  #h = dwo.applyShiftsF(u)
  ##plotSlice(s2,s1,f[50],title="f")
  ##plotSlice(s2,s1,g[50],title="g")
  ##plotSlice(s2,s1,h[50],title="h")
  ##plotSlice(s2,s1,u[50],cmap=jet,cbar="u(tau)")
  wh,wsh,tz,vint,zt,szt,pv=[],[],[],[],[],[],[]
  rsum = 0
  u = mul(dt,u)
  for i in range(nw):
    x2i = inro((x2w[i]-s2.first)/s2.delta)
    x3i = inro((x3w[i]-s3.first)/s3.delta)
    ns = ssy[i].count
    fs = ssy[i].first
    f2 = copy(ns,inro((fs-f1)/dt),f[x3i][x2i])
    u2 = copy(ns,inro((fs-f1)/dt),u[x3i][x2i])
    h0 = DynamicWarpingWT.applyShiftsS(f2,u2,ssy[i])
    wsh0 = Sampling(len(h0),dt,u2[0]+fs)
    #h0,fha = dwc.getWarpedSyntheticSeismogram(ssy[i],u2,sy[i])
    #wsh0  = Sampling(len(h0),dt,fha[0]) 
    tz0,vint0,pv0 = updateTDandVelocity(sz[i],u2,tzl[i],ssy[i])
    zt0,szt0 = getzt(sz[i],tz0)
    tz0 = add(tz0,fs-f1)
    print wsh0.first
    print tz0[0]
    wh.append(h0); wsh.append(wsh0); pv.append(pv0);
    szt.append(szt0);
    tz.append(tz0); vint.append(vint0); zt.append(zt0)
    trc = f[x3i][x2i]
    grc = g[x3i][x2i]
    wn = str(wells[i])
    ##plotCurve(sy[i],ssy[i],h2,s1,title='init syn syn img'+wn)
    ##plotCurve(trc,s1,sy[i],ssy[i],title='init syn '+wn)
    ##plotCurve(grc,s1,h0,wsh0,title='tied syn cut'+wn)
    ##plotCurve(grc,s1,h2,s1,title='tied syn img '+str(wells[i]))
    ##plotCurve(pv0,sz[i],title='% v change '+wn,d="z")
    #pred+=printError(s1,wsh[i],h0,f[x3i][x2i]); szt.append(szt0)
    #r = correlation(h0,trc,wsh0,s1)
    #rsum += r
    #print 'R = ',r,'for well ',wn
  print 'Average R =',rsum/nw
  return h,u,wh,wsh,tz,vint,zt,szt,pv,name

def getNewSyntheticSampling(h,u,ssy,sy=None):
  ni = ssy.count
  fi = ssy.first
  li = ssy.last
  n1 = s1.count
  d1 = s1.delta
  f1 = s1.first
  # find beginning and end of the synthetic seismogram
  fii = inro((fi-f1)/d1)
  lii = inro((li-f1)/d1)
  u = add(floats(s1.getValues()),u)
  fo = u[fii]
  lo = u[lii]
  fio = inro((fo-f1)/d1)
  no = inro((lo-fo)/d1)
  if no+fio>n1: no -= no+fio-n1
  sh = Sampling(no,d1,fo)
  ho = copy(no,fio,h)
  #plotCurve(h,s1,ho,sh,title='b4 and after')
  #plotCurve(ho,sh,title='syn after')
  #if sy: plotCurve(sy,ssy,title='syn b4')
  #plotCurve(u,s1,title='shifts')
  return ho,sh

def getCutU(u2,ssy):
  ni = ssy.count
  fi = ssy.first
  d1 = s1.delta
  f1 = s1.first
  # find beginning and end of the synthetic seismogram
  fii = inro((fi-f1)/d1)
  u = copy(ni,fii,u2)
  u = add(floats(ssy.getValues()),u)
  return u

def getProfile(wells,wellline,lcon,f,gi,h3d,u3d,edg=0.2,hzs=None,plot3d=None,plot1Dtie=None,wln="1",plot1Donly=None):
  imgh = WTutils().getSlice(290,f)
  tops,fmsp,datums = getWellTops(wells)
  ail,vl,dl,gl,sz,x2w,x3w,ids  = getMultipleLogs(wells,datums)
  x2wl,x3wl = getMultipleLogs(wellline,datums,conly=True)
  x2s,x3s,sn,x2wn = getCurveThruWells(x2wl,x3wl,s2.delta,edg=edg)
  plotCoordSlice(s2,s3,x2s,x3s,x2w,x3w,imgh,paper ="seismic coord hslice"+wln)
  ail2,vl2,dl2,sz2,tzl2,tzn2,vintn2,sy2,ssy2,hn2,shn2,pv2,phases \
                                        = reduceWellList(wells,wellline,lcon)
  cf = getImageAlongCurve(f,x2s,x3s)
  gi2 = getImageAlongCurve(gi,x2s,x3s)
  ws = ' '+str(wellline)
  if not plot1Donly:
    h2n = getImageAlongCurve(h3d,x2s,x3s)
    u2n = getImageAlongCurve(u3d,x2s,x3s)
    plotSlice(sn,s1,gi2,title='3D inital image'+ws,paper ='2dpnsynimg'+wln)
    plotSlice(sn,s1,h2n,title='3D warped synthetic image'+ws,paper ='2dpwarpsynimg'+wln)
    #plotSlice(sn,s1,u2n,title='3D vertical time shifts'+ws,cbar='time shifts (s)',cmap=jet)
    plot2DLogs(s1,sn,cf,shn2,hn2,x2wn,title='tied synthetics'+ws,paper ='2dpwarpsyn'+wln)
    plot2DLogs(s1,sn,cf,ssy2,sy2,x2wn,title='init synthetics'+ws,paper ='2dpinitsyn'+wln)
    plot2DLogs(s1,sn,cf,shn2,hn2,x2wn,title='tied synthetics'+ws,paper ='2dpwarpsyncb'+wln,cbar="Amplitude")
    if hzs: 
    	tops,fmsp = hzs
    	gh2 = get2DHorizons(x2s,x3s)
    	plot2DLogs(s1,sn,cf,shn2,hn2,x2wn,gh=gh2,tops=tops,fmsp=fmsp,tz=tzn2,sz=sz2)#,paper='2dpwarpsynhz'+wln)
  cv = [2.75,5.5]
  votli = makeVot(sz2,vl2,tzl2,ssy2)
  plot2DLogs(s1,sn,cf,ssy2,votli,x2wn,title='vel init',cbar='Velocity (km/s)',\
  								lcmap=jet,lclips=cv,paper='2dpveliseis'+wln,velc=True)
  votl = makeVot(sz2,vl2,tzn2,shn2)
  plot2DLogs(s1,sn,cf,shn2,votl,x2wn,title='vel shift3',cbar='Velocity (km/s)',\
        lcmap=jet,lclips=cv,paper='2dpvelseistie'+wln,velc=True)
  bgs = 1.5
  vtensors = makeTensors(cf,p1=1.0)
  vlmg = interpolateLogs2(vl2,sz2,sn,x2wn,tzl2,bgs,ts=vtensors)
  vsmg = interpolateLogs2(vl2,sz2,sn,x2wn,tzn2,bgs,ts=vtensors)
  plotSlice(sn,s1,cf,v=vlmg,title='init',clips2=cv,paper='2dpvelinit'+wln)
  plotSlice(sn,s1,cf,v=vsmg,title='tied',clips2=cv,paper='2dpveltie'+wln)
  if plot1Dtie:
    h,sh,tzi = plot1Dtie
    h2i,sh2i,tzi2 = reduceWellList(wells,wellline,[h,sh,tzi])
    plot2DLogs(s1,sn,cf,sh2i,h2i,x2wn,title='swt'+ws,paper='2dp1dwarpsyn'+wln)
    plot2DLogs(s1,sn,cf,sh2i,h2i,x2wn,title='swt'+ws,paper='2dp1dwarpsyncb'+wln,cbar="Amplitude")
    votl1 = makeVot(sz2,vl2,tzi2,sh2i)
    v1mg = interpolateLogs2(vl2,sz2,sn,x2wn,tzi,bgs,ts=vtensors)
    plot2DLogs(s1,sn,cf,sh2i,votl1,x2wn,title='vel shift1',\
        cbar='Velocity (km/s)',lcmap=jet,lclips=cv,\
        paper='2dpvel1Dseis'+wln,velc=True)
    plotSlice(sn,s1,cf,v=v1mg,title='swt',clips2=cv,paper='2dpswtveltie'+wln)
  if plot3d:
    sy,ssy,hn,shn = lcon[7],lcon[8],lcon[9],lcon[10]
    plot3(f,"untied synthetics",logs=[sy,ssy,x2w,x3w],\
    			plots=True,lmap=[x2w,x3w,x2s,x3s],paper ='ipinitsyn')
    plot3(gi,"initial synthetic image",logs=[sy,ssy,x2w,x3w],\
    			plots=True,lmap=[x2w,x3w,x2s,x3s],paper ='ipinitimg')
    plot3(h3d,"warped synthetic image",logs=[hn,shn,x2w,x3w],\
    			plots=True,lmap=[x2w,x3w,x2s,x3s],paper ='ipwarpimg')
    plot3(f,"warped synthetics",logs=[hn,shn,x2w,x3w],\
  			plots=True,lmap=[x2w,x3w,x2s,x3s],paper ='ipwarpsyn')
    if plot1Dtie:	
  	  plot3(f,"1D tie synthetics",logs=[h,sh,x2w,x3w],\
  					plots=True,lmap=[x2w,x3w,x2s,x3s],paper ='ipindwarpsyn')
  return x2s,x3s

def goSingleWellPlots(f,x2w,x3w,lcon,well,hzs=None,swt=None,vpl=None,syn=None):
  mt = "mwt"
  if swt:
    mt = "swt"
  iw = wells.index(well)
  wid = mt+str(well)
  ttitle = "UWI "+str(ids[iw])
  ail,vl,dl,sz,tzl,tzn,vln,sy,ssy,hn,shn,pvn,phases = lcon
  traces = getTraces(f,x2w[iw],x3w[iw],tonly=True)
  print ttitle
  z = sz[0].delta
  tr = traces[0]
  #plotCurveHorz2(tr,s1,sy[iw],ssy[iw],tr,s1,hn[iw],shn[iw],paper="synspanel"+wid,title=ttitle)
  #plotVelPanel(sz[iw],tzl[iw],vl[iw],pvn[iw],sy12=tzn[iw],sy22=vln[iw],hlim=[-15,15],paper="veltzpanel"+wid,title=ttitle)
  plotCurveHorz2(tr,s1,sy[iw],ssy[iw],tr,s1,hn[iw],shn[iw],slides="slides_synspanel"+wid,title=ttitle)
  plotVelPanel(sz[iw],tzl[iw],vl[iw],pvn[iw],sy12=tzn[iw],sy22=vln[iw],hlim=[-15,15],slides="slides_veltzpanel"+wid,title=ttitle)
  plotLogPanel(sz[iw],vl[iw],dl[iw],reflectivity(ail[iw]),paper="logpanel"+wid)
  #if hzs:
  #	x2wl,x3wl  = getMultipleLogs(wellline2,datums,conly=True)
  #	x2s,x3s,sn,x2wn = getCurveThruWells(x2wl,x3wl,s2.delta,edg=0.1)
  #	cf = getImageAlongCurve(f,x2s,x3s)
  if syn:
    # Test simple syn, syn w and wout mult, w and wout attenuation (5 panels)
    ph = phases[iw]
    tru = getTraces(syn,x2w[iw],x3w[iw],tonly=True)
    trn = normalizeRMSlocal(tru[0],100)
    uwi = ids[iw]
    wset = getLogDataset("d")
    propDir = getPropDir()
    wt = WellTie(uwi,dz,s1,wset)
    #wt.makeNewSeismograms()
    wt.setCutsOff()
    tz0 = wt.tz0
    wt.setNormalizationOff()
    wt.makeSimpleSeismogram(fp)
    sfs = wt.f
    ssfs = wt.sf
    wt.setMultiplesOff()
    wt.makePropagatorSeismogram(propDir,fp,1e6)
    sfp1 = wt.f
    ssfp1 = wt.sf
    wt.setMultiplesOn()
    wt.makePropagatorSeismogram(propDir,fp,1e6)
    sfp2 = wt.f
    ssfp2 = wt.sf
    wt.makePropagatorSeismogram(propDir,fp,100)
    sfp3 = wt.f
    ssfp3 = wt.sf
    # Test before and after normalization
    wt.setCutsOn()
    wt.makePropagatorSeismogram(propDir,fp,100)
    sfp3c = wt.f
    ssfp3c = wt.sf
    wt.setNormalizationOn()
    wt.makeSimpleSeismogram(fp,ph)
    sfsn = wt.f
    ssfsn = wt.sf
    wt.makePropagatorSeismogram(propDir,fp,100,ph)
    sfp3n = wt.f
    ssfp3n = wt.sf
    vlims=[ssfs.first,ssfp3.last]
    hlims=[-0.36,0.36]
    hlimsn=[-3.01,3.01]
    plotSeismogramPanel(s1,s1,ssfp3c,tru[0],trn,sfp3c,ssy4=ssfp3n,sy4=sfp3n,paper="normpanel"+wid,hlim1=hlimsn,hlim2=hlimsn,hlim3=hlims,hlim4=hlimsn)
    plotSeismogramPanel(ssfs,ssfp1,ssfp2,sfs,sfp1,sfp2,ssy4=ssfp3,sy4=sfp3,paper="simpproppanel"+wid,hlim1=hlims,hlim2=hlims,hlim3=hlims,hlim4=hlims,vlim=vlims)
    # Test prop ties vs. simple ties
    wt.makeSimpleSeismogram(fp,ph)
    wt.computeSingleWellTie(traces,-0.12,0.12,-1.0,1.0,0.02)
    u1,h1,sh1,tz1,vl1,pv1 = wt.u,wt.h,wt.sh,wt.tz1,wt.v1,wt.pv
    wt.makePropagatorSeismogram(propDir,fp,100,ph)
    wt.computeSingleWellTie(traces,-0.12,0.12,-1.0,1.0,0.02)
    u2,h2,sh2,tz2,vl2,pv2 = wt.u,wt.h,wt.sh,wt.tz1,wt.v1,wt.pv
    ssfsnt = Sampling(ssfsn.count,ssfsn.delta,ssfsn.first+0.004)
    plotCurve(sfp3n,ssfp3n,sfsn,ssfsnt)
    hlimsn2=[-4.01,4.01]
    plotSeismogramPanel(ssfsnt,ssfp3n,None,sfsn,sfp3n,None,paper="normsimpprop"+wid,hlim1=hlimsn2,hlim2=hlimsn2,rm3=True,size=[343,433])
    #u1,h1,sh1=getWarpingShiftsOld(ssfsn,s1,sfsn,traces,-0.10,0.10,0.02)
    #u2,h2,sh2=getWarpingShiftsOld(ssfp3n,s1,sfp3n,traces,-0.10,0.10,0.02)
    #tz1,vl1,pv1 = updateTDandVelocity(sz[iw],u1,tzl[iw],ssfsn)
    #tz2,vl2,pv2 = updateTDandVelocity(sz[iw],u2,tzl[iw],ssfp3n)
    plotVelPanel(sz[iw],tzl[iw],vl[iw],pv1,sy12=tz1,sy22=vl1,paper="simptiecomp"+wid,title=ttitle+"s",hlim=[-20,20])
    plotVelPanel(sz[iw],tzl[iw],vl[iw],pv2,sy12=tz2,sy22=vl2,paper="proptiecomp"+wid,title=ttitle+"p",hlim=[-20,20])
    plotCurveHorz2(tr,s1,sfsn,ssfsn,tr,s1,h1,sh1,paper="simptiecompsyn"+wid,title=ttitle)
    plotCurveHorz2(tr,s1,sfp3n,ssfp3n,tr,s1,h2,sh2,paper="proptiecompsyn"+wid,title=ttitle)
    if vpl:
      wt.makePropagatorSeismogram(propDir,fp,100,ph)
      wt.computeSingleWellTie(traces,-0.5,0.5,smin,smax,0.50)
      u1,h1,sh1,tz1,vl1,pv1 = wt.u,wt.h,wt.sh,wt.tz1,wt.v1,wt.pv
      wt.computeSingleWellTie(traces,-0.1,0.1,smin,smax,0.10)
      u2,h2,sh2,tz2,vl2,pv2 = wt.u,wt.h,wt.sh,wt.tz1,wt.v1,wt.pv
      wt.computeSingleWellTie(traces,-0.1,0.1,smin,smax,0.02)
      u3,h3,sh3,tz3,vl3,pv3 = wt.u,wt.h,wt.sh,wt.tz1,wt.v1,wt.pv
      #u1,h1,sh1=getWarpingShiftsOld(ssfp3n,s1,sfp3n,traces,-1.0,1.0,1.0) 
      #u2,h2,sh2=getWarpingShiftsOld(ssfp3n,s1,sfp3n,traces,-0.1,0.1,0.1) 
      #u3,h3,sh3=getWarpingShiftsOld(ssfp3n,s1,sfp3n,traces,-0.1,0.1,0.02)
      #tz1,vl1,pv1 = updateTDandVelocity(sz[iw],u1,tzl[iw],ssfp3n)
      #tz2,vl2,pv2 = updateTDandVelocity(sz[iw],u2,tzl[iw],ssfp3n)
      #tz3,vl3,pv3 = updateTDandVelocity(sz[iw],u3,tzl[iw],ssfp3n)

      print 'sh3f=',sh3.first
      plotVelPanel(sz[iw],vl[iw],vl[iw],vl[iw],vl1,vl2,vl3,avel=True,paper='velconsts'+wid,title=ttitle)
      plotVelPanel(sz[iw],pv1,pv2,pv3,pvel=True,paper='pvconsts'+wid,title=ttitle)
      plotVelPanel(sz[iw],tzl[iw],vl[iw],pv1,sy12=tz1,sy22=vl1,paper="comp1panel"+wid,title=ttitle)
      plotVelPanel(sz[iw],tzl[iw],vl[iw],pv2,sy12=tz2,sy22=vl2,paper="comp2panel"+wid,title=ttitle)
      plotVelPanel(sz[iw],tzl[iw],vl[iw],pv3,sy12=tz3,sy22=vl3,paper="comp3panel"+wid,title=ttitle)
      plotCurveHorz2(tr,s1,sfp3n,ssfp3n,tr,s1,h1,sh1,paper="syns1panel"+wid,title=ttitle)
      plotCurveHorz2(tr,s1,sfp3n,ssfp3n,tr,s1,h2,sh2,paper="syns2panel"+wid,title=ttitle)
      plotCurveHorz2(tr,s1,sfp3n,ssfp3n,tr,s1,h3,sh3,paper="syns3panel"+wid,title=ttitle)

   

def readAndMakePlots(nm,f,sy,ssy,tzl,sz,x2w,x3w,bgs,gi,h3d,u3d,hn,shn,tzn,vintn,pvn,name):
  n1 = s1.count; n2 = s2.count; n3 = s3.count;
  bgs = 0.5
  skind = "nsynimg"+nm
  sname = skind+"_"+name+"_"+str(bgs)+"_"+str(n1)+"_"+str(n2)+"_"+str(n3)+".dat"
  dkind = "ntz"+nm
  #dname = dkind+"_"+name+"_"+str(bgs)+"_"+str(n1)+"_"+str(n2)+"_"+str(n3)
  dname = "ntz7_8_-0.06_-1.0_-1.0_0.12_1.0_1.0_0.02_0.5_0.5_-100_100_85_131_576_0.5_576_131_85"
  # vel interps
  #vkinds = ["nvel"+nm,"svel"+nm,"ivel"+nm]
  #vname = vkind+"_"+name+"_"+str(bgs)+".dat"
  fpath = getTpDir()+"csm/interps/"
  sfname = fpath+sname
  dfname = fpath+dname
  gn = like3f(gi)
  gz,tz3,sd = getDepthImage(f,dfname)
  lz = sd.last
  fesz = [sz[0],sz[1],sz[2],sz[3],sz[5],sz[7]]#,sz[4]]#,sz[7]]
  fetzl = [tzl[0],tzl[1],tzl[2],tzl[3],tzl[5],tzl[7]]#,tzl[4]]#,tzl[7]]
  fetzn = [tzn[0],tzn[1],tzn[2],tzn[3],tzn[5],tzn[7]]#,tzn[4]]#,tzn[7]]
  plotAllTDCurves(fesz,fetzl,flz=lz,slides='initialtdcurves')
  plotAllTDCurves(fesz,fetzn,flz=lz,slides='updatedtdcurves')
  lt = max(tzn)
  esz = [sz[3]]
  etz = [tzn[3]]
  plotAllTDCurves(esz,etz,flz=lz,slides='iextrap')
  etz,esz = extraptz(esz[0],etz[0],sd.last,lt)
  plotAllTDCurves([esz],[etz],flz=lz,slides='sextrap')
  vavg = div(mul(2,floats(esz.getValues())),etz)
  vavg[0]=vavg[1]
  plotAllVavgCurves([esz],[vavg],flz=lz,slides='vextrap')
  nw = len(fesz)
  etzn,eszn = [None]*nw,[None]*nw
  for i in range(nw):
    etzn[i],eszn[i] = extraptz(fesz[i],fetzn[i],sd.last,lt)
  plotAllTDCurves(eszn,etzn,flz=lz,slides='aextrap')
  vavgb = []
  for i in range(nw):
    vavgb.append(div(mul(2,floats(eszn[i].getValues())),etzn[i]))
    vavgb[i][0]=vavgb[i][1]
  plotAllVavgCurves(eszn,vavgb,flz=lz,slides='vabextrap')
  vavg = []
  for i in range(nw):
    vavg.append(div(mul(2,floats(fesz[i].getValues())),fetzn[i]))
    vavg[i][0]=vavgb[i][1]
  plotAllVavgCurves(fesz,vavg,flz=lz,slides='vanextrap')
  b0,b1 = getRegTerms(vavg,fesz,100)
  vavgg,vsz = [None]*nw,[None]*nw
  lv = max(vavg)
  for i in range(nw):
    vavgg[i],vsz[i]=extrapolate(vavg[i],fesz[i],sd,lv,sd.last,b1,b0)
  plotAllVavgCurves(vsz,vavgg,flz=lz,slides='vagextrap')
  etzn = []
  for i in range(nw):
    etzn.append(div(mul(2,floats(vsz[i].getValues())),vavgg[i]))
  plotAllTDCurves(vsz,etzn,flz=lz,slides='tzagextrap')
  #gn = readImage(gn,sfname)
  #plot3(f,"untied synthetics",logs=[sy,ssy,x2w,x3w],\
  #		  plots=True,lmap=[x2w,x3w,x2s,x3s],paper ='ipinitsyn')
  #plot3(gi,"initial synthetic image",logs=[sy,ssy,x2w,x3w],\
  #			plots=True,lmap=[x2w,x3w,x2s,x3s],paper ='ipinitimg')
  #plot3(h3d,"warped synthetic image",logs=[hn,shn,x2w,x3w],\
  #			plots=True,lmap=[x2w,x3w,x2s,x3s],paper ='ipwarpimg')
  #plot3(f,"warped synthetics",logs=[hn,shn,x2w,x3w],\
  #			plots=True,lmap=[x2w,x3w,x2s,x3s],paper ='ipwarpsyn')
  s4 = Sampling(inro((max(tz3)-s1.first)/s1.delta),s1.delta,s1.first)
  #plot3(f,"seismic time image",sd=s4,plots=True,slides='gtime')
  #plot3(tz3,"time-depth image",plots=True,cmap=jet,slides='gtime')
  #plot3(gz,"seismic depth image",plots=True,slides ='gdepth')
  plot3P(s4,s2,s3,f,"seismic image",paper='3Pf',slc=[0.85,5,1.65])
  plot3P(sd,s2,s3,gz,"depth image",lims=[0.17,sd.last],z=True,paper='3Pgz',slc=[1.51,5,1.65])
  plot3P(sd,s2,s3,tz3,"time-depth image",paper='3Ptz',z=True,cmap=jet,cbar='Time (s)',slc=[1.51,5,1.65])
  vimg = zerofloat(sd.count,s2.count,s3.count)
  for i3 in range(s3.count):
    for i2 in range(s2.count):
      for i1 in range(1,sd.count):
        vimg[i3][i2][i1] = 2*sd.getValue(i1)/tz3[i3][i2][i1]
  plot3P(sd,s2,s3,vimg,"vavg image",paper='s3Pvavg',z=True,cmap=jet,cbar='Velocity (km/s)',slc=[1.51,5,1.65])
  #plot3P(s1,s2,s3,gi,"synthetic image",paper='3Pgi')
  #plot3P(s1,s2,s3,gn,"new synthetic image",paper='3Pgn')
  #plot3P(s1,s2,s3,h3d,"warped synthetic image",paper='3Ph')

def plotNewCoords(f,wells1,wells2):
  edg = 0.1
  imgh = WTutils().getSlice(290,f)
  tops,fmsp,datums = getWellTops(wells1)
  ail,vl,dl,gl,sz,x2w,x3w,ids  = getMultipleLogs(wells1,datums)
  x2wl,x3wl = getMultipleLogs(wells1,datums,conly=True)
  x2s,x3s,sn,x2wn = getCurveThruWells(x2wl,x3wl,s2.delta,edg=edg)
  tops,fmsp,datums = getWellTops(wells2)
  ail,vl,dl,gl,sz,x2w2,x3w2,ids  = getMultipleLogs(wells2,datums)
  x2wl,x3wl = getMultipleLogs(wells2,datums,conly=True)
  x2s2,x3s2,sn,x2wn = getCurveThruWells(x2wl,x3wl,s2.delta,edg=edg)
  pt2 = [x2s2,x3s2,x2w2,x3w2] 
  plotCoordSlice(s2,s3,x2s,x3s,x2w,x3w,imgh,pt2,paper="seismic coord hslice 0")
  plotCoordSlice(s2,s3,x2s,x3s,x2w,x3w,imgh,pt2,True,paper="seismic coord hslice 1")
  plotCoordSlice(s2,s3,x2s,x3s,x2w,x3w,imgh,pt2,True,True,paper="seismic coord hslice 2")
    
def getDepthImage(f,dfname):
  n1 = s1.count; n2 = s2.count; n3 = s3.count;
  sfname = dfname+"_samplings.dat"
  si = zerodouble(3)
  ais = ArrayInputStream(sfname)
  ais.readDoubles(si)
  ais.close()
  sd = Sampling(int(si[0]),si[1],si[2])
  tz3 = zerofloat(sd.count,n2,n3)
  ais = ArrayInputStream(dfname+".dat")
  ais.readFloats(tz3)
  ais.close()
  #tz3 = smooth23(tz3)
  gz = imageTimeToDepth3(f,tz3)
  return gz,tz3,sd




#------------------------------------------------------------------------------#

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

