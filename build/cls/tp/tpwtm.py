# well ties
# @author Andrew Munoz, Colorado School of Mines
# @version 01.21.2014

from dtw import *
from tputils import *
from wtutils import *
from imports import *

# Params:
fp = 35
q  = 100
dz = getDz()
dz2 = getDz2()
# Image warping
dr1   = 0.02 # 2D smoothness vertically 
dr2   = 0.50 # 2D smoothness vertically 
dr3   = 0.50 # 2D smoothness vertically 
r1min,r1max = -0.06,0.12 # vertical constraint # use -0.06,0.12,0.02,0.5,-1.0,1.0
r2min,r2max = -0.50,0.50 # vertical constraint # use -0.06,0.12,0.02,0.5,-1.0,1.0
r3min,r3max = -0.50,0.50 # vertical constraint # use -0.06,0.12,0.02,0.5,-1.0,1.0
smin,smax = -0.050,0.200 # min and max time shifts for 2D/3D warping

propDir = getPropDir()
csmDir = getCsmDir()
wset = getLogDataset("d")
phase = 44

def main(args):
  global s1,s2,s3
  cut1=[0.10,1.25];cut2=[3.40,6.65];cut3=[0.6,2.7];
  g,s3,s2,s1,gzi,s1z,go = getImage(normalize=True,\
                  cut1=cut1,cut2=cut2,cut3=cut3,rnrm=True,z=True)
  wells3 = [3,12,16,15,7,4,2]; # good wells for 3D (took out 11)
  uwis = deepWellSet()
  wells2a = [3,16,15,7,2]; 
  wells2b = [3,12,4,11]; 
  wells2c = [3,16,15,7]; 
  #goMultipleTies2(g,uwis,wells2a)
  #goMultipleTies2(g,uwis,wells2b)
  #goMultipleTies2(g,uwis,wells2c)
  goAllTest(g,uwis,wells3,wells2a)
  #goAllTest(g,uwis,wells3,wells2b)
  ##goMultipleTies3(g,uwis,wells3,s1z,gzi)
  #goSingleTie(g,uwis[15],15)

def goMultipleTies3(g,uwia,wells,s1z,gzi):
  uwis = getUwis(uwia,wells)
  wid = str(wells)
  # Compute well ties
  mwt = MultipleWellTies(uwis,dz,g,s1,s2,s3,wset,csmDir)
  mwt.setStatic(uwis[1]);
  mwt.setStatic(uwis[7]);
  mwt.makeSyntheticSeismograms(propDir,fp,q,phase)
  bgs = 1.0; p0 = 0.0; p1 = 1.0; p2 = 1.0;
  f = mwt.makeInitialSyntheticImage3(bgs,p0,p1,p2)
  h = mwt.computeMultipleWellTies3(r1min,r1max,r2min,r2max,r3min,r3max,\
                                                dr1,dr2,dr3,smin,smax,f)
  #mwt.makeNewTimeDepthFunction()
  # Compute velocities
  #bgs = 1.0; p0 = 0.0; p1 = 1.0; p2 = 1.0;
  #vi = mwt.interpolateInitialVelocity3(bgs,p0,p1,p2)
  #vn = mwt.interpolateUpdatedVelocity3(bgs,p0,p1,p2)
  #vs = mwt.interpolateInitialVelocityTied3(bgs,p0,p1,p2)
  # Compute time-depths
  bgs = 0.5; p0 = 0.0; p1 = 1.0; p2 = 1.0;
  maxz = mwt.getMaxz()
  sd = Sampling(ince(maxz/dz2),dz2,0)
  #extrapolateTest(mwt,sd)
  tzv = mwt.interpolateAverageVelocity3(bgs,p0,p1,p2,sd)
  tzv = smooth23(tzv)
  #tzv = mwt.maketz2(sd,va)
  tzn = mwt.interpolateUpdatedTimeDepths3(bgs,p0,p1,p2,sd)
  tzn = smooth23(tzn)
  #tzi = mwt.interpolateInitialTimeDepths3(bgs,p0,p1,p2,sd)
  gz1 = mwt.convertTimeToDepth3(g,tzv)
  gz2 = mwt.convertTimeToDepth3(g,tzn)
  #gz3 = mwt.convertTimeToDepth3(g,tzi)
  gzic = copy(sd.count,s2.count,s3.count,0,0,0,gzi)
  dwc = DynamicWarpingWT(-0.20,0.20,sd,s2,s3)
  dwc.setSmoothing(0.1,0.5,0.5)
  uc = dwc.findShifts(sd,gz1,sd,gzic)
  #print 'max uc',max(uc)
  #print 'min uc',min(uc)
  #print 'sum uc',sum(uc)
  #plot3(uc,"difference shifts for depth images",plots=True,sd=sd,cmap=jet,slides="gzdiff3")

  #plot3(g,"seismic time image",plots=True)
  #plot3(f,"initial synthetic image",plots=True)
  #plot3(h,"warped synthetic image",plots=True)
  ##plot3(vi,"initial velocity",plots=True,cmap=jet)
  ##plot3(vn,"updated velocity",plots=True,cmap=jet)
  ##plot3(vs,"shifted velocity",plots=True,cmap=jet)
  #plot3(tzv,"time-depth velocity",plots=True,sd=sd,cmap=jet)
  #plot3(tzn,"time-depth updated",plots=True,sd=sd,cmap=jet)
  ##plot3(tzi,"time-depth initial",plots=True,sd=sd,cmap=jet)
  #plot3(g,"seismic time image",plots=True)
  #plot3(gz1,"seismic depth image v",plots=True,sd=sd,g=tzv,cmap2=jet)
  #plot3(gz2,"seismic depth image n",plots=True,sd=sd,g=tzn,cmap2=jet)
  #plot3(gzic,"seismic depth transform",plots=True,sd=sd)
  ##plot3(gz3,"seismic depth image i",plots=True,sd=sd)
  ds = 1.345
  plot3P(sd,s2,s3,gzic,"depth image",lims=[0.17,sd.last],z=True,slides='gziP',slc=[ds,5,1.65])
  ds = 1.285
  plot3P(sd,s2,s3,uc,"u",z=True,cmap=jet,cbar="Depth shift (km)",slides='gzdP',slc=[ds,5,1.65])
  plot3P(sd,s2,s3,tzn,"tz1",z=True,cmap=jet,cbar="Time (s)",slides='tz1P',slc=[ds,5,1.65])
  plot3P(sd,s2,s3,tzv,"tz2",z=True,cmap=jet,cbar="Time (s)",slides='tz2P',slc=[ds,5,1.65])
  plot3P(sd,s2,s3,gz1,"depth image",lims=[0.17,sd.last],z=True,slides='gz1P',slc=[ds,5,1.65])
  plot3P(sd,s2,s3,gz2,"depth image",lims=[0.17,sd.last],z=True,slides='gz2P',slc=[ds,5,1.65])
  plot3P(s1,s2,s3,g,"time image",lims=[s1.first,1.2],slides='gztP',slc=[0.85,5,1.65])

  return
  f1,h1,sf1,sh1,u1,v0,v1,tz0,tz1,sz,pv,x2,x3 = getWellTies(mwt)
  for i in range(len(wells)):
    ix2 = inro((x2[i]-s2.first)/s2.delta)
    ix3 = inro((x3[i]-s3.first)/s3.delta)
    plotCurve(g[ix3][ix2],s1,h1[i],sh1[i])
    plotVelPanel(sz[i],tz0[i],v0[i],pv[i],sy12=tz1[i],sy22=v1[i],hlim=[-15,15])


def goMultipleTies2(g,uwia,wells):
  wells = sortWells(wells)
  uwis = getUwis(uwia,wells)
  wid = str(wells)
  # Make 2D seismic image with profile through wells
  x2,x3 = MultipleWellTies(uwis,dz,s1,wset,csmDir).getWellCoordinates()
  lx2,lx3,s2n,x2n = getCurveThruWells(x2,x3,s2.delta,edg=0.1)
  g2 = getImageAlongCurve(g,lx2,lx3)
  # Compute well ties
  mwt = MultipleWellTies(uwis,dz,g2,s1,s2n,wset,csmDir)
  mwt.makeSyntheticSeismograms(propDir,fp,q,phase)
  mwt.fixX2(x2n)
  bgs = 2.0; p0 = 0.0; p1 = 1.0;
  f = mwt.makeInitialSyntheticImage2(bgs,p0,p1)
  h = mwt.computeMultipleWellTies2(r1min,r1max,r2min,r2max,dr1,dr2,smin,smax,f)
  # Compute velocities
  bgs = 1.0; p0 = 0.0; p1 = 1.0; 
  vi = mwt.interpolateInitialVelocity2(bgs,p0,p1)
  vn = mwt.interpolateUpdatedVelocity2(bgs,p0,p1)
  vs = mwt.interpolateInitialVelocityTied2(bgs,p0,p1)
  # Compute time-depths
  bgs = 0.5; p0 = 0.0; p1 = 1.0; 
  maxz = mwt.getMaxz()
  sd = Sampling(ince(maxz/dz2),dz2,0)
  #extrapolateTest(mwt,sd)
   
  tzv = mwt.interpolateAverageVelocity2(bgs,p0,p1,sd)
  #tzv = mwt.maketz2(sd,va)
  tzn = mwt.interpolateUpdatedTimeDepths2(bgs,p0,p1,sd)
  tzi = mwt.interpolateInitialTimeDepths2(bgs,p0,p1,sd)
  gz1 = mwt.convertTimeToDepth2(g2,tzv)
  gz2 = mwt.convertTimeToDepth2(g2,tzn)
  gz3 = mwt.convertTimeToDepth2(g2,tzi)
  
  pack2 = [s2n,x2n,f,g2,h,vi,vn,vs,sd,tzv,tzn,tzi,gz1,gz2,gz3,mwt]
  makePlots2(pack2)


def makePlots2(pack2):
  s2n,x2n,f,g2,h,vi,vn,vs,sd,tzv,tzn,tzi,gz1,gz2,gz3,mwt = pack2
  f1,h1,sf1,sh1,u1,v0,v1,tz0,tz1,sz,pv,x2,x3 = getWellTies(mwt)
  #nw = len(x2n)
  #for i in range(nw):
  #  ix2 = inro((x2n[i]-s2n.first)/s2n.delta)
  #  plotCurve(g2[ix2],s1,h1[i],sh1[i])
  #  plotVelPanel(sz[i],tz0[i],v0[i],pv[i],sy12=tz1[i],sy22=v1[i],hlim=[-20,20])
  
  plot2DLogs(s1,s2n,g2,sf1,f1,x2n,title='initial synthetics')
  plot2DLogs(s1,s2n,g2,sh1,h1,x2n,title='warped synthetics')
  plotSlice(s2n,s1,g2,title='seismic image')
  plotSlice(s2n,s1,f,title='initial image')
  plotSlice(s2n,s1,h,title='warped image')
  plotSlice(s2n,s1,vi,cmap=jet,cbar="velocity (km/s)",title='initial v')
  plotSlice(s2n,s1,vn,cmap=jet,cbar="velocity (km/s)",title='updated v')
  plotSlice(s2n,s1,vs,cmap=jet,cbar="velocity (km/s)",title='shifted v')
  plotSlice(s2n,sd,tzv,cmap=jet,cbar="time (s)",title='velocity td')
  plotSlice(s2n,sd,tzn,cmap=jet,cbar="time (s)",title='updated td')
  plotSlice(s2n,sd,tzi,cmap=jet,cbar="time (s)",title='initial td')
  plotSlice(s2n,s1,g2,title='seismic image')
  plotSlice(s2n,sd,gz1,title='seismic depth image1')
  plotSlice(s2n,sd,gz2,title='seismic depth image2')
  plotSlice(s2n,sd,gz3,title='seismic depth image3')

  v0t = makeVot(sz,v0,tz0,sf1)
  v1t = makeVot(sz,v1,tz1,sh1)
  v2t = makeVot(sz,v0,tz1,sh1)
  plot2DLogs(s1,s2n,g2,sf1,v0t,x2n,title='initial v',cbar='Velocity (km/s)',\
									lcmap=jet,lclips=[2.75,5.5],velc=True)
  plot2DLogs(s1,s2n,g2,sh1,v1t,x2n,title='updated v',cbar='Velocity (km/s)',\
									lcmap=jet,lclips=[2.75,5.5],velc=True)
  plot2DLogs(s1,s2n,g2,sh1,v2t,x2n,title='shifted v',cbar='Velocity (km/s)',\
									lcmap=jet,lclips=[2.75,5.5],velc=True)


def goAllTest(g3,uwisa,wells3,wells2):
  wells3 = sortWells(wells3)
  uwis = getUwis(uwisa,wells3)
  wells2 = sortWells(wells2)
  uwis2 = getUwis(uwisa,wells2)
  #goSinglesOnly(g3,uwis,uwis2)
  #goSinglesThenMultiple(g3,uwis,uwis2)
  goMultipleOnly(g3,uwis,uwis2,w2=wells2)

def goSinglesOnly(g3,uwis,uwis2,new=None):
  # Do single wellties only
  mwt = MultipleWellTies(uwis,dz,g3,s1,s2,s3,wset,csmDir)
  #mwt.findPhases()
  mwt.setStatic(uwis[2]);
  if new:
    mwt.makeNewInterpolatedLogs()
  mwt.makeSyntheticSeismograms(propDir,fp,q,phase)
  mwt.computeSingleWellTies(r1min,r1max,-1.0,1.0,dr1)
  title = "swt"
  goVelocitiesAndPlot(g3,mwt,uwis2,title)
  return
  pack = getWellTies(mwt,uwis2)
  f1,h1,sf1,sh1,u1,v0,v1,tz0,tz1,sz,pv,x2,x3 = pack
  lx2,lx3,s2n,x2n = getCurveThruWells(x2,x3,s2.delta,edg=0.1)
  g2 = getImageAlongCurve(g3,lx2,lx3)
  plot2DLogs(s1,s2n,g2,sf1,f1,x2n,title=title+'_li',slides='li_'+title)
  plot2DLogs(s1,s2n,g2,sh1,h1,x2n,title=title+'_lt',slides='lt_'+title)
  plot2DLogs(s1,s2n,g2,sf1,f1,x2n,title=title+'_li',slides='lio_'+title,rmb=True)
  plot2DLogs(s1,s2n,g2,sh1,h1,x2n,title=title+'_lt',slides='lto_'+title,rmb=True)

def goSinglesThenMultiple(g3,uwis,uwis2,new=None):
  # Do single then multiple wellties
  mwt = MultipleWellTies(uwis,dz,g3,s1,s2,s3,wset,csmDir)
  #mwt.findPhases()
  mwt.makeSyntheticSeismograms(propDir,fp,q,phase)
  mwt.computeSingleThenMultiple()
  mwt.computeSingleWellTies(r1min,r1max,-1.0,1.0,dr1)
  bgs = 1.0; p0 = 0.0; p1 = 1.0; p2 = 1.0;
  if new:
    mwt.makeNewSyntheticImage()
    mwt.makeNewMultipleWellTie()
    mwt.makeNewInterpolatedLogs()
  f3 = mwt.makeInitialSyntheticImage3(bgs,p0,p1,p2)
  h3 = mwt.computeMultipleWellTies3(r1min,r1max,r2min,r2max,r3min,r3max,\
                                                dr1,dr2,dr3,smin,smax,f3)
  title = "swttmwt"
  #goPlotSynthetics(g3,f3,h3,mwt,uwis2,title)
  goVelocitiesAndPlot(g3,mwt,uwis2,title)

def goMultipleOnly(g3,uwis,uwis2,new=None,w2=None):
  # Do multiple wellties only
  mwt = MultipleWellTies(uwis,dz,g3,s1,s2,s3,wset,csmDir)
  mwt.setStatic(uwis[1])
  if len(uwis)>7: mwt.setStatic(uwis[7]);
  mwt.makeSyntheticSeismograms(propDir,fp,q,phase)
  bgs = 1.0; p0 = 0.0; p1 = 1.0; p2 = 1.0;
  if new:
    mwt.makeNewSyntheticImage()
    mwt.makeNewMultipleWellTie()
    mwt.makeNewInterpolatedLogs()
  f3 = mwt.makeInitialSyntheticImage3(bgs,p0,p1,p2)
  h3 = mwt.computeMultipleWellTies3(r1min,r1max,r2min,r2max,r3min,r3max,\
                                                dr1,dr2,dr3,smin,smax,f3)
  title = "smwt"
  maxz = mwt.getMaxz()
  sd = Sampling(ince(maxz/dz2),dz2,0)
  #extrapolateTest(mwt,sd)
  goPlotSynthetics(g3,f3,h3,mwt,uwis2,title,w2)
  #goVelocitiesAndPlot(g3,mwt,uwis2,title)

def goVelocitiesAndPlot(g3,mwt,uwis2,title):
  cv = [2.75,5.5]
  # Make 2D seismic image with profile through wells
  maxz = mwt.getMaxz()
  sd = Sampling(ince(maxz/dz2),dz2,0)
  pack = getWellTies(mwt,uwis2)
  f1,h1,sf1,sh1,u1,v0,v1,tz0,tz1,sz,pv,x2,x3 = pack
  lx2,lx3,s2n,x2n = getCurveThruWells(x2,x3,s2.delta,edg=0.1)
  g2 = getImageAlongCurve(g3,lx2,lx3)
  # Make velocities
  bgs = 1.0; p0 = 0.0; p1 = 1.5; p2 = 1.5;
  vi = mwt.interpolateInitialVelocity3(bgs,p0,p1,p2)
  vn = mwt.interpolateInitialVelocityTied3(bgs,p0,p1,p2)
  #va = mwt.interpolateVavg(bgs,p0,p1,p2,sd,None)
  vi2 = getImageAlongCurve(vi,lx2,lx3)
  vn2 = getImageAlongCurve(vn,lx2,lx3)
  #va2 = getImageAlongCurve(va,lx2,lx3)
  vi2l = makeVot(sz,v0,tz0,sf1)
  vn2l = makeVot(sz,v0,tz1,sh1)
  makeVelocityPlots(s1,s2n,x2n,sf1,sh1,g2,vi2,vn2,vi2l,vn2l,title,cv)
  #plotSlice(s2n,s1,g2,v=va2,title='test'+title,clips2=cv,paper='test')

def makeVelocityPlots(s1,s2n,x2n,sf1,sh1,g2,vi2,vn2,vi2l,vn2l,title,cv):
  plot2DLogs(s1,s2n,g2,sf1,vi2l,x2n,title=title+'_linit',cbar='Velocity (km/s)',\
  								lcmap=jet,lclips=cv,paper='linit_'+title,velc=True)
  plot2DLogs(s1,s2n,g2,sh1,vn2l,x2n,title=title+'',cbar='Velocity (km/s)',\
  								lcmap=jet,lclips=cv,paper='ltied_'+title,velc=True)
  plotSlice(s2n,s1,g2,v=vi2,title='init_'+title,clips2=cv,paper='init_'+title)
  plotSlice(s2n,s1,g2,v=vn2,title='tied_'+title,clips2=cv,paper='tied_'+title)


def goPlotSynthetics(g3,f3,h3,mwt,uwis2,title,w2=None):
  pack = getWellTies(mwt,uwis2)
  f1,h1,sf1,sh1,u1,v0,v1,tz0,tz1,sz,pv,x2,x3 = pack
  lx2,lx3,s2n,x2n = getCurveThruWells(x2,x3,s2.delta,edg=0.1)
  g2 = getImageAlongCurve(g3,lx2,lx3)
  f2 = getImageAlongCurve(f3,lx2,lx3)
  h2 = getImageAlongCurve(h3,lx2,lx3)
  plotSlice(s2n,s1,g2,title='seismic_'+title,slides='seismic_'+title)
  plotSlice(s2n,s1,h2,title='wsi_'+title,paper='wsi_'+title)
  plotSlice(s2n,s1,f2,title='si_'+title,paper='si_'+title)
  plot2DLogs(s1,s2n,g2,sf1,f1,x2n,title=title+'_li',paper='li_'+title)
  plot2DLogs(s1,s2n,g2,sh1,h1,x2n,title=title+'_lt',paper='lt_'+title)
  plot2DLogs(s1,s2n,g2,sf1,f1,x2n,title=title+'_li',slides='lio_'+title,rmb=True)
  plot2DLogs(s1,s2n,g2,sh1,h1,x2n,title=title+'_lt',slides='lto_'+title,rmb=True)
  plotMultipleWellTiesIndiv(mwt,h3)
  if w2:
    tops,fmsp,datums = getWellTops(w2)
    gh2 = get2DHorizons(lx2,lx3)
    plot2DLogs(s1,s2n,g2,sh1,h1,x2n,gh=gh2,tops=tops,fmsp=fmsp,tz=tz1,sz=sz,slides='synt_fmshzs')
    plot2DLogs(s1,s2n,g2,sf1,f1,x2n,gh=gh2,tops=tops,fmsp=fmsp,tz=tz0,sz=sz,slides='syni_fmshzs')
  #goTopAnalysis(tops,sf1,f1,tz0,sz,gh2) 

def goTopAnalysis(tops,sf1,f1,tz0,sz,gh2):
  fms = ["DKOT","LKOT","CRMT","RDPK","A Sand","C1 Dolo"]
  ght,ghy,csty = gh2
  

def plotSingleWellTies(mwt):
  f1,h1,sf1,sh1,u1,v0,v1,tz0,tz1,sz,pv,x2,x3 = getWellTies(mwt)
  g3 = mwt.getSeismic3()
  nw = len(sf1)
  for i in range(nw):
    ix2 = inro((x2[i]-s2.first)/s2.delta)
    ix3 = inro((x3[i]-s3.first)/s3.delta)
    plotCurve(g3[ix3][ix2],s1,h1[i],sh1[i])
    plotVelPanel(sz[i],tz0[i],v0[i],pv[i],sy12=tz1[i],sy22=v1[i],hlim=[-20,20])

def plotMultipleWellTiesIndiv(mwt,h3):
  f1,h1,sf1,sh1,u1,v0,v1,tz0,tz1,sz,pv,x2,x3 = getWellTies(mwt)
  g3 = mwt.getSeismic3()
  nw = len(sf1)
  for i in range(nw):
    ix2 = inro((x2[i]-s2.first)/s2.delta)
    ix3 = inro((x3[i]-s3.first)/s3.delta)
    plotCurve(g3[ix3][ix2],s1,h3[ix3][ix2],s1)
    #plotVelPanel(sz[i],tz0[i],v0[i],pv[i],sy12=tz1[i],sy22=v1[i],hlim=[-20,20])


def extrapolateTest(mwt,sd):
  va1,va2,sza,sda,tz0,tz1,tz2,lea,vle = [],[],[],[],[],[],[],[],[]
  #vep,ad = computePredErrors(mwt,sd)
  for well in mwt.wellties:
    mwt.tzeXOn()
    le = mwt.extrapolateLinearX(well.sz,sd,well.tz1)
    #plotCurve(le,sd,well.tz1,well.sz)
    vf = div(mul(2,floats(well.sz.getValues())),well.tz1)
    mwt.tzeXOff()
    ve = mwt.extrapolateLinearX(well.sz,sd,vf)
    #lpf = LocalPredictionFilter(200,0,16)
    #vft = lpf.getTrend(vf)
    #vf0 = sub(vf,vft)
    #a = lpf.calculateCoeff(vft)
    #e = lpf.predictionError(vft,a)
    #p = sub(vft,e)
    #plotCurve(vf,well.sz,title="vavg")
    #plotCurve(vft,well.sz,title="trend")
    #plotCurve(vf0,well.sz,title="detrend")
    #plotCurve(a,well.sz,title="coeff")
    #plotCurve(e,well.sz,title="pred error")
    #plotCurve(p,well.sz,title="pred vavg")
    #plotCurve(vft,well.sz,p,well.sz,title="detrend with pred error")
    vle1 = zerofloat(sd.count)
    for i in range(1,len(vle1)):
      vle1[i] = 2*sd.getValue(i)/le[i]
    vle1[0] = vle1[1]
    lea.append(le)
    vle.append(vle1)#div(mul(2,floats(sd.getValues())),le))
    sza.append(well.sz)
    sda.append(sd)
    va1.append(vf)
    va2.append(ve)
    tz0.append(well.tz0)
    tz1.append(well.tz1)
    tz2.append(div(mul(2,floats(sd.getValues())),ve))
    #plotCurve(ve,sd,vf,well.sz)
  plotAllVavgCurves(sza,va1,flz=0,slides="vf")
  plotAllVavgCurves(sda,va2,flz=0,slides="ve")
  #plotAllVavgCurves(sda,vep,flz=0,slides="vep")
  plotAllTDCurves(sza,tz0,flz=0,slides="tz0")
  plotAllTDCurves(sza,tz1,flz=0,slides="tz1")
  plotAllTDCurves(sda,tz2,flz=0,slides="tz1e")
  plotAllTDCurves(sda,lea,flz=0,slides="tz1el")
  plotAllVavgCurves(sda,vle,flz=0,slides="vel")
  ex = 3
  sda = [sda[ex]]
  sza = [sza[ex]]
  tz0 = [tz1[ex]]
  tz1 = [tz1[ex]]
  tz2 = [tz2[ex]]
  lea = [lea[ex]]
  vle = [vle[ex]]
  ex = str(ex)
  plotAllTDCurves(sza,tz0,flz=0,slides="extz0_"+ex)
  plotAllTDCurves(sza,tz1,flz=0,slides="extz1_"+ex)
  plotAllTDCurves(sda,tz2,flz=0,slides="extz2_"+ex)
  plotAllTDCurves(sda,lea,flz=0,slides="exlea_"+ex)
  plotAllVavgCurves(sda,vle,flz=0,slides="exvle_"+ex)

def computePredErrors(mwt,sd):
  vea = []
  sig = 100
  nz = sd.count
  ad = zerofloat(nz)
  ix = fillfloat(1,nz)
  # compute the coefficients
  lpf = LocalPredictionFilter(sig,2,16)
  for well in mwt.wellties:
    tzr = mwt.resampleLinear(well.sz,sd,well.tz1)
    sz2 = Sampling(len(tzr),sd.delta,well.sz.first)
    ifx = sd.indexOfNearest(well.sz.first)
    nzl = len(tzr)
    vf = div(mul(2,floats(sz2.getValues())),tzr)
    vft = lpf.getTrend(vf)
    vf0 = sub(vf,vft)
    a = lpf.calculateCoeff(vft)
    for i in range(nzl):
      ad[i+ifx] += a[i]
      ix[i+ifx] += 1
  # average the coefficients
  for i in range(nzl):
      if ix[i]>0: 
        ad[i] /= ix[i]
  # apply the coefficients
  for well in mwt.wellties:
    tzr = mwt.resampleLinear(well.sz,sd,well.tz1)
    nzl = len(tzr)
    ifx = sd.indexOfNearest(well.sz.first)
    lfx = sd.indexOfNearest(well.sz.last)
    sz2 = Sampling(len(tzr),sd.delta,well.sz.first)
    vf = div(mul(2,floats(sz2.getValues())),tzr)
    vft = lpf.getTrend(vf)
    vf0 = vf #sub(vf,vft)
    ve = zerofloat(nz)
    ne = nz-ifx-nzl
    for i in range(nzl):
      ve[i+ifx] = vf0[i]*ad[i+ifx]
    for i in range(ne):
      for k in range(nzl):
        ve[i+lfx] += ad[i+lfx]*ve[ifx+nzl+i-k]
    for i in range(ifx):
      for k in range(nzl):
        ve[ifx-i] += ad[ifx-i]*ve[ifx-i+k]
    vea.append(ve)
    SimplePlot.asPoints(ve)
  return vea,ad


def getLogsAlongCurve(v,x2,s2n,sh):
  nw = len(x2)
  f1 = s1.first
  d1 = s1.delta
  f2 = s2n.first
  d2 = s2n.delta
  l = []
  for i in range(nw):
    ix2 = inro((x2[i]-f2)/d2)
    nx1 = sh[i].count
    jx1 = inro((sh[i].first-f1)/d1)
    l.append(copy(nx1,jx1,v[ix2]))
  return l

def getWellTies(mwt,uwis=None):
  f1,h1,sf1,sh1,u1,v0,v1,tz0,tz1,sz,pv,x2,x3 = \
    [],[],[],[],[],[],[],[],[],[],[],[],[]
  i=0
  for well in mwt.wellties:
    if uwis and not uwis[i]==well.id: continue
    f1.append(well.f)
    h1.append(well.h)
    sf1.append(well.sf)
    sh1.append(well.sh)
    u1.append(well.u)
    v0.append(well.v)
    v1.append(well.v1)
    tz0.append(well.tz0)
    tz1.append(well.tz1)
    sz.append(well.sz)
    pv.append(well.pv)
    x2.append(well.x2f)
    x3.append(well.x3f)
    i+=1
    if (i==len(uwis)): break
  return f1,h1,sf1,sh1,u1,v0,v1,tz0,tz1,sz,pv,x2,x3

def goSingleTie(g,uwi,well):
  wid = str(well)
  wt = WellTie(uwi,dz,s1,wset)
  wt.makePropagatorSeismogram(propDir,fp,q,phase)
  #wt.makePropagatorSeismogram(propDir,fp,q,0)
  x3  = wt.x3f; x2  = wt.x2f
  gi,ix2,ix3 = getTraces(g,x2,x3,2)
  #wt.findPhase()
  wt.computeSingleWellTie(gi,r1min,r1max,smin,smax,dr1)
  tz0 = wt.tz0 
  f   = wt.f 
  sf  = wt.sf
  sz  = wt.sz
  tz1 = wt.tz1
  h   = wt.h 
  sh  = wt.sh
  v0  = wt.v
  v1  = wt.v1
  pv  = wt.pv
  u   = wt.u
  php = wt.getPhasePlot()
  gi0 = gi[0]
  sft = Sampling(sf.count,sf.delta,sh.first)
  #plotPhaseRMS(php,len(php),1,paper='phaseplotswt'+wid)
  lim = [s1.first,s1.last]
  plotCurve(gi0,s1,slides="swt_seis"+wid)
  plotCurve(f,sf,slides="swt_untie"+wid,vlim=lim)
  plotCurve(f,sft,slides="swt_shift"+wid,vlim=lim)
  plotCurve(h,sh,slides="swt_tie"+wid,vlim=lim)
  plotCurve(gi0,s1,f,sf,slides="swt_untie_seis_"+wid)
  plotCurve(gi0,s1,h,sh,slides="swt_tie_seis_"+wid)
  plotVelPanel(sz,tz0,v0,pv,sy12=tz1,sy22=v1,hlim=[-15,15],slides="tzp"+wid)
  plotLogPanel(sz,v0,wt.d,wt.r,slides="logpanel"+wid)
  #pack = [f,sf,h,sh,u,sz,tz0,tz1,v0,v1,pv,gi,ix2,ix3]
  #return pack



#########################
# For bug finding only: #
#########################
def goOldShifts(g):
  uwi,set = getWellFromSet(well) 
  ai,x2,x3,z0,v,d,sz = getLogs(uwi,set)
  sf,f,tz0 = getPropagatorSeismogram(sz,v,d,q,phase=phase,normalize=True)
  gi,ix2,ix3 = getTraces(g,x2,x3)
  u,h,sh = getWarpingShiftsOld(sf,s1,f,gi,r1min,r1max,dr1)
  tz1,v1,pv = updateTDandVelocity(sz,u,tz0,sf)
  gi0 = gi[0]
  plotCurve(gi0,s1,f,sf,title='old tie')
  plotCurve(gi0,s1,h,sh,title='old tie')
  plotCurve(tz0,sz,tz1,sz,title='old tie')
  return u,sf,tz1,sz


#------------------------------------------------------------------------------#

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

