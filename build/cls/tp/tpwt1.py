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
# Image warping
dr1   = 0.02 # 2D smoothness vertically 
dr2   = 0.50 # 2D smoothness vertically 
dr3   = 0.50 # 2D smoothness vertically 
r1min,r1max = -0.10,0.10 # vertical constraint # use -0.06,0.12,0.02,0.5,-1.0,1.0
r2min,r2max = -0.50,0.50 # vertical constraint # use -0.06,0.12,0.02,0.5,-1.0,1.0
r3min,r3max = -0.50,0.50 # vertical constraint # use -0.06,0.12,0.02,0.5,-1.0,1.0
smin,smax = -0.100,0.100 # min and max time shifts for 2D/3D warping

propDir = getPropDir()
wset = getLogDataset("d")
phase = 44


def main(args):
  global s1,s2,s3
  cut1=[0.10,1.25];cut2=[3.40,6.65];cut3=[0.6,2.7];
  g,s3,s2,s1,fo = getImage(normalize=True,cut1=cut1,cut2=cut2,cut3=cut3,rnrm=True)
  #wells = [3,12,16,15,7,4,2,11]; nm="7" # good wells for 3D
  well = 15
  uwi = getIDFromDeepSet(well)
  goSingleTie(g,uwi,well)


def goSingleTie(g,uwi,well):
  wid = str(well)
  wt = WellTie(uwi,dz,s1,wset)
  wt.makePropagatorSeismogram(propDir,fp,q,phase)
  x3  = wt.xf; x2  = wt.yf
  gi,ix2,ix3 = getTraces(g,x2,x3)
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
  plotCurve(gi0,s1,f,sf,title='untie')
  plotCurve(gi0,s1,h,sh,title='tied')
  plotCurve(tz0,sz,tz1,sz,title='td-curves')
  #plotPhaseRMS(php,len(php),1,paper='phaseplotswt'+wid)
  pack = [f,sf,h,sh,u,sz,tz0,tz1,v0,v1,pv,gi,ix2,ix3]
  return pack



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


