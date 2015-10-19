"""
Test demo for dynamic warping for well ties
Author: Andrew Munoz, Colorado School of Mines
Date: 02/27/13
"""
import sys

from java.awt import *
from java.io import *
from java.lang import *
from java.util import *
from javax.swing import *
from java.awt.Color import *

from edu.mines.jtk.awt import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util.ArrayMath import *

# CHANGE AS NEEDED:
import DynamicWarpingWells

def main(args):
	n1 = 301 
	n2 = 1001
	f0 = 200 #true f0
	rmsNoise = 0.0
	maxStrain = 0.5
	g,f = makeSin(n2,n1,f0,maxStrain,rmsNoise,C=1.0) 
	
	goTest(n1,n2,f0,f,g)

def goTest(nf,ng,f0,f,g):
	dt = 0.001
	smin = 0
	smax = ng-nf
	strainlim = 1.0
	dw = DynamicWarpingWells(smin,smax)
	dw.setStrainMax(strainlim)
	u = dw.findShifts(f,g)
	h = dw.applyShifts(u,f)
	u = mul(u,dt)

	sh = Sampling(len(h),dt,u[0])	
	sf = Sampling(nf,dt,0.0)
	sg = Sampling(ng,dt,0.0) 

	plotSequences(g,f,sg,sf)
	plotSequences(g,h,sg,sh)

########################################################################
## utils

def makeSin(n1,n2,start,maxStrain,rmsNoise,C=None):
  si = SincInterp()
  si.setExtrapolation(SincInterp.Extrapolation.CONSTANT)
  seedN2 = 10000
  fpeak  = 0.125
  amp = maxStrain*(n2/(2*PI))
  g = zerofloat(n2)
  t = zerofloat(n2)
  f = sig1(n1,rmsNoise,fpeak)
  w = float((2*PI/n2))
  if C:
    w = float((2*PI/n2))*C
  for i in range (0,n2):
    t[i]=start+i-amp*sin(i*w)   
  si.interpolate(n1,1.0,0.0,f,n2,t,g)    
  g = addNoise(rmsNoise,fpeak,g,seedN2)
  return f,g

def sig1(n1,rmsNoise, fpeak):
  seed   = 23232
  seedN1 = 31415 
  f = makeRandomEvents(n1,seed) 
  f = addRickerWavelet(fpeak,f)
  f = mul(1.0/max(abs(f)),f)
  f = addNoise(rmsNoise,fpeak,f,seedN1)
  return f

def makeRandomEvents(n,seed=0):
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  return pow(mul(2.0,sub(randfloat(r,n),0.5)),15.0)

def addRickerWavelet(fpeak,f):
  n = len(f)
  ih = int(3.0/fpeak)
  nh = 1+2*ih
  h = zerofloat(nh)
  for jh in range(nh):
    h[jh] = ricker(fpeak,jh-ih)
  g = zerofloat(n)
  Conv.conv(nh,-ih,h,n,0,f,n,0,g)
  return g

def ricker(fpeak,time):
  x = PI*fpeak*time
  return (1.0-2.0*x*x)*exp(-x*x)

def addNoise(nrms,fpeak,f,seed=0):
  n = len(f)
  if seed!=0:
    r = Random(seed)
  else:
    r = Random()
  nrms *= max(abs(f))
  g = mul(2.0,sub(randfloat(r,n),0.5))
  g = addRickerWavelet(fpeak,g)
  #rgf = RecursiveGaussianFilter(3.0)
  #rgf.apply1(g,g)
  frms = sqrt(sum(mul(f,f))/n)
  grms = sqrt(sum(mul(g,g))/n)
  g = mul(g,nrms*frms/grms)
  return add(f,g)

def plotSequences(f,g,s1=None,s2=None):
  n1 = len(f)
  n2 = len(g)
  if s1 is None:
    s1 = Sampling(n1,1.0,0.0)
  if s2 is None:
    s2 = Sampling(n2,1.0,0.0)
  panel = PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT)
  fv = panel.addPoints(s1,f)
  gv = panel.addPoints(s2,g)
  gv.setLineColor(Color.RED)
  panel.setHLabel("f & g")
  panel.setVLabel("time (s)")
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(WindowConstants.EXIT_ON_CLOSE)
  frame.setFontSize(18)
  frame.setSize(300,800)
  frame.setVisible(True)


#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
