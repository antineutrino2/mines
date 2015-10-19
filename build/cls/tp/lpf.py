# Lab 7
# Local Prediction Filter
# @author Andrew Munoz
# @version 22.11.11

import sys
import time
import java.util.Random
from java.awt import *
from java.lang import *
from java.util import *
from java.nio import *
from javax.swing import *
from java.awt.Color import *
from edu.mines.jtk.awt import *
from edu.mines.jtk.io import *
from edu.mines.jtk.mosaic import *
from edu.mines.jtk.dsp import *
from edu.mines.jtk.util import *
from edu.mines.jtk.util import Stopwatch 
from edu.mines.jtk.util.ArrayMath import *
from edu.mines.jtk.util.RandomFloat import *

from tp import * 

#############################################################################

def main(args):
  pef1()
  #pef2()
  #goBenchmark()

"""
Prediction Error Filter 1D
"""
def pef1():
  sigma = 50 
  niter = 16
  #x = image1D() 
  x = randomSequence() 
  y = like1(x) 
  p = like1(x) 
  for d in [0,1,2]:
    lpf = LocalPredictionFilter(sigma,d,niter) # (sigma,dir,niter)- dir=0,1,2
    t = lpf.getTrend(x) 
    y = sub(x,t)
    a = lpf.calculateCoeff(y)
    e = lpf.predictionError(y,a)
    p = sub(y,e)
    ds = str(d)
    sig = str(sigma)
    plot1D(e,"error:"+ds)#,png='error1D-'+ds+'_sig='+sig)
    plot1D(x,"input:"+ds,p,"predicted:"+ds)#,png='predicted1D-'+ds+'_sig='+sig)
    plot1D(a,"coefficients:"+ds)#,png='coefficients1D-'+ds+'_sig='+sig)
    plot1D(t,"trend:"+ds)#,png='trend1D-'+ds+'_sig='+sig)
    plot1D(y,"x without trend:"+ds)#,png='xminustrend1D-'+ds+'_sig='+sig)
    plot1D(x,"input:"+ds)#,png='input-'+ds+'_sig='+sig)
    plot1D(x,"input:"+ds,y,"input without trend:"+ds,e,"error:"+ds )#,png='x_y_e-'+ds+'_sig='+sig)

"""
Prediction Error Filter 2D
"""
def pef2():
  sigma = 20
  niter = 16
  x = image2D()
  y = like2(x)
  for d in [0,1,2,3,4]:
    lpf = LocalPredictionFilter(sigma,d,niter) # (sigma,dir,niter)- dir=0,1,2,3,4
    t = lpf.getTrend(x)
    y = sub(x,t)
    obja = lpf.calculateCoeff(x)
    a1 = obja[0]
    a2 = obja[1]
    obje = lpf.predictionError(y,a1,a2)
    e1 = obje[0]
    e2 = obje[1]
    e3 = obje[2]
    ds = str(d)
    sig = str(sigma)
    if (d==4): 
      plot2D(e1,"error:"+ds+' sigma='+sig                     )#,png="2Derror-"+ds+'_sig='+sig)
    else:                                                     
      plot2D(e3,"error-all:"+ds+' sigma='+sig                 )#,png="2DerrorA-"+ds+'_sig='+sig)
      plot2D(e2,"error-horizontal:"+ds+' sigma='+sig          )#,png="2DerrorH-"+ds+'_sig='+sig)
      plot2D(e1,"error-vertical:"+  ds+' sigma='+sig          )#,png="2DerrorV-"+ds+'_sig='+sig)
    plot2D(a2,"coefficients horizontal(a2):"+ds+' sigma='+sig )#,png="2Da2-"+ds+'_sig='+sig)
    plot2D(a1,"coefficients vertical(a1):"+  ds+' sigma='+sig )#,png="2Da1-"+ds+'_sig='+sig)
    plot2D(t,"trend:"+ds+' sigma='+sig                        )#,png="2Dtrend-"+ds+'_sig='+sig)
    plot2D(y,"x without trend:"+ds+' sigma='+sig              )#,png="2Dxmtrend-"+ds+'_sig='+sig)
    plot2D(x,"input:"+ds+' sigma='+sig                        )#,png="2Dinput-"+ds+'_sig='+sig)
                                                                
"""                                                             
Benchmark                                                       
"""                                                             
def goBenchmark():
  sigma = 1 
  sw = Stopwatch()
  maxtime = 2.0
  x = image2D()
  n1, n2 = len(x[0]), len(x)
  def benchMethod(method,name,d=0,niter=0):
    nsmooth = 0
    sw.restart()
    while sw.time()<maxtime:
      method(x,sigma,d,niter)
      nsmooth += 1
    sw.stop()
    mflops = int(6.0e-6*n1*n2*nsmooth/sw.time())
    print name + ": mflops- ",mflops
  benchMethod(apply2,"2D LPF- d=0, niter=16",16)
  benchMethod(apply2,"2D LPF- d=4, niter=16",16)

  

#############################################################################

def like1(x):
  return zerofloat(len(x))

def like2(x):
  return zerofloat(len(x[0]),len(x))

def image2D():
  s1 = Sampling(1501,0.002,0.000)
  s2 = Sampling(357,0.025,0.000)
  n1,n2 = s1.count,s2.count
  image = zerofloat(n1,n2)
  ais = ArrayInputStream('./tpslice66.dat')
  ais.readFloats(image)
  ais.close()
  return image

def image1D():
  image = image2D()
  seq = zerofloat(len(image[0]))
  seq = image[150]
  return seq

def randomSequence():
  n = 1001
  x = zerofloat(n)
  for i in range(n):
    rn = Random().nextFloat();
    x[i] = rn
  return x

def randomImage1():
  n1,n2 = 1001,1002
  x = zerofloat(n1,n2)
  for j in range(n2):
    for i in range(n1):
      rn = Random().nextFloat();
      x[j][i] = rn
  return x

def mean(x):
  return ExponentialSmoother.mean(x)

def plot2D(x,title,png=None):
  sp = SimplePlot.asPixels(x)
  sp.setTitle(title)
  sp.addColorBar()
  sp.setSize(1000,1000)
  sp.plotPanel.setColorBarWidthMinimum(100)
  sp.setVLabel('Time (ms)')
  sp.setHLabel('Trace')
  if png:
    sp.paintToPng(360,3.33,"./pics/"+png+'.png') 
 

def plot1D(x1,title1,x2=None,title2=None,x3=None,title3=None,x4=None,title4=None,png=None):
  sp = SimplePlot()
  ln1 = PointsView(x1)
  ln1.setLineColor(BLACK)
  sp.add(ln1)
  sp.setTitle(title1)
  sp.setHLabel('Time (ms)')
  sp.setVLabel('Amplitude')
  if x2:
    ln2 = PointsView(x2)
    ln2.setLineColor(RED)
    sp.add(ln2)
    sp.setTitle(title1+"-black, "+title2+"-red")
  if x3:
    ln3 = PointsView(x3)
    ln3.setLineColor(BLUE)
    sp.add(ln3)
    sp.setTitle(title1+"-black, "+title2+"-red, "+title3+"-blue")
  if x4:
    ln4 = PointsView(x4)
    ln4.setLineColor(ORANGE)
    sp.add(ln4)
    sp.setTitle(title1+"-black, "+title2+"-red, "+title3+"-blue, "+title4+"-orange")
  sp.setSize(1000,500)
  if png:
    sp.paintToPng(360,3.33,"./pics/"+png+'.png') 

"""
For Benchmarking
"""
def apply1(x,sigma,d=None,niter=None):
  lpf = LocalPredictionFilter(sigma,d,niter)
  lpf.apply(x)

def apply2(x,sigma,d=None,niter=None):
  lpf = LocalPredictionFilter(sigma,d,niter)
  lpf.apply(x)

#############################################################################
# Do everything on Swing thread.

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

