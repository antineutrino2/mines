"""
Jython utilities for BP 2D TTI dataset
Author: Andrew Munoz, Colorado School of Mines
Version: 2014.01.30
"""
from imports import *

#############################################################################
# Internal constants

_home = "/Users/amunoz/Home/data/bp/"
_sgydir = _home+"sgy/"
_datdir = _home+"dat/"

#############################################################################
# get data

def getModel(which):
  """
  Reads in a model with the correct samplings 
    samplings s1,s2,s3
  Example: getModel("vp")
  """
  global s1m,s2m
  n1,n2 = 1801,12596
  d1,d2 = 0.00625,0.00625 # km
  f1,f2 = 0.0,0.0 # km
  s1m,s2m = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  modeldir = _datdir+"model/"
  if which=="vp":
    g = readBpImage(modeldir+"vp")
  elif which=="ep":
    g = readBpImage(modeldir+"epsilon")
  elif which=="de":
    g = readBpImage(modeldir+"delta")
  elif which=="th":
    g = readBpImage(modeldir+"theta")
  else:
    print "unrecognized subset:",name
    System.exit
  return g

def getModelSamplings():
  return s1m,s2m

def getStack():
  """
  Reads in a seismic stack with the correct samplings 
    samplings s1,s2,s3
  Example: getModel("vp")
  """
  global s1s,s2s
  n1,n2 = 1151,12596
  d1,d2 = 0.008,0.00625 # s, km
  f1,f2 = 0.0,0.0 # s, km
  s1s,s2s = Sampling(n1,d1,f1),Sampling(n2,d2,f2)
  isodir = _datdir+"iso/"
  g = readBpImage(isodir+"bp_pstm_iso_stk")
  return g

def getStackSamplings():
  return s1s,s2s

#############################################################################
# read/write files

def readBpImage(fname):
  """ 
  Reads an image from a file with specified name.
  name: base name of image file; e.g., "tpsz"
  """
  fileName = fname+".dat"
  n1,n2 = s1.count,s2.count
  image = zerofloat(n1,n2)
  ais = ArrayInputStream(fileName)
  ais.readFloats(image)
  ais.close()
  return image

def writeBpImage(fname,image):
  """ 
  Writes an image to a file with specified name.
  name: base name of image file; e.g., "tpgp"
  image: the image
  """
  fileName = fname+".dat"
  aos = ArrayOutputStream(fileName)
  aos.writeFloats(image)
  aos.close()
  return image

#############################################################################
# graphics





#############################################################################
# other

def addNoise(s,x):
  """ adds s percent random noise to reduce problems with zero traces """
  xdif = max(x)-min(x)
  n1,n2,n3 = s1.count,s2.count,s3.count
  r = randfloat(n1,n2,n3)
  x = add(mul(sub(r,0.5),s*xdif),x)
  return x

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
