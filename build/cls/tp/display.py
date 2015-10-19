"""
Displays Teapot Dome data.
"""
from tputils import *

#setupForSubset("subt_251_4_500")
#setupForSubset("subz_401_4_600")

setupForSubset("subt_1501_2_0")
#setupForSubset("subt_3002_1_0")
#setupForSubset("subz_2762_2_0")

def main(args):
  #displaySlices()
  displaySubset()

def displaySubset():
  world = World()
  x = readImage("tpst")
  addImageToWorld(world,x)
  #addAllHorizonsToWorld(world)
  addHorizonToWorld(world,"KF2F2WC")
  addHorizonToWorld(world,"FallRiverDKOT")
  addHorizonToWorld(world,"FallRiverDKOT")
  addHorizonToWorld(world,"CrowMountainCRMT")
  addHorizonToWorld(world,"TensleepASand")
  addHorizonToWorld(world,"TensleepBbaseC1Dolo")
  addLogsToWorld(world,"d","vd")#,"g")
  makeFrame(world)

def displaySlices():
  displaySlice("tpsz",ColorMap.GRAY)
  displaySlice("tpgv",ColorMap.JET)
  displaySlice("tpgd",ColorMap.JET)
  displaySlice("tpgp",ColorMap.JET)
  displaySlice("tpgg",ColorMap.JET)

def displaySlice(name,cmap):
  x = readSlice3("s3_84/"+name)
  sp = SimplePlot(SimplePlot.Origin.UPPER_LEFT)
  sp.addColorBar()
  pv = sp.addPixels(x)
  pv.setInterpolation(PixelsView.Interpolation.NEAREST)
  if cmap!=None:
    pv.setColorModel(cmap)

#############################################################################
run(main)
