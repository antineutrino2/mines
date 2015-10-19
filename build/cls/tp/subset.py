"""
Makes directories with subsets of seismic time or depth images.
For example, the directory with name "subz_401_4_600" contains a
subset of a seismic depth image with 400 samples, 4 m sampling 
interval, and the depth of the first sample is 600 m.
Author: Dave Hale, Colorado School of Mines
Version: 2009.07.25
"""
from imports import *

#############################################################################
def main(args):
  #makeSubset("st",st,Sampling(876,0.002,0.250))
  #makeSubset("st",st,Sampling(1001,0.002,0.000))
  #makeSubset("st",st,Sampling(1501,0.002,0.000))
  #makeSubset("sz",sz,Sampling(2762,0.002,0.000))
  #makeSubset("sz",sz,Sampling(401,0.004,0.400))
  #makeSubset("sz",sz,Sampling(401,0.004,0.600))
	## Demo Volume
	makeSubset23("st",st,Sampling(301,0.004,0.250),\
										s2,Sampling(101,0.025,2.500),\
										s3,Sampling(101,0.025,1.000))

csmDir = "../../../../../data/tp/csm/"
seismictDir = csmDir+"seismict/"
seismiczDir = csmDir+"seismicz/"
st = Sampling(1501,0.002,0.000) # time sampling of input
sz = Sampling(2762,0.002,0.000) # depth sampling of input
s2 = Sampling(357,0.025,0.000) # x2 sampling (both input and output)
s3 = Sampling(161,0.025,0.000) # x3 sampling (both input and output)

def makeSubset(what,s1i,s1o):
  n1i,d1i,f1i = s1i.count,s1i.delta,s1i.first
  n1o,d1o,f1o = s1o.count,s1o.delta,s1o.first
  n1s,d1s,f1s = str(n1o),str(round(d1o*1000)),str(round(f1o*1000))
  k1i = round(d1o/d1i)
  j1i = round((f1o-f1i)/d1i)
  tz = what[1]
  fileName = "tp"+what+".dat"
  inDir = csmDir+"seismic"+tz+"/"
  outDir = inDir+"sub"+tz+"_"+n1s+"_"+d1s+"_"+f1s+"/"
  #outDir = inDir+"sub"+tz+"_3002_1_"+f1s+"/" # You must give the file a specific name
  inFile = inDir+fileName
  outFile = outDir+fileName
  File(outDir).mkdir()
  xi = zerofloat(n1i)
  xo = zerofloat(n1o)
  ais = ArrayInputStream(inFile)
  aos = ArrayOutputStream(outFile)
  for i3 in range(s3.count):
    for i2 in range(s2.count):
      ais.readFloats(xi)
      copy(n1o,j1i,k1i,xi,0,1,xo)
      #xo = oversample(xo,s1i)
      aos.writeFloats(xo)
  ais.close()
  aos.close()
  makeImageFrame(outFile,s1o,s2,s3)

def makeSubset23(what,s1i,s1o,s2i,s2o,s3i,s3o):
  n1i,d1i,f1i = s1i.count,s1i.delta,s1i.first
  n1o,d1o,f1o = s1o.count,s1o.delta,s1o.first
  n2i,d2i,f2i = s2i.count,s2i.delta,s2i.first
  n2o,d2o,f2o = s2o.count,s2o.delta,s2o.first
  n3i,d3i,f3i = s3i.count,s3i.delta,s3i.first
  n3o,d3o,f3o = s3o.count,s3o.delta,s3o.first
  n1s,d1s,f1s = str(n1o),str(round(d1o*1000)),str(round(f1o*1000))
  n2s,d2s,f2s = str(n2o),str(round(d2o*1000)),str(round(f2o*1000))
  n3s,d3s,f3s = str(n3o),str(round(d3o*1000)),str(round(f3o*1000))
  k1i = round(d1o/d1i)
  k2i = round(d2o/d2i)
  k3i = round(d3o/d3i)
  j1i = round((f1o-f1i)/d1i)
  j2i = round((f2o-f2i)/d2i)
  j3i = round((f3o-f3i)/d3i)
  tz = what[1]
  fileName = "tp"+what+".dat"
  inDir = csmDir+"seismic"+tz+"/"
  outDir = inDir+"sub"+tz+"_"+n1s+"_"+d1s+"_"+f1s+\
													"_"+n2s+"_"+d2s+"_"+f2s+\
													"_"+n3s+"_"+d3s+"_"+f3s+"/"
  #outDir = inDir+"sub"+tz+"_3002_1_"+f1s+"/" # You must give the file a specific name
  inFile = inDir+fileName
  outFile = outDir+fileName
  File(outDir).mkdir()
  xi = zerofloat(n1i,n2i,n3i)
  xo = zerofloat(n1o,n2o,n3o)
  ais = ArrayInputStream(inFile)
  aos = ArrayOutputStream(outFile)
  ais.readFloats(xi)
  copy(n1o,n2o,n3o,j1i,j2i,j3i,k1i,k2i,k3i,xi,0,0,0,1,1,1,xo)
  aos.writeFloats(xo)
  ais.close()
  aos.close()
  makeImageFrame(outFile,s1o,s2o,s3o)

def makeImageFrame(file,s1,s2,s3):
  n1,n2,n3 = s1.count,s2.count,s3.count
  d1,d2,d3 = s1.delta,s2.delta,s3.delta
  x = zerofloat(n1,n2,n3)
  ais = ArrayInputStream(file)
  ais.readFloats(x)
  ais.close()
  frame = SimpleFrame.asImagePanels(s1,s2,s3,x)
  frame.setSize(1200,900)
  view = frame.getOrbitView()
  view.setAxesScale(1.0,1.0,(n2*d2)/(n1*d1))
  view.setScale(2.0)
  view.setAzimuth(-50.0)
  return frame

# Used for sub 2 ms data
def oversample(x1,xs1):
  n1 = len(x1)
  n2 = 3002
  d2 = 0.001
  f2 = 0.0
  x2 = zerofloat(n2)
  si = SincInterpolator()
  si.setUniform(n1,xs1.getDelta(),xs1.getFirst(),x1)
  si.interpolate(n2,d2,f2,x2)
  return x2

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
if __name__=="__main__":
  SwingUtilities.invokeLater(RunMain()) 
