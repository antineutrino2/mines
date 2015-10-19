from tputils import *
from wtutils import *
from imports import *

_tpDir = "/Users/amunoz/Home/data/tp/csm/"
seismicDir = getSeismictDir()
horzDir = _tpDir+"horz/"


def main(args):
	h1 = 'test'
	global n1,n2,n3
	global s1,s2,s3
	u,s1,s2,s3 = readTpstData()
	n1 = s1.count; n2 = s2.count; n3 = s3.count
	#u1,u2,u3 = makeVectors(u)	
	u1,u2,u3 = readVectors()
	pickHorizon(u,u1,u2,u3)

def pickHorizon(u,u1,u2,u3):
	rp = zerofloat(3,2)
	rp[0][0] = 50.0 ; 	rp[1][0] = 100.0 #x2  
	rp[0][1] = 50.0 ; 	rp[1][1] = 100.0 #x3
	rp[0][2] = 440.0; 	rp[1][2] = 440.0 #x1
	# slice
	yt = 50
	ipx = zerofloat(n1,n3)
	for i3 in range(n3):
		for i1 in range(n1):
			ipx[i3][i1] = u[i3][yt][i1]
	se = SurfaceExtractor()
	se.setWeights(50.0,0.0) #smoothing
	surf1 = se.surfaceInitialization(n2,n3,n1-1.0,rp)
	plotSurface(surf1)
	surf1 = se.surfaceUpdateFromSlopes(u1,u2,u3,surf1,rp,n1-1.0,ipx,yt);
	plotSurface(surf1)

def readImage(which):
  u = zerofloat(n1,n2,n3) 
  ais = ArrayInputStream(_tpDir+which+'.dat')
  ais.readFloats(u)
  ais.close()
  return u

def makeVectors(u):
  print 'making normal vectors...'
  u3 = zerofloat(n1,n2,n3)
  u2 = zerofloat(n1,n2,n3)
  u1 = zerofloat(n1,n2,n3)
  lof = LocalOrientFilter(5.0)
  lof.applyForNormal(u,u1,u2,u3)
  writeVectors(u1,u2,u3)
  return u1,u2,u3

def readVectors():
  print 'reading normal vectors...'
  u1 = zerofloat(n1,n2,n3)
  u2 = zerofloat(n1,n2,n3)
  u3 = zerofloat(n1,n2,n3)
  ais1 = ArrayInputStream(horzDir+"nvec/u1.dat")
  ais2 = ArrayInputStream(horzDir+"nvec/u2.dat")
  ais3 = ArrayInputStream(horzDir+"nvec/u3.dat")
  ais1.readFloats(u1)
  ais2.readFloats(u2)
  ais3.readFloats(u3)
  ais1.close()
  ais2.close()
  ais3.close()
  return u1,u2,u3

def writeVectors(u1,u2,u3):
  print 'writing normal vectors...'
  aos1 = ArrayOutputStream(horzDir+"nvec/u1.dat")
  aos2 = ArrayOutputStream(horzDir+"nvec/u2.dat")
  aos3 = ArrayOutputStream(horzDir+"nvec/u3.dat")
  aos1.writeFloats(u1)
  aos2.writeFloats(u2)
  aos3.writeFloats(u3)
  aos1.close()
  aos2.close()
  aos3.close()

def plotSurface(u):
	hs = SimpleFrame(AxesOrientation.XOUT_YRIGHT_ZUP)
	hs.asTriangles(True,s3,s2,u)

#------------------------------------------------------------------------------#

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())


