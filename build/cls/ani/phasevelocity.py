from imports import *

ac = FLT_PI/180.0
fa = 0.0
da = 1.0
na = 90
sa = Sampling(na,da,fa)

def main(args):
	#q1()
	q2()

def q1():
	goParams1()
	vpe = goExact()
	vpw = goWeak()
	title = "exact vs. weak (dashed) delta="+str(dlt)
	plot(vpe,vpw,title,png="q1-dlt="+str(dlt))

def q2():
	tim = 2
	goParams2(tim)
	fv = 1.4
	lv = 4.0
	dv = 0.001
	nv = int(round((lv-fv)/dv))+1
	sv = Sampling(nv,dv,fv)
	vpe = zerofloat(na,nv)
	for iv in range(nv):
		vpvs = fv+dv*iv
		f = 1.0-(1.0/(vpvs*vpvs))
		vpe[iv] = goExact(f)
	title = "exact: delta="+str(dlt)+" eps="+str(eps)
	plot2(sa,sv,vpe,title=title,png="q2-tim="+str(tim)+"-vpvs="+str(vpvs))
	
def goExact(fi=None):
	if fi: fe = fi
	else: fe = f
	vp = zerofloat(na)
	for ia in range(na):
		t = ac*(ia*da+fa)
		sin2 = sin(t)*sin(t)
		si2n = sin(2.0*t)*sin(2.0*t)
		vp[ia] = sqrt(1.0+eps*sin2-fe/2.0+(fe/2.0)*\
			sqrt(pow(1.0+2.0*eps*sin2/fe,2.0)-2.0*(eps-dlt)*si2n/fe))
	return vp
		
def goWeak(fi=None):
	if fi: fe = fi
	else: fe = f
	vp = zerofloat(na)
	for ia in range(na):
		t = ac*(ia*da+fa)
		sin2 = sin(t)*sin(t)
		sin4 = sin2*sin2
		cos2 = cos(t)*cos(t)
		vp[ia] = sqrt(1.0+2.0*dlt*sin2*cos2+2.0*eps*sin4+(4.0/fe)*(eps-dlt)*\
							sin4*cos2*(eps*sin2+dlt*cos2))
	return vp


def goParams1():
	global vpvs,eps,dlt,f,da,na
	global fa,da
	global na
	vpvs = 2.0
	f = 1.0-(1.0/(vpvs*vpvs))
	eps  = 0.6
	dlt = 0.1
	#dlt = 0.3
	#dlt = 0.5

def goParams2(tim):
	global vpvs,eps,dlt,f,da,na
	global fa,da
	global na
	if tim==1:
		eps  = 0.6
		dlt = 0.1
	if tim==2:
		dlt = 0.5
		eps  = 0.6


def plot(x1,x2=None,title=None,png=None):
	pp = PlotPanel(PlotPanel.Orientation.X1RIGHT_X2UP)
	l1 = pp.addPoints(sa,x1)
	if x2:
		l2 = pp.addPoints(sa,x2)
		l2.setStyle('r--')
	if title:
		pp.setTitle(title)
	pp.setVLabel("velocity")
	pp.setHLabel("angle (degrees)")
	frame = PlotFrame(pp)
	frame.setFontSize(18)
	frame.setSize(900,500)
	frame.setVisible(True)
	if png:
		frame.paintToPng(720,3.33,'png/'+png+'.png')

def plot2(s1,s2,x1,title=None,png=None):
	pp = PlotPanel(PlotPanel.Orientation.X1RIGHT_X2UP)
	p1 = pp.addPixels(s1,s2,x1)
	p1.setColorModel(ColorMap.PRISM)
	if title:
		pp.setTitle(title)
	pp.setVLimits(s2.first,s2.last)
	pp.addColorBar("phase velocity")
	pp.setVLabel("vp/vs")
	pp.setHLabel("angle (degrees)")
	frame = PlotFrame(pp)
	frame.setFontSize(18)
	frame.setSize(900,500)
	frame.setVisible(True)
	if png:
		frame.paintToPng(720,3.33,'png/'+png+'.png')


	

#############################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
if __name__=="__main__":
  SwingUtilities.invokeLater(RunMain()) 
