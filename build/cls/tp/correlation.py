#test




from imports import *

def main(args):
	x,y = makeData()
	mx = mean(x)
	my = mean(y)
	#r = num(x,mx,y,my)/den(x,mx,y,my)
	r = correlation(x,y)
	print 'r = ',r
	plot(x,y)

def makeData():
	sc = 0.2
	n = 30
	on = 1.0/n
	x = zerofloat(n)
	y = zerofloat(n)
	r = RandomFloat(43294)
	for i in range(n):
		x[i] = on*i+sc*r.normal()
		y[i] = on*i+sc*r.normal()
	return x,y

def mean(a):
	return sum(a)/len(a)

def num(a,ma,b,mb):
	return sum(mul(sub(a,ma),sub(b,mb)))

def den(a,ma,b,mb):
	return sqrt(sum(pow(sub(a,ma),2.0))*sum(pow(sub(b,mb),2.0)))

def correlation(x,y):
	mx = mean(x)
	my = mean(y)
	sx = sub(x,mx)
	sy = sub(y,my)
	sxsx = mul(sx,sx)
	sysy = mul(sy,sy)
	r = sum(mul(sx,sy))/sqrt(sum(sxsx)*sum(sysy))
	return r



def plot(x,y):
	si = SimplePlot()
	p = si.addPoints(x,y)
	p.setStyle('kO')
	p.setMarkSize(5)
	

#------------------------------------------------------------------------------#

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
