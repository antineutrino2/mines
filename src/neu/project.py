from imports import *
from neu import *

prints = True
slides = False
papers = True 
pngDir = "/Users/amunoz/Documents/Neuroscience/pngs/"
pfx = "msp"
#pfx = "ssp"
#pfx = ""

dt = 0.001
tmin = 0
tmax = 100
nt = round((tmax-tmin)/dt)
st = Sampling(nt,dt,tmin)

def main(args):
  #goCaConcentrationTest()
  computeNMDACa()
  #modelInmda()

################################################################################
# Use parameters from Hoonuraiah (2013) to replicate Shouval (2002)
# Goldman-Hodgkin-Katz convention for current

def computeNMDACa():
  minf = 0
  maxf = 25
  df = 1
  sf = Sampling((maxf-minf)/df,df,minf)
  # spike
  c0 = 100e-6 # initial inside concentration of ca
  tauc = 30
  vr = 30.77
  inj = 50.0
  nb = 1
  #v = hodgkinHuxely(st,inj)
  v = Utils.pyrimidalCellFiring(st,inj,nb)
  bv = bov(v)
  #bovm(v,vr)
  plot1(st,v,"Time (ms)","V (mV)",png="vot")
  plot1(v,bv,"V (mV)","B(V)",png="bov")
  ica = mul(bv,sub(vr,v))
  plot1(v,ica,"V (mV)","Y",png="Y")
  ca = rk4dca(st,c0,v,tauc,ica)
  plot1(st,v,"Time (ms)","Voltage (mV)",png="HH")
  plot1(st,ca,"Time (ms)","Ca Concentration",png="HH ca")
  plot1(ca,v,"Ca Concentration","Voltage (mV)",png="HH cav")
  w = goCaConcentration(ca,v)


"""
RK4 solver for dc/dt
"""
def rk4dc(st,c0,v,tau):
  nt = st.count
  dt = st.delta
  c = zerofloat(nt)
  c[0] = c0
  inmda = zerofloat(nt)
  for i in range(nt-1):
    t = st.getValue(i)
    v1 = dt*dcdt(t,c[i],Inmda(v[i],t,c[i]),tau) 
    v2 = dt*dcdt(t+dt/2.0,c[i]+v1/2.0,Inmda(v[i],t,c[i]),tau)
    v3 = dt*dcdt(t+dt/2.0,c[i]+v2/2.0,Inmda(v[i],t,c[i]),tau)
    v4 = dt*dcdt(t+dt,c[i]+v3,Inmda(v[i],t,c[i]),tau)
    c[i+1] = c[i] + (v1+2*v2+2*v3+v4)/6.0
    inmda[i] = Inmda(v[i],t,c[i])
    #print inmda[i]
  plot1(st,inmda,"Time (ms)","Current (A)")
  return div(c,max(c))

# Inmda at a certian time and voltage, and Ca at a certian time
# Honnuraiah (2013)
def dcdt(t,ca,inmda,tau):
  dpt = 0.1 #um depth of shell
  caf = 100 #nM steady state value of Ca
  F = 1.0364268820905e-8 # mAs
  Ac = 5e-4 #mol/C*cm Ca accumulation factor
  dcdt = -Ac*inmda + (caf-ca)/tau # Sheynikhovich (2013)
  #dcdt = -10000*inmda/(3.6*dpt*F) + (caf-ca)/tau
  return dcdt

# Honnuraiah (2013)
def Inmda(v,t,ca):
  # constants
  pca = 10.6 # relative permeability ratio
  pnmda = 1 # maximum permeability of NMDA receptor
  taud = 50 #ms decay time constant
  taut = 5  #ms rise time
  mgo = 2 #mM Mg concentration
  z = 2 # valence of Ca
  cao = 2 # outside concentration of Ca
  #ca = 100e-6 # initial inside concentration of ca
  R = 8.314 #J*mol-1*K-1 is the gas constant
  T = 308 #K is the average body temperature
  F = 96485 #C/mol is the Faraday constant
  ENa = +55 #mV equilibrium Na
  EK = -90 #mV equilibrium K
  # functions
  a = getA(taud,taut,t)
  sti = sot(a,taud,taut,t)
  MgBV = 1/(1+mgo*exp(-0.062*v)/3.57)
  zF = z*F
  fort = (zF/(R*T))/(1000.0*1000.0)
  inmda = pnmda*pca*sti*MgBV*zF*fort*(ca-cao*exp(-v*fort)/(1-exp(-v*fort)))
  return inmda

# s(t) Governs the kinetics of the NMDA current
def sot(a,taud,taut,t):
  return a*(exp(-t/taud)-exp(-t/taut))

# normalization constant for s(t) such that 0<=s(t)<=1
def getA(taud,taut,t):
  c = 1
  # if taut>taud, find minimum instead of maximum
  if taut>taud:
    tauc = taud
    taud = taut
    taut = tauc
    c = -1
  b = taud/taut
  tmax = taud*log(b)/(b-1)
  a = sot(1,taud,taut,tmax)
  return c*a

############################################
# Hodgkin-Huxely Solution with euler method

# Compute the post-synaptic action potential
def hodgkinHuxely(st,inj):
  dt = st.delta
  nt = st.count
  v = zerofloat(nt)
  # initial conditions
  Vm = -65.0
  n0 = 0.0003
  m0 = 0.0011
  h0 = 0.9998
  #n0 = 0.9
  #m0 = 0.9
  #h0 = 0.0
  v = eulerHH(st,Vm,n0,m0,h0,inj)
  return v

def eulerHH(st,Vm,n0,m0,h0,inj):
  dt = st.delta
  nt = st.count
  f = zerofloat(nt)
  f[0] = Vm
  for i in range(nt-1):
      t = st.getValue(i)
      f[i+1] = f[i] + dt*vprime(f[i],m0,n0,h0,inj)
      dn1 = nprime(t,n0,-f[i])
      dm1 = mprime(t,m0,-f[i])
      dh1 = hprime(t,h0,-f[i])
      n1 = n0 + dt*dn1
      m1 = m0 + dt*dm1
      h1 = h0 + dt*dh1
      n0 = n1
      m0 = m1
      h0 = h1
  return f

# FROM HH Lab 
# Activation of Na+ channels
def mprime(t,m,V):
  # Calculates the derivative of m and a point t given Membrane potential V
  k1p = 0.1*(V+25)/(exp((V+25)/10)-1)
  k1m = 4*exp(V/18)
  p = (1-m)*k1p-m*k1m
  return p
# Activation of K+ channels
def nprime(t,n,V):
  # Calculates the derivative of n and a point t given Membrane potential V
  k1p = 0.01*(V+10)/(exp((V+10)/10)-1)
  k1m = 0.125*exp(V/80)
  p = (1-n)*k1p-n*k1m
  return p
# Deactivation of Na+ channels
def hprime(t,h,V):
  # Calculates the derivative of n and a point t given Membrane potential V
  k1p = 0.07*exp(V/20)
  k1m = 1/(exp((V+30)/10)+1)
  p = (1-h)*k1p-h*k1m
  return p
def vprime(Vm,m,n,h,inj):
  C=1
  gk=36
  gna=200#120
  gl= 0.3
  gca = 0.13
  Ek= -80#-12
  Ena= 50#115
  El= -65#10.613
  Eca = 75
  dV = (1/C)*(-gk*(n*n*n*n)*(Vm - Ek) - gna*(m*m*m)*h*(Vm-Ena) - gl*(Vm-El)+inj)
  return dV


# FROM Warman (1994) 
# Activation of Na+ channels
def mprimeW(t,m,V):
  # Calculates the derivative of m and a point t given Membrane potential V
  am = -1.74*(V-11)/(exp((V-11)/-12.94)-1.0)
  bm = 0.06*(V-5.9)/(exp((V-5.9)/4.47)-1.0)
  p = (1-m)*am-m*bm
  return p
# Activation of K+ channels
def nprimeW(t,n,V):
  # Calculates the derivative of n and a point t given Membrane potential V
  an = -18.0*V/((exp(V/-25))-1.0)
  bn = 3.6*(V-10)/(exp((V-10)/12)-1.0)
  p = (1-n)*an-n*bn
  return p
# Deactivation of Na+ channels
def hprimeW(t,h,V):
  # Calculates the derivative of n and a point t given Membrane potential V
  ah = 3.0/(exp((V+80)/10))
  bh = 12.0/(exp((V-77)/-27)-1.0)
  p = (1-h)*ah-h*bh
  return p
def vprimeW(Vm,m,n,h,inj):
  C=1
  gk=23
  gna=450
  gl=0.3
  Ek= -85 
  Ena= 45
  El= -67 
  dV = (1/C)*(-gk*(n*n*n*n)*(Vm - Ek) - gna*(m*m*m)*h*(Vm-Ena) - gl*(Vm-El)+inj)
  return dV


# FROM Grobler (1994) 
# Activation of Na+ channels
def mprimeG(t,m,V):
  # Calculates the derivative of m and a point t given Membrane potential V
  am = 1/(exp((-45.0-V)/4.0)+1.0)
  bm = 1/(exp(( 30.0+V)/4.0)+1.0)
  p = (1-m)*am-m*bm
  return p
# Activation of K+ channels
def nprimeG(t,n,V):
  # Calculates the derivative of n and a point t given Membrane potential V
  an = -18.0*V/((exp(V/-25))-1.0)
  bn = 3.6*(V-10)/(exp((V-10)/12)-1.0)
  p = (1-n)*an-n*bn
  return p
# Deactivation of Na+ channels
def hprimeG(t,h,V):
  # Calculates the derivative of n and a point t given Membrane potential V
  ah = 3.0/(exp((V+80)/10))
  bh = 12.0/(exp((V-77)/-27)-1.0)
  p = (1-h)*ah-h*bh
  return p
# Deactivation of Na+ channels
def qprimeG(t,V,Ca):
  ah = 3.0/(exp((V+80)/10))
  bh = 12.0/(exp((V-77)/-27)-1.0)
  p = (1-h)*ah-h*bh
  return p
def vprimeG(Vm,m,n,h,inj):
  C=1
  gk=  0.15
  gna= 0.03
  gl=  0.015
  gca= 0.100
  gkc= 0.15
  Ek= -95.0
  Ena= 50.0
  El= -67.0
  Eca= 75.0
  dV = (1/C)*(-gk*(n*n*n*n)*(Vm - Ek) - gna*(m*m*m)*h*(Vm-Ena) - gl*(Vm-El)+inj)
  return dV



def vprime2(Vm,m,n,h,inj):
  C=1
  gk=5
  gna=42
  gl=0.35
  Ek=-90
  Ena=55
  El=-30
  dV = (1/C)*(-gk*(n*n*n*n)*(Vm - Ek) - gna*(m*m*m)*h*(Vm-Ena) - gl*(Vm-El)+inj)
  return dV

################################################################################
# Use parameters from Hoonuraiah (2013) to replicate Shouval (2002)

def goCaConcentrationTest():
  # Calcium Sampling for concentration
  ca = getRampCa(0.0,1.0)
  #ca = getLTP()
  #ca = getLDP()
  #ca = getLDPLTP()
  #ca = getLTPLDP2()
  #ca = getLTPLDP()
  #ca = slowLTPFastLDP()
  c0 = 100e-6
  tau = 30
  vr = 30.77
  v = getRampCa(-90.0,55.0)
  bv = bov(v)
  #bovm(v,vr)
  plot1(v,bv,"V (mV)","B(V)",png="bov")
  ica = mul(bv,sub(vr,v))
  plot1(v,ica,"V (mV)","Y",png="Y")
  ca = rk4dca(st,c0,v,tau,ica)
  goCaConcentration(ca,v)
  bovm(v,vr,st,c0,tau,ica)

def bov(v,mgo=1):
  nv = len(v)
  bv = zerofloat(nv)
  for i in range(len(v)):
    bv[i] = 1/(1+mgo*exp(-0.062*v[i])/3.57)
  return bv

def bovm(v,vr,st,c0,tau,ica):
  mgi = 0.0
  mge = 5.0
  dmg = 1.0
  nmg = round((mge-mgi)/dmg)
  bvm,cam = [],[]
  for i in range(nmg):
    bvm.append(bov(v,mgi+dmg*i))
    ica = mul(bvm[i],sub(vr,v))
    ca = rk4dca(st,c0,v,tau,ica)
    cam.append(goCaConcentration(ca,v))
  plot1m(v,bvm,"V (mV)","B(V)",png="bovm")
  plot1m(st,cam,"Time (ms)","W",png="wam")
    
"""
RK4 solver for dc/dt
"""
def rk4dca(st,c0,v,tau,ica):
  return Utils.rk4dca(st,c0,tau,ica)
  nt = st.count
  dt = st.delta
  c = zerofloat(nt)
  c[0] = c0
  for i in range(nt-1):
    t = st.getValue(i)
    v1 = dt*dcdt2(t,c[i],ica[i],tau) 
    v2 = dt*dcdt2(t+dt/2.0,c[i]+v1/2.0,ica[i],tau)
    v3 = dt*dcdt2(t+dt/2.0,c[i]+v2/2.0,ica[i],tau)
    v4 = dt*dcdt2(t+dt,c[i]+v3,ica[i],tau)
    c[i+1] = c[i] + (v1+2*v2+2*v3+v4)/6.0
  return div(c,max(c))

def dcdt2(t,ca,ica,tau):
  return ica - ca/tau


def slowLTPFastLDP():
  tau1 = 50
  tau2 = 10
  return getLTPLDP(tau1,tau2)
   
def getLTPLDP2(tau1=30,tau2=-10):
  t = rampfloat(tmin,-dt,nt)
  ca = sub(exp(div(t,tau1)),exp(div(t,tau2)))
  return div(ca,-max(ca))

def getLTPLDP(tau1=30,tau2=30):
  nt2 = nt/2
  t1 = rampfloat(tmin,-dt,nt2)
  t2 = rampfloat(tmin, dt,nt2)
  ca1 = exp(div(t1,tau1))
  ca2 = exp(div(t2,tau2))
  ca3 = zerofloat(nt)
  mca1 = max(ca1)
  mca2 = max(ca2)
  for i in range(nt2):
    ca3[i+nt2] = -0.5*ca1[i]/mca1+0.5
    ca3[i    ] =  0.5*ca2[i]/mca2+0.5
  return ca3

def getLDPLTP(tau1=30,tau2=30):
  nt2 = nt/2
  t1 = rampfloat(tmin, dt,nt2)
  t2 = rampfloat(tmin,-dt,nt2)
  ca1 = exp(div(t1,tau1))
  ca2 = exp(div(t2,tau2))
  ca3 = zerofloat(nt)
  mca1 = max(ca1)
  mca2 = max(ca2)
  for i in range(nt2):
    ca3[i    ] = -0.5*ca1[i]/mca1+0.5
    ca3[i+nt2] =  0.5*ca2[i]/mca2+0.5
  return ca3

def getLDP():
  t = rampfloat(tmin,-dt,nt)
  ca = exp(div(t,30))
  return div(ca,-max(ca))

def getLTP():
  t = rampfloat(tmin,-dt,nt)
  ca = exp(div(t,30))
  return div(ca,max(ca))

def getRampCa(fc,lc):
  dc = (lc-fc)/nt
  return rampfloat(fc,dc,nt)

def goCaConcentration(ca,v):
  # Calcium Sampling for concentration
  w0 = 0.25
  can = div(ca,max(ca))
  eta = getLearningRate(can)
  omega = getOmega(can)
  dw = dwdta(st,len(can),w0,eta,omega)
  w = rk4dw(st,can,eta,omega,w0)
  hlb = "Ca Concentration"
  plot1(st,ca,"Time (ms)","Ca concentration",png="ca")
  plot1(v,ca,"V (mV)","Ca concentration",png="cav")
  plot1(ca,dw,hlb,"dw/dt",png="dwdt")
  plot1(ca,eta,hlb,"Learning Rate",png="eta")
  plot1(ca,omega,hlb,"Omega",vlim=[0.0,1.1],png="omega")
  plot1(st,w,"Time (ms)","W",png="w")
  return w


def getLearningRate(ca):
  p1 = 1 #s
  p2 = 0.1 #s
  p3 = p2*1e-4
  p4 = 3
  tau = add(p1,div(p2,add(p3,pow(ca,p4))))
  return div(1,tau)

def getOmega(ca):
  a1 = 0.35
  a2 = 0.55
  b1 = 80
  b2 = 80
  omega = add(0.25,add(div(1,add(1,exp(mul(-b2,add(ca,-a2))))),\
            mul(-0.25,div(1,add(1,exp(mul(-b1,add(ca,-a1))))))))
  return omega

"""
RK4 solver for dw/dt
"""
def rk4dw(st,ca,eta,omega,w0):
  return Utils.rk4dw(st,ca,eta,omega,w0)
  nt = st.count
  dt = st.delta
  w = zerofloat(nt)
  w[0] = w0
  l = 1
  for i in range(nt-1):
    v1 = dt*dwdt(st.getValue(i),w[i],eta[i],omega[i],l) 
    v2 = dt*dwdt(st.getValue(i)+dt/2.0,w[i]+v1/2.0,eta[i],omega[i],l)
    v3 = dt*dwdt(st.getValue(i)+dt/2.0,w[i]+v2/2.0,eta[i],omega[i],l)
    v4 = dt*dwdt(st.getValue(i)+dt,w[i]+v3,eta[i],omega[i],l)
    w[i+1] = w[i] + (v1+2*v2+2*v3+v4)/6.0
  return w

"""
Calcium control hypothesis ODE
@param t     time
@param eta   Ca dependent learning rate (constant)
@param omega Ca dependent sign magnitude (constant)
@param w     synaptic plasticity weight 
@param l     decay constant 
"""
def dwdt(t,w,eta,omega,l=1):
  return eta*(omega-w*l)

def dwdta(st,na,w0,eta,omega,l=1):
  dw = zerofloat(na)
  for i in range(na):
    dw[i] = dwdt(st,w0,eta[i],omega[i])
  return dw

################################################################################
# plots

def plot1(s1,x1,hlabel=None,vlabel=None,vlim=None,hlim=None,png=None):
  pp = PlotPanel(PlotPanel.Orientation.X1RIGHT_X2UP)
  l1 = pp.addPoints(s1,x1)
  if vlabel:
    pp.setVLabel(vlabel)
  if hlabel:
    pp.setHLabel(hlabel)
  if vlim:
    pp.setVLimits(vlim[0],vlim[1])
  if hlim:
    pp.setHLimits(hlim[0],hlim[1])
  frame = PlotFrame(pp)
  frame.setFontSize(24)
  frame.setSize(700,600)
  frame.setVisible(True)
  if png and prints: 
    printSlide(frame,png,0.8,0.8,4.0/3.0)
    printPaper(frame,png)

def plot1m(s1,x1,hlabel=None,vlabel=None,vlim=None,hlim=None,png=None):
  pp = PlotPanel(PlotPanel.Orientation.X1RIGHT_X2UP)
  for i in range(len(x1)):
    l1 = pp.addPoints(s1,x1[i])
  if vlabel:
    pp.setVLabel(vlabel)
  if hlabel:
    pp.setHLabel(hlabel)
  if vlim:
    pp.setVLimits(vlim[0],vlim[1])
  if hlim:
    pp.setHLimits(hlim[0],hlim[1])
  frame = PlotFrame(pp)
  frame.setFontSize(24)
  frame.setSize(700,600)
  frame.setVisible(True)
  if png and prints: 
    printSlide(frame,png,0.8,0.8,4.0/3.0)
    printPaper(frame,png)


def plot(x1,x2,x3=None,x4=None,title=None,vlabel=None,hlabel=None,png=None):
  pp = PlotPanel(PlotPanel.Orientation.X1RIGHT_X2UP)
  l1 = pp.addPoints(x1,x2)
  pp.setVLabel("time (s)")
  pp.setHLabel("offset (km)")
  if x3:
    l2 = pp.addPoints(x1,x3)
    l2.setLineColor(Color.RED)
  if x4:
    l3 = pp.addPoints(x1,x4)
    l3.setLineColor(Color.BLUE)
  if title:
    pp.setTitle(title)
  if titles:
    pp.setTitle(png)
  if vlabel:
    pp.setVLabel(vlabel)
  if hlabel:
    pp.setHLabel(hlabel)
  frame = PlotFrame(pp)
  frame.setFontSize(24)
  frame.setSize(600,600)
  frame.setVisible(True)
  if png and prints: 
    printSlide(frame,png,0.8,0.8,4.0/3.0)
    printPaper(frame,png)

"""
Method to print slides.
@param frame a PlotFrame class
@param fname the png file name
@param fw    the fraction of the slide width that the figure occupies
@param fh    the fraction of the slide height that the figure occupies
@param ar    the aspect ratio of the slide (use keynote default pixel widths)
@param sc    scalar for frame size that is dependent on screen resolution 
              (my MBP 15" needs sc=2)
"""
def printSlide(frame,fname,fw,fh,ar,sc=2):
  if not slides:
    return
  swp = 1920 # 16/9 keynote slide default 
  if ar==4.0/3.0: swp = 1024 # 4/3 keynote slide default
  fwi = int(swp*fw/sc)+1
  fhi = int(swp/ar*fh/sc)+1
  frame.setSize(fwi,fhi)
  frame.setFontSizeForSlide(fw,fh,ar) 
  frame.paintToPng(swp*fw,1.0,pngDir+"s_"+fname+".png")

# Dimensions for Geophysics paper
def printPaper(frame,fname,sc=8,colm=1,dpi=720):
  if not papers:
    return
  fs = 11.0
  if colm==1: 
    swp = 2400.0
  if colm==2: 
    swp = 3120.0
  if colm==3: 
    swp = 5040.0
  ar = 4.0/3.0;
  fwi = int(swp/sc)+1
  fhi = int((swp/ar)/sc)+1
  frame.setSize(fwi,fhi)
  frame.setFontSizeForPrint(fs,swp/10.0) 
  frame.paintToPng(dpi,swp/dpi,pngDir+"p_"+pfx+fname+".png")


################################################################################
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
if __name__=="__main__":
  SwingUtilities.invokeLater(RunMain()) 

