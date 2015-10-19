"""
Plot methods for well ties in Teapot Dome.
Author: Andrew Munoz, Colorado School of Mines
Version: 2013.04.09
"""
from dtw import *
from imports import *
from wtutils import *

_pngDir = "./tppics/"
wtu = WTutils()

gray = ColorMap.GRAY
jet = ColorMap.JET
rwb = ColorMap.RED_WHITE_BLUE
xrxu = PlotPanel.Orientation.X1RIGHT_X2UP
xdxr = PlotPanel.Orientation.X1DOWN_X2RIGHT
aplt = PlotPanel.AxesPlacement.LEFT_TOP
inear = PixelsView.Interpolation.NEAREST
eoc = PlotFrame.EXIT_ON_CLOSE

def plotMatrix(c,sf,slag,u=None,cp=None,lim=None,png=None,paper=None,slides=None):
  n1,nlag = len(c[0]),len(c)
  #slag = Sampling(nlag,1.0,-(nlag-1)/2)
  panel = PlotPanel(xrxu)
  cv = panel.addPixels(sf,slag,c)
  cv.setInterpolation(inear)
  cv.setColorModel(jet)
  panel.setVLimits(slag.first,slag.last)
  if u:
    uv = panel.addPoints(sf,u)
    uv.setLineColor(Color.WHITE)
    uv.setLineWidth(3)
    panel.setVLimits(min(u)*0.5,max(u)*1.5)
  if cp:
    cv.setPercentiles(0,cp)
  if lim:
  	panel.setVLimits(lim[0],lim[1])
  panel.setHLimits(sf.first,sf.last)
  panel.setHLabel("time (s)")
  panel.setVLabel("time lag (s)")
  panel.addColorBar()
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(eoc)
  frame.setFontSize(18)
  frame.setSize(1000,800)
  frame.setVisible(True)
  if png:
    frame.paintToPng(720,3.33,_pngDir+png+'.png')
  if paper:
    frame.setSize(732,641)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')
  if slides:
    frame.setSize(550,500)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,3.50,_pngDir+slides+'.png')
  
def plotSequences(f,g,s1=None,s2=None,u=None,title=None):
  n1 = len(f)
  n2 = len(g)
  if s1 is None:
    s1 = Sampling(n1,1.0,0.0)
  if s2 is None:
    s2 = Sampling(n2,1.0,0.0)
  panel = PlotPanel(xdxr)
  fv = panel.addPoints(s1,f)
  gv = panel.addPoints(s2,g)
  if u:
    gv = panel.addPoints(u,g)
  gv.setLineColor(Color.RED)
  panel.setVLabel("time (s)")
  if title:
    panel.setTitle(title)
  frame = PlotFrame(panel)
  frame.setDefaultCloseOperation(eoc)
  frame.setFontSize(18)
  frame.setSize(300,800)
  frame.setVisible(True)

"""
Plots a panel of well logs
"""
def plotLogPanel(sz,v,d,rf,paper=None,slides=None):
  p1 = PlotPanel(1,3,xdxr)
  p1.setVLabel("Depth (km)")
  vp = p1.addPoints(0,0,sz,v)
  dp = p1.addPoints(0,1,sz,d)
  rf = p1.addPoints(0,2,sz,rf)
  p1.setHLabel(0,"Velocity (km/s)")
  p1.setHLabel(1,"Density (g/cc)")
  p1.setHLabel(2,"Reflectivity")
  p1.setHLimits(0,min(v),max(v))
  p1.setHLimits(1,1.51,max(d))
  p1.setHInterval(1,0.5)
  #p1.setHLimits(2,-0.13,0.13)
  p1.setHInterval(2,0.1)
  #dp.setLineColor(BLUE)
  #gp.setLineColor(GREEN)
  frame = PlotFrame(p1)
  frame.setSize(800,1500)
  frame.setFontSize(28)
  frame.setDefaultCloseOperation(eoc);
  if paper:
    frame.setSize(443,643)
    frame.setFontSizeForPrint(8.0,222.0)
    frame.paintToPng(720,3.08,_pngDir+paper+'.png')
  if slides:
    frame.setSize(690,561)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)

"""
Plots a panel of 3 seismograms
"""
def plotSeismogramPanel(ssy1,ssy2,ssy3,sy1,sy2,sy3,paper=None,slides=None):
  p1 = PlotPanel(1,3,xdxr)
  p1.setVLabel("Time (s)")
  vp = p1.addPoints(0,0,ssy1,sy1)
  dp = p1.addPoints(0,1,ssy2,sy2)
  rf = p1.addPoints(0,2,ssy3,sy3)
  p1.setHLabel(0,"simple")
  p1.setHLabel(1,"multiples")
  p1.setHLabel(2,"multiples and Q")
  #p1.setHLimits(0,-0.51,0.51)
  p1.setHLimits(1,-0.29,0.20)
  #p1.setHLimits(2,-0.51,0.51)
  #dp.setLineColor(BLUE)
  #gp.setLineColor(GREEN)
  frame = PlotFrame(p1)
  frame.setSize(800,1200)
  frame.setFontSize(28)
  frame.setDefaultCloseOperation(eoc);
  if paper:
    frame.setSize(443,643)
    frame.setFontSizeForPrint(8.0,222.0)
    frame.paintToPng(720,3.08,_pngDir+paper+'.png')
  if slides:
    frame.setSize(690,561)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)

"""
Plots a panel of 3 velocities
"""
def plotVelPanel(sz,sy1,sy2,sy3,paper=None,slides=None):
  p1 = PlotPanel(1,3,xdxr)
  p1.setVLabel("Depth (km)")
  vp = p1.addPoints(0,0,sz,sy1)
  dp = p1.addPoints(0,1,sz,sy2)
  rf = p1.addPoints(0,2,sz,sy3)
  p1.setHLabel(0,"velocity (km/s)")
  p1.setHLabel(1,"velocity (km/s)")
  p1.setHLabel(2,"velocity (km/s)")
  #p1.setHLimits(0,-0.51,0.51)
  #p1.setHLimits(1,-0.29,0.20)
  #p1.setHLimits(2,-0.51,0.51)
  #dp.setLineColor(BLUE)
  #gp.setLineColor(GREEN)
  frame = PlotFrame(p1)
  frame.setSize(800,1200)
  frame.setFontSize(28)
  frame.setDefaultCloseOperation(eoc);
  if paper:
    frame.setSize(443,643)
    frame.setFontSizeForPrint(8.0,222.0)
    frame.paintToPng(720,3.08,_pngDir+paper+'.png')
  if slides:
    frame.setSize(690,561)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)


"""
Plots 3 curves on top of eachother
"""
def plotCurvePanel3(sf,sg,f,g,sh=None,h=None,tlim=None,paper=None,slides=None):
	if h:
		p1 = PlotPanel(3,1,xrxu,aplt)
		pf = p1.addPoints(0,0,sf,f)
		pg = p1.addPoints(1,0,sg,g)
		pg2f = p1.addPoints(0,0,sh,h)
		pg2g = p1.addPoints(2,0,sh,h)
		pg2f.setLineColor(RED)
		pg2g.setLineColor(RED)
		#p1.setVLabel(0,"f & gw")
		#p1.setVLabel(1,"g & gw")
	else: 
		p1 = PlotPanel(2,1,xrxu,aplt)
		pf = p1.addPoints(0,0,sf,f)
		pg = p1.addPoints(1,0,sg,g)
	p1.setHLabel(0,"Time (s)")
	#p1.setVLabel(0,"f")
	#p1.setVLabel(1,"g")
	if tlim:
		p1.setHLimits(0,tlim)
	else: p1.setHLimits(0,sf.last)
	frame = PlotFrame(p1)
	frame.setDefaultCloseOperation(eoc)
	frame.setFontSize(28)
	frame.setSize(1500,350)
	if paper:
		frame.setSize(1100,413)
		frame.setFontSizeForPrint(8.0,240.0)
		frame.paintToPng(720,3.33,_pngDir+paper+'.png')
	if slides:
		#frame.setSize(1000,500)
		frame.setSize(1000,350)
		frame.setFontSizeForSlide(1.0,1.0)
		frame.paintToPng(720,7.00,_pngDir+slides+'.png')
		#if g2:
		#  frame.setSize(965,455)
		#else: frame.setSize(965,349)
		#frame.setFontSizeForSlide(1.0,1.0)
		#frame.paintToPng(720,7.00,_pngDir+slides+'.png')
	frame.setVisible(True)

def plotCurvePanel2(sf,sg,sh,f,g,h,lims=None,paper=None,slides=None):
	p1 = PlotPanel(2,1,xrxu,aplt)
	p1.addPoints(0,0,sg,g)
	p1.addPoints(1,0,sg,g)
	pf = p1.addPoints(0,0,sf,f)
	ph = p1.addPoints(1,0,sh,h)
	pf.setLineColor(RED)
	ph.setLineColor(RED)
	#p1.setVLabel(0,"f & g")
	#p1.setVLabel(1,"f & h")
	p1.setHLabel(0,"Time (s)")
	if lims:
		p1.setHLimits(lims[0],lims[1])
	frame = PlotFrame(p1)
	frame.setDefaultCloseOperation(eoc)
	frame.setFontSize(28)
	frame.setSize(1500,350)
	if paper:
		frame.setSize(712,418)
		frame.setFontSizeForPrint(8.0,222.0)
		frame.paintToPng(720,3.08,_pngDir+paper+'.png')
	if slides:
		#frame.setSize(1000,500)
		frame.setSize(1000,350)
		frame.setFontSizeForSlide(1.0,1.0)
		frame.paintToPng(720,7.00,_pngDir+slides+'.png')
		#if g2:
		#  frame.setSize(965,455)
		#else: frame.setSize(965,349)
		#frame.setFontSizeForSlide(1.0,1.0)
		#frame.paintToPng(720,7.00,_pngDir+slides+'.png')
	frame.setVisible(True)

"""
Plots comparison of reflectvity to its resulting synthetic seismogram
"""
def plotCurveW(sz,rf,ss,sy,paper=None):
  p1 = PlotPanel(1,1,xdxr)
  p2 = PlotPanel(1,1,xdxr)
  p1.addPoints(0,0,sz,rf)
  p2.addPoints(0,0,ss,sy)
  p1.setVLabel(0,"Depth (km)")
  p2.setVLabel(0,"t (s)")
  p1.setHLimits(0,-0.13,0.13)
  p2.setHLimits(0,-0.5,0.5)
  frame = PlotFrame(p1,p2,PlotFrame.Split.HORIZONTAL)
  frame.setSize(300,1500)
  frame.setDefaultCloseOperation(eoc);
  frame.setVisible(True)
  if paper:
    frame.setSize(584,757)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')

"""
Plots comparison of synthetic seismogram before and after normalization
"""
def plotCurveN(st,tr,trn,ss,sy,syn,lim=None,paper=None):
  p1 = PlotPanel(1,2,xdxr)
  p2 = PlotPanel(1,2,xdxr)
  p1.addPoints(0,0,st,tr)
  p1.addPoints(0,1,st,trn)
  p2.addPoints(0,0,ss,sy)
  p2.addPoints(0,1,ss,syn)
  p1.setVLabel(0,"Time (s)")
  p2.setVLabel(0,"t (s)")
  p2.setHLimits(0,-1.1,1.1)
  frame = PlotFrame(p1,p2,PlotFrame.Split.HORIZONTAL)
  frame.setSize(300,1500)
  frame.setDefaultCloseOperation(eoc);
  frame.setVisible(True)
  if paper:
    frame.setSize(584,757)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper+'.png')

"""
Plots time-depth curves before and after warping
"""
def plotTDCurves(sz,tauz,tz,tzn=None,paper=None,slides=None):
  p1 = PlotPanel(xdxr)
  l1 = p1.addPoints(sz,tz)
  l2 = p1.addPoints(sz,tauz)
  l1.setLineColor(RED)
  if tzn:
    l3 = p1.addPoints(sz,tzn)
    l3.setLineColor(RED)
    l1.setLineColor(BLUE)
  p1.setHLabel('Time (s)')
  p1.setVLabel('Depth (km)')
  frame = PlotFrame(p1)
  frame.setSize(700,500)
  frame.setFontSizeForPrint(8.0,240.0)
  frame.setDefaultCloseOperation(eoc);
  if paper:
		#colm = 222.0
    colm = 148.0
    frame.setSize(513,379) 
    frame.setFontSizeForPrint(8.0,colm)
    frame.paintToPng(720,colm/72.0,_pngDir+paper+'.png')
  if slides:
    frame.setSize(1000,500)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
    #frame.setSize(513,379)
    #frame.setFontSizeForSlide(1.0,1.0)
    #frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)

def plotTDandV(sz,v1,v2,tauz,tz,h2label=None,paper=None,slides=None):
	p1 = PlotPanel(1,2,xdxr)  
	p1.mosaic.setWidthElastic(0,67)
	p1.mosaic.setWidthElastic(1,33)
	l1 = p1.addPoints(0,0,sz,tz)                              
	l2 = p1.addPoints(0,0,sz,tauz)                          	
	l1.setLineColor(RED)                                      
	p1.setHLabel(0,'Time (s)')                                
	p1.setVLabel(0,'Depth (km)')                            	
	if h2label:                                             	
		p1.setHLabel(1,h2label)                               	
	l1 = p1.addPoints(0,1,sz,v1)                            	
	l1.setLineColor(BLACK)                                  	
	l2 = p1.addPoints(0,1,sz,v2)                            	
	l2.setLineColor(RED)                                    	
	frame = PlotFrame(p1)                                   	
	frame.setSize(1000,1000)	                              	
	if slides:                                              	
		frame.setSize(800,800)                                	
		frame.setFontSizeForSlide(1.0,1.0)                    	
		frame.paintToPng(720,7.00,_pngDir+slides+'.png')      	
	if paper:                                               	
		colm = 222.0                                          	
		frame.setSize(712,524)                                	
		frame.setFontSizeForPrint(8.0,colm)                   	
		frame.paintToPng(720,colm/72.0,_pngDir+paper+'.png')  	
	frame.setVisible(True)                                  	
	frame.setDefaultCloseOperation(eoc)

"""
Plots a curve horizontally and optionally another as a red in the same panel
"""
def plotCurveHorz(c1,s1,c2=None,s2=None,title=None,lim=None,hlabel=None,png=None,d=None,slides=None,paper=None):
	p1 = PlotPanel(xrxu)
	l1 = p1.addPoints(s1,c1)
	l1.setLineColor(BLACK)
	if c2 and s2:
		l2 = p1.addPoints(s2,c2)
		l2.setLineColor(RED)
	if d=="z":
		p1.setHLabel("Depth (km)")
	if d=="t":
		p1.setHLabel("Time (s)")
	if lim:
		p1.setVLimits(lim[0],lim[1])
	if title: 
		p1.setTitle(title)
	if hlabel:
		p1.setVLabel(hlabel)
	frame = PlotFrame(p1)
	frame.setSize(1500,300)
	if png:
		frame.paintToPng(720,3.33,_pngDir+png+'.png') 
	if slides:
		frame.setSize(800,300)
		frame.setFontSizeForSlide(1.0,1.0)
		frame.paintToPng(720,7.00,_pngDir+slides+'.png')
	if paper:
		#colm = 90.0
		colm=222.0
		frame.setSize(581,266)
		frame.setFontSizeForPrint(8.0,colm)
		frame.paintToPng(720,colm/72.0,_pngDir+paper+'.png')
	frame.setVisible(True)
	frame.setDefaultCloseOperation(eoc);

"""
Plots a curve vertically and optionally another as a red in the same panel
"""
def plotCurve(c1,s1,c2=None,s2=None,title=None,lim=None,hlabel=None,png=None,d=None,slides=None,paper=None,psize=None):
	p1 = PlotPanel(xdxr)
	l1 = p1.addPoints(s1,c1)
	l1.setLineColor(BLACK)
	if c2 and s2:
		l2 = p1.addPoints(s2,c2)
		l2.setLineColor(RED)
	p1.setVLabel("samples")
	if d=="z":
		p1.setVLabel("Depth (km)")
	if d=="t":
		p1.setVLabel("Time (s)")
	if lim:
		p1.setHLimits(-lim,lim)
	if title: 
		p1.setTitle(title)
	if hlabel:
		p1.setHLabel(hlabel)
	frame = PlotFrame(p1)
	frame.setSize(300,1500)
	if png:
		frame.paintToPng(720,3.33,_pngDir+png+'.png') 
	if slides:
		frame.setSize(300,800)
		frame.setFontSizeForSlide(1.0,1.0)
		frame.paintToPng(720,7.00,_pngDir+slides+'.png')
	if paper:
		colm = 90.0
		frame.setSize(200,379)
		frame.setFontSizeForPrint(8.0,colm)
		frame.paintToPng(720,colm/72.0,_pngDir+paper+'.png')
	frame.setVisible(True)
	frame.setDefaultCloseOperation(eoc);


"""
Plots a 2-D seismic slice, and well logs through that slice, and seismic horizons on the slice
"""
def plot2DLogs(s1,s2,img,sh,h,x2w,wwidth=1,tops=None,fmsp=None,tz=None,sz=None,gh=None,\
               title=None,nolog=None,paper=None,slides=None,cmap=gray,lims=None,cbar=None,\
								psize=None,lcmap=gray,clips=None):
	nw = len(sh)
	g = copy(img)
	#if isEven(wwidth): wwidth -= 1
	for i in range(nw):
		x2i = inro((x2w[i]-s2.first)/s2.delta)
		shi = infl((sh[i].first-s1.first)/s1.delta)
		nj = len(h[i])
		for j in range(nj):
			g[x2i][j+shi] = h[i][j]
	sp = PlotPanel(xdxr) 
	im1 = sp.addPixels(s1,s2,g)
	im1.setColorModel(cmap)
	im1.setInterpolation(inear)
	if clips:
		im1.setClips(clips[0],clips[1])
	if gh:
		# Add Horizons
		for i in range(len(gh[0])):
			pv = sp.addPoints(gh[0][i],gh[1][i])
			pv.setStyle(gh[2][i])
			pv.setLineWidth(2)
	# Add tops
	if tops:
		for iw in range(nw):
			nfms = len(fmsp)
			x2i = int(round((x2w[i]-s2.first)/s2.delta))
			for ic in range(nfms):
				tpd = tops[iw][ic]
				if tpd>=sz[iw].first:
					iz = sz[iw].indexOfNearest(tpd)
					ttpd = tz[iw][iz]
					ttpd = fillfloat(ttpd,1)
					tcolor = fmsp[ic][1]
					pvt = sp.addPoints(ttpd,xvs)
					pvt.setStyle(tcolor); pvt.setLineWidth(3)
	sp.setHLabel("Distance (km)")
	sp.setVLabel("Time (s)")
	if lims:
		sp.setLimits(lims[0],lims[1],lims[2],lims[3])
	else: 
		sp.setLimits(s2.first,s1.first,s2.last,s1.last)
	if title:
		sp.setTitle(title)
	if paper or slides:
		sp.removeTitle()
	if cbar:
		sp.addColorBar()
	frame = PlotFrame(sp)
	frame.setSize(1000,1000)
	frame.setFontSize(24)
	if paper:
		frame.setSize(712,524)
		if psize:
			frame.setSize(psize[0],psize[1])
		frame.setFontSizeForPrint(8.0,222.0)
		frame.paintToPng(720,3.08,_pngDir+paper+'.png')
	if slides:
		frame.setSize(1000,700)
		frame.setFontSizeForSlide(1.0,1.0)
		frame.paintToPng(720,7.00,_pngDir+slides+'.png')
	frame.setVisible(True)
	frame.setDefaultCloseOperation(eoc);

def plotLogPanelsTops(st,sz,a,g,wells,fmsp,tops,zt,tz,gh,x2w):
	nw = len(a); 
	alim = [0,20]
	#alim = [1,8]
	glim = [0,200]
	p1 = PlotPanel(1,nw*2,xdxr)
	p1.setVLabel("Depth (km)")
	nfms = len(fmsp)
	fmtps,fmhrz = [],[]
	for ic in range(nfms):
		for clr in fmsp[ic][1]:
			fmtps.append((fmsp[ic][0],clr+"O"))
			fmhrz.append((fmsp[ic][0],clr+"S"))
			break
	hc = zeroint(nfms,nw)
	for iw in range(nw):
		x2c = x2w[iw]
		for ic in range(nfms):
			ghtm = gh[1][ic][0]
			for j in range(1,len(gh[1][ic])):
				ghtp = gh[1][ic][j]
				if (ghtm<=x2c and x2c<ghtp):
					hc[iw][ic] = j
				ghtm = ghtp
	iwp = 0
	def addTops(x,iwp,iw,mrk):
		dist = 1.0
		for ic in range(nfms):
			fmpt = fillfloat(-99.9,sz[iw].count)
			tpd = tops[iw][ic]
			if tpd>sz[iw].first:
				tcolor = fmtps[ic][1]
				iz = sz[iw].indexOfNearest(tpd)
				fmpt[iz] = x[iz]
				pvt = p1.addPoints(0,iwp,sz[iw],fmpt)
				pvt.setStyle(tcolor); pvt.setMarkSize(mrk)
			fmpt = fillfloat(-99.9,sz[iw].count)
			hrt = int((gh[0][ic][hc[iw][ic]]-tz[iw][0])/st.delta)
			if hrt<len(zt[iw]) and hrt>0:
				hpd = zt[iw][hrt]
				iz = sz[iw].indexOfNearest(hpd)
				fmpt[iz] = x[iz]
				pvh = p1.addPoints(0,iwp,sz[iw],fmpt)
				tcolorh = fmhrz[ic][1]
				pvh.setStyle(tcolorh); pvh.setMarkSize(mrk)
				#if tpd>0.0 and ic==0 and iw==3: 
				#	dist = abs(hpd-tpd)
	for iw in range(nw):
		ap = p1.addPoints(0,iwp,sz[iw],a[iw])
		p1.setHLabel(iwp,"I"+str(wells[iw]))
		p1.setHLimits(iwp,alim[0],alim[1])
		dist = addTops(a[iw],iwp,iw,10)
		iwp+=1
	for iw in range(nw):
		gp = p1.addPoints(0,iwp,g[iw][1],g[iw][0])
		p1.setHLabel(iwp,"G"+str(wells[iw]))
		p1.setHLimits(iwp,glim[0],glim[1])
		addTops(g[iw][0],iwp,iw,10)
		iwp+=1
	frame = PlotFrame(p1)
	frame.setSize(1500,1500)
	frame.setFontSize(28)
	frame.setDefaultCloseOperation(eoc);
	frame.setVisible(True)



"""
Plots a 2-D seismic slice
"""
def plotSlice(sy,st,img,title=None,lims=None,png=None,slides=None,paper=None,cmap=gray,tens=None,cbar=None,psize=None,clips=None,pts=None):
	p1 = PlotPanel(xdxr)
	px = p1.addPixels(st,sy,img)
	if clips:
		px.setClips(clips[0],clips[1])
	p1.setHLabel("Distance (km)")
	p1.setVLabel("Time (s)")
	if title:
		p1.setTitle(title)
	px.setColorModel(cmap)
	if cbar:
		p1.addColorBar(cbar)
		p1.setColorBarWidthMinimum(100)
	if lims:
		p1.setLimits(lims[0],lims[1],lims[2],lims[3])
	if tens:
		tv = TensorsView(st,sy,tens)
		tv.setOrientation(TensorsView.Orientation.X1DOWN_X2RIGHT)
		tv.setLineColor(Color.YELLOW)
		tv.setLineWidth(2)
		#tv.setEllipsesDisplayed(40)
		fac1 = 40; fac2 = 20
		e1 = Sampling(st.count/fac1+1,st.delta*fac1,st.first)
		e2 = Sampling(sy.count/fac2+1,sy.delta*fac2,sy.first)
		tv.setEllipsesDisplayed(e1,e2)
		tv.setScale(0.55)
		#tv.setScale(10.00)
		tile = p1.getTile(0,0)
		tile.addTiledView(tv)
	if pts:
		ptv = p1.addPoints(pts[0],pts[1])
		ptv.setStyle("rO")
		ptv.setMarkSize(8)
	if slides: p1.removeTitle()
	if paper: p1.removeTitle()
	frame = PlotFrame(p1)
	frame.setSize(1000,1000)
	frame.setFontSize(24)
	frame.setDefaultCloseOperation(eoc);
	frame.setVisible(True)
	if png:
		frame.paintToPng(1000,5.00,_pngDir+png+'.png')  
	if slides:
		frame.setSize(1000,700)
		frame.setFontSizeForSlide(1.0,1.0)
		frame.paintToPng(720,7.00,_pngDir+slides+'.png')
	if paper:
		frame.setSize(712,524)
		if psize:
			frame.setSize(psize[0],psize[1])
		frame.setFontSizeForPrint(8.0,222.0)
		frame.paintToPng(720,3.08,_pngDir+paper+'.png')
	frame.setVisible(True)

def plotCoordSlice(s2,s3,x2s,x3s,x2w,x3w,tpimg,slides=None,paper=None):
	pp = PlotPanel()
	pp.setHLabel("Crossline (km)")
	pp.setVLabel("Inline (km)")
	pp.addPixels(s2,s3,tpimg)
	ln = pp.addPoints(x2s,x3s)
	ln.setStyle("y-")
	ln.setLineWidth(2)
	pt = pp.addPoints(x2w,x3w)
	pt.setStyle("yO")
	frame = PlotFrame(pp)
	frame.setSize(1000,1000)
	frame.setFontSize(24)
	frame.setDefaultCloseOperation(eoc)
	if slides:
		frame.setSize(1000,700)
		frame.setFontSizeForSlide(1.0,1.0)
		frame.paintToPng(720,7.00,_pngDir+slides+'.png')
	if paper:
		frame.setSize(556,314)
		frame.setFontSizeForPrint(8.0,222.0)
		frame.paintToPng(720,3.08,_pngDir+paper+'.png')
	frame.setVisible(True)

def plotHSlice(s2,s3,tpimg,x2w=None,x3w=None,slides=None,paper=None):
	pp = PlotPanel()
	pp.setHLabel("Distance (km)")
	pp.setVLabel("Distance (km)")
	pp.addPixels(s2,s3,tpimg)
	if x2w and x3w:
		pt = pp.addPoints(x2w,x3w)
		pt.setStyle("yO")
		pt.setMarkSize(5)
	frame = PlotFrame(pp)
	frame.setSize(1000,1000)
	frame.setFontSize(24)
	frame.setDefaultCloseOperation(eoc)
	if slides:
		frame.setSize(1000,700)
		frame.setFontSizeForSlide(1.0,1.0)
		frame.paintToPng(720,7.00,_pngDir+slides+'.png')
	if paper:
		frame.setSize(556,314)
		frame.setFontSizeForPrint(8.0,222.0)
		frame.paintToPng(720,3.08,_pngDir+paper+'.png')
	frame.setVisible(True)

def plotPanel3(s1,s2,s3,x,clrpix=None,maps=None,mpoints=None,\
		label1=None,label2=None,label3=None,cbar=gray,coord=None,paper=None,slides=None):
	svp = Viewer3P(s1,s2,s3,x)
	if maps:	
		svp.getPP().pixelsView23.tile.addTiledView(maps)
	if mpoints:
		svp.getPP().pixelsView23.tile.addTiledView(mpoints)
	svp.setColorModel1(cbar)
	if label1:
		svp.setLabel1(label1)
	if label2:
		svp.setLabel2(label2)
	if label3:
		svp.setLabel3(label3)

	svp.show()
	



"""
Plots phase rotation errors
"""
def plotPhaseRMS(mrms,nph,phf,png=None,paper=None,slides=None):
  sp = SimplePlot()
  sph = Sampling(nph/phf,phf,0.0)
  sp.addPoints(mrms)
  sp.setHLabel("Phase (degrees)")
  sp.setVLabel("Minimum Dn")
  sp.setTitle("Minumum RMS vs. Phase")
  sp.setSize(1500,1000)
  if png:
    sp.paintToPng(1200,6.00,_pngDir+png+'.png')
  if paper:
    sp.removeTitle()
    sp.setSize(1000,800)
    sp.setFontSizeForPrint(8.0,240.0)
    sp.paintToPng(720,3.33,_pngDir+paper+'.png')  
  if slides:
    sp.removeTitle()
    sp.setSize(950,800)
    sp.setFontSizeForSlide(1.0,1.0)
    sp.paintToPng(720,7.00,_pngDir+slides+'.png')

"""
Plots u0 vs dmin
"""
def plotuo(us,sl,png=None,paper=None,slides=None):
  sp = SimplePlot()
  sp.addPoints(sl,us)
  sp.setHLabel("u0 time shift (s)")
  sp.setVLabel("minimum D")
  sp.setSize(1200,1000)
  if png:
    sp.paintToPng(1200,6.00,_pngDir+png+'.png')
  if paper:
    sp.removeTitle()
    sp.setSize(1000,800)
    sp.setFontSizeForPrint(8.0,240.0)
    sp.paintToPng(720,3.33,_pngDir+paper+'.png')  
  if slides:
    sp.removeTitle()
    sp.setSize(950,800)
    sp.setFontSizeForSlide(1.0,1.0)
    sp.paintToPng(720,7.00,_pngDir+slides+'.png')

def plotCompare(s1,s2,x1,x2,title):
	pp = PlotPanel(xdxr)
	pp.setTitle(title)
	pp.setVLabel('time (s)')
	ln1 = pp.addPoints(s1,x1)
	ln2 = pp.addPoints(s2,x2)
	ln2.setLineColor(Color.RED)
	#pp.setVLimits(max(s1.first,s2.first),max(s1.last,s2.last))
	frame = PlotFrame(pp)
	frame.setDefaultCloseOperation(eoc)
	frame.setSize(500,1000)
	frame.setVisible(True)

"""
Plots just a gaussian
"""
def gaussianPlot(sig):
  nf = 801
  f1 = zerofloat(nf)
  f2 = zerofloat(nf)
  f1[int(nf/2)] = 100000.0
  RecursiveGaussianFilter(sig).apply2(f1,f2)
  f2 = mul(-1,f2)
  si = SimplePlot()
  si.addPoints(f2)
  si.setSize(700,500)
  si.paintToPng(900,3.33,_pngDir+"RICKER"+".png")

"""
Histogram plot 
@author Farhad Bazarghani
"""
def goHistogram(data,xlabel,png=None):
  hg = Histogram(data)
  h = hg.getCounts()
  hf = zerofloat(len(h))
  for i in range (1,len(h)):
    hf[i]=float(h[i])
  sam = hg.getBinSampling()
  sp = SimplePlot()
  sp1 = sp.addSequence(sam,hf)
  sp.setTitle("histogram"+" "+xlabel)
  sp.setVLabel("# of points")
  sp.setHLabel(xlabel)
  if png:
    sp.paintToPng(360,3.33,_pngDir+png+'.png') 


"""
Plots a wavelet impulse response
"""
def waveletTest():
  fp = 35 
  nt = 200; nz = 700;
  st = Sampling(nt,0.001,0.0)
  sz = Sampling(nz,0.0001524,0.0)
  rf = zerofloat(nz)
  #rf[0] = 1.0
  #rf[(nz-1)/3] = 1.0
  rf[(nz-1)/2] = 1.0
  #rf[2*(nz-1)/3] = 1.0
  #rf[(nz-2)] = 1.0
  vl = fillfloat(2.0,nz)
  tzl = WTutils().tzVlog(vl,sz)
  #sy = WTutils().syntheticSeismogram(fp,rf,tzl,st,sz,"ricker")
  sy = WTutils().syntheticSeismogram(fp,rf,tzl,st,sz,"morlet")
  #sy = applyPhaseRot(92,sy)
  ss = Sampling(len(sy),st.delta,tzl[0])
  pp = PlotPanel(xdxr)
  pp.addPoints(ss,sy)
  pp.setVLabel("Time (s)")
  #pp.setVLimits(0,0.1);
  frame = PlotFrame(pp)
  frame.setDefaultCloseOperation(eoc)
  frame.setSize(350,600)
  frame.setFontSizeForSlide(1.0,1.0)
  #frame.paintToPng(720,4.00,_pngDir+'swavelet.png')
  #frame.setFontSizeForPrint(8.0,240.0)
  #frame.paintToPng(720,7.00,_pngDir+'pwavelet.png')
  frame.setVisible(True)


"""
Plots teaser figure for CWP report 2012 of synthetic and trace before and after alignment
"""
def plotTeaser(tr,st,y,ys,sst,ss,tm=None,lims=None,paper=None,slides=None,paper2=None,paper3=None):
  if tm:
    p1 = PlotPanel(1,3,xdxr)
    p1.setHLabel(2,"tau")
    stau = Sampling(ss.count,ss.delta,tm[0])
    p1.addPoints(0,2,stau,tm)
  else: 
    p1 = PlotPanel(2,1)#,xdxr)
  l1 = p1.addPoints(0,0,st,tr)
  l1.setLineWidth(2)
  p1.setHLabel("Time (s)")
  #ss2 = Sampling(ss.count,ss.delta,sst.first)
  l2 = p1.addPoints(0,0,ss,y)
  l2.setStyle('r-')
  l2.setLineWidth(2)
  l3 = p1.addPoints(1,0,st,tr)
  l3.setLineWidth(2)
  l4 = p1.addPoints(1,0,sst,ys)
  l4.setStyle('r-')
  l4.setLineWidth(2)
  if lims:
    p1.setVLimits(0,-lims,lims)
    p1.setVLimits(1,-lims,lims)
    if tm: p1.setVLimits(2,-lims,lims)
  p1.setHLimits(0,sst.last+0.5)
  frame = PlotFrame(p1)
  frame.setSize(900,1500)
  frame.setDefaultCloseOperation(eoc);
  frame.setVisible(True)
  if paper:
    frame.setSize(1360,354)
    frame.setFontSizeForPrint(8.0,504.0)
    frame.paintToPng(720,7.00,_pngDir+paper+'.png')
  if paper2:
    frame.setSize(904,379)
    frame.setFontSizeForPrint(8.0,240.0)
    frame.paintToPng(720,3.33,_pngDir+paper2+'.png')
  if paper3:
    frame.setSize(904,379)
    frame.setFontSizeForPrint(8.0,312.0)
    frame.paintToPng(720,4.33,_pngDir+paper3+'.png')
  if slides:
    frame.setSize(500,900)
    frame.setFontSizeForSlide(1.0,1.0)
    frame.paintToPng(720,7.00,_pngDir+slides+'.png')
  frame.setVisible(True)


def plot3P(s1,s2,s3,x,title,cmap=gray,p=None,m=None,d=None,v=None,\
						clips1=None,clips2=None,ss=None,lims=None,slc=None,con=None,paper=None):
	svp = Viewer3P(s1,s2,s3,x)
	if v:
		svp.addPixels(v)
		svp.addColorBar("Velocity (km/s)")
		svp.setColorModel2(ColorMap.getJet(0.75))
		if clips2:
			svp.setClips2(clips2[0],clips2[1])
	if d:
		svp.addPixels(d)
		svp.addColorBar("Depth (km)")
		svp.setColorModel2(ColorMap.getJet(0.75))
		if clips2:
			svp.setClips2(clips2[0],clips2[1])
		if con:
			svp.addContours(s1,s2,s3,d)
	if m:
		svp.getPP().pixelsView23.tile.addTiledView(m)
	if p:
		#svp.addPnts23(p[0],p[1])
		svp.addPts(p)
	#svp.setTitle(title)
	#svp.setSize(1200,1000)
	if cmap:
		svp.setColorModel1(cmap)
	if lims:
		svp.setLimits1(lims[0],lims[1])
	if clips1:
		svp.setClips1(clips1[0],clips1[1])
	svp.setLabel1("Time (s)")
	svp.setLabel2("Crossline (km)")
	svp.setLabel3("Inline (km)")
	if slc:
		svp.setSlices(inro((slc[0]-s1.first)/s1.delta),inro(slc[1]/s2.delta),inro(slc[2]/s3.delta))
	#svp.setSlices(80,160,343)
	#svp.setSlices(343,160,80)
	#svp.setSlices(563,174,84)
	#svp.setSlices(575,228,84)
	svp.setSize(600,600)
	if paper:
		svp.setFontSizeForPrint(8.0,222.0)
		svp.paintToPng(720,3.08,_pngDir+paper+'.png')
	svp.show()

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
