package ani;

import java.util.*;
import java.awt.Color;
import edu.mines.jtk.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;

/** 
 * Anisotropic 2D ray tracing code in horizontal layers for a single source.
 * Uses Christoffel equation to solve for phase and group velocities and 
 * angles and traveltimes for P-SV rays in a layered anisotropic medium. 
 *
 * @author Andrew Munoz, CSM
 * @version 11.09.13
 */

public class RayTrace {

	/**
	 * Adds a layer to the model.
	 * @param z layer depth
	 * @param vp p-wave velocity
	 * @param vs s-wave velocity
	 * @param ep epsilon
	 * @param de delta
	 */
	public void addLayer(
		double z, double vp, double vs, double ep, double de, boolean reflect) 
	{
		addLayer((float)z,(float)vp,(float)vs,(float)ep,(float)de,reflect);
	}
	public void addLayer(
		float z, float vp, float vs, float ep, float de, boolean reflect) 
	{
		Layer layer = _layers.get(z);
		if (layer!=null) {
			layer.vp = vp;
			layer.vs = vs;
			layer.ep = ep;
			layer.de = de;
			layer.reflect = reflect;
	} else {
			layer = new Layer(z,vp,vs,ep,de,reflect);
			_layers.put(z,layer);
		}
		if (reflect) _reflector=true;
	}

	/**
	 * Compute rays.
	 * Compute rays with one ray per angle. Rays symmetric around vertical axis
	 * @param x  maximum source-receiver offset
	 * @param pa initial sampling of phase angles (in degrees)
	 */
	public void computeRays(Sampling pa) {
		int npa   = pa.getCount();
		float dpa = (float)pa.getDelta();
		float fpa = (float)pa.getFirst();
		int nr = npa; // number of rays equal to number of angles 
		addRays(nr); // make empty rays
		int nl = _layers.size();
		nl = (_reflector)?nl*2:nl;
		for (int ia=0; ia<npa; ++ia) {
			
			float a = toRadians(fpa+dpa*ia);

			// Get the ray correponding to a phase takeoff angle
			Ray pray,sray;
			pray = _prays.get(ia);
			sray = _srays.get(ia);

			// Initialize first layer
			Float zs = _layers.firstKey();
			Layer ly = _layers.get(zs);

			float vp,vs;
			float ap = a;
			float as = a;
		
			float tp=0f; float ts=0f; float xp=0f; float xs=0f;
			// Compute rays until all reflect and reach top
			for (int l=0; l<nl; ++l) {
				
				vp = ly.vp;
				vs = ly.vs;
				float[]   pv = phaseV(ap,as,vp,vs,ly.ep,ly.de);
				float[][] gv = groupV(ap,as,vp,vs,ly.ep,ly.de);
				xp += 2*tan(gv[0][1])*ly.z;
				xs += 2*tan(gv[1][1])*ly.z;
				tp += 2*sqrt(tan(gv[0][1])*ly.z*tan(gv[0][1])*ly.z+ly.z*ly.z)/gv[0][0];
				ts += 2*sqrt(tan(gv[1][1])*ly.z*tan(gv[1][1])*ly.z+ly.z*ly.z)/gv[1][0];
				pray.set(l,
				   tp,					  			// traveltime
					 pv[0], 	 	 	  			// phase velocity
					 gv[0][0], 						// group velocity
					 toDegrees(ap),			  // phase angle
					 toDegrees(gv[0][1]),	// group angle
					 xp);									// ray offset
				sray.set(l,
				   ts,					  			// traveltime
					 pv[1], 	 	 	  			// phase velocity
					 gv[1][0], 						// group velocity
					 toDegrees(as),			  // phase angle
					 toDegrees(gv[1][1]),	// group angle
					 xs);									// ray offset

				if (ly.reflect) {
					int lo=l;
					for (int rl=l+1; rl<nl*2 && lo>=0; ++rl) {
						tp+=pray.t[lo];
						ts+=sray.t[lo];
						pray.set(rl,tp,pray.pv[lo],pray.gv[lo],pray.pa[lo],pray.ga[lo],pray.x[lo]);
						sray.set(rl,ts,sray.pv[lo],sray.gv[lo],sray.pa[lo],sray.ga[lo],sray.x[lo]);
						--lo;
					}
					break;
				} else {
					if (l==nl-1) break;
				}
				float pr1 = sin(ap)/pv[0];
				float sr1 = sin(as)/pv[1];
				zs = _layers.higherKey(zs);
				ly = _layers.get(zs);
				// Compute theta2 by trying v2 for a variety of angles 
				// and minimize to solve snell's law
				float dpm = Float.MAX_VALUE; float dsm = Float.MAX_VALUE;
				int m = 10; float map2 = ap; float mas2 = as; 
				int npad2 = -npa/2; int npaf = npa;
				for (int i2=1; i2<=10*m; i2+=m) {
					float da = 1.0f/i2;
					float tmap2 = map2;
					for (int i1=npad2; i1<npaf; ++i1) {
						float ap2 = i1*da+toDegrees(map2);
						if (ap2>90.0 || ap2<-90.0) continue;
						ap2 = toRadians(ap2);
						float pv2 = phaseVP(ap2,ly.vp,ly.vs,ly.ep,ly.de);
						float pr2 = sin(ap2)/pv2;
						float dp = abs(pr1-pr2);
						if (dp<dpm) {
							dpm = dp;
							tmap2 = ap2;
						}
					}
					map2 = tmap2;
					//dpm = Float.MAX_VALUE;
				}
				for (int i2=1; i2<=10*m; i2+=m) {
					float da = 1.0f/i2;
					float tmas2 = mas2;
					for (int i1=npad2; i1<npaf; ++i1) {
						float as2 = i1*da+toDegrees(mas2);
						if (as2>90.0 || as2<-90.0) continue;
						as2 = toRadians(as2);
						float sv2 = phaseVS(as2,ly.vp,ly.vs,ly.ep,ly.de);
						float sr2 = sin(as2)/sv2;
						float ds = abs(sr1-sr2);
						if (ds<dsm) {
							dsm = ds;
							tmas2 = as2;
						}
					}
					mas2 = tmas2;
					//dsm = Float.MAX_VALUE;
				}
				ap = map2;
				as = mas2;
			}
		}
	}

	/**
	 * @return t[wave][ray][layer]
	 */
	public float[][][] getTravelTimes() {
		int nr = _prays.size();
		float[][][] t = new float[2][nr][];
		for (int ir=0; ir<nr; ++ir) {
			t[0][ir] = _prays.get(ir).t;
			t[1][ir] = _srays.get(ir).t;
		}
		return t;
	}

	/**
	 * @return pv[wave][ray][layer]
	 */
	public float[][][] getPhaseVelocities() {
		int nr = _prays.size();
		float[][][] pv = new float[2][nr][];
		for (int ir=0; ir<nr; ++ir) {
			pv[0][ir] = _prays.get(ir).pv;
			pv[1][ir] = _srays.get(ir).pv;
		}
		return pv;
	}

	/**
	 * @return gv[wave][ray][layer]
	 */
	public float[][][] getGroupVelocities() {
		int nr = _prays.size();
		float[][][] gv = new float[2][nr][];
		for (int ir=0; ir<nr; ++ir) {
			gv[0][ir] = _prays.get(ir).gv;
			gv[1][ir] = _srays.get(ir).gv;
		}
		return gv;
	}

	/**
	 * @return ga[wave][ray][layer]
	 */
	public float[][][] getGroupAngles() {
		int nr = _prays.size();
		float[][][] ga = new float[2][nr][];
		for (int ir=0; ir<nr; ++ir) {
			ga[0][ir] = _prays.get(ir).ga;
			ga[1][ir] = _srays.get(ir).ga;
		}
		return ga;
	}
	/**
	 * @return pa[wave][ray][layer]
	 */
	public float[][][] getPhaseAngles() {
		int nr = _prays.size();
		float[][][] pa = new float[2][nr][];
		for (int ir=0; ir<nr; ++ir) {
			pa[0][ir] = _prays.get(ir).pa;
			pa[1][ir] = _srays.get(ir).pa;
		}
		return pa;
	}
	/**
	 * @return t[wave][layer][ray]
	 */
	public float[][][] getTravelTimesRays() {
		int nr = _prays.size();
		int nz = _layers.size();
		nz = (_reflector)?nz*2:nz;
		float[][][] t = new float[2][nz][nr];
		for (int iz=0; iz<nz; ++iz) {
			for (int ir=0; ir<nr; ++ir) {
				t[0][iz][ir] = _prays.get(ir).t[iz];
				t[1][iz][ir] = _srays.get(ir).t[iz];
		}}
		return t;
	}
	/**
	 * @return ga[wave][layer][ray]
	 */
	public float[][][] getGroupAnglesRays() {
		int nr = _prays.size();
		int nz = _layers.size();
		nz = (_reflector)?nz*2:nz;
		float[][][] ga = new float[2][nz][nr];
		for (int iz=0; iz<nz; ++iz) {
			for (int ir=0; ir<nr; ++ir) {
				ga[0][iz][ir] = _prays.get(ir).ga[iz];
				ga[1][iz][ir] = _srays.get(ir).ga[iz];
		}}
		return ga;
	}
	/**
	 * @return pa[wave][layer][ray]
	 */
	public float[][][] getPhaseAnglesRays() {
		int nr = _prays.size();
		int nz = _layers.size();
		nz = (_reflector)?nz*2:nz;
		float[][][] pa = new float[2][nz][nr];
		for (int iz=0; iz<nz; ++iz) {
			for (int ir=0; ir<nr; ++ir) {
				pa[0][iz][ir] = _prays.get(ir).pa[iz];
				pa[1][iz][ir] = _srays.get(ir).pa[iz];
		}}
		return pa;
	}

	// Returns the total, x, z, and angle components of the group velocity
	private float[][] groupV(float ap, float as, float vp, float vs, float ep, float de) {
		return new float[][]{groupVP(ap,vp,vs,ep,de),groupVS(as,vp,vs,ep,de)};
	}
	private float[] groupVP(float a, float vp, float vs, float ep, float de) {
		return groupVC(a,vp,vs,ep,de,1f);
	}
	private float[] groupVS(float a, float vp, float vs, float ep, float de) {
		return groupVC(a,vp,vs,ep,de,-1f);
	}
	private float[] groupVC(float a, float vp, float vs, float ep, float de, float w) {
		float f = 1f-(vs*vs)/(vp*vp);
		float v = phaseVC(a,vp,vs,ep,de,w);
		float dv = derivativePhaseVelocity(a,vp,ep,de,w,f,v);
		float vg  = v*sqrt(1f+(dv/v)*(dv/v));
		float vgx = v*sin(a)+dv*cos(a);
		float vgz = v*cos(a)-dv*sin(a);
		return new float[]{vg,atan(vgx/vgz)};
	}

	// Returns the total, x, and z component of the phase velocity
	private float[] phaseV(float ap, float as, float vp, float vs, float ep, float de) {
		return new float[]{phaseVP(ap,vp,vs,ep,de),phaseVS(as,vp,vs,ep,de)};
	}
	private float phaseVP(float a, float vp, float vs, float ep, float de) {
		return phaseVC(a,vp,vs,ep,de,1f); 
	}
	private float phaseVS(float a, float vp, float vs, float ep, float de) {
		return phaseVC(a,vp,vs,ep,de,-1f); 
	}
	private float phaseVC(float a, float vp, float vs, float ep, float de, float w) {
		float f = 1f-(vs*vs)/(vp*vp);
		float sins = sin(a)*sin(a);
		float coss = cos(a)*cos(a);
		float v = vp*sqrt(1f+ep*sins-f/2f+
							w*f/2f*sqrt(1f + 4f*ep*ep*sins*sins/(f*f)+(4f*sins/f)*(2f*de*coss-ep*cos(2*a))));
		return v; 
	}
	// Analytic derivative of the phase velocity 
	private float derivativePhaseVelocity(
		float a, float vp, float ep, float de, float w, float f, float v) 
	{
		float sins = sin(a)*sin(a);
		float coss = cos(a)*cos(a);
		//float den1 = sqrt(pow(1f+2f*ep*sins/f,2f)-2f*(ep-de)*sin(2f*a)*sin(2f*a)/f);
		//float num = w*(8f*ep*sin(a)*cos(a)*(2f*ep*sins/f+1f) - 
		//						   8f*(ep-de)*sin(2*a)*cos(2*a))/4f*den1 + 2f*ep*sin(a)*cos(a);
		//float den = 2f*sqrt(w*0.5f*f*den1 - 0.5f*f + ep*sins+1f);
		float p1 = (vp*vp*0.5f)*ep*sin(2f*a)/v; 
		float p2 = (vp*vp/v)*(0.5f*sin(2f*a)*(2f*de*coss-ep*cos(2f*a)) 
								+ sins*(ep*sin(2f*a) - de*sin(2f*a)) 
								+ (2f/f)*ep*ep*sins*sin(a)*cos(a));
		float p3 = sqrt(1f+(4f/f)*sins*(2f*de*coss-ep*cos(2f*a)) 
										+ (2f*ep/f)*(2f*ep/f)*sins*sins);
		return p1 + w*p2/p3;
	}

////////////////////////////////////////////////////////////////////////////////
// private

	private TreeMap<Integer,Ray> _prays = new TreeMap<Integer,Ray>(); // P-Rays
	private TreeMap<Integer,Ray> _srays = new TreeMap<Integer,Ray>(); // SV-Rays
	private TreeMap<Float,Layer> _layers = new TreeMap<Float,Layer>(); // Layers
	private float _nr; // number of rays
	private boolean _reflector=false;

	private class Ray {
		// All quantities are per layer traveled and for a P or SV waves
		int     r;   // ray number
		int     w;   // type of wave (1 for p, -1 for sv)
		float[] t;   // traveltimes at layer interface
		float[] x; 	 // ray offset
		float[] pa;  // Phase Angles
		float[] ga;  // Group Angles
		float[] pv;  // Phase Velocity Vector magnitude
		float[] gv;  // Group Velocity Vector magnitude

		public Ray(
			int r, int w, float[] t, float[] pv, 
			float[] gv, float[] pa, float[] ga, float[] x)
		{
			this.r  = r; 		
			this.w  = w;
			this.t  = t; 		
			this.x  = x; 		
			this.pv = pv; 
			this.gv = gv; 
			this.pa = pa; 	
			this.ga = ga; 	
		}		
		public void set(
			int l, float t, float pv, float gv, float pa, float ga, float x)
		{
			this.t[l]  = t; 
			this.x[l]  = x;
			this.pv[l] = pv;
			this.gv[l] = gv;
			this.pa[l] = pa; 
			this.ga[l] = ga;
		}
	}

	private class Layer {
		float  z; // Thickness
		float vp; // P-wave Velocity
		float vs; // S-wave Velocity
		float ep; // Epsilon
		float de; // Delta
		boolean reflect; // Reflect from bottom of layer if true

		public Layer(float z, float vp, float vs, float ep, float de) {
			this(z,vp,vs,ep,de,false);	
		}
		public Layer(float z, float vp, float vs, float ep, float de, boolean reflect) {
			this.z = z;
			this.vp = vp;
			this.vs = vs;
			this.ep = ep;
			this.de = de;
			this.reflect = reflect;
		}
	}

	private void addRays(int nr) {
		int m = (_reflector)?2:1;
		int nl = _layers.size()*m; // include upgoing and downgoing 
		for (int ir=0; ir<nr; ++ir) {
			float[] t  = new float[nl];
			float[] x  = new float[nl];
			float[] pa = new float[nl];
			float[] ga = new float[nl];
			float[] pv = new float[nl];
			float[] gv = new float[nl];
			Ray pray = _prays.get(ir);
			pray = new Ray(ir, 1,t,pv,gv,pa,ga,x);
			_prays.put(ir,pray);
		}
		for (int ir=0; ir<nr; ++ir) {
			float[] t  = new float[nl];
			float[] x  = new float[nl];
			float[] pa = new float[nl];
			float[] ga = new float[nl];
			float[] pv = new float[nl];
			float[] gv = new float[nl];
			Ray sray = _srays.get(ir);
			sray = new Ray(ir,-1,t,pv,gv,pa,ga,x);
			_srays.put(ir,sray);
		}
	}

	public void addTraveltimesP(PlotPanel pp, float xmax) {
		int nr = _prays.size();	
		float[][] c = new float[2][nr];
		for (int ir=0; ir<nr; ++ir) {
			float[] t = _prays.get(ir).t;
			float[] x = _prays.get(ir).x;
			makeTTCoordinates(ir,t,x,c);
		}
		PointsView pv = pp.addPoints(c[1],c[0]);
		pp.setHLimits(-xmax,xmax);	
		pp.setVLimits(min(c[1]),1.2);	
	}	
	public void addTraveltimesS(PlotPanel pp, float xmax) {
		int nr = _prays.size();	
		float[][] c = new float[2][nr];
		for (int ir=0; ir<nr; ++ir) {
			float[] t = _srays.get(ir).t;
			float[] x = _srays.get(ir).x;
			makeTTCoordinates(ir,t,x,c);
		}
		PointsView pv = pp.addPoints(c[1],c[0]);
		pp.setHLimits(-xmax,xmax);	
		pp.setVLimits(min(c[1]),3.5);	
	}
	private void makeTTCoordinates(
		int i, float[] t, float[] x, float[][] c) 
	{
		int nz = t.length;
		int nzm = nz-1;
		c[1][i] = t[nzm];
		c[0][i] = x[nzm];
	}
	public void addRaysPlotP(PlotPanel pp, float xmax) {
		int nr = _prays.size();	
		float[] z = addLayersPlot(pp,xmax);
		for (int ir=0; ir<nr; ++ir) {
			float[]   x = _prays.get(ir).x;
			float[][] c = makeRayCoordinates(z,x);
			PointsView pv = pp.addPoints(c[0],c[1]);
		}
		pp.setHLimits(-xmax,xmax);
	}
	public void addRaysPlotS(PlotPanel pp, float xmax) {
		int nr = _srays.size();	
		float[] z = addLayersPlot(pp,xmax);
		for (int ir=0; ir<nr; ++ir) {
			float[]   x = _srays.get(ir).x;
			float[][] c = makeRayCoordinates(z,x);
			PointsView pv = pp.addPoints(c[0],c[1]);
		}
		pp.setHLimits(-xmax,xmax);
	}
	private float[][] makeRayCoordinates(float[] z, float[] x) {
		int nz = z.length;	
		float[] zn = new float[nz];
		float[] xn = new float[nz];
		int ne = (_reflector)?nz-1:nz;
		for (int iz=1; iz<ne; ++iz) {
			zn[iz] = z[iz];
			xn[iz] = x[iz-1];
		}
		return new float[][]{zn,xn};
	}
	public float[] addLayersPlot(PlotPanel pp, float x) {
		int nz = _layers.size();	
		float z=0; int i=1;
		int nz2 = (_reflector)?nz*2+2:nz+1;
		float[] za = new float[nz2];
		for(Map.Entry<Float,Layer> k : _layers.entrySet()) {
			z += k.getKey().floatValue();	
			PointsView pv = pp.addPoints(new float[]{z,z},new float[]{-x,x});
			pv.setLineColor(Color.RED);
			pp.setVLimits(0f,z+0.1f);
			za[i] = z; ++i;
		}
		if (_reflector) {
			for (int iz=nz2/2,ip=0; iz<nz2; ++iz,++ip)
				za[iz] = za[ip];
		}
		return za;
	}
		
////////////////////////////////////////////////////////////////////////////////
// testing only

	public static void main(String[] args) {
    test1();
	}

  public static void test1() {
		float xmax = 2.0f;
		float da = 1.0f;
		float fa = -90f;
		float la =  90f;
		// Layer properties
		float z1  = 0.5f; float z2  = 1.0f; // layer thickness 
		float vp1 = 3.0f; float vp2 = 3.0f; 
		float vs1 = 1.8f; float vs2 = 1.5f; 
		float ep1 = 0.2f; float ep2 = -0.2f; 
		float de1 = 0.1f; float de2 = 0.3f; 
		Sampling sa = new Sampling((int)((la-fa)/da)+1,da,fa);
		RayTrace rt = new RayTrace();
		rt.addLayer(z1,vp1,vs1,ep1,de1,false);
		rt.addLayer(z2,vp2,vs2,ep2,de2,false);
		rt.computeRays(sa);
		float[][][] t  = rt.getTravelTimes();
		float[][][] ga = rt.getGroupAngles();
		float[][][] pa = rt.getPhaseAngles();
		float[][][] pv = rt.getPhaseVelocities();
		float[][][] gv = rt.getGroupVelocities();
		float[][][] tr  = rt.getTravelTimesRays();
		float[][][] gar = rt.getGroupAnglesRays();
		float[][][] par = rt.getPhaseAnglesRays();
		// Plot
		int l = 0;
		plotSlownessSurface(par[0][l],getSlownesses(pv[0],l),"p-phase slowness (s/km)");
		plotSlownessSurface(par[1][l],getSlownesses(pv[1],l),"s-phase slowness (s/km)");
		plotSlownessSurface(gar[0][l],getSlownesses(gv[0],l),"p-group slowness (s/km)");
		plotSlownessSurface(gar[1][l],getSlownesses(gv[1],l),"s-group slowness (s/km)");
		plotRaysP(rt,xmax,"p-rays-nep");
		plotRaysS(rt,xmax,"s-rays-nep");
		plotTraveltimesP(rt,xmax,"p-tt-nep");
		plotTraveltimesS(rt,xmax,"s-tt-nep");
		//plot(floats(sz.getValues()),t[0][0],"depth (km)","phase angle (deg)");
	}

	private static void plotTraveltimesP(RayTrace rt, float xmax, String png) {
		PlotPanel pp = new PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT);
		rt.addTraveltimesP(pp,xmax);
		pp.setTitle("P-rays");
		pp.setVLabel("traveltime (s)");
		pp.setHLabel("half-offset (km)");
		PlotFrame pf = new PlotFrame(pp);
		pf.setVisible(true);
		pf.paintToPng(720,3.33,png+".png");
	}

	private static void plotTraveltimesS(RayTrace rt, float xmax, String png) {
		PlotPanel pp = new PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT);
		rt.addTraveltimesS(pp,xmax);
		pp.setTitle("SV-rays");
		pp.setVLabel("traveltime (s)");
		pp.setHLabel("half-offset (km)");
		PlotFrame pf = new PlotFrame(pp);
		pf.setVisible(true);
		pf.paintToPng(720,3.33,png+".png");
	}

	private static void plotRaysP(RayTrace rt, float xmax, String png) {
		PlotPanel pp = new PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT);
		rt.addRaysPlotP(pp,xmax);
		pp.setTitle("P-rays");
		pp.setVLabel("depth (km)");
		pp.setHLabel("distance (km)");
		PlotFrame pf = new PlotFrame(pp);
		pf.setVisible(true);
		pf.paintToPng(720,3.33,png+".png");
	}
	private static void plotRaysS(RayTrace rt, float xmax, String png) {
		PlotPanel pp = new PlotPanel(PlotPanel.Orientation.X1DOWN_X2RIGHT);
		rt.addRaysPlotS(pp,xmax);
		pp.setTitle("SV-rays");
		pp.setVLabel("depth (km)");
		pp.setHLabel("distance (km)");
		PlotFrame pf = new PlotFrame(pp);
		pf.setVisible(true);
		pf.paintToPng(720,3.33,png+".png");
	}


	private static void plot(float[] x, float[] y, String vaxis, String haxis) {
	  SimplePlot si = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
		si.addPoints(x,y);
		si.setSize(500,800);
		si.setHLabel(haxis);
		si.setVLabel(vaxis);
	}

	private static void plotSlownessSurface(float[] a, float[] p, String haxis) {
	  SimplePlot si = new SimplePlot();
		si.addPoints(a,p);
		si.setSize(500,500);
		si.setHLabel("angle (deg)");
		si.setVLabel(haxis);
	}

	private static float[] floats(double[] x) {
		int n = x.length;
		float[] y = new float[n];
		for (int i=0; i<n; ++i) 
			y[i] = (float)x[i];
		return y;
	}

	private static float[] getSlownesses(float[][] pv, int il) {
		int nr = pv.length;
		float[] p = new float[nr];
		for(int ir=0; ir<nr; ++ir) 
			p[ir] = 1f/pv[ir][il];
		return p;
	}


};
