package tp;

import java.io.File;
import wt.SeismicModel1D;
import dtw.DynamicWarpingWT;
import dtw.DynamicWarpingSWT;
import edu.mines.jtk.io.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.SimplePlot;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.SincInterp;
import edu.mines.jtk.dsp.SincInterp.Extrapolation;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Compute a single well tie for a Teapot Dome well.
 *
 * TODO
 *
 * @author Andrew Munoz, Colorado School of Mines
 * @version 12.19.13
 */

public class WellTie {
  
  public long id; // 12-digit log API
  public float x3f,x2f,x1f; // first samples for x,y,z; x,y relative to seismic 
  public float[] v,d,g,a,r; // null-sample corrected logs
  public float[] f; // intial synthetic seismogram
  public float[] h; // tied synthetic seismogram
  public float[] u; // single well-tie shifts
  public float[] v1; // updated velocity
  public float[] pv; // percentage difference in velocities in depth
  public float[] tz0,tz1; // tau(z) time-depth functions
  public Sampling sz; // sampling of logs
  public Sampling s1; // seismic time sampling
  public Sampling sf; // sampling of intial synthetic
  public Sampling sh; // sampling of tied synthetic
  public String wset; // well set location
  public float phase = 0.0f; // constant phase of synthetic seismogram
  
  // This value represents missing values in all well curve data.
  public static final float NULL_VALUE = -999.2500f;

  // Static correction velocity for the seismic survey (km/s)
  public static final float STATIC_VEL = 2.7432f; // 9000 ft/s
  
  //TODO list
  // 1. get logs and edit
  // 2. simple synthetic seismogram computation and fixing
  // 2. prop synthetic seismograms computation and fixing
  // 3. single-well ties and input parameters
  // 4. associate with multiple-well ties
  // 5. link to time-depth computations 

  /**
   * Reads in logs and corrects sample values
   * @param id log API number
   * @param dz log sampling interval
   * @param s1 sampling of the seismic 
   * @param wset file path to binary well log set
   */
  public WellTie(long id, double dz, Sampling s1, String wset) {
    this.id = id;
    this.s1 = s1;
    this.wset = wset;
    WellLog.Data wdata = WellLog.Data.readBinary(wset);
    WellLog logs = wdata.get(id);
    LogSamples ls = new LogSamples();
    ls.fixLogsForSynthetics(logs);
    v = ls.vl;
    d = ls.dl;
    g = ls.gl;
    despike(3);
    a = mul(v,d); // acoustic impedance
    r = reflectivity(a);
    x3f = ls.xl;
    x2f = ls.yl;
    x1f = ls.zl;
    sz = new Sampling(a.length,dz,x1f);
    _nz = a.length;
    _dz = (float)dz;
    _fz = x1f;
    _n1 = s1.getCount();
    _d1 = (float)s1.getDelta();
    _f1 = (float)s1.getFirst();
    _ds = _d1;
    maketz0(v);
  } 

  public WellTie() {}

  /**
   * Computes a single-well tie. 
   * Use smooth dynamic time warping to compute a well tie. Refer 
   * to DynamicWarpingWT class for extra details.
   * @param g  adjacent seismic trace(s) to the well
   * @param r1min lower bound on strain in time.
   * @param r1max upper bound on strain in time.
   * @param smin minimum time shift between synthetic seismogram and seismic 
   * @param smax maximum time shift between synthetic seismogram and seismic
	 * @param dr1 time smoothing parameter
   */
  public void computeSingleWellTie(
    float[] g, float r1min, float r1max, float smin, float smax) 
  {
    computeSingleWellTie(new float[][]{g},r1min,r1max,smin,smax,1.0f);
  }
  public void computeSingleWellTie(
    float[] g, float r1min, float r1max, float smin, float smax, float dr1) 
  {
    computeSingleWellTie(new float[][]{g},r1min,r1max,smin,smax,dr1);
  }
  public void computeSingleWellTie(
    float[][] g, float r1min, float r1max, float smin, float smax) 
  {
    computeSingleWellTie(g,r1min,r1max,smin,smax,1.0f);
  }
  public void computeSingleWellTie(
    float[][] g, float r1min, float r1max, 
    float smin, float smax, float dr1) 
  {
    Check.argument(f!=null,"no synthetic seismogram exists");
    int nl = g[0].length-_ns; 
    DynamicWarpingWT dw = new DynamicWarpingWT(smin,smax,s1);
	  dw.setStrainLimits(r1min,r1max);
	  dw.setSmoothing(dr1);
	  dw.setSyntheticSampling(sf);
    //dw.setInterpolationLinear();
    if (_fphase) {
      float[][] ps = findPhase(dw,f,g); 
      _php = ps[0];
      phase = ps[1][0];
      applyPhaseRotation(phase);
    }
    u = dw.findShifts(sf,f,s1,g);
    h = dw.applyShifts(f,u);
    sh = new Sampling(h.length,_ds,u[0]+_f1);
    update();
  }

  /**
   * Computes a synthetic seismogram.
   * This method uses vertically propagating plane waves to compute
   * synthetic seismograms in a attenuative medium (constant q). 
   * This method computes internal and free surface multiples.
   * @param dataDir the directory for storing computed synthetic seismograms
   * @param fp peak frequency of synthetic seismogram
   * @param q  constant attenuation factor
   * @param ph constant phase rotation
   * @param snorm the sigma for normalization
   */
  public void makePropagatorSeismogram(
    String dataDir, double fp, double q)
  {
    makePropagatorSeismogram(dataDir,fp,q,phase);
  }
  public void makePropagatorSeismogram(
    String dataDir, double fp, double q, double ph) 
  {
    makePropagatorSeismogram(dataDir,fp,q,ph,_snorm);
  }
  public void makePropagatorSeismogram(
    String dataDir, double fp, double q, double ph, double snorm) 
  {
    phase = (float)ph;
    _fp = fp;
    _snorm = snorm;
    double d1 = s1.getDelta();
    int nz = _nz;
	  int n1 = _n1;
    int ns = ince((tz0[nz-1])/d1);
    float dz = _dz;
    float ds = _ds;
    float fz = _fz;
    float f1 = _f1;
    float fs = _fs;
    float fref = (float)(10000.0*d1); // reference frequency (Hz)
    float r1 = 1.0f; // surface reflection coefficient
    float amp = 1.0f; // source amplitude
    _ns = ns; 
    sf = new Sampling(ns,d1,0);
    String fname = getfname(dataDir,id,q,fp,ns,fs,fref,r1,amp,_nomults,_norms);
    if (exists(fname) && !_mseis && dataDir!=null) {
      System.out.println("Reading seismogram...");
      f = new float[ns];
      readSeismograms(fname);
  } else {
      System.out.println("Making seismogram...");
      SeismicModel1D ssm = new SeismicModel1D();
      ssm.turnOffMultiples(_nomults);
      // Land Dynamite
		  //ssm.setSourceType(SeismicModel1D.SourceType.ANISOTROPIC); 
      // Airgun, Vib, Marine dyn
			ssm.setSourceType(SeismicModel1D.SourceType.ISOTROPIC);
      // Hydrophone
			//ssm.setSensorType(SeismicModel1D.SensorType.PRESSURE);
      // Geophone
			ssm.setSensorType(SeismicModel1D.SensorType.VELOCITY);
		  ssm.setSurfaceReflectionCoefficient(r1);
		  ssm.setLayers(sz,v,d,(float)q);
		  ssm.addSensor(0.0);
		  ssm.addSource(0.0,amp);
      ssm.setRickerWavelet(fp);
			ssm.removeSourceSignature(true);
      float[][] rs = ssm.makeSeismograms(ns,d1,fref);
      f = copy(rs[0]);
      writeSeismograms(fname);
    }
    applyPhaseRotation(phase);
    cutSeismogram();
    normalize();
  }

  /**
   * Computes a simple synthetic seismogram.
   * The method computes a sum of shifted wavelets 
   * (ricker with peak frequency fp) scaled by the 
   * reflection coefficients r found from the well logs. 
   * @param ph constant phase rotation
   * @param fp peak frequency of synthetic seismogram
   */
  public void makeSimpleSeismogram(double fp) {
    makeSimpleSeismogram(fp,phase);
  }
  public void makeSimpleSeismogram(double fp, double ph) {
    makeSimpleSeismogram(fp,ph,_snorm);
  }
  public void makeSimpleSeismogram(double fp, double ph, double snorm) {
    phase = (float)ph;
    _fp = fp;
    _snorm = snorm;
    double d1 = s1.getDelta();
    int nz = _nz;
	  int n1 = _n1;
    int ns = ince((tz0[nz-1]-_fs)/d1);
    float dz = _dz;
    float ds = _ds;
    float fz = _fz;
    float f1 = _f1;
    float fs = _fs;
    _ns = ns;
    sf = new Sampling(ns,d1,fs);
    float[] ti = rampfloat(fs,ds,ns);
    f = new float[ns];
    double sig = 10.0/(fp*PI);
    for (int iz=0; iz<nz; ++iz) {
      for (int it=0; it<ns; ++it) {
        if ((tz0[iz]-sig)<=ti[it] && ti[it]<=(tz0[iz]+sig)) {
          f[it] += r[iz]*ricker(fp,tz0[iz]-ti[it]);
        } 
      }
    }
    applyPhaseRotation(phase);
    normalize();
  }

  /**
   * Integrates a velocity log to get a time-depth function tau_0(z).
   * @param v The velocity log
   * @param sz The log sampling
   * @return toz The time-depth function t(z)
   */
  public void maketz0(float[] v) {
    int nz = _nz;
    float fz = _fz;
    float dz = _dz;
		int nm = 50;
		MedianFinder mf = new MedianFinder(nm);
		float vm = mf.findMedian(copy(nm,v));
    float to = 2f*fz*(1f/vm);//+1f/STATIC_VEL)*0.5f;//-1f/pf;
    if (_stat)
      to = 2f*fz*(1f/vm+1f/STATIC_VEL)*0.5f;//-1f/pf;
    //float to = 2f*fz*(1f/v[0]); // for simp comparison
    tz0 = new float[nz];
    float dz2 = dz*2f;
    float tc = 0.0f;
		tz0[0] = tc = to;
    for (int i=1; i<nz; ++i) 
      tz0[i] = tc = tc+dz2/v[i];
    _fs = to;
  }

  public void setNormalizationOff() {
    _norms=false;
  }
  public void setNormalizationOn() {
    // default
    _norms=true;
  }
  public void makeNewSeismograms() {
    _mseis=true;
  }
  public void setMultiplesOn() {
    // default
    _nomults=false;
  }
  public void setMultiplesOff() {
    _nomults=true;
  }
  public void setCutsOff() {
    _cuts=false;
  }
  public void setCutsOn() {
    // default
    _cuts=true;
  }
  public void findPhase() {
    _fphase=true;
  }
  public float[] getPhasePlot() {
    return _php;
  }

  public void setStatic(boolean stat) {
    _stat  = stat; 
  }


////////////////////////////////////////////////////////////////////////////////
// private

  private int _n1,_nz,_ns;
  private float _dz,_d1,_ds,_fz,_f1,_fs;
  private double _fp; // peak frequency
  private int _ifi=0; // initial synthetic sample
  private double _snorm=100; // sigma for normalization
  private boolean _cuts=true; // cut the synthetic seismogram
  private boolean _norms=true; // normalize the synthetic seismograms
  private boolean _mseis=false; // make new synthetic seismograms 
  private boolean _nomults=false; // computes multiples
  private boolean _fphase=false; // computes optimal phase rotation
  private boolean _stat=false;
  private float[] _php; // phase rotation plot


  public void update(float[] u1, float fs) {
  // Updates time-depths, velocity, and percentage difference in velocities.
    u = u1;
    int nu = u1.length;
    int nz = _nz;
    float[] x1 = rampfloat(fs,_d1,nu);
    float[] x2 = rampfloat(_fs,_d1,nu);
    float[] ttau = add(x1,u1);
    // compute tz1 from t(tau=tau(z))
    tz1 = new float[nz];
    pv = new float[nz];
	  CubicInterpolator ci1 = new CubicInterpolator(
	  	CubicInterpolator.Method.LINEAR,nu,x2,ttau);
	  ci1.interpolate(nz,tz0,tz1);
    v1 = div(2f,bderivative(tz1,_dz));
    pv = mul(div(sub(v,v1),v1),-100f);
  }
  private void update() {
    update(u,_f1);
  }

  /**
   * Performs a local RMS normalization on the synthetic seismogram.
   * This method applies an exponential smoother on a 
   * window of samples, s, around a sample in the sequence
   * and normalizes that value by the smoothed samples.
   * @param x The input sequence
   * @param sig The window size, sigma, to smooth samples
   * @return y The normalized sequence
   */
   private void normalize() {
     if (!_norms) return;
     int n = f.length;
     float[] x = copy(f);
     float[] y = new float[n];
     float[] yy = new float[n];
     float[] xx = mul(x,x);
		 RecursiveExponentialFilter es = new RecursiveExponentialFilter(_snorm);
		 //es.setEdges(RecursiveExponentialFilter.Edges.OUTPUT_ZERO_VALUE);
     es.apply(xx,yy);
     for (int i=0; i<n; ++i)
       y[i] = x[i]/sqrt(yy[i]);
     f = copy(y);
   }

  private void applyPhaseRotation(float ph) {
    int n = f.length;
    float phr = (FLT_PI/180f)*ph;
    float[] fr = copy(f); 
    float[] fi = new float[n];
    HilbertTransformFilter hb = new HilbertTransformFilter();
    hb.apply(n,fr,fi);
    f = add(mul(cos(phr),fr),mul(sin(phr),fi));
  }
  
  private void cutSeismogram() {
     if (!_cuts) return;
     // remove zero data above first reflection point
     // and update seismogram and sampling
     //int ifs = infl(((2f*_fz/v[0])-1f/_fp)/_ds); 
     double d1 = s1.getDelta();
     int ifs = infl(_fs/d1); 
     int ns = _ns-ifs;
     float[] f0 = copy(f);
     f = copy(ns,ifs,f0);
     _ns = ns;
     _ifi = 0;
     sf = new Sampling(ns,s1.getDelta(),_fs);
  }

  private static float[][] findPhase(
    DynamicWarpingWT dw, float[] f, float[][] g) 
  {
    int nf = f.length;
    int ng = g[0].length;
    int nh = 361;
    float[] mdm = new float[nh]; 
    HilbertTransformFilter hb = new HilbertTransformFilter();
    float pd = 1; int ihm=0;
    for (int ih=0; ih<nh; ++ih) {
      float[] fr = applyPhaseRot(hb,ih,f);
      float dm = (float)dw.findPhaseError(fr,g);
      mdm[ih] = dm;
      if (ih==0) pd=dm;
      if (dm<pd) {
        ihm = ih;
        pd  = dm;
      }
    }
    float fihm = (float)ihm;
    System.out.println("Optimal phase is "+Float.toString(fihm)+" degrees");
    return new float[][]{mdm,new float[]{fihm}};
  }
  private static float[] applyPhaseRot(
    HilbertTransformFilter hb, int ph, float[] fr) 
  {
    int n = fr.length;
    float phr = (FLT_PI/180f)*ph;
    float[] fi = new float[n];
    hb.apply(n,fr,fi);
    return add(mul(cos(phr),fr),mul(sin(phr),fi));
  }


  private static float[] reflectivity(float[] a) {
    int n = a.length; 
    float[] r = new float[n];
    for (int i=1; i<n; ++i) 
      r[i] = (a[i-1]-a[i])/(a[i-1]+a[i]);
    r[0] = r[1];
    return r;
  }

  private static float[] bderivative(float[] x, float dx) {
    // numerical derivative using a backwards difference 
    int n = x.length;
    float[] y = new float[n];
    for (int i=1; i<n; ++i)
      y[i] = (x[i]-x[i-1])/(dx);
    y[0] = y[1];
    return y;
  }
  private static float[] fderivative(float[] x, float dx) {
    // numerical derivative using a forward difference 
    int n = x.length;
		int nm = n-1;
    float[] y = new float[n];
    for (int i=0; i<nm; ++i)
      y[i] = (x[i+1]-x[i])/(dx);
    y[nm] = y[nm-1];
    return y;
  }
  private static float[] cderivative(float[] x, float dx) {
    // numerical derivative using a centered difference 
    int n = x.length;
		int nm = n-1;
    float[] y = new float[n];
    for (int i=1; i<nm; ++i)
      y[i] = (x[i+1]-x[i-1])/(2*dx);
    y[nm] = y[nm-1];
		y[0] = y[1];
    return y;
  }

  private void readSeismograms(String fname) {
    try {
      ArrayInputStream ais = new ArrayInputStream(fname);
      ais.readFloats(f);
      ais.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }
  private void writeSeismograms(String fname) {
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fname);
      aos.writeFloats(f);
      aos.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  private float ricker(double fpeak, double time) {
    double x = PI*fpeak*time;
    double xx = x*x;
    return (float)((1.0-2.0*xx)*exp(-xx));
  }
  private float morlet(double fpeak, double x) {
    double c = 1/sqrt(2*PI);
    return (float)(c*exp(-x*x/2)*cos(2*PI*fpeak*x));
  }

	private int inro(float x) {
		return round(x);
	}
	private int inro(double x) {
		return (int)round(x);
	}
	private int infl(float x) {
		return (int)floor(x);
	}
	private int infl(double x) {
		return (int)floor(x);
	}
	private int ince(float x) {
		return (int)ceil(x);
	}
	private int ince(double x) {
		return (int)ceil(x);
	}
  private static float[] floats(double[] x) {
    int n1 = x.length;
    float[] y = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      y[i1] = (float)x[i1];
    }
    return y;
  }
  private static float[] floats(int[] x) {
    int n1 = x.length;
    float[] y = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      y[i1] = (float)x[i1];
    }
    return y;
  }


  private static String getfname(
    String dataDir, long id, double q, double fp, int ns, 
    float fs, float fref, float r1, float amp, boolean mlt, boolean nrm) 
  {
    String s1 = Long.toString(id);
    String s2 = Double.toString(q);
    String s3 = Double.toString(fp);
    String s4 = Float.toString((float)ns);
    String s5 = Float.toString(fs);
    String s6 = Float.toString(fref);
    String s7 = Float.toString(r1);
    String s8 = Float.toString(amp);
    String s9 = Boolean.toString(mlt);
    String s10 = Boolean.toString(nrm);
    String d = "_";
    String fname = dataDir+s1+d+s2+d+s3+d+s4+d+s5+d+
                           s6+d+s7+d+s8+d+s9+d+s10+".dat";
    return fname;
  }

  private static boolean exists(String fname) {
    return new File(fname).exists();
  }

  // Thanks to Dave Hale for this median despike filter
  public void despike(int nmed) {
    for (int imed=0; imed<nmed; ++imed) {
      despike(v);
      despike(d);
      despike(g);
    }
  }
  private void despike(float[] f) {
    if (f==null) 
      return;
    int n = f.length;
    for (int i=1; i<n-1; ++i) {
      float fa = f[i-1];
      float fb = f[i  ];
      float fc = f[i+1];
      if (fa!=NULL_VALUE && fb!=NULL_VALUE && fc!=NULL_VALUE)
        f[i] = med3(fa,fb,fc);
    }
  }
  private float med3(float a, float b, float c) {
    return a<b ? 
           (b<c ? b : (a<c ? c : a)) : 
           (b>c ? b : (a>c ? c : a));
  }


};
