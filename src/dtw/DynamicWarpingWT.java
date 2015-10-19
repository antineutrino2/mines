package dtw;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.awt.ColorMap;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Dynamic warping for single-well ties and multiple-well ties. 
 * This computes automatic well ties using smooth dynamic time warping and 
 * smooth dynamic image warping. For single-well ties, a seismic trace is 
 * warped to the synthetic seismogram because the time shifts can only be
 * measured with respect to the range of times in the synthetic seismogram.
 * For multiple well ties, a seismic image is warped to a synthetic image,
 * and again, the time shifts are computed using inverse linear interpolation.
 * <p>
 * Smooth dynamic image warping description written by Dr. Dave Hale:
 * <p>
 * Constraints are placed on strains, the rates at which shifts change in any
 * dimension. For example, strain r1 = du/d1 is the derivative of shift
 * u(x1,x2,...) with respect to the x1 coordinate, and is constrained to lie
 * between lower and upper bounds r1min and r1max. Default bounds are 
 * r1min = -1.0 and r1max = 1.0.
 * <p>
 * In many applications of warping, strains derived from estimated shifts may
 * be an important by-product. However, when computing shifts, only a finite
 * number of strains are permitted, and this quantization may yield
 * unrealistic estimates of strain. To address this problem, this dynamic
 * warping provides control over the sampling of strains, in the form of a
 * smoothness value. The number of strains sampled is proportional to this
 * value.
 * <p>
 * Smoothness values represent approximate intervals for a quasi-uniform
 * subsampling grid on which shifts are computed. For example, for a time
 * sampling interval of 4 ms, one might specify a smoothness of 200 ms, so
 * that shifts are computed on a grid that is 50 times more coarse than the 4
 * ms sampling grid. However, strains on this coarse grid would be sampled 50
 * times more finely. After initially computing shifts on the coarse grid,
 * shifts are interpolated onto the finer grid.
 * <p>
 * Thanks to Stefan Compton for help and design of many of these methods.
 * @author Andrew Munoz, Colorado School of Mines
 * @version 11.16.13
 */

public class DynamicWarpingWT {

	/**
   * Constructs a 1D dynamic warping.
   * @param smin minimum time shift between synthetic seismogram and seismic 
   * @param smax maximum time shift between synthetic seismogram and seismic
   * @param s1 time sampling of seismic trace.
   */
  public DynamicWarpingWT(float smin, float smax, Sampling s1) {
		this(smin,smax,s1,null,null);	
	}

	/**
   * Constructs a 2D dynamic warping.
   * @param smin minimum time shift between synthetic image and seismic 
   * @param smax maximum time shift between synthetic image and seismic
   * @param s1 time sampling of seismic image
   * @param s2 distance sampling of seismic image
   */
  public DynamicWarpingWT(float smin, float smax, Sampling s1, Sampling s2) {
		this(smin,smax,s1,s2,null);
  }

 	/**
   * Constructs a 3D dynamic warping.
   * @param smin minimum time shift between synthetic image and seismic 
   * @param smax maximum time shift between synthetic image and seismic
   * @param s1 time sampling of seismic image
   * @param s2 distance sampling of seismic image
   * @param s3 distance sampling of seismic image
   */
  public DynamicWarpingWT(
		float smin, float smax, Sampling s1, Sampling s2, Sampling s3) 
	{
		double ds = s1.getDelta();
		int ismin = (int) ceil(smin/ds);
    int ismax = (int)floor(smax/ds);
		int nl = 1+ismax-ismin;
    _ismin = ismin;
    _ismax = ismax;
    _su = new Sampling(nl,ds,ismin*ds);
		_nl = nl;
		_ss = s1;
   	_s1 = s1; 
   	_s2 = s2; 
   	_s3 = s3; 
		_r1min = -1.0f;
    _r1max =  1.0f;
    _r2min = -1.0f;
    _r2max =  1.0f;
    _r3min = -1.0f;
    _r3max =  1.0f;
		_dg1 = 1;
		_dg2 = 1;
		_dg3 = 1;
		_si = new SincInterp();
		_si.setExtrapolation(SincInterp.Extrapolation.CONSTANT);
		_ui = new InterpShifts();
    _ui.setMethod(InterpShifts.Method.SPLINE);
  }

	/**
	 * Sets a new time sampling for synthetics.
	 * Initially assumed to be the same as the seismic trace.
	 * @param ss time sampling of synthetic seismogram.
	 */
 	public void setSyntheticSampling(Sampling ss) {
		double ds = _su.getDelta();
		double fs = _su.getFirst();
		Check.argument(ss.getDelta()==ds,
			"synthetic sampling interval must equal seismic sampling interval");
		_ss = ss;
		_nl = _s1.getCount()-_ss.getCount()+1;
		_su = new Sampling(_nl,ds,fs);
	}

  /**
   * Sets bounds on strains.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1min lower bound on strain in time.
   * @param r1max upper bound on strain in time.
   */
  public void setStrainLimits(double r1min, double r1max) {
    setStrainLimits(r1min,r1max,-1.0,1.0,-1.0,1.0);
  }

  /**
   * Sets bounds on strains. 
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1min lower bound on strain in time.
   * @param r1max upper bound on strain in time.
   * @param r2min lower bound on strain in distance.
   * @param r2max upper bound on strain in distance.
   */
  public void setStrainLimits(
    double r1min, double r1max,
    double r2min, double r2max)
  {
    setStrainLimits(r1min,r1max,r2min,r2max,-1.0,1.0);
  }

  /**
   * Sets bounds on strains.
   * Default lower and upper bounds are -1.0 and 1.0, respectively.
   * @param r1min lower bound on strain in time.
   * @param r1max upper bound on strain in time.
   * @param r2min lower bound on strain in distance.
   * @param r2max upper bound on strain in distance.
   * @param r3min lower bound on strain in distance.
   * @param r3max upper bound on strain in distance.
   */
  public void setStrainLimits(
    double r1min, double r1max,
    double r2min, double r2max,
    double r3min, double r3max)
  {
    _r1min = (float)r1min; _r1max = (float)r1max;
    _r2min = (float)r2min; _r2max = (float)r2max;
    _r3min = (float)r3min; _r3max = (float)r3max;
  }

	/**
	 * Sets the smoothing parameter.
	 * These parameters determine the time grid spacing in 
	 * smooth dynamic time warping. Default is 1.0.
	 * @param dr1 time smoothing parameter
	 */
	public void setSmoothing(double dr1) {
		setSmoothing(dr1,1.0,1.0);
	}

	/**
	 * Sets the smoothing parameters.
	 * These parameters determine the grid spacings for each dimension in 
	 * smooth dynamic image warping. Defaults are 1.0.
	 * @param dr1 time smoothing parameter
	 * @param dr2 distance smoothing parameter
	 */
	public void setSmoothing(double dr1, double dr2) {
		setSmoothing(dr1,dr2,1.0);
	}

	/**
	 * Sets the smoothing parameters.
	 * These parameters determine the grid spacings for each dimension in 
	 * smooth dynamic image warping. Defaults are 1.0.
	 * @param dr1 time smoothing parameter
	 * @param dr2 distance smoothing parameter
	 * @param dr3 distance smoothing parameter
	 */
	public void setSmoothing(double dr1, double dr2, double dr3) {
    Check.argument(dr1<=1.0,"dr1<=1.0");
    Check.argument(dr2<=1.0,"dr2<=1.0");
    Check.argument(dr3<=1.0,"dr3<=1.0");
		_dg1 = (int)ceil(1/dr1);
		_dg2 = (int)ceil(1/dr2);
		_dg3 = (int)ceil(1/dr3);
	}

  /**
   * Set the method used for interpolating shifts. 
	 * LINEAR interpolation ensures the constraints are followed exactly, but
	 * results in blocky time shifts. SPLINE is the best option for well ties
	 * even though it may slightly violate the constraints because the time 
	 * shifts and the derivative of the time shifts are smoothly varying, which
	 * is the most resonable assumption.
   */
  public void setInterpolationLinear() {
    _ui.setMethod(InterpShifts.Method.LINEAR);
  }
  public void setInterpolationMonotonic() {
    _ui.setMethod(InterpShifts.Method.SPLINE);
  }
  public void setInterpolationSpline() { // default
    _ui.setMethod(InterpShifts.Method.SPLINE);
  }


	/**
   * Get the shift sampling.
   * @return the shift sampling.
   */
  public Sampling getShiftSampling() {
    return _su;
  }

	/**
	 * Find shifts for synthetic seismograms and seismic traces.
	 * Shifts are computed to align the seismic trace g(g1) with the
	 * synthetic seismogram f(tau) where shifts g1 are sparse samples of g(t). 
	 * The sparse shifts are interpolated back onto the fine grid of t.
	 * Then, since the synthetic seismogram is typically aligned to the 
	 * seismic trace, the shifts need to be inverse linear interpolated 
	 * to get the desired alignment.  
	 * <p>
	 * The minimum interval of g1 samples is dg1.
	 * @param f synthetic seismogram.
	 * @param g nearby seismic trace(s)
   * @return shifts that align the seismic trace to the synthetic seismogram
	 */  
	public float[] findShifts(
    Sampling sf, float[] f, 
    Sampling sg, float[][] g) 
  {
		int n1 = f.length;
		int n2 = g[0].length;
		int nt = g.length;
		int[] g1 = Subsample.subsample(n1,_dg1);
		int nc = g1.length;
		//float[][] e = computeErrors1(f,g);
		//Sampling su = _su;
		//Sampling se = _s1;
    //int ns = su.getCount();
    //int ne = se.getCount();
    //float[][] e = new float[ne][ns];
		//computeErrors(sf,f,sg,g[0],e);
		float[][] e = computeErrors1(f,g);
    float[][] es = smoothErrors(e,_r1min,_r1max,g1);
    float[][][] dm = accumulateForward(es,g1,_r1min,_r1max);
		float[] uc = backtrackReverse(dm[0],dm[1]);
    float fs = (float)_su.getFirst();
		return sub(_ui.interpolate(_ss,g1,uc),fs);
	}
	public float[] findShifts(Sampling sf, float[] f, Sampling sg, float[] g) {
		return findShifts(sf,f,sg,new float[][]{g}); // one trace
	}
	
  /**
   * Find shifts for 2D images. 
	 * Shifts are computed for the alignemnt of seismic image g(g1,g2) to 
	 * synthetic image f(x,tau), where g1 and g2 are sparse samples of g(x,t). 
	 * The sparse shifts are then interpolated back to the fine grid. The 
   * sparse samples g1 are constant for each g2. 
   * <p>
   * The minimum interval of g1 and g2 samples is dg1 and
   * dg2 respectively. These intervals are variable in 
   * order to include the first and last indices of the fine grid.
   * @param f synthetic image
   * @param g seismic image
   * @return shifts that align the seismic image to the synthetic image
   */
  public float[][] findShifts(
		final Sampling sf, final float[][] f, 
    final Sampling sg, final float[][] g)
  {
    final int dg1 = _dg1;
    final int dg2 = _dg2;
		final int n1 = f[0].length;
		final int n2 = f.length;
    final int[] g1 = Subsample.subsample(n1,dg1);
		final int[] g2 = Subsample.subsample(n2,dg2);
    print("g1:"); dump(g1); print("g2"); dump(g2);
		_g1 = g1; _g2 = g2;
    
    Stopwatch s = new Stopwatch();
    s.start();
    print("Smoothing 1st dimension...");
    final float[][][] es1 = smoothErrors1(sf,f,sg,g,_r1min,_r1max,g1);
    print("Finished 1st dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es1);
    
    s.restart();
    print("Smoothing 2nd dimension...");
    final float[][][] es = smoothErrors2(es1,_r2min,_r2max,g2);
    print("Finished 2nd dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es);
    
    final int ng2 = es.length;
    final float[][] u = new float[ng2][];
    Parallel.loop(ng2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][][] dm = accumulateForward(es[i2],g1,_r1min,_r1max);
      u[i2] = backtrackReverse(dm[0],dm[1]);
    }});
    return _ui.interpolate(_s1,_s2,g1,g2,u);
  }
  
  /**
   * Find shifts for 3D images. 
	 * Shifts are computed for the alignemnt of seismic image 
	 * g(g1,g2,g3) to synthetic image f(x,y,tau), where g1 and 
	 * g2 are sparse samples of g(x,y,t). 
	 * The sparse shifts are then interpolated back to the fine grid. The 
   * sparse samples g1 are constant for each g2 and g3. 
   * <p>
   * The minimum interval of g1, g2, and g3 samples is dg1, dg2, and 
   * dg3 respectively. These intervals are variable in 
   * order to include the first and last indices of the fine grid.
   * @param f synthetic image
   * @param g seismic image
   * @return shifts that align the seismic image to the synthetic image
   */
  public float[][][] findShifts(
		final Sampling sf, final float[][][] f, 
    final Sampling sg, final float[][][] g)
  {
    final int dg1 = _dg1;
    final int dg2 = _dg2;
    final int dg3 = _dg3;
		final int n1 = f[0][0].length;
		final int n2 = f[0].length;
		final int n3 = f.length;
    final int[] g1 = Subsample.subsample(n1,dg1);
    final int[] g2 = Subsample.subsample(n2,dg2);
    final int[] g3 = Subsample.subsample(n3,dg3);
    print("g1:"); dump(g1); print("g2"); dump(g2); print("g3"); dump(g3);
		_g1 = g1; _g2 = g2; _g3 = g3;
    
    Stopwatch s = new Stopwatch();
    s.start();
    print("Smoothing 1st dimension...");
    final float[][][][] es1 = smoothErrors1(sf,f,sg,g,_r1min,_r1max,g1);
    print("Finished 1st dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es1);
    
    s.restart();
    print("Smoothing 2nd dimension...");
    final float[][][][] es12 = smoothErrors2(es1,_r2min,_r2max,g2);
    print("Finished 2nd dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es12);
    
    s.restart();
    print("Smoothing 3rd dimension...");
    final float[][][][] es = smoothErrors3(es12,_r3min,_r3max,g3);
    print("Finished 3rd dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es);
    
    final int ng3 = es.length;
    final int ng2 = es[0].length;
    final float[][][] u = new float[ng3][ng2][];
    Parallel.loop(ng3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<ng2; i2++) {
        float[][][] dm = accumulateForward(es[i3][i2],g1,_r1min,_r1max);
        u[i3][i2] = backtrackReverse(dm[0],dm[1]);
      }
    }});
		return _ui.interpolate(_s1,_s2,_s3,g1,g2,g3,u);
  }

	/**
	 * Applies shifts to synthetic seismogram.
	 * Uses inverse linear interpolation to compute the shifts needed to
	 * warp the synthetic seismogram to the seismic trace.
	 * @param f synthetic seismogram 
	 * @param u shifts that align the seismic trace to synthetic seismogram 
	 * @return warped synthetic seismogram 
	 */
  public float[] applyShifts(float[] f, float[] u, Sampling ss) {
    int ns = f.length;
		float ds = (float)ss.getDelta();
		float fs = (float)ss.getFirst();
		float[] s = inverseInterpolation(ss,u);
		int nn = s.length;
    float[] h = new float[nn];
    for (int i1=0; i1<nn; ++i1) 
      h[i1] = _si.interpolate(ns,ds,fs,f,s[i1]);
    return h;
  }
  public float[] applyShifts(float[] f, float[] u) {
    return applyShifts(f,u,_ss);
  }

	/**
	 * Applies shifts to 2D synthetic image.
	 * Uses inverse linear interpolation to compute the shifts needed to
	 * warp the synthetic image to the seismic image.
	 * @param f synthetic image
	 * @param u shifts that align the seismic image to synthetic image
	 * @return warped synthetic image
	 */
  public float[][] applyShifts(final float[][] f, final float[][] u) {
		final int n1 = f[0].length;
		final int n2 = f.length;
		final float d1 = (float)_s1.getDelta();
		final float f1 = (float)_s1.getFirst();
    final float[][] hf = new float[n2][n1];
		final float[][] s = inverseInterpolation(_s1,u);
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        hf[i2][i1] = _si.interpolate(n1,d1,f1,f[i2],s[i2][i1]);
      }
    }});
    return hf;
  }

	/**
	 * Applies shifts to 3D synthetic image.
	 * Uses inverse linear interpolation to compute the shifts needed to
	 * warp the synthetic image to the seismic image.
	 * @param f synthetic image
	 * @param u shifts that align the seismic image to synthetic image
	 * @return warped synthetic image
	 */
  public float[][][] applyShifts(final float[][][] f, final float[][][] u) {
		final int n1 = f[0][0].length;
		final int n2 = f[0].length;
		final int n3 = f.length;
		final float d1 = (float)_s1.getDelta();
		final float f1 = (float)_s1.getFirst();
    final float[][][] hf = new float[n3][n2][n1];
		final float[][][] s = inverseInterpolation(_s1,u);
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; i2++) {
        for (int i1=0; i1<n1; i1++) {
          hf[i3][i2][i1] = _si.interpolate(n1,d1,f1,f[i3][i2],s[i3][i2][i1]);
        }
      }
    }});
    return hf;
  }

		/**
	 * Applies shifts to 3D synthetic image.
	 * Uses inverse linear interpolation to compute the shifts needed to
	 * warp the synthetic image to the seismic image.
	 * @param f synthetic image
	 * @param u shifts that align the seismic image to synthetic image
	 * @return warped synthetic image
	 */
  public float[][][] applyShiftsTest(final float[][][] g, final float[][][] u) {
		final int n1 = g[0][0].length;
		final int n2 = g[0].length;
		final int n3 = g.length;
		final float d1 = (float)_s1.getDelta();
		final float f1 = (float)_s1.getFirst();
    final float[][][] hf = new float[n3][n2][n1];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; i2++) {
        for (int i1=0; i1<n1; i1++) {
          hf[i3][i2][i1] = _si.interpolate(n1,d1,f1,g[i3][i2],
              						_s1.getValue(i1)+u[i3][i2][i1]);
        }
      }
    }});
    return hf;
  }

	/**
	 * Applies shifts to 2D synthetic image.
	 * Uses inverse linear interpolation to compute the shifts needed to
	 * warp the synthetic image to the seismic image.
	 * @param f synthetic image
	 * @param u shifts that align the seismic image to synthetic image
	 * @return warped synthetic image
	 */
  public float[][] applyShiftsX(final float[][] f, final float[][] u) {
		final int n1 = f[0].length;
		final int n2 = f.length;
		final float d1 = (float)_s1.getDelta();
		final float f1 = (float)_s1.getFirst();
    final float[][] hf = new float[n2][n1];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        hf[i2][i1] = _si.interpolate(n1,d1,f1,f[i2],_s1.getValue(i1)+u[i2][i1]);
      }
    }});
    return hf;
  }

	/**
	 * Applies shifts to synthetic seismogram.
	 * Uses inverse linear interpolation to compute the shifts needed to
	 * warp the synthetic seismogram to the seismic trace.
	 * @param f synthetic seismogram 
	 * @param u shifts that align the seismic trace to synthetic seismogram 
	 * @return warped synthetic seismogram 
	 */
  public float[] applyShiftsX(float[] f, float[] u, Sampling ss) {
    int ns = f.length;
    int nsm = ns-1;
		float ds = (float)_s1.getDelta();
		float fs = (float)_s1.getFirst();
    float[] h = new float[ns];
    for (int i1=0; i1<ns; ++i1) 
      h[i1] = _si.interpolate(ns,ds,fs,f,_s1.getValue(i1)+u[i1]);
    return h;
  }

	 /** 
		* Computes warped synthetic seismograms from trace in synthetic image.
		* @param u  time shifts at the synthetic (in samples)
		* @param ss synthetic seismogram samplings
		* @param x  initial synthetic seismogram
		* @return h[ns] synthetic seismogram, and new first sampling fs
		*/
	public float[][] getWarpedSyntheticSeismogram(
		Sampling ss, float[] u, float[] x) 
	{
		int nt = _s1.getCount();
		int ns = ss.getCount();
		float ft = (float)_s1.getFirst();
		float dt = (float)_s1.getDelta();
		float fs = (float)ss.getFirst();
		float ls = (float)ss.getLast();
		int fio  = round((fs-ft)/dt);
		int lio  = round((ls-ft)/dt);
		int no   = lio-fio+1;
		if (no+fio>nt) no -= no+fio-nt;
		float[] us = copy(no,fio,u);
		Sampling sn = new Sampling(no,dt,fs);
		float[] s = inverseInterpolation(sn,us);
		int n2 = s.length;
		float[] h = new float[n2];
		for (int i=0; i<n2; ++i)
			h[i] = _si.interpolate(nt,dt,fs+s[0],x,s[i]);
		return new float[][]{h,new float[]{s[0]}};
	}

  public float[][][] accumulateForward(
      float[][] e, int[] g, float rmin, float rmax)
  {
    float[][] d = like(e);
    float[][] m = like(e);
    accumulateFromSparse(1,e,d,m,g,rmin,rmax);
    return new float[][][]{d,m};
  }
  
  public static float[][][] accumulateForwardSparse(
      float[][] e, float rmin, float rmax, int[] g)
  {
    int nl = e[0].length;
    int ng = g.length;
    float[][] d = new float[ng][nl];
    float[][] m = new float[ng][nl];
    accumulateSparse(1,rmin,rmax,g,e,d,m);
    return new float[][][]{d,m};
  }
  
  public static float[][][] accumulateReverseSparse(
      float[][] e, float rmin, float rmax, int[] g)
  {
    int nl = e[0].length;
    int ng = g.length;
    float[][] d = new float[ng][nl];
    float[][] m = new float[ng][nl];
    accumulateSparse(-1,rmin,rmax,g,e,d,m);
    return new float[][][]{d,m};
  }
  
  public float[] backtrackReverse(float[][] d, float[][] m) {
    float[] u = new float[d.length];
    backtrack(-1,_su,d,m,u);
    return u;
  }
  
  public float[] backtrackForward(float[][] d, float[][] m) {
    float[] u = new float[d.length];
    backtrack(1,_su,d,m,u);
    return u;
  }

	/**
   * Returns normalized alignment errors for single well ties.
   * @param f synthetic seismogram
   * @param g seismic traces
   * @return alignment errors.
   */	
	public float[][] computeErrors1(float[] f, float[][] g) {
    int n1 = f.length;
    int n2 = g[0].length;
		int nt = g.length;
    int nl = _nl;
    //nl -= n1+1;
    //nl = (nl>1)?nl:1;
		//if (nl==1) System.out.println("nl = 1");
    int tmin = _ismin+
      (int)round((_ss.getFirst()-_s1.getFirst())/_ss.getDelta())-1;
    tmin = (tmin>0)?tmin:0;
    int nm = n2-tmin-n1+1; // maximum number of lags
		nl = (nl<nm)?nl:nm;
    updatenl(nl);
    float[][] e = new float[n1][nl];
		for (int it=0; it<nt; ++it) 
    	for (int i1=0; i1<n1; ++i1) 
     	  for (int il=0,i2=tmin+i1+il; il<nl; ++il,++i2)
        	e[i1][il] += error(f[i1],g[it][i2]);
    normalizeErrors(e);
    return e;
  }

  /**
   * Returns normalized alignment errors for single well ties.
   * @param f synthetic seismogram
   * @param g seismic traces
   * @return alignment errors.
   */	
  public float[][] computeErrors1(
    Sampling sf, float[] f, 
    Sampling sg, float[][] g)
  {
    int nf = f.length;
    int ng = g[0].length;
    int nt = g.length;
		Sampling su = _su;
		Sampling se = _s1;
    int ns = su.getCount();
    int ne = se.getCount();
		float ff = (float)sf.getFirst();
		float fg = (float)sg.getFirst();
		float dg = (float)sg.getDelta();
		float fe = (float)se.getFirst();
		float de = (float)se.getDelta();
		float fu = (float)su.getFirst();
		float du = (float)su.getDelta();
    float[][] e = new float[ne][ns];
    for (int it=0; it<nt; ++it) {
      float[] fi = new float[ne];
      float[] gi = new float[ne];
      _si.interpolate(sf,f,se,fi);
      for (int is=0; is<ns; ++is) {
        _si.interpolate(ng,dg,fg,g[it],ne,de,fe+su.getValue(is),gi);
        for (int ie=0; ie<ne; ++ie)
          e[ie][is] += error(fi[ie],gi[ie]);
      }
    }
    int nl = _nl;
    int nm = ng-nf+1; // maximum number of lags
		nl = (nl<nm)?nl:nm;
    int tmin = _ismin+round((ff-fg)/dg);
    tmin = (tmin>0)?tmin:0;
    float[][] ec = copy(nl,ng,tmin,0,e);
    return ec;
  }

  /**
   * Normalizes values to be in range [0,1].
   * @param e input/output array.
   */
  public static void normalizeErrors(float[][] e) {
    int nl = e[0].length;
    int n1 = e.length;
    float emin = e[0][0];
    float emax = e[0][0];
    for (int i1=0; i1<n1; ++i1) {
      for (int il=0; il<nl; ++il) {
        float ei = e[i1][il];
        if (ei<emin) emin = ei;
        if (ei>emax) emax = ei;
      }
    }
    shiftAndScale(emin,emax,e);
  }

  /**
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
   */
  public static void normalizeErrors(float[][][] e) {
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final float[][][] ef = e;
    MinMax mm = Parallel.reduce(n2,new Parallel.ReduceInt<MinMax>() {
    public MinMax compute(int i2) {
      float emin =  Float.MAX_VALUE;
      float emax = -Float.MAX_VALUE;
      for (int i1=0; i1<n1; ++i1) {
        for (int il=0; il<nl; ++il) {
          float ei = ef[i2][i1][il];
          if (ei<emin) emin = ei;
          if (ei>emax) emax = ei;
        }
      }
      return new MinMax(emin,emax);
    }
    public MinMax combine(MinMax mm1, MinMax mm2) {
      return new MinMax(min(mm1.emin,mm2.emin),max(mm1.emax,mm2.emax));
    }});
    shiftAndScale(mm.emin,mm.emax,e);
  }
  
  /**
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
   */
  public static void normalizeErrors(float[][][][] e) {
    final int nl = e[0][0][0].length;
    final int n1 = e[0][0].length;
    final int n2 = e[0].length;
    final int n3 = e.length;
    final float[][][][] ef = e;
    MinMax mm = Parallel.reduce(n3,new Parallel.ReduceInt<MinMax>() {
      public MinMax compute(int i3) {
        float emin =  Float.MAX_VALUE;
        float emax = -Float.MAX_VALUE;
        for (int i2=0; i2<n2; ++i2) {
          for (int i1=0; i1<n1; ++i1) {
            for (int il=0; il<nl; ++il) {
              float ei = ef[i3][i2][i1][il];
              if (ei<emin) emin = ei;
              if (ei>emax) emax = ei;
            }
          }  
        }
        return new MinMax(emin,emax);
      }
      public MinMax combine(MinMax mm1, MinMax mm2) {
        return new MinMax(min(mm1.emin,mm2.emin),max(mm1.emax,mm2.emax));
      }});
    shiftAndScale(mm.emin,mm.emax,e);
  }

	
	public int[] getG1() {
		return _g1;
	}

	public int[] getG2() {
		return _g2;
	}

	public int[] getG3() {
		return _g3;
	}

  public double findPhaseError(float[] f, float[][] g) {
    int n1 = f.length;
  	int[] g1 = Subsample.subsample(n1,_dg1);
		int nc = g1.length;
    float[][] e = computeErrors1(f,g);
    float[][] es = smoothErrors(e,_r1min,_r1max,g1);
    float[][][] dm = accumulateForward(es,g1,_r1min,_r1max);
    int ld = dm[0].length-1;
    double dmin = min(dm[0][ld]);
    return dmin;
  }



  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _nl; // number of lags
	private int[] _g1; // x1 grid samples
	private int[] _g2; // x2 grid samples
	private int[] _g3; // x3 grid samples
	private Sampling _su; // shift samplings
	private Sampling _ss; // synthetic seismogram samplings
	private Sampling _s1; // seismic samplings for x1
	private Sampling _s2; // seismic samplings for x2
	private Sampling _s3; // seismic samplings for x3
	private float _r1min,_r1max; // strain limits in x1 direction
	private float _r2min,_r2max; // strain limits in x2 direction
	private float _r3min,_r3max; // strain limits in x3 direction
	private int _dg1; // grid interval in x1 direction
	private int _dg2; // grid interval in x2 direction
	private int _dg3; // grid interval in x3 direction
  private int _ismin,_ismax; // min and max time shifts in samples
  private float _epow = 2f; // exponent used for alignment errors |f-g|^e
  private SincInterp _si; // for warping with non-integer shifts
  private InterpShifts _ui; // for interpolating non-uniform shifts

  private void updatenl(int nl) {
    _nl = nl;
    _su = new Sampling(nl,_su.getDelta(),_su.getFirst());
  }

  private float error(float f, float g) {
    return pow(abs(f-g),_epow);
  }

  private void computeErrors(
		Sampling sf, float[] f, 
		Sampling sg, float[] g, float[][] e)
  {
    int ns = e[0].length;
    int ne = e.length;
    int nf = f.length;
    int ng = g.length;
		Sampling su = _su;
		Sampling se = _s1;
		float fg = (float)sg.getFirst();
		float dg = (float)sg.getDelta();
		float fe = (float)se.getFirst();
		float de = (float)se.getDelta();
		float fu = (float)su.getFirst();
		float du = (float)su.getDelta();
    float[] fi = new float[ne];
    float[] gi = new float[ne];
    _si.interpolate(sf,f,se,fi);
    for (int is=0; is<ns; ++is) {
      _si.interpolate(ng,dg,fg,g,ne,de,fe+su.getValue(is),gi);
      for (int ie=0; ie<ne; ++ie)
        e[ie][is] = error(fi[ie],gi[ie]);
    }
  }  

  private void computeErrors(
		Sampling sf, float[] f, 
		Sampling sg, float[][] g, float[][] e)
  {
    int ns = e[0].length;
    int ne = e.length;
    int nf = f.length;
    int ng = g[0].length;
    int nt = g.length;
		Sampling su = _su;
		Sampling se = _s1;
		float fg = (float)sg.getFirst();
		float dg = (float)sg.getDelta();
		float fe = (float)se.getFirst();
		float de = (float)se.getDelta();
		float fu = (float)su.getFirst();
		float du = (float)su.getDelta();
    for (int it=0; it<nt; ++it) {
      float[] fi = new float[ne];
      float[] gi = new float[ne];
      _si.interpolate(sf,f,se,fi);
      for (int is=0; is<ns; ++is) {
        _si.interpolate(ng,dg,fg,g[it],ne,de,fe+su.getValue(is),gi);
        for (int ie=0; ie<ne; ++ie)
          e[ie][is] += error(fi[ie],gi[ie]);
      }
    }
  }  

  /**
   * Returns smooth alignment errors on the sparse grid defined
   * by the indices of g. 
   * @param e 2D array of alignment errors.
   * @param rmin minimum slope.
   * @param rmax maximum slope.
   * @param g sparse grid indices.
   * @return smoothed alignment errors with size 
   *  [g.length][e[0].length].
   */
  private static float[][] smoothErrors(
      float[][] e, float rmin, float rmax, int[] g)      
  {
    int ng = g.length;
    int nel = e[0].length;
    float[][] ef = new float[ng][nel];
    float[][] er = new float[ng][nel];
    float[][] es = new float[ng][nel];
    accumulateSparse( 1,rmin,rmax,g,e,ef,null);
    accumulateSparse(-1,rmin,rmax,g,e,er,null);
    float scale = 1.0f/e.length;
    for (int i1=0; i1<ng; i1++) 
      for (int il=0; il<nel; il++) 
        es[i1][il] = scale*(ef[i1][il]+er[i1][il]-e[g[i1]][il]);
    return es;
  }
  
  /**
   * Returns alignment errors smoothed in the first dimension.
   * Returned errors are sparse in the first dimension, and
   * unchanged in the second dimension. Alignment errors are
   * computed on the fly.
   * @param f the synthetic image.
   * @param g the seismic image.
   * @param g1 first dimension sparse grid indices.
   * @return smoothed alignment errors with size 
   *  [n2][g1.length][_nl].
   */
  private float[][][] smoothErrors1(
		final Sampling sf, final float[][] f, 
		final Sampling sg, final float[][] g,
		final float r1min, final float r1max, final int[] g1)
  {
    final int ng1 = g1.length;
		final int n2 = f.length;
		final int n1 = f[0].length;
    final float[][][] es1 = new float[n2][ng1][_nl];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] e = new float[n1][_nl];
      computeErrors(sf,f[i2],sg,g[i2],e);
      es1[i2] = smoothErrors(e,r1min,r1max,g1);
    }});
    return es1;
  }
  
  /**
   * Returns alignment errors smoothed in the first dimension.
   * Returned errors are sparse in the first dimension, and
   * unchanged in the second dimension. Alignment errors are
   * computed on the fly.
   * @param f the synthetic image.
   * @param g the seismic image.
   * @param g1 first dimension sparse grid indices specified
   *  for all n2 (size[n2][])
   * @return smoothed alignment errors with size 
   *  [n2][g1.length][_nl].
   */
  private float[][][] smoothErrors1(
      final Sampling sf, final float[][] f, 
			final Sampling sg, final float[][] g,
			final float r1min, final float r1max, final int[][] g1)
  {
    final int ng1 = g1[0].length;
		final int n2 = f.length;
		final int n1 = f[0].length;
    final float[][][] es1 = new float[n2][ng1][_nl];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] e = new float[n1][_nl];
      computeErrors(sf,f[i2],sg,g[i2],e);
      es1[i2] = smoothErrors(e,r1min,r1max,g1[i2]);
    }});
    return es1;
  }
  
  /**
   * Returns alignment errors smoothed in the first dimension.
   * Returned errors are sparse in the first dimension, and
   * unchanged in the second and third dimension. Alignment 
   * errors are computed on the fly.
   * @param f the synthetic image.
   * @param g the seismic image.
   * @param g1 first dimension sparse grid indices.
   * @return smoothed alignment errors with size 
   *  [n3][n2][g1.length][_nl].
   */
  private float[][][][] smoothErrors1(
      final Sampling sf, final float[][][] f, 
			final Sampling sg, final float[][][] g,
			final float r1min, final float r1max, final int[] g1)
  {
    final int ng1 = g1.length;
		final int n3 = f.length;
		final int n2 = f[0].length;
		final int n1 = f[0][0].length;
    final float[][][][] es1 = new float[n3][n2][ng1][_nl];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; i2++) {
        float[][] e = new float[n1][_nl];
        computeErrors(sf,f[i3][i2],sg,g[i3][i2],e);
        es1[i3][i2] = smoothErrors(e,r1min,r1max,g1);  
      }
    }});
    return es1;
  }
  
  /**
   * Returns alignment errors smoothed in the first dimension.
   * Returned errors are sparse in the first dimension, and
   * unchanged in the second and third dimension. Alignment
   * errors are computed on the fly.
   * @param f the synthetic image.
   * @param g the seismic image.
   * @param g1 first dimension sparse grid indices specified
   *  for all n3 and n2 (size[n3][n2][]).
   * @return smoothed alignment errors with size 
   *  [n3][n2][g1.length][_nl].
   */
  private float[][][][] smoothErrors1(
      final Sampling sf, final float[][][] f, 
			final Sampling sg, final float[][][] g, 
			final float r1min, final float r1max, final int[][][] g1)
  {
    final int ng1 = g1[0][0].length;
		final int n3 = f.length;
		final int n2 = f[0].length;
		final int n1 = f[0][0].length;
    final float[][][][] es1 = new float[n3][n2][ng1][_nl];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; i2++) {
        float[][] e = new float[n1][_nl];
        computeErrors(sf,f[i3][i2],sg,g[i3][i2],e);
        es1[i3][i2] = smoothErrors(e,r1min,r1max,g1[i3][i2]);  
      }
    }});
    return es1;
  }
  
  /**
   * Returns alignment errors smoothed in the second dimension.
   * Returned errors are sparse in the second dimension, and
   * unchanged in the first dimension.
   * @param e alignment errors.
   * @param r2min minimum slope in the second dimension.
   * @param r2max maximum slope in the second dimension.
   * @param g2 second dimension sparse grid indices.
   * @return smoothed alignment errors with size 
   *  [g2.length][e[0].length][e[0][0].length].
   */
  private static float[][][] smoothErrors2(
      final float[][][] e, 
      final float r2min, final float r2max, final int[] g2)
  {
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final int ng2 = g2.length;
    final float[][][] es = new float[ng2][n1][nl]; // smoothed errors
    Parallel.loop(n1,new Parallel.LoopInt() {
    public void compute(int i1) {
      float[][]  e2 = new float[n2][nl]; // errors at index i1
      for (int i2=0; i2<n2; ++i2)
        e2[i2] = e[i2][i1];
      float[][] es2 = smoothErrors(e2,r2min,r2max,g2);
      for (int i2=0; i2<ng2; i2++)
        es[i2][i1] = es2[i2];
    }});
    return es;
  }

  /**
   * Returns alignment errors smoothed in the second dimension.
   * Returned errors are sparse in the second dimension, and
   * unchanged in the first dimension and third dimension.
   * @param e alignment errors.
   * @param r2min minimum slope in the second dimension.
   * @param r2max maximum slope in the second dimension.
   * @param g2 second dimension sparse grid indices.
   * @return smoothed alignment errors with size 
   *  [e.length][g2.length][e[0][0].length][e[0][0][0].length].
   */
  private static float[][][][] smoothErrors2(
      final float[][][][] e, 
      final float r2min, final float r2max, final int[] g2)
  {
    final int n3 = e.length;
    final float[][][][] es = new float[n3][][][]; // smoothed errors
    for (int i3=0; i3<n3; i3++)
      es[i3] = smoothErrors2(e[i3],r2min,r2max,g2);
    return es;
  }

  /**
   * Returns alignment errors smoothed in the third dimension.
   * Returned errors are sparse in the third dimension, and
   * unchanged in the first and second dimension.
   * @param e alignment errors.
   * @param r3min minimum slope in the third dimension.
   * @param r3max maximum slope in the third dimension.
   * @param g3 third dimension sparse grid indices.
   * @return smoothed alignment errors with size 
   *  [g3.length][e[0].length][e[0][0].length][e[0][0][0].length].
   */
  private static float[][][][] smoothErrors3(
      final float[][][][] e, 
      final float r3min, final float r3max, final int[] g3)
  {
    final int nl = e[0][0][0].length;
    final int n1 = e[0][0].length;
    final int n2 = e[0].length;
    final int n3 = e.length;
    final int ng3 = g3.length;
    final float[][][][] es = new float[ng3][n2][n1][nl]; // smoothed errors
    Parallel.loop(n1,new Parallel.LoopInt() {
    public void compute(int i1) {
      for (int i2=0; i2<n2; i2++) {
        float[][]  e3 = new float[n3][nl]; // smooth errors at index i1,i2
        for (int i3=0; i3<n3; i3++)
          e3[i3] = e[i3][i2][i1];
        float[][] es3 = smoothErrors(e3,r3min,r3max,g3);
        for (int i3=0; i3<ng3; i3++)
          es[i3][i2][i1] = es3[i3];
      }
    }});
    return es;
  }
  
  private static void accumulateFromSparse(
      int dir, float[][] e, float[][] d, float[][] m, int[] g,
      float rmin, float rmax)
  {
    int nl = e[0].length;
    int nlm1 = nl-1;
    int ni = e.length;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1; // beginning index
    int ie = (dir>0)?ni:-1;  // end index
    int is = (dir>0)?1:-1;   // stride
    int ic = (dir>0)?-1:1;   // contraint dir, forward=-lag, reverse=+lag
    float dmax = ni;
    int ii=ib;
    for (int il=0; il<nl; ++il)
      d[ii][il] = e[ii][il];
    ii+=is;
    for (; ii!=ie; ii+=is) {
      int ji = ii-is;
      int iemax = g[ii];
      int iemin = g[ji];
      int dg = abs(iemax-iemin); // sparse grid delta
      int kmin = (int)ceil( rmin*dg);
      int kmax = (int)floor(rmax*dg);
      for (int il=0; il<nl; ++il) {
        float dm = dmax;
        int mi = 0;
        for (int k=kmin; k<=kmax; k++) {
          int rk = k*ic;
          int ik = il+rk;
          if (ik<0 || ik>nlm1)
            continue;
          float dc = d[ji][ik];
          if (dc<dm) {
            dm = dc;
            mi = rk;
          }
        }
        m[ji][il] = mi;
        d[ii][il] = dm+e[ii][il];
      }
    }
  }
  
  private static void accumulateSparse(
      int dir, float rmin, float rmax, int[] g,
      float[][] e, float[][] d, float[][] m)
  {
    int nl   = e[0].length;
    int ng   = g.length;
    int ngm1 = ng-1;
    int ibg = (dir>0)?0:ngm1; // beginning index
    int ieg = (dir>0)?ng:-1;  // end index
    int is = (dir>0)?1:-1;   // stride
    int isp = ibg; // sparse grid index
    int ie = g[isp]; // error index
    
    // Initialize accumulation values
    for (int il=0; il<nl; ++il)
      d[isp][il] = e[ie][il];
    isp += is;
    // Loop over all sparse grid points.
    for (; isp!=ieg; isp+=is) {
      int ispm1 = isp-is; // previous sparse grid index
      ie = g[isp]; // new error index
      int je = g[ispm1]; // min error index for interpolation.
      int dg = ie-je; // sparse grid delta
      int kmin, kmax;
      if (dg>0) {
        kmin = (int) ceil(-rmax*dg);
        kmax = (int)floor(-rmin*dg);
      } else {
        kmin = (int) ceil(-rmin*dg);
        kmax = (int)floor(-rmax*dg);
      }
      assert kmin<=kmax : "kmin="+kmin+", kmax="+kmax;
      float[] dm = new float[nl];
      fill(Float.MAX_VALUE,d[isp]);
      // loop over all slope indices
      for (int k=kmin; k<=kmax; k++) {
        int ils = max(0,-k);
        int ile = min(nl,nl-k);
        for (int il=ils; il<ile; il++)
          dm[il] = d[ispm1][il+k] + e[ie][il];
        float r = (float)k/(float)dg; // slope
        if (r==0) { // zero slope, no interpolation necessary
          for (int x=je+is; x!=ie; x+=is)
            for (int il=ils; il<ile; il++)
              dm[il] += e[x][il]; 
        } else { // linearly interpolate
          for (int x=je+is; x!=ie; x+=is) {
            float ky = r*(ie-x);
            int k1 = (int)ky;
            if (ky<0.0f) --k1;
            int k2 = k1+1;
            float w1 = k2-ky;
            float w2 = 1.0f-w1;
            for (int il=ils; il<ile; il++)
              dm[il] += w1*e[x][k1+il]+w2*e[x][k2+il];
          }  
        }
        // update previous errors and record moves.
        for (int il=ils; il<ile; il++) {
          if (dm[il]<d[isp][il]) {
            d[isp][il] = dm[il];
            if (m!=null)
              m[ispm1][il] = k;
          }
        }
      }
    }
  }
  
  private static void backtrack(
      int dir, Sampling shifts, float[][] d, float[][] m, float[] u) 
    {
      int nl = d[0].length;
      int ni = d.length;
      int nlm1 = nl-1;
      int nim1 = ni-1;
      int ib = (dir>0)?0:nim1; // begining index
      int ie = (dir>0)?nim1:0; // end index
      int is = (dir>0)?1:-1;   // stride
      int ii = ib;
      // Set initial lag for the case that all errors at ii are equal.
      int il = (dir>0)?0:nlm1; 
      float dl = d[ii][il]; // Current accumulated error value.
          
      // Find minimum lag value(dl) and index(il) at trace index ii.
      for (int jl=0; jl<nl; ++jl) {
        if (d[ii][jl]<dl) {
          dl = d[ii][jl];
          il = jl;
        }
      }
      u[ii] = (float)shifts.getValue(il);
      while (ii!=ie) {
        ii += is;
        il += (int)m[ii][il];
        u[ii] = (float)shifts.getValue(il);
      }
    }
  
  /**
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private static void shiftAndScale(float emin, float emax, float[][] e) {
//    System.out.println("shiftAndScale: emin="+emin+" emax="+emax);
    int nl = e[0].length;
    int n1 = e.length;
    float eshift = emin;
    float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    for (int i1=0; i1<n1; ++i1) {
      for (int il=0; il<nl; ++il) {
        e[i1][il] = (e[i1][il]-eshift)*escale;
      }
    }
  }

  /**
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private static void shiftAndScale(float emin, float emax, float[][][] e) {
//    System.out.println("shiftAndScale: emin="+emin+" emax="+emax);
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final float eshift = emin;
    final float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    final float[][][] ef = e;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        for (int il=0; il<nl; ++il) {
          ef[i2][i1][il] = (ef[i2][i1][il]-eshift)*escale;
        }
      }
    }});
  }
  
  /**
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private static void shiftAndScale(float emin, float emax, float[][][][] e) {
//    System.out.println("shiftAndScale: emin="+emin+" emax="+emax);
    final int nl = e[0][0][0].length;
    final int n1 = e[0][0].length;
    final int n2 = e[0].length;
    final int n3 = e.length;
    final float eshift = emin;
    final float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    final float[][][][] ef = e;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          for (int il=0; il<nl; ++il) {
            ef[i3][i2][i1][il] = (ef[i3][i2][i1][il]-eshift)*escale;
          }
        }  
      }
    }});
  }

  private static float min3(float a, float b, float c) {
    return b<=a?(b<=c?b:c):(a<=c?a:c); // if equal, choose b
  }
  
  private static float min4(float a, float b, float c, float d) {
    float min = a;
    if (b<=min)
      min = b;
    if (c<=min)
      min = c;
    if (d<=min)
      min = d;
    return min;
  }

  private static float[] like(float[] a) {
    return new float[a.length];
  }
  private static float[][] like(float[][] a) {
    return new float[a.length][a[0].length];
  }
  private static float[][][] like(float[][][] a) {
    return new float[a.length][a[0].length][a[0][0].length];
  }
	private static float[] floats(double[] x) {
		int n = x.length;
		float[] y = new float[n];
		for (int i=0; i<n; ++i)
			y[i] = (float)x[i];
		return y;
	}
	private static float[] floats(int[] x) {
		int n = x.length;
		float[] y = new float[n];
		for (int i=0; i<n; ++i)
			y[i] = (float)x[i];
		return y;
	}


  private static class MinMax {
    float emin,emax;
    MinMax(float emin, float emax) {
      this.emin = emin;
      this.emax = emax;
    }
  }

  private static void print(String s) {
    System.out.println(s);
  }
	private static void print(int s) {
    System.out.println(s);
  }

  public static float[] applyShiftsS(float[] f, float[] u, Sampling ss) {
    int ns = f.length;
		float ds = (float)ss.getDelta();
		float fs = (float)ss.getFirst();
		float[] s = inverseInterpolation(ss,u);
		int nn = s.length;
    float[] h = new float[nn];
    SincInterp si = new SincInterp();
    for (int i1=0; i1<nn; ++i1) 
      h[i1] = si.interpolate(ns,ds,fs,f,s[i1]);
    return h;
  }

  private static float[] inverseInterpolation(
		Sampling ss, float[] u) 
	{
    int n1 = u.length;
    int n1m = n1-1;
		float ds = (float)ss.getDelta();
		float fs = (float)ss.getFirst();
		float[] r = add(u,floats(ss.getValues()));
    int n2 = (int)(ceil(r[n1m]/ds)-floor(r[0]/ds));
		float[] s = new float[n2];
   	inverseLinearInterpolation(n1,ds,fs,r,n2,ds,u[0]+fs,s,r[0],r[n1m]);
		return s;
	}

  private static float[][] inverseInterpolation(
		final Sampling ss, final float[][] u) 
	{
    final int n1 = u[0].length;
    final int n2 = u.length;
    final int n1m = n1-1;
		final float ds = (float)ss.getDelta();
		final float fs = (float)ss.getFirst();
		final float[][] s = new float[n2][n1];
		final float[] t = floats(ss.getValues()); 
    Parallel.loop(n2,new Parallel.LoopInt() {
      public void compute(int i2) {
			  float[] r = add(u[i2],t);
    	  int nn = (int)(ceil(r[n1m]/ds)-floor(r[0]/ds));
			  nn = (nn<n1)?n1:nn;
			  float[] st = new float[nn];
    	  inverseLinearInterpolation(n1,ds,fs,r,nn,ds,r[0],st,r[0],r[n1m]);
        int j1 = round(u[i2][0]/ds);
        j1 = (j1>0)?j1:0;
			  copy(n1-j1,0,st,j1,s[i2]);
		}});
		return s;
	}

	private static float[][][] inverseInterpolation(
		final Sampling ss, final float[][][] u) 
	{
    final int n1 = u[0][0].length;
    final int n2 = u[0].length;
    final int n3 = u.length;
    final int n1m = n1-1;
		final float ds = (float)ss.getDelta();
		final float fs = (float)ss.getFirst();
		final float[][][] s = new float[n3][n2][n1];
		final float[][][] rs = new float[n3][n2][n1];
		final float[] t = floats(ss.getValues());
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; i2++) {
				float[] r = add(u[i3][i2],t);
    		int nn = (int)(ceil(r[n1m]/ds)-floor(r[0]/ds));
			  nn = (nn<n1)?n1:nn;
			  float[] st = new float[nn];
    		inverseLinearInterpolation(n1,ds,fs,r,nn,ds,r[0],st,r[0],r[n1m]);
        int j1 = round(u[i3][i2][0]/ds);
        j1 = (j1>0)?j1:0;
			  copy(n1-j1,0,st,j1,s[i3][i2]);
			}
		}});
		return s;
	}
	
	/**
   * Inverse linear interpolation.
   * @author Dave Hale, Colorado School of Mines
   * @version 1989.06.02
   * @author Translated to java by Andrew Munoz, Colorado School of Mines
   * @version 2012.09.15
   *
   * Compute a regularly-sampled, monotonically increasing function x(y) from a 
   * regularly-sampled, monotonically increasing function y(x) by inverse linear 
   * interpolation.
   * @param nx	number of samples of y(x)
   * @param dx	x sampling interval; dx &gt 0.0 is required
   * @param fx	first x
   * @param y	array[nx] of y(x) values; y[0] &lt y[1] &lt ... &lt y[nx-1] required
   * @param ny	number of samples of x(y)
   * @param dy	y sampling interval; dy &gt 0.0 is required
   * @param fy	first y
   * @param xylo  x value assigned to x(y) when y is less than smallest y(x)
   * @param xyhi  x value assigned to x(y) when y is greater than largest y(x)
   * @return x  array[ny] of x(y) values
   */
  private static void inverseLinearInterpolation(
    int nx, float dx, float fx, float[] y, 
    int ny, float dy, float fy, float[] x, float xylo, float xyhi) 
  { 
    int nxi,nyo,jxi1,jxi2,jyo;
    float dxi,fxi,dyo,fyo,fyi,yo,xi1,yi1,yi2,yid,q; 
    nxi = nx; dxi = dx; fxi = fx;
    nyo = ny; dyo = dy; fyo = fy;
    fyi = y[0];
    // loop over output y less than smallest input y
    for (jyo=0,yo=fyo; jyo<nyo; jyo++,yo+=dyo) {
      if (yo>=fyi) break;
      x[jyo] = xylo;
    }
    // loop over output y between smallest and largest input y
    if (jyo==nyo-1 && yo==fyi) {
      x[jyo++] = fxi;
      yo += dyo;
    }
    jxi1 = 0;
    jxi2 = 1;
    xi1 = fxi;
    while (jxi2<nxi && jyo<nyo) {
      yi1 = y[jxi1];
      yi2 = y[jxi2];
      if (yi1<=yo && yo<=yi2) {
        yid = abs(yi2-yi1);
	x[jyo++] = (yid>0.0f)?xi1+dxi*(yo-yi1)/(yi2-yi1):xi1;
        yo += dyo;
    } else {
        jxi1++;
        jxi2++;
        xi1 += dxi;
      }
    }
    // loop over output y greater than largest input y
    while (jyo<nyo) x[jyo++] = xyhi;
  }

	public static class InterpShifts {

	  public enum Method {
    	LINEAR,
   	  MONOTONIC,
   	  SPLINE
  	}

		public InterpShifts() {
			_doS = true; // default SPLINE
		}

  	public void setMethod(Method method) {
   	 _method = method;
		 set();
  	}

  	public float[] interpolate(Sampling s1, int[] g1i, float[] u) {
  	  int ng1 = g1i.length;
  	  int n1 = s1.getCount();
			float[] g1 = makeGrid(s1,g1i);
  	  float[] ui = new float[n1];
			if (_doL) {
  	  	ci1 = new CubicInterpolator(
					CubicInterpolator.Method.LINEAR,ng1,g1,u);
			} else if (_doM) {
  	  	ci1 = new CubicInterpolator(
					CubicInterpolator.Method.MONOTONIC,ng1,g1,u);
			} else if (_doS) {
  	  	ci1 = new CubicInterpolator(
					CubicInterpolator.Method.SPLINE,ng1,g1,u);
			}
  	  ci1.interpolate(floats(s1.getValues()),ui);
  	  return ui;
  	}

  	public float[][] interpolate(
  	  Sampling s1, Sampling s2, int[] g1i, int[] g2i, float[][] u)
		{
  	  int ng1 = g1i.length;
  	  int ng2 = g2i.length;	
  	  int n1 = s1.getCount();
  	  int n2 = s2.getCount();
			float[] g1 = makeGrid(s1,g1i);
			float[] g2 = makeGrid(s2,g2i);
  	  float[][] ui = new float[n2][n1];
			if (_doL) {
  	  	BilinearInterpolator2 li = new BilinearInterpolator2(ng1,ng2,g1,g2,u); 
  	   	li.interpolate(s1,s2,ui);
  	  	return ui;
			}
			if (_doM) {
  	  	ci2 = new BicubicInterpolator2(BicubicInterpolator2.Method.MONOTONIC,
															 				 BicubicInterpolator2.Method.MONOTONIC,
															 				 ng1,ng2,g1,g2,u); 
			} else if (_doS) {
  	  	ci2 = new BicubicInterpolator2(BicubicInterpolator2.Method.SPLINE,
															 				 BicubicInterpolator2.Method.SPLINE,
															 				 ng1,ng2,g1,g2,u); 
			}
  	  ci2.interpolate(s1,s2,ui);
  	  return ui;
		}

  	public float[][][] interpolate(
  	  Sampling s1, Sampling s2, Sampling s3,
  	  int[] g1i, int[] g2i, int[] g3i, float[][][] u)
		{
  	  int ng1 = g1i.length;
  	  int ng2 = g2i.length;	
  	  int ng3 = g3i.length;	
  		int n1 = s1.getCount();
  		int n2 = s2.getCount();
  		int n3 = s3.getCount();
			float[] g1 = makeGrid(s1,g1i);
			float[] g2 = makeGrid(s2,g2i);
			float[] g3 = makeGrid(s3,g3i);
  	  float[][][] ui = new float[n3][n2][n1];
			if (_doL) {
  	  	TrilinearInterpolator3 li = new 
					TrilinearInterpolator3(ng1,ng2,ng3,g1,g2,g3,u); 
  	  	li.interpolate(s1,s2,s3,ui);
				return ui;
			}
			if (_doM) {
  	  	ci3 = new TricubicInterpolator3(TricubicInterpolator3.Method.MONOTONIC,
															 				  TricubicInterpolator3.Method.MONOTONIC,
															 				  TricubicInterpolator3.Method.MONOTONIC,
															 				  ng1,ng2,ng3,g1,g2,g3,u); 
			} else if (_doS) {
  	  	ci3 = new TricubicInterpolator3(TricubicInterpolator3.Method.SPLINE,
															 				  TricubicInterpolator3.Method.SPLINE,
															 				  TricubicInterpolator3.Method.SPLINE,
															 				  ng1,ng2,ng3,g1,g2,g3,u); 
			}
  	  ci3.interpolate(s1,s2,s3,ui);
			return ui;
		}

  	private Method _method;
		private boolean _doL,_doM,_doS;
		private	CubicInterpolator     ci1;
  	private BicubicInterpolator2  ci2;
  	private TricubicInterpolator3 ci3;

		private void set() {
			_doL = false; _doM = false; _doS = false;
  	  switch (_method) {
  	    case LINEAR:    _doL = true; break;
  	    case MONOTONIC: _doM = true; break;
  	    case SPLINE:    _doS = true; break;
  	    default: throw new IllegalArgumentException(
  	      _method.toString()+" is not a recognized interpolation method.");
  	  }
  	}
		private static float[] makeGrid(
			Sampling s1, int[] g) 
		{
			int n = g.length;
			float f = (float)s1.getFirst();
			float d = (float)s1.getDelta();
			float[] gf = new float[n];
			for (int i=0; i<n; ++i) 
				gf[i] = f+d*g[i];
			return gf;
		}
		private static float[][] makeGrid(
			Sampling s1, int n2, int[] g) 
		{
			float[] gr1 = makeGrid(s1,g);
			float[][] gr2 = new float[n2][gr1.length];
			for (int i=0; i<n2; ++i) 
				gr2[i] = gr1;
			return gr2;
		}
		private static float[][][] makeGrid(
			Sampling s1, int n2, int n3, int[] g) 
		{
			float[][] gr2 = makeGrid(s1,n2,g);
			float[][][] gr3 = new float[n3][n2][gr2[0].length];
			for (int i=0; i<n3; ++i) 
				gr3[i] = gr2;
			return gr3;
		}
		private static float[] floats(double[] x) {
			int n = x.length;
			float[] y = new float[n];
			for (int i=0; i<n; ++i)
				y[i] = (float)x[i];
			return y;
	}

	}
  
  /* testing only */
	public static class InterpShiftsX {

	  public enum Method {
    	LINEAR,
   	  MONOTONIC,
   	  SPLINE
  	}

		public InterpShiftsX() {
			_doS = true; // default SPLINE
		}

  	public void setMethod(Method method) {
   	 _method = method;
		 set();
  	}

  	private Method _method;
		private boolean _doL,_doM,_doS;
		private	CubicInterpolator     ci1;
  	private BicubicInterpolator2  ci2;
  	private TricubicInterpolator3 ci3;

		private void set() {
			_doL = false; _doM = false; _doS = false;
  	  switch (_method) {
  	    case LINEAR:    _doL = true; break;
  	    case MONOTONIC: _doM = true; break;
  	    case SPLINE:    _doS = true; break;
  	    default: throw new IllegalArgumentException(
  	      _method.toString()+" is not a recognized interpolation method.");
  	  }
  	}
    private static CubicInterpolator makeInterpolator1(
      float[] x, float[] y) 
    {
      return new CubicInterpolator(CubicInterpolator.Method.SPLINE,x,y);
    }
    private static CubicInterpolator makeInterpolator2(
      float[] x, float[] y) 
    {
      return new CubicInterpolator(CubicInterpolator.Method.SPLINE,x,y);
    }
    private static CubicInterpolator makeCubicInterpolator(
      float[] x, float[] y, float[] yd) 
    {
      return new CubicInterpolator(x,y,yd);
    }
    
   public static float[] interpolate(
      Sampling s1, int[] k1s, float[] uk) 
    {
      int n1 = s1.getCount();
      int nk1 = k1s.length;
      float[] xk1 = new float[nk1];
      for (int jk1=0; jk1<nk1; ++jk1)
        xk1[jk1] = (float)s1.getValue(k1s[jk1]);
      CubicInterpolator ci = makeInterpolator1(xk1,uk);
      float[] u = new float[n1];
      for (int j1=0; j1<n1; ++j1) {
        float x1 = (float)s1.getValue(j1);
        u[j1] = ci.interpolate(x1);
      }
      return u;
    }
    public static float[][] interpolate(
      Sampling s1, Sampling s2, int[] k1s, int[] k2s, float[][] ukk) 
    {
      //trace("ukk:"); dump(ukk);
      int n1 = s1.getCount();
      int n2 = s2.getCount();
      int nk1 = k1s.length;
      int nk2 = k2s.length;

      // Coarse sampling of 1st and 2nd dimensions.
      float[] xk1 = new float[nk1];
      for (int jk1=0; jk1<nk1; ++jk1)
        xk1[jk1] = (float)s1.getValue(k1s[jk1]);
      float[] xk2 = new float[nk2];
      for (int jk2=0; jk2<nk2; ++jk2)
        xk2[jk2] = (float)s2.getValue(k2s[jk2]);

      // Compute 1st derivatives in 1st dimension.
      float[][] vkk = new float[nk2][nk1];
      for (int jk2=0; jk2<nk2; ++jk2) {
        CubicInterpolator ci = makeInterpolator1(xk1,ukk[jk2]);
        ci.interpolate1(xk1,vkk[jk2]);
      }
        
      // Interpolate in 2nd dimension.
      float[] uk2 = new float[nk2];
      float[] vk2 = new float[nk2];
      float[][] uk = new float[n2][nk1];
      float[][] vk = new float[n2][nk1];
      for (int jk1=0; jk1<nk1; ++jk1) {
        for (int jk2=0; jk2<nk2; ++jk2) {
          uk2[jk2] = ukk[jk2][jk1];
          vk2[jk2] = vkk[jk2][jk1];
        }
        CubicInterpolator ciu = makeInterpolator2(xk2,uk2);
        CubicInterpolator civ = makeInterpolator2(xk2,vk2);
        for (int j2=0; j2<n2; ++j2) {
          float x2 = (float)s2.getValue(j2);
          uk[j2][jk1] = ciu.interpolate(x2);
          vk[j2][jk1] = civ.interpolate(x2);
        }
      }

      // Interpolate 1st dimension.
      float[][] u = new float[n2][n1];
      for (int j2=0; j2<n2; ++j2) {
        CubicInterpolator ci = makeCubicInterpolator(xk1,uk[j2],vk[j2]);
        for (int j1=0; j1<n1; ++j1) {
          float x1 = (float)s1.getValue(j1);
          u[j2][j1] = ci.interpolate(x1);
        }
      }
      return u;
    }
 }


};

