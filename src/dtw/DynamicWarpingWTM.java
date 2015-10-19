package dtw;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.SincInterp;
import edu.mines.jtk.interp.BicubicInterpolator2;
import edu.mines.jtk.interp.BilinearInterpolator2;
import edu.mines.jtk.interp.CubicInterpolator;
import edu.mines.jtk.interp.CubicInterpolator.Method;
import edu.mines.jtk.interp.TrilinearInterpolator3;
import edu.mines.jtk.interp.TricubicInterpolator3;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Dynamic warping for multiple wells. 
 * Thanks to Stefan Compton for authoring many of these methods
 * @author Andrew Munoz and Stefan Compton, CSM
 * @version 11.16.13
 */
public class DynamicWarpingWTM {
  
  public DynamicWarpingWTM(float[][] f, float[][] g, int smin, int smax) {
    Check.argument(f.length==g.length,"f.length==g.length");
    int n1f = f[0].length;
    int n1g = g[0].length;
    _fr = 1;
		int nl = 1+smax-smin;
    int[] el = new int[]{n1f,nl};
    _ss = new Sampling(nl,1.0,smin);
    _nel = el[1];
    _ne1 = el[0];
    _n1f = n1f;
    _n1g = n1g;
    _n2 = f.length;
    _ff2 = f;
    _gg2 = g;
    _si = new SincInterp();
  }
  
  public DynamicWarpingWTM(float[][][] f, float[][][] g, int smin, int smax) {
    Check.argument(f.length==g.length,"f.length==g.length");
    Check.argument(f[0].length==g[0].length,"f[0].length==g[0].length");
    int n1f = f[0][0].length;
    int n1g = g[0][0].length;
    _fr = 1;
		int nl = 1+smax-smin;
    int[] el = new int[]{n1f,nl};
    _ss = new Sampling(nl,1.0,smin);
    _nel = el[1];
    _ne1 = el[0];
    _n1f = n1f;
    _n1g = n1g;
    _n2 = f[0].length;
    _n3 = f.length;
    _ff3 = f;
    _gg3 = g;
    _si = new SincInterp();
  }
  
  /**
   * Returns the number of PP samples used for alignment errors.
   * This length is computed from the average Vp/Vs given in the 
   * constructor. 
   * @return the number of PP samples used for alignment errors.
   */
  public int getPPErrorLength() {
    return _ne1;
  }
  
  /**
   * Returns the number of lags. The number of lags is computed from
   * the average Vp/Vs given in the constructor. It is the difference
   * between the scaled PP trace and PS trace length. 
   * @return the number of lags.
   */
  public int getNumberOfLags() {
    return _nel;
  }
  
  /**
   * Find shifts for 2D images. Shifts are computed for ff(g1,g2),
   * where g1 and g2 are sparse samples of ff(x1,x2). The sparse 
   * shifts are then interpolated back to the fine grid. The 
   * sparse samples g1 are constant for each g2. 
   * <p>
   * The minimum interval of g1 and g2 samples is ceil(1/dr1) and
   * ceil(1/dr2) respectively. These intervals are variable in 
   * order to include the first and last indices of the fine grid.
   * @param r1min minimum slope in first dimension.
   * @param r1max maximum slope in first dimension. If |r1max|&lt
   *  dr1. Then the maximum slope is zero.
   * @param dr1 controls the sparse grid interval in the first 
   *  dimension. A smaller value results in a larger interval and 
   *  smoother shifts.
   * @param r2min minimum slope in second dimension.
   * @param r2max maximum slope in second dimension. If |r2max|&lt
   *  dr2. Then the maximum slope is zero.
   * @param dr2 controls the sparse grid interval in the second 
   *  dimension. A smaller value results in a larger interval and 
   *  smoother shifts. 
   * @return shifts for 2D images.
   */
  public float[][] findShifts(
      final float r1min, final float r1max, final float dr1,
      final float r2min, final float r2max, final float dr2,
			final float[] x2map)
  {
    final int dg1 = (int)ceil(1.0f/dr1);
    final int dg2 = (int)ceil(1.0f/dr2);
    final int[] g1 = Subsample.subsample(_ne1,dg1);
		final int[] g2 = Subsample.subsample(_n2,dg2,x2map);
    print("g1:"); dump(g1); print("g2"); dump(g2);
		_g1 = g1;
		_g2 = g2;
    
    Stopwatch s = new Stopwatch();
    s.start();
    print("Smoothing 1st dimension...");
    final float[][][] es1 = smoothErrors1(_ff2,_gg2,r1min,r1max,g1);
    print("Finished 1st dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es1);
    
    s.restart();
    print("Smoothing 2nd dimension...");
    final float[][][] es = smoothErrors2(es1,r2min,r2max,g2);
    print("Finished 2nd dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es);
    
    final int ng2 = es.length;
    final float[][] u = new float[ng2][];
    Parallel.loop(ng2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][][] dm = accumulateForward(es[i2],g1,r1min,r1max);
      u[i2] = backtrackReverse(dm[0],dm[1]);
    }});
    return interpolate(_ne1,_n2,g1,g2,u,false,(float)_ss.getFirst());
  }
	public float[][] findShifts(float r1min, float r1max, float dr1,
      												float r2min, float r2max, float dr2) {
	return findShifts(r1min,r1max,dr1,r2min,r2max,dr2,new float[1]);
	}

  
  /**
   * Find shifts for 3D images. Shifts are computed for ff(g1,g2,g3),
   * where g1, g2, and g3 are sparse samples of ff(x1,x2,x3). The
   * sparse shifts are then interpolated back to the fine grid. The 
   * sparse samples g1 are constant for every g2 and g3. 
   * <p>
   * The minimum interval of g1, g2, and g3 samples is ceil(1/dr1),
   * ceil(1/dr2), and ceil(1/dr3) respectively. These intervals are
   * variable in order to include the first and last indices of the
   * fine grid.
   * @param r1min minimum slope in first dimension.
   * @param r1max maximum slope in first dimension. If |r1max|&lt
   *  dr1. Then the maximum slope is zero.
   * @param dr1 controls the sparse grid interval in the first 
   *  dimension. A smaller value results in a larger interval and 
   *  smoother shifts.
   * @param r2min minimum slope in second dimension.
   * @param r2max maximum slope in second dimension. If |r2max|&lt
   *  dr2. Then the maximum slope is zero.
   * @param dr2 controls the sparse grid interval in the second 
   *  dimension. A smaller value results in a larger interval and 
   *  smoother shifts. 
   * @param r3min minimum slope in third dimension.
   * @param r3max maximum slope in third dimension. If |r3max|&lt
   *  dr3. Then the maximum slope is zero.
   * @param dr3 controls the sparse grid interval in the third 
   *  dimension. A smaller value results in a larger interval and 
   *  smoother shifts. 
   * @return shifts for 3D images.
   */
  public float[][][] findShifts(
      final float r1min, final float r1max, final float dr1, 
      final float r2min, final float r2max, final float dr2,
      final float r3min, final float r3max, final float dr3)
  {
    final int dg1 = (int)ceil(1.0f/dr1);
    final int dg2 = (int)ceil(1.0f/dr2);
    final int dg3 = (int)ceil(1.0f/dr3);
    final int[] g1 = Subsample.subsample(_ne1,dg1);
    final int[] g2 = Subsample.subsample( _n2,dg2);
    final int[] g3 = Subsample.subsample( _n3,dg3);
    print("g1:"); dump(g1); print("g2"); dump(g2); print("g3"); dump(g3);
    
    Stopwatch s = new Stopwatch();
    s.start();
    print("Smoothing 1st dimension...");
    final float[][][][] es1 = smoothErrors1(_ff3,_gg3,r1min,r1max,g1);
    print("Finished 1st dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es1);
    
    s.restart();
    print("Smoothing 2nd dimension...");
    final float[][][][] es12 = smoothErrors2(es1,r2min,r2max,g2);
    print("Finished 2nd dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es12);
    
    s.restart();
    print("Smoothing 3rd dimension...");
    final float[][][][] es = smoothErrors3(es12,r3min,r3max,g3);
    print("Finished 3rd dimension smoothing in "+s.time()+" seconds");
    normalizeErrors(es);
    
    final int ng3 = es.length;
    final int ng2 = es[0].length;
    final float[][][] u = new float[ng3][ng2][];
    Parallel.loop(ng3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<ng2; i2++) {
        float[][][] dm = accumulateForward(es[i3][i2],g1,r1min,r1max);
        u[i3][i2] = backtrackReverse(dm[0],dm[1]);
      }
    }});
    return interpolate(_ne1,_n2,_n3,g1,g2,g3,u,false,(float)_ss.getFirst());
  }

  public float[][][] applyShifts(float[][] u) {
    final int n1u = u[0].length;
    final int n1um = n1u-1;
		final float[][] uf = u;
    final float[][] hf = new float[_n2][_n1f];
		final float[][] s = inverseInterpolation(uf);
    Parallel.loop(_n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1u; ++i1) {
        hf[i2][i1] = _si.interpolate(_n1f,1.0,0.0,_ff2[i2],s[i2][i1]);
      }
    }});
    return new float[][][]{hf,s};
  }

  public float[][][][] applyShifts(float[][][] u) {
    final int n1u = u[0][0].length;
    final int n1um = n1u-1;
		final float[][][] uf = u;
    final float[][][] hf = new float[_n3][_n2][_n1f];
		final float[][][] s = inverseInterpolation(uf);
    Parallel.loop(_n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<_n2; i2++) {
        for (int i1=0; i1<n1u; i1++) {
          hf[i3][i2][i1] = _si.interpolate(_n1f,1.0,0.0,_ff3[i3][i2],
              						s[i3][i2][i1]);
        }
      }
    }});
    return new float[][][][]{hf,s};
  }

  public float[][] applyShiftsF(float[][] u) {
    final int n1u = u[0].length;
    final int n1um = n1u-1;
		final float[][] uf = u;
    final float[][] hf = new float[_n2][_n1f];
    Parallel.loop(_n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1u; ++i1) {
        hf[i2][i1] = _si.interpolate(_n1g,1.0,0.0,_gg2[i2],i1+uf[i2][i1]);
      }
    }});
    return hf;
  }

  public float[][][] applyShiftsF(float[][][] u) {
    final int n1u = u[0][0].length;
    final int n1um = n1u-1;
		final float[][][] uf = u;
    final float[][][] hf = new float[_n3][_n2][_n1f];
    Parallel.loop(_n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<_n2; i2++) {
        for (int i1=0; i1<n1u; i1++) {
          hf[i3][i2][i1] = _si.interpolate(_n1g,1.0,0.0,_gg3[i3][i2],
              						i1+uf[i3][i2][i1]);
        }
      }
    }});
    return hf;
  }

	 /** 
		* Computes warped synthetic seismograms from time shifts
		* @param u  time shifts at the synthetic (in samples)
		* @param ss synthetic seismogram samplings
		* @param st seismic sampling
		* @param x  initial synthetic seismogram
		* @return h[ns] synthetic seismogram, and new first sampling fs
		*/
	public float[][] getWarpedSyntheticSeismogram(
		Sampling st, Sampling ss, float[] u, float[] x) 
	{
		int nt = st.getCount();
		int ns = ss.getCount();
		float ft = (float)st.getFirst();
		float dt = (float)st.getDelta();
		float fs = (float)ss.getFirst();
		float ls = (float)ss.getLast();
		int fio  = round((fs-ft)/dt);
		int lio  = round((ls-ft)/dt);
		int no   = round((ls-fs)/dt);
		if (no+fio>nt) no -= no+fio-nt;
		float[] us = copy(no,fio,u);
		float[] tt = rampfloat((float)fio,1f,no);
		float[] s = inverseInterpolation(us,tt);
		int n2 = s.length;
		float[] h = new float[n2];
		for (int i=0; i<n2; ++i)
			h[i] = _si.interpolate(nt,1.0,0.0,x,s[i]);	
		return new float[][]{h,new float[]{s[0]*dt}};
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
    backtrack(-1,_ss,d,m,u);
    return u;
  }
  
  public float[] backtrackForward(float[][] d, float[][] m) {
    float[] u = new float[d.length];
    backtrack(1,_ss,d,m,u);
    return u;
  }

  /**
   * Interpolates u, sampled on sparse grid g, to a uniformly 
   * sampled array of length n.
   * @param n length of the output array.
   * @param g sparse grid indices.
   * @param u sparse shifts.
   * @param doLinear {@code true} for linear interpolation,
   *  {@code false} for cubic.
   * @return the interpolated shifts.
   */
  public static float[] interpolate(
      int n, int[] g, float[] u, boolean doLinear) 
  {
    int ng = g.length;
    float[] gf = new float[ng];
    float[] ui = new float[n ];
    for (int ig=0; ig<ng; ig++)
      gf[ig] = (float)g[ig];
    Method m = (doLinear)?Method.LINEAR:CIM;
    CubicInterpolator ci = new CubicInterpolator(m,ng,gf,u);
    for (int i=0; i<n; i++)
      ui[i] = ci.interpolate(i);
    return ui;
  }

  /**
   * Interpolates u, sampled on sparse grid g, to a uniformly 
   * sampled array of size [n2][n1].
   * @param n1 length of first dimension of output array.
   * @param n2 length of second dimension of output array.
   * @param g1 first dimension sparse grid indices.
   * @param g2 second dimension sparse grid indices.
   * @param u sparse shifts.
   * @param doLinear {@code true} for bilinear interpolation,
   *  {@code false} for bicubic.
   * @return the interpolated shifts.
   */
  public static float[][] interpolate(
      int n1, int n2, int[] g1, int[] g2, float[][] u, boolean doLinear, float ds)
  {
    int ng1 = g1.length;
    int ng2 = g2.length;
    Check.argument(ng1==u[0].length,"ng1==u[0].length: "+ng1+"=="+u[0].length);
    Check.argument(ng2==u.length,"ng2==u.length: "+ng2+"=="+u.length);
    float[] g1f = new float[ng1];
    float[] g2f = new float[ng2];
    for (int ig1=0; ig1<ng1; ig1++)
      g1f[ig1] = (float)g1[ig1];
    for (int ig2=0; ig2<ng2; ig2++)
      g2f[ig2] = (float)g2[ig2];
    Sampling s1 = new Sampling(n1,1.0,0.0);
    Sampling s2 = new Sampling(n2,1.0,0.0);
    if (doLinear) {
      BilinearInterpolator2 bli = new BilinearInterpolator2(g1f,g2f,u);
      return add(bli.interpolate(s1,s2),ds); 
    } else {
      BicubicInterpolator2 bci = 
          new BicubicInterpolator2(BCIM1,BCIM2,g1f,g2f,u);
      return add(bci.interpolate(s1,s2),ds); 
    }
  }
  
  /**
   * Interpolates u, sampled on sparse grid g, to a uniformly 
   * sampled array of size [n3][n2][n1].
   * @param n1 length of first dimension of output array.
   * @param n2 length of second dimension of output array.
   * @param n3 length of third dimension of output array.
   * @param g1 first dimension sparse grid indices.
   * @param g2 second dimension sparse grid indices.
   * @param g3 third dimension sparse grid indices.
   * @param u sparse shifts.
   * @param doLinear {@code true} for trilinear interpolation,
   *  {@code false} for tricubic.
   * @return the interpolated shifts.
   */
  public static float[][][] interpolate(
      int n1, int n2, int n3, int[] g1, int[] g2, int[] g3, float[][][] u,
      boolean doLinear, float ds)
  {
    int ng1 = g1.length;
    int ng2 = g2.length;
    int ng3 = g3.length;
    Check.argument(
        ng1==u[0][0].length,"ng1==u[0][0].length: "+ng1+"=="+u[0][0].length);
    Check.argument(ng2==u[0].length,"ng2==u[0].length: "+ng2+"=="+u[0].length);
    Check.argument(ng3==u.length,"ng3==u.length: "+ng3+"=="+u.length);
    float[] g1f = new float[ng1];
    float[] g2f = new float[ng2];
    float[] g3f = new float[ng3];
    for (int ig1=0; ig1<ng1; ig1++)
      g1f[ig1] = (float)g1[ig1];
    for (int ig2=0; ig2<ng2; ig2++)
      g2f[ig2] = (float)g2[ig2];
    for (int ig3=0; ig3<ng3; ig3++)
      g3f[ig3] = (float)g3[ig3];
    Sampling s1 = new Sampling(n1,1.0,0.0);
    Sampling s2 = new Sampling(n2,1.0,0.0);
    Sampling s3 = new Sampling(n3,1.0,0.0);
    if (doLinear) {
      TrilinearInterpolator3 tli = new TrilinearInterpolator3(g1f,g2f,g3f,u);
      return add(tli.interpolate(s1,s2,s3),ds);
    } else {
      TricubicInterpolator3 tci = 
          new TricubicInterpolator3(TCIM1,TCIM2,TCIM3,g1f,g2f,g3f,u);
      return add(tci.interpolate(s1,s2,s3),ds);
    }
  }
    
  public static float[][][] shiftVolume(
      float[] u1, float[] uS, int fr, float[][][] e) {
    int n1 = e.length;
    int nl1 = e[0].length;
    int nlS = e[0][0].length;
    Check.argument(n1==u1.length, "n1==u1.length");
    Check.argument(n1==uS.length, "n1==uS.length");
    float[][][] sv = zerofloat(nlS,nl1,n1);
    for (int i1=0; i1<n1; ++i1) {
      int il1 = (int)u1[i1];
      int ilS = (int)uS[i1]*fr;
      sv[i1][il1][ilS] = 1.0f;
    }
    return sv;
  }
  
  public static float[] extrapolateShifts(float[] u, int n1) {
    int nu = u.length;
    float v = u[nu-1];
    float[] eu = new float[n1];
    for (int iu=0; iu<nu; ++iu) {
      eu[iu] = u[iu];
    }
    for (int iu=nu; iu<n1; ++iu) {
      eu[iu] = v;
    }
    return eu;
  }
  
  public static float[][] extrapolateShifts(float[][] u, int n1) {
    int nu = u[0].length;
    int n2 = u.length;
    float[][] eu = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
      float v = u[i2][nu-1];
      for (int iu=0; iu<nu; ++iu) {
        eu[i2][iu] = u[i2][iu];
      }
      for (int iu=nu; iu<n1; ++iu) {
        eu[i2][iu] = v;
      }
    }
    return eu;
  }
  
  public static float[][][] extrapolateShifts(float[][][] u, int n1) {
    int nu = u[0][0].length;
    int n2 = u[0].length;
    int n3 = u.length;
    float[][][] eu = new float[n3][n2][n1];
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        float v = u[i3][i2][nu-1];
        for (int iu=0; iu<nu; iu++) {
          eu[i3][i2][iu] = u[i3][i2][iu];
        }
        for (int iu=nu; iu<n1; ++iu) {
          eu[i3][i2][iu] = v;
        }
      }  
    }
    return eu;
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

  /**
   * Returns errors in an array with lag the slowest dimension.
   * Useful only for visualization of errors. Other methods in this
   * class assume that lag is the fastest dimension in arrays of errors.
   * @param e array[n1][nl] of errors.
   * @return transposed array[nl][n1] of errors.
   */
  public static float[][] transposeLag(float[][] e) {
    int nl = e[0].length;
    int n1 = e.length;
    float[][] t = new float[nl][n1];
    for (int il=0; il<nl; ++il) {
      for (int i1=0; i1<n1; ++i1) {
        t[il][i1] = e[i1][il];
      }
    }
    return t;
  }

  /**
   * Returns errors in an array with lag the slowest dimension.
   * Useful only for visualization of errors. Other methods in this
   * class assume that lag is the fastest dimension in arrays of errors.
   * @param e array[n2][n1][nl] of errors.
   * @return transposed array[nl][n2][n1] of errors.
   */
  public static float[][][] transposeLag(float[][][] e) {
    int nl = e[0][0].length;
    int n1 = e[0].length;
    int n2 = e.length;
    float[][][] t = new float[nl][n2][n1];
    for (int il=0; il<nl; ++il) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          t[il][i2][i1] = e[i2][i1][il];
        }
      }
    }
    return t;
  }
  
  public static float[][][] transposeLag12(float[][][] e) {
    int nl = e[0][0].length;
    int n1 = e[0].length;
    int n2 = e.length;
    float[][][] t = new float[n2][nl][n1];
    for (int i2=0; i2<n2; ++i2) {
      for (int il=0; il<nl; ++il) {
        for (int i1=0; i1<n1; ++i1) {
          t[i2][il][i1] = e[i2][i1][il];
        }
      }
    }
    return t;
  }
  
  public static float[][][] transposeLag23(float[][][] e) {
    int nl = e[0][0].length;
    int n1 = e[0].length;
    int n2 = e.length;
    float[][][] t = new float[n1][n2][nl];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2) {
        t[i1][i2] = e[i2][i1];
      }
    }
    return t;
  }
  
  public static float[][][][] transposeLag12(float[][][][] e) {
    int nl = e[0][0][0].length;
    int n1 = e[0][0].length;
    int n2 = e[0].length;
    int n3 = e.length;
    float[][][][] t = new float[n3][n2][nl][n1];
    for (int i3=0; i3<n3; ++i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int il=0; il<nl; ++il) {
          for (int i1=0; i1<n1; ++i1) {
            t[i3][i2][il][i1] = e[i3][i2][i1][il];
          }
        }
      }  
    }
    return t;
  }
  
  public static float[][][][] getEmptySparseErrors(int n3) {
    return new float[n3][][][];
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




  ///////////////////////////////////////////////////////////////////////////
  // private

  private int _nel; // number of lags
  private int _ne1; // number of alignment error samples
  private int _n1f; // number of pp samples before scaling
  private int _n1g; // number of pp samples before scaling
  private int _n2; // number of traces
  private int _n3; // number of ensembles
  private float[] _pp1; // pp trace
  private float[] _ps1; // ps trace
  private float[][] _ff2; // pp traces
  private float[][] _gg2; // ps traces
  private float[][][] _ff3; // pp traces
  private float[][][] _gg3; // ps traces
	private int[] _g1;
	private int[] _g2;
	private int[] _g3;
  private int _fr; // fractional shift factor (1.0/_fr is the shift interval)
	private Sampling _ss; // lag samplings
  private float _epow = 2; // exponent used for alignment errors |f-g|^e
  private SincInterp _si; // for warping with non-integer shifts
	private static final CubicInterpolator.Method CIM = 
      CubicInterpolator.Method.SPLINE;
  private static final BicubicInterpolator2.Method BCIM1 =
      BicubicInterpolator2.Method.SPLINE;
  private static final BicubicInterpolator2.Method BCIM2 =
      BicubicInterpolator2.Method.SPLINE;
  private static final TricubicInterpolator3.Method TCIM1 = 
      TricubicInterpolator3.Method.SPLINE;
  private static final TricubicInterpolator3.Method TCIM2 = 
      TricubicInterpolator3.Method.SPLINE;
  private static final TricubicInterpolator3.Method TCIM3 = 
      TricubicInterpolator3.Method.SPLINE;

  private float error(float f, float g) {
    return pow(abs(f-g),_epow);
  }

  //private void computeErrors(float[] f, float[] g, float[][] e) {
  //  int n1max = e.length;
  //  int nl1 = e[0].length;
  //  for (int i1=0; i1<n1max; ++i1) {
  //    for (int il1=0,j1=i1*_fr; il1<nl1; ++il1,++j1) {
  //      e[i1][il1] = error(f[i1],g[j1]);
  //    }
  //  }
  //}
  private void computeErrors(
		float[] f, float[] g, float[][] e)
  {
    int ns = e[0].length;
    int ne = e.length;
    int nf = f.length;
    int ng = g.length;
		Sampling ss = _ss;
    float[] fi = new float[ne];
    float[] gi = new float[ne];
    _si.interpolate(nf,1.0,0.0,f,ne,1.0,0.0,fi);
    for (int is=0; is<ns; ++is) {
      _si.interpolate(ng,1.0,0.0,g,ne,1.0,ss.getValue(is),gi);
      for (int ie=0; ie<ne; ++ie)
        e[ie][is] = error(fi[ie],gi[ie]);
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
    for (int i1=0; i1<ng; i1++) {
      for (int il=0; il<nel; il++) {
        es[i1][il] = scale*(ef[i1][il]+er[i1][il]-e[g[i1]][il]);
      }
    }
    return es;
  }
  
  /**
   * Returns alignment errors smoothed in the first dimension.
   * Returned errors are sparse in the first dimension, and
   * unchanged in the second dimension. Alignment errors are
   * computed on the fly.
   * @param pp the PP image.
   * @param ps the PS image.
   * @param r1min minimum slope in the first dimension.
   * @param r1max maximum slope in the first dimension.
   * @param g1 first dimension sparse grid indices.
   * @return smoothed alignment errors with size 
   *  [_n2][g1.length][_nel].
   */
  private float[][][] smoothErrors1(final float[][] pp, final float[][] ps,
      final float r1min, final float r1max, final int[] g1)
  {
    final int ng1 = g1.length;
    final float[][][] es1 = new float[_n2][ng1][_nel];
    Parallel.loop(_n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] e = new float[_ne1][_nel];
      computeErrors(pp[i2],ps[i2],e);
      es1[i2] = smoothErrors(e,r1min,r1max,g1);
    }});
    return es1;
  }
  
  /**
   * Returns alignment errors smoothed in the first dimension.
   * Returned errors are sparse in the first dimension, and
   * unchanged in the second dimension. Alignment errors are
   * computed on the fly.
   * @param pp the PP image.
   * @param ps the PS image.
   * @param r1min minimum slope in the first dimension.
   * @param r1max maximum slope in the first dimension.
   * @param g1 first dimension sparse grid indices specified
   *  for all _n2 (size[_n2][])
   * @return smoothed alignment errors with size 
   *  [_n2][g1.length][_nel].
   */
  private float[][][] smoothErrors1(
      final float[][] pp, final float[][] ps,
      final float r1min, final float r1max, final int[][] g1)
  {
    final int ng1 = g1[0].length;
    final float[][][] es1 = new float[_n2][ng1][_nel];
    Parallel.loop(_n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      float[][] e = new float[_ne1][_nel];
      computeErrors(pp[i2],ps[i2],e);
      es1[i2] = smoothErrors(e,r1min,r1max,g1[i2]);
    }});
    return es1;
  }
  
  /**
   * Returns alignment errors smoothed in the first dimension.
   * Returned errors are sparse in the first dimension, and
   * unchanged in the second and third dimension. Alignment 
   * errors are computed on the fly.
   * @param pp the PP image.
   * @param ps the PS image.
   * @param r1min minimum slope in the first dimension.
   * @param r1max maximum slope in the first dimension.
   * @param g1 first dimension sparse grid indices.
   * @return smoothed alignment errors with size 
   *  [_n3][_n2][g1.length][_nel].
   */
  private float[][][][] smoothErrors1(
      final float[][][] pp, final float[][][] ps,
      final float r1min, final float r1max, final int[] g1)
  {
    final int ng1 = g1.length;
    final float[][][][] es1 = new float[_n3][_n2][ng1][_nel];
    Parallel.loop(_n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<_n2; i2++) {
        float[][] e = new float[_ne1][_nel];
        computeErrors(pp[i3][i2],ps[i3][i2],e);
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
   * @param pp the PP image.
   * @param ps the PS image.
   * @param r1min minimum slope in the first dimension.
   * @param r1max maximum slope in the first dimension.
   * @param g1 first dimension sparse grid indices specified
   *  for all _n3 and _n2 (size[_n3][_n2][]).
   * @return smoothed alignment errors with size 
   *  [_n3][_n2][g1.length][_nel].
   */
  private float[][][][] smoothErrors1(
      final float[][][] pp, final float[][][] ps,
      final float r1min, final float r1max, final int[][][] g1)
  {
    final int ng1 = g1[0][0].length;
    final float[][][][] es1 = new float[_n3][_n2][ng1][_nel];
    Parallel.loop(_n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<_n2; i2++) {
        float[][] e = new float[_ne1][_nel];
        computeErrors(pp[i3][i2],ps[i3][i2],e);
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
  
  /**
   * Accumulation for 3D alignment Errors
   * @param dir
   * @param b1
   * @param bS
   * @param e
   * @param d
   */
  private static void accumulate(
      int dir, int b1, int bS, float[][][] e, float[][][] d)
  {
    int nl1 = e[0].length;
    int nl1m1 = nl1-1;
    int nlS = e[0][0].length;
    int nlSm1 = nlS-1;
    int ni = e.length;
    int nim1 = ni-1;
    int ib = (dir>0)?0:nim1; // beginning index
    int ie = (dir>0)?ni:-1;  // end index
    int is = (dir>0)?1:-1;   // stride
    int ic = (dir>0)?-1:1;   // contraint dir, forward=-lag, reverse=+lag
    for (int il1=0; il1<nl1; ++il1)
      for (int ilS=0; ilS<nlS; ++ilS)
        d[ib][il1][ilS] = 0.0f;
    for (int ii=ib; ii!=ie; ii+=is) {
      int ji = max(0,min(nim1,ii-is));
      int jb1 = max(0,min(nim1,ii-is*b1));
      int jbS = max(0,min(nim1,ii-is*bS));
      for (int il1=0; il1<nl1; ++il1) {
        int il1c = il1+ic;
        il1c = (il1c==-1)?0:(il1c==nl1)?nl1m1:il1c;
        for (int ilS=0; ilS<nlS; ++ilS) {
          int ilSc = ilS+ic;
          ilSc = (ilSc==-1)?0:(ilSc==nlS)?nlSm1:ilSc;
          float dc1 = d[jb1][il1c][ilS ];
          float dcS = d[jbS][il1 ][ilSc];
          float dc1S= d[jbS][il1c][ilSc];
          float di  = d[ji ][il1 ][ilS ];
          for (int kb1=ji; kb1!=jb1; kb1-=is) {
            dc1 += e[kb1][il1c][ilS ];
          }
          for (int kbS=ji; kbS!=jbS; kbS-=is) {
            dcS += e[kbS][il1 ][ilSc];
            dc1S+= e[kbS][il1c][ilSc];
          }
          d[ii][il1][ilS] = min4(dc1,dcS,dc1S,di)+e[ii][il1][ilS];
        }
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
      u[ii] = (float)il;//shifts.getValue(il);
      while (ii!=ie) {
        ii += is;
        il += (int)m[ii][il];
        u[ii] = (float)il;//shifts.getValue(il);
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

  private static float[] inverseInterpolation(float[] u, float[] t) {
    int n1u = u.length;
    int n1um = n1u-1;
		float[] r = add(u,t);
    int n2 = (int)(ceil(r[n1um])-floor(r[0]));
		float[] s = new float[n2];
   	inverseLinearInterpolation(n1u,1f,0f,r,n2,1f,u[0],s,r[0],r[n1um]);
		return s;
	}

  private static float[][] inverseInterpolation(final float[][] u) {
    final int n1u = u[0].length;
    final int n2u = u.length;
    final int n1um = n1u-1;
		final float[][] s = new float[n2u][n1u];
		final float[] t = rampfloat(0f,1f,n1u);
    Parallel.loop(n2u,new Parallel.LoopInt() {
    public void compute(int i2) {
    	//int n2 = (int)(ceil(r[i2][n1um])-floor(r[i2][0]));
			float[] r = add(u[i2],t);
    	inverseLinearInterpolation(n1u,1f,0f,r,n1u,1f,u[i2][0],s[i2],r[0],r[n1um]);
		}});
		return s;
	}

	private static float[][][] inverseInterpolation(final float[][][] u) {
    final int n1u = u[0][0].length;
    final int n2u = u[0].length;
    final int n3u = u.length;
    final int n1um = n1u-1;
		final float[][][] s = new float[n3u][n2u][n1u];
		final float[] t = rampfloat(0f,1f,n1u);
    Parallel.loop(n3u,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2u; i2++) {
    		//int n2 = (int)(ceil(r[i3][i2][n1um])-floor(r[i3][i2][0]));
				float[] r = add(t,u[i3][i2]);
    		inverseLinearInterpolation(n1u,1f,0f,r,n1u,1f,u[i3][i2][0],s[i3][i2],r[0],r[n1um]);
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

};

