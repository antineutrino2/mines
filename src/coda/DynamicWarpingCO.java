package coda;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.interp.CubicInterpolator.Method;
import static edu.mines.jtk.util.ArrayMath.*;


// for testing only:
import java.awt.Color;
import edu.mines.jtk.awt.ColorMap;

public class DynamicWarpingCO {

	public enum Interpolation {
		LINEAR,
		MONOTONIC,
		SPLINE,
	};

  /**
   * Constructs a dynamic warping for specified bounds on shifts.
   * Shifts vary on the coarse interval given by 1/ds
   * @param nl number of lags where nl = (n2-n1)+1.
   * @param smin lower bound on shift derivative u'
   * @param smax upper bound on shift derivative u'
	 * @param scgs coarse grid sampling
   */
  public DynamicWarpingCO(
		int smin, int smax, double rmin, double rmax, double dr) 
	{
    Check.argument(rmax>=0,"rmax>=0");
    Check.argument(rmin<=0,"rmin<=0");
    Check.argument(dr<=1.0,"dr<=1.0");
		int nl = 1+smax-smin;
		_scgs = (int)ceil(1/dr);
		_nl = nl;
		_rmin = rmin;
		_rmax = rmax;
		_dr = dr;
		_ss = new Sampling(nl,1.0,smin);
		_si = new SincInterp(); 
  }

	public void setStrainLimits(double rmin, double rmax) {
    Check.argument(rmax>=0,"rmax>=0");
    Check.argument(rmin<=0,"rmin<=0");
		_rmin = rmin;
		_rmax = rmax;
	}

	public void setSmoothing(double dr) {
    Check.argument(dr<=1.0,"dr<=1.0");
		_scgs = (int)ceil(1/dr);
		_dr = dr;
	}

	public void setInterpolation(Interpolation interp) {
		if (interp==Interpolation.SPLINE) {
			interp = Interpolation.SPLINE;
		} else if (interp==Interpolation.MONOTONIC) {
			interp = Interpolation.MONOTONIC;
		} else if (interp==Interpolation.LINEAR) {
			interp = Interpolation.LINEAR;
		}
    _interp = interp;
  }

  /**
   * Sets the exponent used to compute alignment errors |f-g|^e.
   * The default exponent is 2.
   * @param e the exponent.
   */
  public void setErrorExponent(double e) {
    _epow = (float)e;
  }

 /**
   * Returns normalized alignment errors for all samples and fractional lags.
   * @param f array[n1] for the sequence f[i1].
   * @param g array[n2] for the sequence g[i2].
   * @return array[n1][nl] of alignment errors.
   */
  public float[][] computeErrors(float[] f, float[] g) {
    int n1 = f.length;
    int n2 = g.length;
    int nl = _nl;
    float[][] e = new float[n1][nl];
		computeErrors(f,g,e);
    normalizeErrors(e);
    return e;
  }


	/**
   * Computes shifts for specified sequences.
   * @param e array of alignment errors.
   * @param d array of accumulated alignment errors.
   * @param u output array of shifts u.
   */
  public float[] findShiftsR(float[][] e) {
		int nl = e[0].length;
		int nt = e.length;
		int dc = _scgs;
		int[] cgs = subsample(nt,dc); _g1 = cgs;
		int nc = cgs.length;
		float[][] d = new float[nc][nl];
		int[][] m = new int[nc][nl];
    float[] uc = new float[nc];
		accumulateR(_rmin,_rmax,cgs,e,d,m,1,-1);
		backtrack(d,m,uc,-1);
		float[] u = interpolateSparseShifts(nt,cgs,uc);
		return u;
  }
  public float[] findShiftsR(float[] f, float[] g, float[][] e) {
		int nl = e[0].length;
		int nt = e.length;
		int dc = _scgs;
		HilbertTransformFilter ht = new HilbertTransformFilter();
		float[] ft = new float[nt];
		ht.apply(nt,f,ft);
		float[] h = sqrt(add(pow(f,2.0f),pow(ft,2.0f)));
		int lc = findlc1(f,g,dc);
		int[] cgs = subsample(h,dc,(float)lc); _g1 = cgs;
		int k = findlc2(lc,cgs);
		int nc = cgs.length;
		float[][] d = new float[nc][nl];
		int[][] m = new int[nc][nl];
    float[] uc = new float[nc];
    int zl = _ss.indexOfNearest(0);
		accumulateR(_rmin,_rmax,cgs,e,d,m,k,zl);
		//accumulateR(_rmin,_rmax,cgs,e,d,m,zl,zl);
    //int[] ps = findlc3(f,g);
		//accumulateR(_rmin,_rmax,cgs,e,d,m,ps[0],ps[1]);
		backtrack(d,m,uc,zl);
		float[] u = interpolateSparseShifts(nt,cgs,uc);
		return u;
  }
  public float[] findShiftsR(float[] f, float[] g, float[][] e, int ngg) {
		int nl = e[0].length;
		int nt = e.length;
		int dc = _scgs;
		HilbertTransformFilter ht = new HilbertTransformFilter();
		float[] ft = new float[nt];
		ht.apply(nt,f,ft);
		float[] h = sqrt(add(pow(f,2.0f),pow(ft,2.0f)));
		int lc = findlc1(f,g,dc);
    lc = 0;
		int[] cgs = subsample(h,dc,ngg,(float)lc); _g1 = cgs;
		int k = findlc2(lc,cgs);
		int nc = cgs.length;
		float[][] d = new float[nc][nl];
		int[][] m = new int[nc][nl];
    float[] uc = new float[nc];
    int zl = _ss.indexOfNearest(0);
    zl = -1;
		accumulateR(_rmin,_rmax,cgs,e,d,m,zl,zl);
		backtrack(d,m,uc,zl);
		float[] u = interpolateSparseShifts(nt,cgs,uc);
		return u;
  }
	private int findlc1(float[] f, float[] g, int dc) {
    int nf = f.length;
    float mf = sum(abs(f))/nf;
    int lc = 0;
		for (int i=0; i<nf; ++i) {
      if (f[i]<0.1*mf) lc = i-dc;
      else if (g[i]<0.1*mf) lc = i-dc;
      else break;
    }
    return lc;
  }
	private int findlc2(int lc, int[] cgs) {
		int nc = cgs.length;
		int k=1;
		for (int i=0; i<nc; ++i) {
			int j = cgs[i];
			if (j>=lc) {
				k=i;
				break;
			}
		}
		return k;
	}
	private int[] findlc3(float[] f, float[] g) {
    int nf = f.length;
    float mf = sum(abs(f))/nf;
    int ng = g.length;
    float mg = sum(abs(g))/ng;
    int lc1 = 0;
    int lc2 = 0;
		for (int i=0; i<nf; ++i) {
      if (f[i]<0.01*mf) lc1 = i;
      else break;
    }
		for (int i=0; i<ng; ++i) {
      if (g[i]<0.01*mg) lc2 = i;
      else break;
    }
    return new int[]{lc1,lc1-lc2};
  }


	public float[][] accumulate(float[] f, float[] g, float[][] e) {
		int nl = e[0].length;
		int nt = e.length;
		int dc = _scgs;
		HilbertTransformFilter ht = new HilbertTransformFilter();
		float[] ft = new float[nt];
		ht.apply(nt,f,ft);
		float[] h = sqrt(add(pow(f,2.0f),pow(ft,2.0f)));
		int lc = findlc1(f,g,dc);
		int[] cgs = subsample(h,dc,(float)lc); _g1 = cgs;
		int k = findlc2(lc,cgs);
		int nc = cgs.length;
		float[][] d = new float[nc][nl];
		int[][] m = new int[nc][nl];
    float[] uc = new float[nc];
    int zl = (int)abs(_ss.getFirst());
		accumulateR(_rmin,_rmax,cgs,e,d,m,k,zl);
		return d;
	}

  /**
   * Computes a sequence warped by applying specified shifts.
   * @param u input array of shifts.
   * @param f input array for the sequence to be warped.
   * @param h output array for the warped sequence.
   */
  public float[] applyShifts(float[] f, float[] u) {
    int n1 = u.length;
    int nf = f.length;
		double ssf = _ss.getFirst();
    float[] h = new float[n1];
    for (int i1=0; i1<n1; ++i1)
      h[i1] = _si.interpolate(nf,1.0,0.0,f,i1+u[i1]);
    return h;
  }

	/**
   * Returns an approximately uniformly-sampled subset of indices in [0,n).
   * Indices in the subset are chosen to be approximately uniform, with the
   * difference between consecutive indices not less than the specified
   * minimum increment kmin. Because the first and last indices 0 and n-1 are
   * included in the subset, n must be greater than the minimum increment
   * kmin.
   * @param n number of indices in the set {0,1,2,...,n-1}.
   * @param kmin minimum increment between indices in the subset.
   * @return array of indices in the subset.
	 * @author Dave Hale, CSM
   */
  private static int[] subsample(int n, int kmin) {
    if (kmin>=n)
      kmin = n-1;
    int m = 1+(n-1)/kmin;
    double d = (double)(n-1)/(double)(m-1);
    int[] j = new int[m];
    for (int i=0; i<m; ++i)
      j[i] = (int)(i*d+0.5);
    return j;
  }

  /**
   * Uses cubic interpolation to interpolate u, sampled on grid g,
   * to a uniformly sampled array of length n.
   * @param n
   * @param g
   * @param u
   * @return the interpolated shifts.
	 * @author Stefan Compton, CSM
   */
  public float[] interpolateSparseShifts(
      int n, int[] g, float[] u) 
  {
    int ng = g.length;
    float fs = (float)_ss.getFirst();
    float[] gf = new float[ng];
    float[] ui = new float[n ];
    for (int ig=0; ig<ng; ig++)
      gf[ig] = (float)g[ig];
		Method m = CubicInterpolator.Method.SPLINE;
		if (_interp==Interpolation.SPLINE) {
			m = CubicInterpolator.Method.SPLINE;
		} else if (_interp==Interpolation.MONOTONIC) {
			m = CubicInterpolator.Method.MONOTONIC;
		} else if (_interp==Interpolation.LINEAR) {
			m = CubicInterpolator.Method.LINEAR;
		}
    CubicInterpolator ci = new CubicInterpolator(m,ng,gf,u);
    int fg = g[0];
    for (int i=fg; i<n; i++)
      ui[i] = ci.interpolate(i)+fs;
    return ui;
  }

  public static int[] subsample(float[] f, int d, float lcf) {
    int lc = (int)lcf;
    int n = f.length;
    int nm = n-1;
    int[] i = rampint(0,1,n);
    float[] fa = abs(f);
    quickIndexSort(fa,i);
//    StringBuffer sb = new StringBuffer();
//    for (int ii=0; ii<10; ii++) {
//      sb.append(i[nm-ii]+",");
//    }
//    System.out.println("Top 10: "+sb.toString());
    
    // Get maximum amplitute indices as long as they
    // are greater than or equal to d.
    List<Integer> gl = new ArrayList<Integer>();
    gl.add(0); gl.add(nm);
    if (lc>0) gl.add(lc);
    for (int ii=0; ii<n; ii++) {
      int im = i[nm-ii]; // next highest amplitude index
      if (im>lc || lc==0) {
        int s = gl.size();
        boolean add = true;
        for (int ig=0; ig<s; ig++) {
          if (abs(im-gl.get(ig))<d) {
            add = false;
            break;
          }
        }
        if (add)
          gl.add(im);
      }
    }
    Collections.sort(gl);
    int nl = gl.size();
    int[] g = new int[nl];
    for (int il=0; il<nl; il++) {
      g[il] = gl.get(il);
    }
    return g;
  }

 public static int[] subsample(float[] f, int d, int ng, float lcf) {
    int lc = (int)lcf;
    int n = f.length;
    int nm = n-1;
    int[] i = rampint(0,1,n);
    float[] fa = abs(f);
    quickIndexSort(fa,i);
//    StringBuffer sb = new StringBuffer();
//    for (int ii=0; ii<10; ii++) {
//      sb.append(i[nm-ii]+",");
//    }
//    System.out.println("Top 10: "+sb.toString());
    
    // Get maximum amplitude indices, adjusting the interval
    // d so that the subsampled array contains ng samples.
    List<Integer> gl = new ArrayList<Integer>();
    for (; d>=1; d--) {
      gl.clear(); 
      gl.add(0);  gl.add(nm); 
      if (lc>0) gl.add(lc); 
      for (int ii=0; ii<n; ii++) {
        int im = i[nm-ii]; // next highest amplitude index
        if (im>lc || lc==0) {
          int s = gl.size();
          if (s==ng) // size matches desired sparse array size, quit
            break;
          boolean add = true;      
          for (int ig=0; ig<s; ig++) {
            if (abs(im-gl.get(ig))<d) {
              add = false;
              break;
            }
          }
          if (add)
            gl.add(im);
        }
      }
      if (gl.size()==ng) { // size matches desired sparse array size, quit
        //System.out.println("final d="+d);
        break;
      }
    }
    Collections.sort(gl);
    int nl = gl.size();
    int[] g = new int[nl];
    for (int il=0; il<nl; il++) {
      g[il] = gl.get(il);
    }
    return g;
  }
  
	public int[] getG1() {
		return _g1;
	}

  //////////////////////////////////////////////////////////////////////////////
  // private

  private int _nl; // number of lags
  private float _epow = 2; // exponent used for alignment errors |f-g|^e
	private int _scgs; // coarse grid sampling
  private double _rmin; // upper int strain limit of u'
  private double _rmax; // lower int strain limit of u'
	private double _dr; // one over the grid sampling
	private int[] _g1;
	private SincInterp _si;
	private Sampling _ss;
  private Interpolation _interp = Interpolation.SPLINE;

  private float error(float f, float g) {
    return pow(abs(f-g),_epow);
  }
	
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

//  private static void accumulateR(
//    double rmin, double rmax, int[] cgs,
//    float[][] e, float[][] d, int[][] m, int lc) 
//  {
//    accumulateR(rmin,rmax,cgs,e,d,m,lc,-1);
//  }
//   /**
//	 * Accumulation of errors on a variably coarse grid.
//	 * Determined by the array cgs. The max and minimum 
//	 * accumulation slopes correspond to rmin and rmax.
//	 */
//  private static void accumulateR(
//    double rmin, double rmax, int[] cgs,
//    float[][] e, float[][] d, int[][] m, int lc, int zl) 
//  {
//    int nl = e[0].length;
//    int nt = e.length;
//		int nc = cgs.length;
//    int nlm1 = nl-1;
//		// Compute d[0]
//    for (int il=0; il<nl; ++il)
//			d[0][il] = e[0][il];
//    float mx = max(e); 
//    if (zl>=0) {
//      for (int it=0; it<lc; ++it)
//        for (int i1=0; i1<nl; ++i1)
//          e[it][i1] = mx;
//      for (int it=0; it<lc; ++it)
//        e[it][zl] = 0f;
//    }
//    // Compute d[1,nlm1]
//		// Loop over sparse grid points
//		float me = max(e)*nt; // default max value for d
//    for (int ic=1; ic<nc; ++ic) {
//      int icm1 = ic-1;
//			int cgse = cgs[ic];
//			int cgsb = cgs[icm1];
//			int cgsd = cgse-cgsb;
//			int kmax = (int)ceil( rmax*cgsd);
//			int kmin = (int)floor(rmin*cgsd);
//			if (ic<=lc || ic==nc-1) {
//				kmin = 0;
//				kmax = 0;
//			}
//			// Loop over all lags
//      for (int il=0; il<nl; ++il) {
//				int mi = 0;
//				float dm = me;
//				// Loop over all slopes
//        for (int k=kmin; k<=kmax; ++k) {
//	  			int ik = il+k; 
//	  			if (0<=ik && ik<nl) {
//						float dp = d[icm1][ik];
//	    			float dv = interpSlope(k,cgsd,il,cgsb,cgse,dp,e);
//        		if (dv<dm) {
//	      		  dm = dv;
//	      		  mi = k;
//						}
//	      	}
//	  		}
//				m[ic][il] = mi;
//				d[ic][il] = e[cgse][il] + dm;  
//			}
//    }
//  }
//	/**
//	 * @author Stefan Compton, CSM
//	 */
//	private static float interpSlope(   
//    int k, int dg, int il, int cgsb, int cgse, float dp, float[][] e)
//  {
//    float dc = dp;
//    if (k==0) { 
//      for (int x=cgsb+1; x<cgse; ++x) {
//        dc += e[x][il]; 
//      }
//    } else {
//			double slope = -(double)k/dg;
//			// Loop over all points along slope
//      for (int x=cgsb+1; x<cgse; ++x) {
//        double y = il+(slope*(x-cgsb)+k);
//        int y1 = (int)y;
//        int y2 = y1+1;
//        dc += (y2-y)*e[x][y1]+(y-y1)*e[x][y2];
//      }
//    }
//    return dc;
//  }
//  
//  /**
//   * Finds shifts by backtracking in accumulated alignment errors (fractional).
//   */
//  private static void backtrack(
//    float[][] d, int[][] m, float[] u) 
//  {
//    int mi = 0;
//    int nl = d[0].length;
//    int ni = d.length;
//    int nlm1 = nl-1;
//    int il = 0;
//    int ii = ni-1;
//    float dl = d[ii][il];
//    for (int jl=1; jl<nl; ++jl) {
//      if (d[ii][jl]<dl) {
//        dl = d[ii][jl];
//        il = jl;
//      }
//    }
//    u[ii] = il;
//    while (ii>0) { 
//      mi = m[ii][il];
//      il += mi;
//      --ii;
//      u[ii] = il;
//      if (mi!=0) {
//        float du = u[ii]-u[ii+1];
//        u[ii] = u[ii+1]+du;
//      }
//    }
//    u[0]=u[1];
//  }

  /**
	 * Accumulation of errors on a variably coarse grid.
	 * Determined by the array cgs. The max and minimum 
	 * accumulation slopes correspond to rmin and rmax.
	 * Accumulate maximum values.
	 */
  private static void accumulateR(
    double rmin, double rmax, int[] cgs,
    float[][] e, float[][] d, int[][] m, int p1, int p2) 
  {
    int dir = -1;
    int nl = e[0].length;
    int nt = e.length;
		int nc = cgs.length;
    int nlm1 = nl-1;
		float mx = max(e); 
    if (p1>=0) {
      for (int i1=0; i1<nl; ++i1)
        e[p1][i1] = mx;
      e[p1][p2] = 0f;
    }
    int ibg = (dir>0)?0:nc-1; // beginning index
    int ieg = (dir>0)?nc:0;  // end index
    int is = (dir>0)?1:-1;   // stride
    int ic = ibg; // sparse grid index
    int ie = cgs[ic]; // error index
		// Compute d[0]
    for (int il=0; il<nl; ++il)
			d[ic][il] = e[ie][il];
    // Compute d[1,nlm1]
		// Loop over sparse grid points
		float me = max(e)*nt; // default min value for d
    ic += is;
    // Loop over all sparse grid points.
    for (; ic!=ieg; ic+=is) {
      int icm1 = ic-is; // previous sparse grid index
			int cgse = cgs[ic]; // new error index
			int cgsb = cgs[icm1]; // min error index for interpolation.
			int cgsd = cgse-cgsb; // sparse grid delta, possibly negative
      int kmax,kmin;
      if (cgsd>0) {
			  kmax = (int)ceil( -rmin*cgsd);
			  kmin = (int)floor(-rmax*cgsd);
      } else {
			  kmax = (int)ceil( -rmax*cgsd);
			  kmin = (int)floor(-rmin*cgsd);
      }
			// Loop over all lags
      for (int il=0; il<nl; ++il) {
				int mi = 0;
				float dm = me;
				// Loop over all slopes
        for (int k=kmin; k<=kmax; ++k) {
	  			int ik = il+k; 
	  			if (0<=ik && ik<nl) {
						float dp = d[icm1][ik];
	    			float dv = interpSlope(k,cgsd,il,cgsb,cgse,dp,e);
        		if (dv<dm) {
	      		  dm = dv;
	      		  mi = k;
						}
	      	}
	  		}
				m[ic][il] = mi;
				d[ic][il] = e[cgse][il] + dm;  
			}
    }
  }
	/**
	 * @author Stefan Compton, CSM
	 */
	private static float interpSlope(   
    int k, int dg, int il, int cgsb, int cgse, float dp, float[][] e)
  {
    float dc = dp;
    if (k==0) { 
      for (int x=cgsb+1; x<cgse; ++x) {
        dc += e[x][il]; 
      }
    } else {
			double slope = -(double)k/dg;
			// Loop over all points along slope
      for (int x=cgsb+1; x<cgse; ++x) {
        double y = il+(slope*(x-cgsb)+k);
        int y1 = (int)y;
        int y2 = y1+1;
        dc += (y2-y)*e[x][y1]+(y-y1)*e[x][y2];
      }
    }
    return dc;
  }
  
  /**
   * Finds path by backtracking in accumulated alignment errors.
   */
  private static void backtrack(
    float[][] d, int[][] m, float[] u, int p1) 
  {
    int dir = -1;
    int mi = 0;
    int nl = d[0].length;
    int ni = d.length;
    int nim1 = ni-1;
    int nlm1 = nl-1;
    int il = 0;
    int ib = (dir>0)?nim1:0;
    int ie = (dir>0)?1:nim1;
    int is = (dir>0)?1:-1;
    int ii = ib;
    float dl = d[ii][il];
    for (int jl=1; jl<nl; ++jl) {
      if (d[ii][jl]>dl) {
        dl = d[ii][jl];
        il = jl;
      }
    }
    if (dir<0 && p1>=0) il = p1;
    u[ii] = il;
    if (dir<0) {
      ii-=is;
      u[ii] = il;
    }
    while (ii!=ie) { 
      mi = m[ii][il];
      il += mi;
      ii-=is;
      u[ii] = il;
      if (mi!=0) {
        float du = u[ii]-u[ii+is];
        u[ii] = u[ii+is]+du;
      }
    }
    //if (dir>0) u[0] = u[1];
  }



/////////////////////////////////////////////////////////
// static utilities

  private static float errors(float f, float g, float epow) {
    return pow(abs(f-g),epow);
  }

  private static float[][] castf(int[][] x) {
    int n1 = x.length;
    int n2 = x[0].length;
    float[][] y = new float[n1][n2];
    for (int i1=0; i1<n1; ++i1) {
      for (int i2=0; i2<n2; ++i2)
        y[i1][i2] = (float)x[i1][i2];
    }
    return y;
  }

  /**
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
   */
  private static void normalizeErrors(float[][] e) {
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
   * Shifts and scales alignment errors to be in range [0,1].
   * @param emin minimum alignment error before normalizing.
   * @param emax maximum alignment error before normalizing.
   * @param e input/output array of alignment errors.
   */
  private static void shiftAndScale(float emin, float emax, float[][] e) {
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

  private static float min3(float a, float b, float c) {
    return b<=a?(b<=c?b:c):(a<=c?a:c); // if equal, choose b
  }
  
  private static float minf(float[] a, float b, float[] c) {
    float ma = min(a); float mc = min(c);
    return b<=ma?(b<=mc?b:mc):(ma<=mc?ma:mc); // if equal, choose b
  }
   
  private static int mini(float[] a, float b, float[] c) {
    int na = a.length; 
    int nc = c.length; 
    int mi=0,mia=0,mic=0;
    float m=b,ma=a[0],mc=c[0];
    for (int i=0; i<na; ++i) {
      if (ma>a[i]) {
        ma = a[i];
	mia = i;
      }
    }
    for (int i=0; i<nc; ++i) {
      if (mc>c[i]) {
        mc = c[i];
	mic = i;
      }
    }
    if (mc>b  && ma>b  && m>b ) {m=b ; mi= 0    ;} 
    if (b >ma && mc>ma && m>ma) {m=ma; mi=-mia-1;}
    if (b >mc && ma>mc && m>mc) {m=mc; mi= mic+1;}
    return mi; // if equal, choose b
  }

  private static float[] like(float[] a) {
    return new float[a.length];
  }
  private static float[][] like(float[][] a) {
    return new float[a.length][a[0].length];
  }

  private float[] floats(int[] x) {
    int n1 = x.length;
    float[] y = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      y[i1] = (float)x[i1];
    }
    return y;
  }

	private float[][] floats(int[][] x) {
    int n1 = x.length;
    int n2 = x[0].length;
    float[][] y = new float[n1][n2];
    for (int i1=0; i1<n1; ++i1) {
    	for (int i2=0; i2<n1; ++i2) 
      	y[i1][i2] = (float)x[i1][i2];
    }
    return y;
  }


};

