package sem;
import java.util.List;
import java.util.ArrayList;
import java.util.Collections;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;

// for testing only:
import java.awt.Color;
import edu.mines.jtk.awt.ColorMap;

public class DynamicWarpingSolver {

  /**
   * Constructs a dynamic solver for specified bounds on shifts.
   * Shifts vary on the coarse interval given by 1/ds
   * @param nl number of lags 
   * @param rmin lower bound on shift derivative u'
   * @param rmax upper bound on shift derivative u'
	 * @param scgs coarse grid sampling
   */
  public DynamicWarpingSolver(
    int nl, double rmin, double rmax, double dr, int wb) 
  {
    Check.argument(rmax>=0,"rmax>=0");
    Check.argument(rmin<=0,"rmin<=0");
    Check.argument(dr<=1.0,"dr<=1.0");
		if (dr> min(rmax,-rmin)) dr = min(rmax,-rmin);
    _nl1 = nl;
		_scgs = (int)ceil(1/dr);
		_rmin1 = rmin;
		_rmax1 = rmax;
		_dr = dr;
    _wb = wb;
  }
  public DynamicWarpingSolver(int nl, double rmin, double rmax, double dr) {
		this(nl,rmin,rmax,1.0,0);
  }
  public DynamicWarpingSolver(int nl, double rmin, double rmax) {
		this(nl,rmin,rmax,1.0,0);
	}
  public DynamicWarpingSolver(int nl) {
		this(nl,-1.0,1.0,1.0,0);
	}
	/**
   * Constructs a dynamic warping for specified bounds on shifts.
   * Shifts vary on the coarse interval given by 1/ds
   * @param nl1,nl2 number of lags 
   * @param rmin1,rmin2 lower bound on shift derivative u1' and u2'
   * @param rmax1,rmax2 upper bound on shift derivative u1' and u2'
	 * @param scgs coarse grid sampling
   */
  public DynamicWarpingSolver(
		int nl1, int nl2, double rmin1, double rmax1, 
		double rmin2, double rmax2, double dr, int wb) 
	{
    Check.argument(rmax1>=0,"rmax1>=0");
    Check.argument(rmin1<=0,"rmin1<=0");
    Check.argument(rmax2>=0,"rmax2>=0");
    Check.argument(rmin2<=0,"rmin2<=0");
    Check.argument(dr<=1.0,"dr<=1.0");
		if (dr> min(rmax1,-rmin1)) dr = min(rmax1,-rmin1);
		if (dr> min(rmax2,-rmin2)) dr = min(rmax2,-rmin2);
    _nl1 = nl1;
    _nl2 = nl2;
		_scgs = (int)ceil(1/dr);
		_rmin1 = rmin1;
		_rmax1 = rmax1;
		_rmin2 = rmin2;
		_rmax2 = rmax2;
		_dr = dr;
		_wb = wb;
  }
  public DynamicWarpingSolver(
		int nl1, int nl2, double rmin1, double rmax1, double rmin2, double rmax2)
	{
		this(nl1,nl2,rmin1,rmax1,rmin2,rmax2,1.0,0);
	}
  public DynamicWarpingSolver(int nl1, int nl2) {
		this(nl1,nl2,-1.0,1.0,-1.0,1.0,1.0,0);
	}


  ///**
  // * Sets bounds on strains.
  // * Default lower and upper bounds are -1.0 and 1.0, respectively.
  // * @param r1min lower bound on strain in units of s1.
  // * @param r1max upper bound on strain in units of s1.
  // */
  //public void setStrainLimits(double r1min, double r1max) {
  //  setStrainLimits(r1min,r1max,-1.0,1.0,-1.0,1.0);
  //}

  ///**
  // * Sets bounds on strains.
  // * Default lower and upper bounds are -1.0 and 1.0, respectively.
  // * @param r1min lower bound on strain in units of s1.
  // * @param r1max upper bound on strain in units of s1.
  // * @param r2min lower bound on strain in units of s2.
  // * @param r2max upper bound on strain in units of s2.
  // */
  //public void setStrainLimits(
  //  double r1min, double r1max,
  //  double r2min, double r2max)
  //{
  //  _r1min = (float)r1min; _r1max = (float)r1max;
  //  _r2min = (float)r2min; _r2max = (float)r2max;
  //}

	///**
	// * Sets the smoothing parameter.
	// * These parameters determine the time grid spacing in 
	// * smooth dynamic time warping. Default is 1.0.
	// * @param dr1 time smoothing parameter
	// */
	//public void setSmoothing(double dr1) {
	//	setSmoothing(dr1,1.0);
	//}

	///**
	// * Sets the smoothing parameters.
	// * These parameters determine the grid spacings for each dimension in 
	// * smooth dynamic image warping. Defaults are 1.0.
	// * @param dr1 time smoothing parameter
	// * @param dr2 distance smoothing parameter
	// * @param dr3 distance smoothing parameter
	// */
	//public void setSmoothing(double dr1, double dr2) {
  //  Check.argument(dr1<=1.0,"dr1<=1.0");
  //  Check.argument(dr2<=1.0,"dr2<=1.0");
	//	_dg1 = (int)ceil(1/dr1);
	//	_dg2 = (int)ceil(1/dr2);
	//}

  //public setWaterBottom(double wb) {
  //  _wb = (float)wb;
  //}
  public void setWbParam(int p1) {
    setWbParam(p1,0);
  }
  public void setWbParam(int p1, int p2) {
    _p1 = p1;
    _p2 = p1;
  }

	/**
   * Computes shifts for specified sequences.
   * @param e array of alignment errors.
   * @param d array of accumulated alignment errors.
   * @param u output array of shifts u.
   */
  public float[] findSolution(float[][] e) {
		int nl = e[0].length;
		int nt = e.length;
		int dc = _scgs;
		int[] cgs = subsample(nt,dc,(float)_wb);
		int nc = cgs.length;
		float[][] d = new float[nc][nl];
		int[][] m = new int[nc][nl];
    float[] uc = new float[nc];
		accumulate(_rmin1,_rmax1,cgs,e,d,m,_p1,_wb);
		backtrack(d,m,uc);
		float[] u = interpolateSparseShifts(nt,cgs,uc);
		return u;
  }

  public float[] findSolution(float[] f, float[][] e) {
		int nl = e[0].length;
		int nt = e.length;
		int dc = _scgs;
		//HilbertTransformFilter ht = new HilbertTransformFilter();
		//float[] ft = new float[nt];
		//ht.apply(nt,f,ft);
		//float[] h = sqrt(add(pow(f,2.0f),pow(ft,2.0f)));
		//int[] cgs = subsample(f,dc,15,(float)_wb);
		int[] cgs = subsample(f,dc,(float)_wb);
		int nc = cgs.length;
		float[][] d = new float[nc][nl];
		int[][] m = new int[nc][nl];
    float[] uc = new float[nc];
		accumulate(_rmin1,_rmax1,cgs,e,d,m,_p1,_wb);
		backtrack(d,m,uc);
		float[] u = interpolateSparseShifts(nt,cgs,uc);
		return u;
  }

	public float[][] accumulate(float[][] e) {
		int nl = e[0].length;
		int nt = e.length;
		int dc = _scgs;
		int[] cgs = subsample(nt,dc,(float)_wb);
		int nc = cgs.length;
		float[][] d = new float[nc][nl];
		int[][] m = new int[nc][nl];
    float[] uc = new float[nc];
		accumulate(_rmin1,_rmax1,cgs,e,d,m,_p1,_wb);
		return d;
	}

	/**
   * Computes shifts for specified sequences.
   * @param e array of alignment errors.
   * @param u output array of shifts u.
   */
  public float[][] findSolution(float[][][] e) {
		int nl1 = e[0][0].length;
		int nl2 = e[0].length;
		int nt = e.length;
		int dc = _scgs;
		int[] cgs = subsample(nt,dc,(float)_wb);
    //SimplePlot.asPoints(floats(cgs));
		int nc = cgs.length;
		float[][][] d = new float[nc][nl2][nl1];
		int[][][][] m = new int[nc][nl2][nl1][2];
    float[] uc1 = new float[nc];
    float[] uc2 = new float[nc];
		accumulate(_rmin1,_rmax1,_rmin2,_rmax2,cgs,e,d,m,_p1,_p2,_wb);
		backtrack(d[nc-1],m,uc1,uc2);
		// Check slopes
    //float lastLag = nl1;
    //for (int ig1=1; ig1<nc; ig1++) {
    //  int i1 = cgs[ig1];
    //  int i1m1 = cgs[ig1-1];
    //  float na = uc1[ig1] - uc1[ig1-1];
    //  float da = i1-i1m1;
    //  float ra = na/da;
    //  assert (ra+0.0001>=_rmin1 && ra-0.0001<=_rmax1) || uc1[ig1]==lastLag :
    //    "u1: "+"n="+na+", d="+da+", r="+ra;
    //}
    //lastLag = nl2;
    //for (int ig1=1; ig1<nc; ig1++) {
    //  int i1 = cgs[ig1];
    //  int i1m1 = cgs[ig1-1];
    //  float na = uc2[ig1] - uc2[ig1-1];
    //  float da = i1-i1m1;
    //  float ra = na/da;
    //  assert (ra+0.001>=_rmin2 && ra-0.001<=_rmax2) || uc2[ig1]==lastLag :
    //    "u2: "+"n="+na+", d="+da+", r="+ra;
    //}
		float[] u1 = interpolateSparseShifts(nt,cgs,uc1);
		float[] u2 = interpolateSparseShifts(nt,cgs,uc2);
		return new float[][]{u1,u2};
  }

  public float[][] findSolution(float[] f, float[][][] e) {
    int nl1 = e[0][0].length;
    int nl2 = e[0].length;
    int nt = e.length;
    int dc = _scgs;
    int[] cgs = subsample(f,nt,dc,(float)_wb);
    SimplePlot.asPoints(floats(cgs));
    int nc = cgs.length;
    float[][][] d = new float[nc][nl2][nl1];
    int[][][][] m = new int[nc][nl2][nl1][2];
    float[] uc1 = new float[nc];
    float[] uc2 = new float[nc];
    accumulate(_rmin1,_rmax1,_rmin2,_rmax2,cgs,e,d,m,_p1,_p2,_wb);
    backtrack(d[nc-1],m,uc1,uc2);
    float[] u1 = interpolateSparseShifts(nt,cgs,uc1);
    float[] u2 = interpolateSparseShifts(nt,cgs,uc2);
    return new float[][]{u1,u2};
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
  private static int[] subsample(int n, int kmin, float wbf) {
    int wb = (int)wbf;
    n -= wb;
    if (kmin>=n)
      kmin = n-1;
    int m = 1+(n-1)/kmin;
    double d = (double)(n-1)/(double)(m-1);
    int[] j = new int[m];
    wb -= (wb>0)?(int)d:0;
    for (int i=0; i<m; ++i)
      j[i] = wb+(int)(i*d+0.5);
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
    float[] gf = new float[ng];
    float[] ui = new float[n ];
    for (int ig=0; ig<ng; ig++)
      gf[ig] = (float)g[ig];
    CubicInterpolator ci = new CubicInterpolator(
			CubicInterpolator.Method.LINEAR,ng,gf,u);
			//CubicInterpolator.Method.SPLINE,ng,gf,u);
			//CubicInterpolator.Method.MONOTONIC,ng,gf,u);
    for (int i=0; i<n; i++)
      ui[i] = ci.interpolate(i);
    return ui;
  }

  public static int[] subsample(float[] f, int d, float wbf) {
    int wb = (int)wbf;
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
    gl.add(0); gl.add(nm); gl.add(wb);
    for (int ii=0; ii<n; ii++) {
      int im = i[nm-ii]; // next highest amplitude index
      if (im>wb) {
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

 public static int[] subsample(float[] f, int d, int ng, float wbf) {
    int wb = (int)wbf;
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
      gl.add(0); gl.add(nm); gl.add(wb);
      for (int ii=0; ii<n; ii++) {
        int im = i[nm-ii]; // next highest amplitude index
        if (im>wb) {
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
        System.out.println("final d="+d);
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
   
  //////////////////////////////////////////////////////////////////////////////
  // private

  private int _nl1,_nl2; // number of lags
	private int _scgs; // coarse grid sampling
  private double _rmin1,_rmin2; // upper int strain limit of u'
  private double _rmax1,_rmax2; // lower int strain limit of u'
	private double _dr; // one over the grid sampling
	private int _wb=0; // wb depth
	private int _p1=-1; // wb parameter 1
	private int _p2=-1; // wb parameter 2

  /**
	 * Accumulation of errors on a variably coarse grid.
	 * Determined by the array cgs. The max and minimum 
	 * accumulation slopes correspond to rmin and rmax.
	 * Accumulate maximum values.
	 */
  private static void accumulate(
    double rmin, double rmax, int[] cgs,
    float[][] e, float[][] d, int[][] m, int p1, int wb) 
  {
    int nl = e[0].length;
    int nt = e.length;
		int nc = cgs.length;
    int nlm1 = nl-1;
    if (p1>=0) {
      for (int it=0; it<wb; ++it)
       for (int i1=0; i1<nl; ++i1)
         e[it][i1] = 0f;
      for (int it=0; it<wb; ++it)
        e[it][p1] = Float.MAX_VALUE;    
    }
		// Compute d[0]
    for (int il=0; il<nl; ++il)
			d[0][il] = e[0][il];
    // Compute d[1,nlm1]
		// Loop over sparse grid points
		float me = min(e)*nt; // default min value for d
    for (int ic=1; ic<nc; ++ic) {
      int icm1 = ic-1;
			int cgse = cgs[ic];
			int cgsb = cgs[icm1];
			int cgsd = cgse-cgsb;
			int kmax = (int)ceil( rmax*cgsd);
			int kmin = (int)floor(rmin*cgsd);
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
        		if (dv>dm) {
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
    float[][] d, int[][] m, float[] u) 
  {
    int mi = 0;
    int nl = d[0].length;
    int ni = d.length;
    int nlm1 = nl-1;
    int il = 0;
    int ii = ni-1;
    float dl = d[ii][il];
    for (int jl=1; jl<nl; ++jl) {
      if (d[ii][jl]>dl) {
        dl = d[ii][jl];
        il = jl;
      }
    }
    u[ii] = il;
    while (ii>1) { 
      mi = m[ii][il];
      il += mi;
      --ii;
      u[ii] = il;
      if (mi!=0) {
        float du = u[ii]-u[ii+1];
        u[ii] = u[ii+1]+du;
      }
    }
    u[0] = u[1];
  }

	/**
	 * Accumulation of errors on a variably coarse grid.
	 * Determined by the array cgs. The max and minimum 
	 * accumulation slopes correspond to rmin and rmax.
	 * Accumulate maximum values.
	 */
  private static void accumulate(
    double rmin1, double rmax1, double rmin2, double rmax2, 
		int[] cgs, float[][][] e, float[][][] d, int[][][][] m, 
    int p1, int p2, int wb) 
  {
    int nl1 = e[0][0].length;
    int nl2 = e[0].length;
    int nt = e.length;
		int nc = cgs.length;
    int nl1m1 = nl1-1;
    int nl2m1 = nl2-1;
		float mx = max(e)*nt;
    // Fix the water bottom
    if (p1>=0 || p2>=0) {
      for (int it=0; it<wb; ++it)
        for (int i2=0; i2<nl2; ++i2)
          for (int i1=0; i1<nl1; ++i1)
            e[it][i2][i1] = 0f;
    }
    if (p1>=0 && p2>=0) {
      for (int it=0; it<wb; ++it)
        e[it][p2][p1] = Float.MAX_VALUE;    
    }
    else if (p1<0 && p2>=0) {
      for (int it=0; it<wb; ++it)
        for (int ip=0; ip<nl1; ++ip)
          e[it][p2][ip] = Float.MAX_VALUE;    
    }
    else if (p1>=0 && p2<0) {
      for (int it=0; it<wb; ++it)
        for (int ip=0; ip<nl2; ++ip)
          e[it][ip][p1] = Float.MAX_VALUE;    
    }
		// Compute d[0]
    for (int il2=0; il2<nl2; ++il2) 
    	for (int il1=0; il1<nl1; ++il1) 
				d[0][il2][il1] = e[0][il2][il1];
    // Compute d[1,nlm1]
		// Loop over sparse grid points
		float me = min(e)*nt; // default min value for d
    for (int ic=1; ic<nc; ++ic) {
      int icm1 = ic-1;
			int cgse = cgs[ic];
			int cgsb = cgs[icm1];
			int cgsd = cgse-cgsb;
			int kmax1 = (int)ceil( rmax1*cgsd);
			int kmin1 = (int)floor(rmin1*cgsd);
			int kmax2 = (int)ceil( rmax2*cgsd);
			int kmin2 = (int)floor(rmin2*cgsd);
			float[][] dm = new float[nl2][nl1]; 
			// Loop over all slopes
      for (int k2=kmin2; k2<=kmax2; ++k2) {
        int il2s = max(0,-k2);      // lag 2 start
        int il2e = min(nl2,nl2-k2); // lag 2 end
        float r2 = (float)k2/(float)cgsd; // slope
     	 	for (int k1=kmin1; k1<=kmax1; ++k1) {
				  int il1s = max(0,-k1);      // lag 1 start
        	int il1e = min(nl1,nl1-k1); // lag 1 end
        	float r1 = (float)k1/(float)cgsd; // slope
					// Loop over all lags
      		for (int il2=il2s; il2<il2e; ++il2) {
         	 	int il2k2 = il2+k2;
      			for (int il1=il1s; il1<il1e; ++il1) {
            	dm[il2][il1] = d[icm1][il2k2][il1+k1] + e[ic][il2][il1];
          	}
					}
					if (r1==0 && r2==0) { // zero slopes, no interpolation necessary
            for (int x=cgsb+1; x<cgse; ++x) 
              for (int il2=il2s; il2<il2e; il2++) 
                for (int il1=il1s; il1<il1e; il1++) 
                  dm[il2][il1] += e[x][il2][il1];
          } else if (r1==0 && r2!=0) {
            linearInterp(cgse,cgsb,1,r2,il2s,il2e,il1s,il1e,dm,e,true);
          } else if (r1!=0 && r2==0) {
            linearInterp(cgse,cgsb,1,r1,il2s,il2e,il1s,il1e,dm,e,false);
          } else {
            bilinearInterp(cgse,cgsb,1,r2,r1,il2s,il2e,il1s,il1e,dm,e);
          }
          // update previous errors and record moves.
          for (int il2=il2s; il2<il2e; il2++) {
            for (int il1=il1s; il1<il1e; il1++) {
              if (dm[il2][il1]>d[ic][il2][il1]) {
                d[ic][il2][il1] = dm[il2][il1];
                m[icm1][il2][il1][0] = k1;
                m[icm1][il2][il1][1] = k2;
              }
            }
          }
				} // end k1
			} // end k2
    } // end coarse grid
  }
	/**
	 * @author Stefan Compton, CSM
	 */
	private static void linearInterp(
      int ie, int je, int is, float r,
      int il1s, int il1e, int ilSs, int ilSe,
      float[][] dm, float[][][] e, boolean interp1)
  {
    for (int x=je+is; x!=ie; x+=is) {
      float ky = r*(ie-x);
      int k1 = (int)ky;
      if (ky<0.0f) --k1;
      int k2 = k1+1;
      float w1 = k2-ky;
      float w2 = 1.0f-w1;
      if (interp1) {
        for (int il1=il1s; il1<il1e; il1++) {
          int l1 = il1+k1;
          int l2 = il1+k2;
          for (int ilS=ilSs; ilS<ilSe; ilS++)
            dm[il1][ilS] += w1*e[x][l1][ilS]+w2*e[x][l2][ilS];
        }
      } else {
        for (int il1=il1s; il1<il1e; il1++)
          for (int ilS=ilSs; ilS<ilSe; ilS++)
            dm[il1][ilS] += w1*e[x][il1][ilS+k1]+w2*e[x][il1][ilS+k2];
      }
    }
  }
	/**
	 * @author Stefan Compton, CSM
	 */
  private static void bilinearInterp(
      int ie, int je, int is, float r1, float rS,
      int il1s, int il1e, int ilSs, int ilSe,
      float[][] dm, float[][][] e)
  {
    for (int x=je+is; x!=ie; x+=is) {
      float k1y = r1*(ie-x);
      float kSy = rS*(ie-x);
      int k1i1 = (int)k1y;
      int kSi1 = (int)kSy;
      if (k1y<0.0f) --k1i1;
      if (kSy<0.0f) --kSi1;
      int k1i2 = k1i1+1;
      int kSi2 = kSi1+1;
      float d1 = k1i2-k1y;
      float dS = kSi2-kSy;
      float d1S = d1*dS;
      for (int il1=il1s; il1<il1e; il1++) {
        int l1 = il1+k1i1;
        int l2 = il1+k1i2;
        for (int ilS=ilSs; ilS<ilSe; ilS++) {
          float e11 = e[x][l1][kSi1+ilS];
          float e21 = e[x][l1][kSi2+ilS];
          float e12 = e[x][l2][kSi1+ilS];
          float e22 = e[x][l2][kSi2+ilS];
          dm[il1][ilS] += 
            e11 + (e21-e11)*dS + (e12-e11)*d1 + (e11-e21-e12+e22)*d1S;
        }
      }
    }
  }
 private static void backtrack(
    float[][] d, int[][][][] m, float[] u1, float[] u2)
  {
    int n1  = m.length;
    int nl2 = m[0].length;
    int nl1 = m[0][0].length;
    int i1 = n1-1;
    float dm = d[0][0];
    int il1 = 0;
    int il2 = 0;
    // Search for minimum lag indices.
    for (int jl2=0; jl2<nl2; jl2++) {
      for (int jl1=0; jl1<nl1; jl1++) {
        if (d[jl2][jl1]>dm) {
          dm = d[jl2][jl1];
          il1 = jl1;
          il2 = jl2;
        }
      }
    }
    u1[i1] = il1;
    u2[i1] = il2;
    while (i1>0) {
      int mi1 = m[i1][il2][il1][0];
      int mi2 = m[i1][il2][il1][1];
      il1 += mi1;
      il2 += mi2;
      --i1;
      u1[i1] = il1;
      u2[i1] = il2;
      if (mi1!=0) {
        float du1 = u1[i1]-u1[i1+1];
        u1[i1] = u1[i1+1]+du1;
      }
      if (mi2!=0) {
        float du2 = u2[i1]-u2[i1+1];
        u2[i1] = u2[i1+1]+du2;
      }
    }
    u1[0] = u1[1];
    u2[0] = u2[1];
  }
	
/////////////////////////////////////////////////////////
// static utilities

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
  
	private static float[] floats(int[] x) {
    int n1 = x.length;
    float[] y = new float[n1];
    for (int i1=0; i1<n1; ++i1) {
      y[i1] = (float)x[i1];
    }
    return y;
  }

	private static float[][] floats(int[][] x) {
    int n2 = x.length;
    int n1 = x[0].length;
    float[][] y = new float[n2][n1];
    for (int i2=0; i2<n2; ++i2) {
    	for (int i1=0; i1<n1; ++i1) 
      	y[i2][i1] = (float)x[i2][i1];
    }
    return y;
  }


};


