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

public class DynamicSolver {

  /**
   * Constructs a dynamic solver for specified bounds on shifts.
   * Shifts vary on the coarse interval given by 1/ds
   * @param nl number of lags 
   * @param rmin lower bound on shift derivative u'
   * @param rmax upper bound on shift derivative u'
	 * @param scgs coarse grid sampling
   */
  public DynamicSolver(int nl, double rmin, double rmax, double dr, int wb) {
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
  public DynamicSolver(int nl, double rmin, double rmax, double dr) {
		this(nl,rmin,rmax,1.0,0);
  }
  public DynamicSolver(int nl, double rmin, double rmax) {
		this(nl,rmin,rmax,1.0,0);
	}
  public DynamicSolver(int nl) {
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
  public DynamicSolver(
		int nl1, int nl2, double rmin1, double rmax1, 
		double rmin2, double rmax2, double dr, int wb) 
	{
    Check.argument(rmax1>=0,"rmax1>=0");
    Check.argument(rmin1<=0,"rmin1<=0");
    Check.argument(rmax2>=0,"rmax2>=0");
    Check.argument(rmin2<=0,"rmin2<=0");
    Check.argument(dr<=1.0,"dr<=1.0");
		//if (dr> min(rmax1,-rmin1)) dr = min(rmax1,-rmin1);
		//if (dr> min(rmax2,-rmin2)) dr = min(rmax2,-rmin2);
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
  public DynamicSolver(
		int nl1, int nl2, double rmin1, double rmax1, 
    double rmin2, double rmax2, double dr)
	{
		this(nl1,nl2,rmin1,rmax1,rmin2,rmax2,dr,0);
	}
  public DynamicSolver(int nl1, int nl2) {
		this(nl1,nl2,-1.0,1.0,-1.0,1.0,1.0,0);
	}

	public void setStrainLimits1(double rmin, double rmax) {
    Check.argument(rmax>=0,"rmax>=0");
    Check.argument(rmin<=0,"rmin<=0");
		_rmin1 = rmin;
		_rmax1 = rmax;
	}

	public void setStrainLimits2(double rmin, double rmax) {
    Check.argument(rmax>=0,"rmax>=0");
    Check.argument(rmin<=0,"rmin<=0");
		_rmin2 = rmin;
		_rmax2 = rmax;
	}

	public void setTimeVaryingStrainLimits(
    double[] r1min, double[] r1max, 
    double[] r2min, double[] r2max) 
  {
		_rmin1a = r1min;
		_rmax1a = r1max;
		_rmin2a = r2min;
		_rmax2a = r2max;
	}

	public void setSmoothing(double dr) {
    Check.argument(dr<=1.0,"dr<=1.0");
		_scgs = (int)ceil(1/dr);
		_dr = dr;
	}

  public void setWbParam(int p1) {
    setWbParam(p1,_p2);
  }
  public void setWbParam(int p1, int p2) {
    _p1 = p1;
    _p2 = p2;
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
		backtrack(d,m,uc,_p1);
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
    _cgs = cgs;
		int nc = cgs.length;
		float[][] d = new float[nc][nl];
		int[][] m = new int[nc][nl];
    float[] uc = new float[nc];
		accumulate(_rmin1,_rmax1,cgs,e,d,m,_p1,_wb);
		backtrack(d,m,uc,_p1);
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
    int nt = e.length;
		int dc = _scgs;
		int[] cgs = subsample(nt,dc,(float)_wb);
    return findSolution(e,cgs);
  }
  public float[][] findSolution(float[] f, float[][][] e) {
    int dc = _scgs;
    int[] cgs = subsample(f,dc,(float)_wb);
    return findSolution(e,cgs);
  }
  public float[][] findSolution(float[] f, float[][][] e, int ng) {
    int dc = _scgs;
    int[] cgs = subsample(f,dc,ng,(float)_wb);
    return findSolution(e,cgs);
  }
  public float[][] findSolution(float[][][] e, int[] cgs) {
    int nl1 = e[0][0].length;
    int nl2 = e[0].length;
    int nt = e.length;
    //SimplePlot.asPoints(floats(cgs));
    _cgs = cgs;
    int nc = cgs.length;
    float[][][] d = new float[nc][nl2][nl1];
    int[][][][] m = new int[nc][nl2][nl1][2];
    float[] uc1 = new float[nc];
    float[] uc2 = new float[nc];
    if (_rmin1a!=null && _rmax1a!=null && _rmin2a!=null && _rmax2a!=null)
      accumulate(_rmin1a,_rmax1a,_rmin2a,_rmax2a,cgs,e,d,m,_p1,_p2,_wb);
    else
      accumulate(_rmin1,_rmax1,_rmin2,_rmax2,cgs,e,d,m,_p1,_p2,_wb);
    backtrack(d,m,uc1,uc2,_p1,_p2);
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
    gl.add(0); gl.add(nm); 
    if (wb>0) gl.add(wb);
    for (int ii=0; ii<n; ii++) {
      int im = i[nm-ii]; // next highest amplitude index
      if (im>wb || wb==0) {
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
    System.out.println("final ng="+nl);
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
      gl.add(0); gl.add(nm); 
      if (wb>0) gl.add(wb);
      for (int ii=0; ii<n; ii++) {
        int im = i[nm-ii]; // next highest amplitude index
        if (im>wb || wb==0) {
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
  
  public int[] getCgs() {
    return _cgs;
  }
   
  //////////////////////////////////////////////////////////////////////////////
  // private

  private int _nl1,_nl2; // number of lags
	private int _scgs; // coarse grid sampling
  private double _rmin1,_rmin2; // upper int strain limit of u'
  private double _rmax1,_rmax2; // lower int strain limit of u'
  private double[] _rmin1a,_rmin2a; // upper int strain limit of u'
  private double[] _rmax1a,_rmax2a; // lower int strain limit of u'
	private double _dr; // one over the grid sampling
	private int _wb=0; // wb depth
	private int _p1=-1; // wb parameter 1
	private int _p2=-1; // wb parameter 2
  private int[] _cgs; // variable grid coarse sampling

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
    int dir = 1;
    int nl = e[0].length;
    int nt = e.length;
		int nc = cgs.length;
    int nlm1 = nl-1;
		float mx = max(e); 
    if (p1>=0) {
      for (int it=0; it<wb; ++it)
        for (int i1=0; i1<nl; ++i1)
          e[it][i1] = 0f;
      for (int it=0; it<wb; ++it)
        e[it][p1] = mx;
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
		float me = min(e)*nt; // default min value for d
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
    float[][] d, int[][] m, float[] u, int p1) 
  {
    int dir = 1;
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
    if (dir>0) u[0] = u[1];
  }

	/**
	 * Accumulation of errors on a variably coarse grid.
	 * Determined by the array cgs. The max and minimum 
	 * accumulation slopes correspond to rmin and rmax.
	 * Accumulate maximum values.
	 */
  private static void accumulate(
    double rmin1, double rmax1, double rmin2, double rmax2, 
		int[] cgs, float[][][] ei, float[][][] d, int[][][][] m,
    int p1, int p2, int wb) 
  {
    int nt = ei.length;
    accumulate(filldouble(rmin1,nt),filldouble(rmax1,nt),
               filldouble(rmin2,nt),filldouble(rmax2,nt),cgs,ei,d,m,p1,p2,wb);
  }
  private static void accumulate(
    double[] r1min, double[] r1max, double[] r2min, double[] r2max, 
		int[] cgs, float[][][] e, float[][][] d, int[][][][] m,
    int p1, int p2, int wb)
  {
    int dir = -1;
    //float[][][] e = copy(ei);
    int nl1 = e[0][0].length;
    int nl2 = e[0].length;
    int nt = e.length;
		int nc = cgs.length;
    int nl1m1 = nl1-1;
    int nl2m1 = nl2-1;
		float mx = max(e);
    if (p1>=0 || p2>=0) {
      for (int it=0; it<wb; ++it)
        for (int i2=0; i2<nl2; ++i2)
          for (int i1=0; i1<nl1; ++i1)
            e[it][i2][i1] = -mx;
    }
    if (p1>=0 && p2>=0) {
      for (int it=0; it<wb; ++it)
        e[it][p2][p1] = mx;
    }
    else if (p1<0 && p2>=0) {
      for (int it=0; it<wb; ++it)
        for (int ip=0; ip<nl1; ++ip)
          e[it][p2][ip] = mx;
    }
    else if (p1>=0 && p2<0) {
      for (int it=0; it<wb; ++it)
        for (int ip=0; ip<nl2; ++ip)
          e[it][ip][p1] = mx;
    }
    int ibg = (dir>0)?0:nc-1; // beginning index
    int ieg = (dir>0)?nc:0;  // end index
    int is = (dir>0)?1:-1;   // stride
    int ic = ibg; // sparse grid index
    int ie = cgs[ic]; // error index
		// Compute d[0]
    for (int il2=0; il2<nl2; ++il2) {
    	for (int il1=0; il1<nl1; ++il1) {
				d[ic][il2][il1] = e[ie][il2][il1];
			}
		}
    // Compute d[1,nlm1]
		// Loop over sparse grid points
    ic += is;
    // Loop over all sparse grid points.
    for (; ic!=ieg; ic+=is) {
      int icm1 = ic-is; // previous sparse grid index
			int cgse = cgs[ic];
			int cgsb = cgs[icm1];
			int cgsd = cgse-cgsb;
      int k1min,k1max,k2min,k2max;
      if (cgsd>0) { // forward
        k1min = (int) ceil(-r1max[ie]*cgsd);
        k1max = (int)floor(-r1min[ie]*cgsd);
        k2min = (int) ceil(-r2max[ie]*cgsd);
        k2max = (int)floor(-r2min[ie]*cgsd);
      } else { // reverse
        k1min = (int) ceil(-r1min[ie]*cgsd);
        k1max = (int)floor(-r1max[ie]*cgsd);
        k2min = (int) ceil(-r2min[ie]*cgsd);
        k2max = (int)floor(-r2max[ie]*cgsd);
      }
      k1min = k1min>k1max ? k1max : k1min;
      k2min = k2min>k2max ? k2max : k2min;
			float[][] dm = new float[nl2][nl1]; 
      fill(Float.MIN_VALUE,d[ic]);
			// Loop over all slopes
      for (int k2=k2min; k2<=k2max; ++k2) {
        int il2s = max(0,-k2);      // lag 2 start
        int il2e = min(nl2,nl2-k2); // lag 2 end
        float r2 = (float)k2/(float)cgsd; // slope
     	 	for (int k1=k1min; k1<=k1max; ++k1) {
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
            for (int x=cgsb+is; x!=cgse; x+=is) 
              for (int il2=il2s; il2<il2e; il2++) 
                for (int il1=il1s; il1<il1e; il1++) 
                  dm[il2][il1] += e[x][il2][il1];
          } else if (r1==0 && r2!=0) {
            linearInterp(cgse,cgsb,is,r2,il2s,il2e,il1s,il1e,dm,e,true);
          } else if (r1!=0 && r2==0) {
            linearInterp(cgse,cgsb,is,r1,il2s,il2e,il1s,il1e,dm,e,false);
          } else {
            bilinearInterp(cgse,cgsb,is,r2,r1,il2s,il2e,il1s,il1e,dm,e);
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
    float[][][] d, int[][][][] m, float[] u1, float[] u2, int p1, int p2)
  {
    int dir = -1;
    int n1  = m.length;
    int nl2 = m[0].length;
    int nl1 = m[0][0].length;
    int n1m1 = n1-1;
    int il1 = 0;//(dir>0)?0:nl1-1;
    int il2 = 0;//(dir>0)?0:nl2-1;
    int ie = (dir>0)?1:n1m1;
    int ib = (dir>0)?n1m1:0;
    int is = (dir>0)?1:-1;
    int i1 = ib;
    float dm = d[i1][0][0];
    for (int jl2=0; jl2<nl2; jl2++) {
      for (int jl1=0; jl1<nl1; jl1++) {
        if (d[i1][jl2][jl1]>dm) {
          dm = d[i1][jl2][jl1];
          il1 = jl1;
          il2 = jl2;
        }
      }
    }
    //if (dir<0 && p1>=0) il1 = p1;
    //if (dir<0 && p2>=0) il2 = p2;
    u1[i1] = il1;
    u2[i1] = il2;
    //if (dir<0) {
    //  i1-=is;
    //  u1[i1] = il1;
    //  u2[i1] = il2;
    //}
    while (i1!=ie) {
      i1-=is;
      int mi1 = m[i1][il2][il1][0];
      int mi2 = m[i1][il2][il1][1];
      il1 += mi1;
      il2 += mi2;
      u1[i1] = il1;
      u2[i1] = il2;
      if (mi1!=0) {
        float du1 = u1[i1]-u1[i1+is];
        u1[i1] = u1[i1+is]+du1;
      }
      if (mi2!=0) {
        float du2 = u2[i1]-u2[i1+is];
        u2[i1] = u2[i1+is]+du2;
      }
    }
    //u1[0] = u1[1];
    //u2[0] = u2[1];
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


  ///////////////////////////////////////////////////////////////////////////////////////////
  // For Testing only:
/*
  public void interpTest( 
    float[] u, float[] g) 
  {
    int n1 = g.length;
    float[] r = add(rampfloat(0f,1f,n1),u);
    monotonicityTest(r,u);
    int n2 = (int)(ceil(r[n1-1])-floor(r[0]));
    float[] s = new float[n2];
    float[] h = new float[n2];
    inverseLinearInterpolation(n1,1f,0f,r,n2,1f,u[0],s,r[0],r[n1-1]);
    _si.setUniform(n1,1.0,0.0,g);
    for (int i2=0; i2<n2; ++i2) 
      h[i2] = _si.interpolate(s[i2]);
     
    int n3 = (int)(ceil(s[n2-1])-floor(s[0]));
    float[] s2 = new float[n3];
    float[] h2 = new float[n3];
    inverseLinearInterpolation(n2,1f,u[0],s,n3,1f,0f,s2,s[0],s[n2-1]);
    _si.setUniform(n2,1.0,u[0],h);
    for (int i3=0; i3<n3; ++i3) 
      h2[i3] = _si.interpolate(s2[i3]);
    System.out.println("n1="+n1);
    System.out.println("n2="+n2);
    System.out.println("n3="+n3);
    System.out.println("s20="+s2[0]);
    System.out.println("s2[n2m]="+s2[n3-1]);
    SimplePlot sp = new SimplePlot();
    sp.addPoints(g);
    PointsView l1 = sp.addPoints(h2);
    l1.setLineColor(Color.RED);
  }

  /**
   * Returns errors accumulated alignment errors.
   * @param e array of alignment errors.
   * @return array of accumulated errors.
   *//*
  public float[][] accumulate(float[][] e) {
    float[][] d = like(e);
    int[][] m = new int[e.length][e[0].length];
    accumulate(e,d,m);
    return d;
  }

  /**
   * Accumulates alignment errors.
   * @param e input array of alignment errors.
   * @param d output array of accumulated errors.
   *//*
  public void accumulate(float[][] e, float[][] d, int[][] m) {
    accumulate(_kmin,_kmax,_jmin,_jmax,e,d,m);
  }

  /**
   * Computes shifts for specified sequences.
   * @param e array of alignment errors.
   * @param d array of accumulated alignment errors.
   * @param u output array of shifts u.
   *//*
  public float[] findShifts(float[][] d) {
    float[] u = new float[d.length];
    backtrack(d,u);
    return u;
  }

  /**
   * Computes shifts by backtracking in reverse direction.
   * @param d input array of accumulated errors.
   * @param e input array of alignment errors.
   * @param u output array of shifts.
   *//*
  public void backtrack(float[][] d, float[] u) {
    backtrack(_shifts,d,m,u);
  }
  /**
   * Computes shifts for specified sequences efficiently.
   * Saves memory and speed. 
   * @param f array[n1] for the sequence f[i1].
   * @param g array[n2] for the sequence g[i2].
   * @param u output array of shifts u.
   *//*
  public float[] findShiftsFast(float[] f, float[] g) {
    return findShiftsFast(f,g,false);
  }
  public float[] findShiftsFast(float[] f, float[] g, boolean phase) {
    float[] u = like(f);
    float[] dmin = new float[1];
    findShiftsFast(_fr,_kmin,_kmax,_jmin,_jmax,_epow,_shifts,dmin,f,g,u);
    if (phase) return dmin;
    return u;
  }

  //Testing only//////////////////////////////////////////////////
  public float[][] finduo(float[] f, float[] g) {
    int fr = _fr;
    int nl = _shifts.getCount();
    int[] uoi = rampint(0,1,nl);
    float[] ds = new float[nl];
    finduo(_kmin,_kmax,_fr,_shifts,f,g,ds,uoi);
    float[] uoif = div(floats(uoi),fr);
    return new float[][]{uoif,ds};
  }
  public float[] findShiftsX(float[][] e, Sampling sr, float uo) {
    return findShiftsX(e,sr,uo,false);
  }
  public float[] findShiftsX(float[][] e, Sampling sr, float uo, boolean mind) {
    int n1 = e.length;
    int nr = sr.getCount();
    // accumulate
    float[][] d = new float[2][nr];
    int[][] m = new int[n1][nr];
    float[][] s = new float[n1][nr];
    accumulateX(uo,_shifts,sr,e,s,d,m);
    if (mind) return new float[]{min(d[0])};
    // backtrack
    float[] dmin = new float[1];
    float[] u = new float[n1];
    backtrackX(d[0],m,s,u);
    return u;
  }
  ////////////////////////////////////////////////////////////////


 /**
   * Non-linear accumulation of fractional alignment errors.
   *//*
  private static void accumulate(
    int kmin, int kmax, int jmin, int jmax,
    float[][] e, float[][] d, int[][] m) 
  {
    int nl = e[0].length;
    int ni = e.length;
    int nlm1 = nl-1;
    // Compute d[0]
    for (int il=0; il<nl; ++il)
      d[0][il] = e[0][il];
    // Compute d[1] 
    for (int il=0; il<nl; ++il) {
      float dm = d[0][il];
      for (int k=kmin; k<=kmax; ++k) {
        int ik = il+k; 
	if (0<=ik && ik<nl) {
          float dv = d[0][ik];
          if (dv<dm) {
            dm = dv;
            m[1][il] = k;
          }
	}
      }
      d[1][il] = dm + e[1][il];
    }
    // Compute d[2,nlm1]
    float em = max(e)*ni;
    for (int ii=2; ii<ni; ++ii) {
      int im1 = ii-1;
      for (int il=0; il<nl; ++il) {
        float dm = em;
        for (int k=kmin; k<=kmax; ++k) {
	  int ik = il+k; 
	  if (0<=ik && ik<nl) {
	    float dv = d[im1][ik];
            if (dv<dm) {
	      int ktest = k-m[im1][ik];
	      if (jmin<=ktest && ktest<=jmax) {
	        dm = dv;
	        m[ii][il] = k;
	      }
	    }
	  }
	}
	d[ii][il] = dm + e[ii][il];
      }
    }
  }
*/
  /**
   * Non-linear accumulation of fractional alignment errors efficiently in memory.
   *//*
  private static void findShiftsFast(
    int fr, int kmin, int kmax, int jmin, int jmax, float p, 
    Sampling shifts, float[] dmin, float[] f, float[] g, float[] u) 
  {
    int ni = f.length;
    int n2 = g.length;
    int n2i = (n2-1)*fr+1;
    int nl = 1+(n2-ni)*fr;
    int nlm1 = nl-1;
    float[] gi = new float[n2i];
    float[] d0 = new float[nl];
    float[] d1 = new float[nl];
    int[][] m = new int[ni][nl];
    SincInterpolator si2 = new SincInterpolator();
    si2.setUniform(n2,1f,0f,g);
    si2.interpolate(n2i,shifts.getDelta(),0f,gi);
    // Initalize d
    for (int il=0; il<nl; ++il)
      d0[il] = errors(f[0],gi[il],p);
    // Compute m[1] 
    for (int il=0,i2=fr+il; il<nl; ++il,++i2) {
      float dm = d0[il];
      for (int k=kmin; k<=kmax; ++k) {
        int ik = il+k; 
	if (0<=ik && ik<nl) {
          float dv = d0[ik];
          if (dv<dm) {
            dm = dv;
            m[1][il] = k;
          }
	}
      }
      float es = errors(f[1],gi[i2],p);
      d1[il] = dm + es;
    }
    // Compute m[2,nlm1]
    for (int ii=2; ii<ni; ++ii) {
      int im1 = ii-1;
      float em = max(d1)*ni;
      float[] d2 = d1;
      d1 = d0;
      d0 = d2;
      for (int il=0,i2=fr*ii+il; il<nl; ++il,++i2) {
        float dm = em; 
        for (int k=kmin; k<=kmax; ++k) {
	  int ik = il+k; 
	  if (0<=ik && ik<nl) {
	    float dv = d2[ik];
            if (dv<dm) {
	      int ktest = k-m[im1][ik];
	      if (jmin<=ktest && ktest<=jmax) {
	        dm = dv;
	        m[ii][il] = k;
	      }
	    }
	  }
	}
        float es = errors(f[ii],gi[i2],p);
	d1[il] = dm + es;
      }
    }
    // Backtrack
    int il = 0;
    int ii = ni-1;
    float dl = d1[il];
    for (int jl=1; jl<nl; ++jl) {
      if (d1[jl]<dl) {
        dl = d1[jl];
        il = jl;
      }
    }
    dmin[0] = dl;
    u[ii] = (float)shifts.getValue(il);
    while (ii>0) { 
      int mi = m[ii][il];
      il += mi;
      --ii;
      u[ii] = (float)shifts.getValue(il);
      if (mi!=0) {
        float du = u[ii]-u[ii+1];
        u[ii] = u[ii+1]+du;
      }
    }
  }

  private static void accumulateX(
    float uo, Sampling ss, Sampling sr, 
    float[][] e, float[][] s, float[][] d, int[][] m) 
  {
    int n1 = s.length;
    int nr = s[0].length;
    int nl = e[0].length;
    int n1m = n1-1;
    int nrm = nr-1;
    double ds = ss.getDelta();
    int lmax = (int)(nl*ds);
    SincInterpolator si = new SincInterpolator();
    si.setUniform(nl,ds,0f,e[0]);
    for (int ir=0; ir<nr; ++ir) {
      s[0][ir] = uo;
      d[0][ir] = si.interpolate(uo);
    }
    for (int i1=1; i1<n1; ++i1) {
      int im1 = i1-1;
      si.setUniform(nl,ds,0f,e[i1]);
      float[] d0 = d[0]; d[0] = d[1]; d[1] = d0;
      for (int ir=0; ir<nr; ++ir) {
	int mm = 0;
        int irp = ir+1; if (irp>nrm) irp=nrm;
        int irm = ir-1; if (irm<0) irm=0;
	float r = (float)sr.getValue(ir); 
        float dl = d[1][ir ]+abs(si.interpolate(s[im1][ir ]+r));
        float dp = d[1][irp]+abs(si.interpolate(s[im1][irp]+r));
        float dm = d[1][irm]+abs(si.interpolate(s[im1][irm]+r));
	if (dp<dl) {
	  dl = dp;
	  mm = 1;
	}
	if (dm<dl) {
	  dl = dm;
	  mm = -1;
	}
	float sl = s[im1][ir+mm] + r;
        s[i1][ir] = (sl<lmax)?((sl>=0)?sl:0):lmax;
	d[0][ir] = d[1][ir+mm] + abs(si.interpolate(s[i1][ir]));
	m[i1][ir] = mm;
      }
    }
  }


  private static void backtrackX(
    float[] d, int[][] m, float[][] s, float[] u) 
  {
    int mi = 0;
    int nl = d.length;
    int ni = u.length;
    int nlm1 = nl-1;
    int ii = ni-1;
    int il = 0;
    float dl = d[il];
    for (int jl=1; jl<nl; ++jl) {
      if (d[jl]<dl) {
        dl = d[jl];
        il = jl;
      }
    }
    u[ii] = s[ii][il];
    while (ii>0) {
      mi = m[ii][il];
      il += mi;
      --ii;
      u[ii] = s[ii][il];
    }
  }
*/
  /**
   * Finds uo by backwards accumulation and sorts them.
   * int[] uoi = rampint(0,1,nl);
   */
  /*private static void finduo(
    int kmin, int kmax, int fr, Sampling shifts,
    float[] f, float[] g, float[] ds, int[] uoi)
  {
    int n1 = f.length;
    int n2 = g.length;
    int nl = uoi.length;
    int n1m = n1-1;
    int nlm = nl-1;
    int n1m2 = n1-2;
    int n2i = (n2-1)*fr+1;
    float[] d0 = new float[nl];
    float[] d1 = new float[nl];
    float[] gi = new float[n2i];
    SincInterpolator si2 = new SincInterpolator();
    si2.setUniform(n2,1f,0f,g);
    si2.interpolate(n2i,shifts.getDelta(),0f,gi);
    for (int il=0,i2= n1m+il; il<nl; ++il,++i2)
      d1[il] = errors(f[n1m],gi[i2],2);
    for (int i1=n1m2; i1>=0; --i1) {
      int im1 = i1-1;
      float[] d2 = d1; d1 = d0; d0 = d2;
      for (int il=0,i2=il+i1; il<nl; ++il,++i2) {
        float e = errors(f[i1],gi[i2],2);
	float dm = d2[il];
        for (int k=kmin; k<=kmax; ++k) {
          int ik = il+k; 
          if (0<=ik && ik<nl) {
            float dv = d2[ik];
            if (dv<dm) {
              dm = dv;
	    }
	  }
	}
	d1[il] = dm + e;
      }
    }
    // Find uo
    copy(d1,ds);
    quickIndexSort(ds,uoi);
    quickSort(ds);
  }


*/

};

