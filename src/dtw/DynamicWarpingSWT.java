package dtw;
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

public class DynamicWarpingSWT {

  /**
   * Constructs a dynamic warping for specified bounds on shifts.
   * Shifts vary on the coarse interval given by 1/ds
   * @param nl number of lags where nl = (n2-n1)+1.
   * @param smin lower bound on shift derivative u'
   * @param smax upper bound on shift derivative u'
	 * @param scgs coarse grid sampling
   */
  public DynamicWarpingSWT(int nl, double rmin, double rmax, double dr) {
    Check.argument(rmax>=0,"rmax>=0");
    Check.argument(rmin<=0,"rmin<=0");
    Check.argument(dr<=1.0,"dr<=1.0");
		int nr = (int)ceil((rmax-rmin)/dr);
    _nl = nl;
		_scgs = (int)ceil(1/dr);
		_rmin = rmin;
		_rmax = rmax;
		_dr = dr;
  }
  public DynamicWarpingSWT(int nl, double rmin, double rmax) {
	this(nl,rmin,rmax,1.0);
	}
  public DynamicWarpingSWT(int nl) {
	this(nl,-1.0,1.0,1.0);
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
    for (int i1=0; i1<n1; ++i1) {
      for (int il=0,i2=i1+il; il<nl; ++il,++i2)
        e[i1][il] = error(f[i1],g[i2]);
    }
    normalizeErrors(e);
    return e;
  }

	/**
   * Returns normalized alignment errors for all samples and fractional lags.
   * @param f array[n1] for the sequence f[i1].
   * @param g array[nt][n2] for the sequences g[i2].
   * @return array[n1][nl] of alignment errors.
   */
  public float[][] computeErrorsMulti(float[] f, float[][] g, int tmin, int tmax) {
		//Check.argument(smin<=0,"smin<=0");
		//Check.argument(smax>=0,"smax>=0");
    int n1 = f.length;
    int n2 = g[0].length;
		int nt = g.length;
		//int nl = smax-smin+1;
		int ngn = tmax-tmin+1;
		ngn = (ngn<n2)?ngn:n2;
		int nl = ngn-n1+1;
		nl = (nl>1)?nl:1;
		if (nl==1) System.out.println("nl = 1");
    //int nl = _nl;
    float[][] e = new float[n1][nl];
		for (int it=0; it<nt; ++it) {
    	for (int i1=0; i1<n1; ++i1) {
      	for (int il=0,i2=tmin+i1+il; il<nl; ++il,++i2) 
        	e[i1][il] += error(f[i1],g[it][i2]);
			  //int illo = min(nl-1,max(0,-smin-i1)); 
      	//int ilhi = max(0,min(nl,n1-smin-i1)); 
      	//for (int il=illo,i2=smin+i1+il; il<ilhi; ++il,++i2) 
        //	e[i1][il] += error(f[i1],g[it][i2]);
			}
    }
    normalizeErrors(e);
    return e;
  }
	public float[][] computeErrorsMulti(float[] f, float[][] g) {
		return computeErrorsMulti(f,g,0,g[0].length);
	}
	// Compute errors with a top/horizon constraint
  public float[][] computeErrorsMulti(float[] f, float[][] g, int[] wp) {
		int hzi = wp[0];
		int tpi = wp[1];
		int ncl = wp[2];
    int n1 = f.length;
    int n2 = g[0].length;
		int nt = g.length;
    int nl = _nl;
    float[][] e = new float[n1][nl];
		for (int it=0; it<nt; ++it) {
    	for (int i1=0; i1<n1; ++i1) {
      	for (int il=0,i2=i1+il; il<nl; ++il,++i2) 
        	e[i1][il] += error(f[i1],g[it][i2]);
			}
    }
    normalizeErrors(e);
		if (tpi>0) {
			int hil = hzi-tpi;
			float[] tpe = new float[ncl];
			int ib = -(int)floor(ncl/2);
			int ie =  (int)ceil(ncl/2); 
			for (int ij=ib,ii=0;ij<ie;++ii,++ij) 
				tpe[ii] = e[tpi][hil+ij];
			for (int il=0,i2=tpi+il; il<nl; ++il,++i2) 
				e[tpi][il] = FLT_MAX;
			for (int ij=ib,ii=0; ij<ie;++ii,++ij)
				e[tpi][hil+ij] = tpe[ii];
		}
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
		int[] cgs = subsample(nt,dc);
		int nc = cgs.length;
		float[][] d = new float[nc][nl];
		int[][] m = new int[nc][nl];
    float[] uc = new float[nc];
		accumulateR(_rmin,_rmax,cgs,e,d,m);
		backtrack(d,m,uc);
		float[] u = interpolateSparseShifts(nt,cgs,uc);
		return u;
  }
  public float[] findShiftsR(float[] f, float[][] e) {
		int nl = e[0].length;
		int nt = e.length;
		int dc = _scgs;
		HilbertTransformFilter ht = new HilbertTransformFilter();
		float[] ft = new float[nt];
		ht.apply(nt,f,ft);
		float[] h = sqrt(add(pow(f,2.0f),pow(ft,2.0f)));
		int[] cgs = subsample(h,dc,10);
		int nc = cgs.length;
		float[][] d = new float[nc][nl];
		int[][] m = new int[nc][nl];
    float[] uc = new float[nc];
		accumulateR(_rmin,_rmax,cgs,e,d,m);
		backtrack(d,m,uc);
		float[] u = interpolateSparseShifts(nt,cgs,uc);
		return u;
  }


	public float[][] accumulate(float[][] e) {
		int nl = e[0].length;
		int nt = e.length;
		int dc = _scgs;
		int[] cgs = subsample(nt,dc);
		int nc = cgs.length;
		float[][] d = new float[nc][nl];
		int[][] m = new int[nc][nl];
    float[] uc = new float[nc];
		accumulateR(_rmin,_rmax,cgs,e,d,m);
		return d;
	}

  public double findDmin(float[] f, float[][] g, int tmin, int tmax) {
		float[][] e = computeErrorsMulti(f,g,tmin,tmax);
		int nl = e[0].length;
		int nt = e.length;
		int dc = _scgs;
		int[] cgs = subsample(nt,dc);
		int nc = cgs.length;
		float[][] d = new float[nc][nl];
		int[][] m = new int[nc][nl];
		accumulateR(_rmin,_rmax,cgs,e,d,m);
		double dmin = min(d[nc-1]);
		return dmin;
  }
	public double findDmin(float[] f, float[][] g, int[] wp) {
		float[][] e = computeErrorsMulti(f,g,wp);
		int nl = e[0].length;
		int nt = e.length;
		int dc = _scgs;
		int[] cgs = subsample(nt,dc);
		int nc = cgs.length;
		float[][] d = new float[nc][nl];
		int[][] m = new int[nc][nl];
		accumulateR(_rmin,_rmax,cgs,e,d,m);
		double dmin = min(d[nc-1]);
		return dmin;
  }

  /**
   * Computes a sequence warped by applying specified shifts.
   * @param u input array of shifts.
   * @param g input array for the sequence to be warped.
   * @param h output array for the warped sequence.
   */
  public float[] applyShifts(float[] u, float[] g) {
    int n1 = g.length;
    float[] r = add(rampfloat(0f,1f,n1),u);
    //monotonicityTest(r,u);
    int n2 = (int)(ceil(r[n1-1])-floor(r[0]));
    float[] s = new float[n2];
    float[] h = new float[n2];
		SincInterp si = new SincInterp();
    inverseLinearInterpolation(n1,1f,0f,r,n2,1f,u[0],s,r[0],r[n1-1]);
    for (int i2=0; i2<n2; ++i2) 
      h[i2] = si.interpolate(n1,1.0,0.0,g,s[i2]);
    return h;
  }
  public float[][] applyShifts(float[] u, float[] g, boolean shft) {
    int n1 = g.length;
    float[] r = add(rampfloat(0f,1f,n1),u);
    //monotonicityTest(r,u);
    int n2 = (int)(ceil(r[n1-1])-floor(r[0]));
    float[] s = new float[n2];
    float[] h = new float[n2];
		SincInterp si = new SincInterp();
    inverseLinearInterpolation(n1,1f,0f,r,n2,1f,u[0],s,r[0],r[n1-1]);
    for (int i2=0; i2<n2; ++i2) 
      h[i2] = si.interpolate(n1,1.0,0.0,g,s[i2]);
    return new float[][]{h,s};
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
    float[] gf = new float[ng];
    float[] ui = new float[n ];
    for (int ig=0; ig<ng; ig++)
      gf[ig] = (float)g[ig];
    CubicInterpolator ci = new CubicInterpolator(
			CubicInterpolator.Method.SPLINE,ng,gf,u);
			//CubicInterpolator.Method.LINEAR,ng,gf,u);
    for (int i=0; i<n; i++)
      ui[i] = ci.interpolate(i);
    return ui;
  }

  public static int[] subsample(float[] f, int d) {
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
    for (int ii=0; ii<n; ii++) {
      int im = i[nm-ii]; // next highest amplitude index
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
    Collections.sort(gl);
    int nl = gl.size();
    int[] g = new int[nl];
    for (int il=0; il<nl; il++) {
      g[il] = gl.get(il);
    }
    return g;
  }

 public static int[] subsample(float[] f, int d, int ng) {
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
      for (int ii=0; ii<n; ii++) {
        int im = i[nm-ii]; // next highest amplitude index
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

  private int _nl; // number of lags
  private float _epow = 2; // exponent used for alignment errors |f-g|^e
	private int _scgs; // coarse grid sampling
  private double _rmin; // upper int strain limit of u'
  private double _rmax; // lower int strain limit of u'
	private double _dr; // one over the grid sampling

  private float error(float f, float g) {
    return pow(abs(f-g),_epow);
  }

   /**
	 * Accumulation of errors on a variably coarse grid.
	 * Determined by the array cgs. The max and minimum 
	 * accumulation slopes correspond to rmin and rmax.
	 */
  private static void accumulateR(
    double rmin, double rmax, int[] cgs,
    float[][] e, float[][] d, int[][] m) 
  {
    int nl = e[0].length;
    int nt = e.length;
		int nc = cgs.length;
    int nlm1 = nl-1;
		// Compute d[0]
    for (int il=0; il<nl; ++il)
			d[0][il] = e[0][il];
    // Compute d[1,nlm1]
		// Loop over sparse grid points
		float me = max(e)*nt; // default max value for d
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
   * Finds shifts by backtracking in accumulated alignment errors. 
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
      if (d[ii][jl]<dl) {
        dl = d[ii][jl];
        il = jl;
      }
    }
    u[ii] = il;
    while (ii>0) { 
      mi = m[ii][il];
      il += mi;
      --ii;
      u[ii] = il;
      if (mi!=0) {
        float du = u[ii]-u[ii+1];
        u[ii] = u[ii+1]+du;
      }
    }
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
  private void inverseLinearInterpolation(
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
  private void monotonicityTest(float[] x, float[] u) {
    // Prints location of non-monotonic points
    int n = x.length;
    int im1;
    float b = 1.0f;
    for (int i=1; i<n; ++i) {
      im1 = i-1;
      if (x[i]<x[im1]) 
        System.out.format(
	  "r not montonic between %d and %d with a val of %f and %f \n",
	                          im1,    i,            x[im1], x[i]);
      if (abs(u[i]-u[im1])>b) 
        System.out.format(
	  "u not montonic between %d and %d with a val of %f b is %f \n",
	                          im1,    i,        abs(u[i]-u[im1]),b);
    }
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

