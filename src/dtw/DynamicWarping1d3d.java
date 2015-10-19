package dtw;

import edu.mines.jtk.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;

// for testing only:
import java.awt.Color;
import edu.mines.jtk.awt.ColorMap;

public class DynamicWarping1d3d {

  /**
   * Constructs a dynamic warping for specified bounds on shifts.
   * Shifts vary at the interval given by 1.0/ds.
   * @param nl number of lags where nl1= (ng1-nf)*ods+1.
   * @param ds fractional sampling integer interval
   */
  public DynamicWarping1d3d(
		int nl1, int nl2, int nl3,
	  int ods1, int ods2, int ods3) 
	{
    Check.argument(ods1>0,"ods1>0");
    Check.argument(ods2>0,"ods2>0");
    Check.argument(ods3>0,"ods3>0");
    _nl1 = nl1; _nl2 = nl2; _nl3 = nl3;
    _ods1 = ods1; _ods2 = ods2; _ods3 = ods3;
    _k1min = -ods1; _k2min = -ods2; _k3min = -ods3;
    _k1max =  ods1; _k2max =  ods2; _k3max =  ods3;
    _j1min = -ods1; _j2min = -ods2; _j3min = -ods3;
    _j1max =  ods1; _j2max =  ods2; _j3max =  ods3;
    _shifts1 = new Sampling(_nl1,1.0/_ods1,0.0);
    _shifts2 = new Sampling(_nl2,1.0/_ods2,(_nl2+1)/2.0); // centered
    _shifts3 = new Sampling(_nl3,1.0/_ods3,(_nl3+1)/2.0); // centered
  }
  public DynamicWarping1d3d(int nl1, int nl2, int ods1, int ods2) {
	  this(nl1,nl2,nl2,ods1,ods2,ods2);
	}
  public DynamicWarping1d3d(int nl1, int ods1) {
	  this(nl1,nl1,nl1,ods1,ods1,ods1);
	}

  /** 
	 * Sets the constraints on the vertical time shift lag
   * @param k1min lower bound on shift derivative u'
   * @param k1max upper bound on shift derivative u'
   * @param j1min lower bound on shift second derivative u''
   * @param j1max upper bound on shift second derivative u''
	 */
	public void setl1(int k1min, int k1max, int j1min, int j1max) {
    Check.argument(k1max>=0,"k1max>=0");
    Check.argument(k1min<=0,"k1min<=0");
    Check.argument(j1max>=0,"j1max>=0");
    Check.argument(j1min<=0,"j1min<=0");
    _k1min = k1min;
    _k1max = k1max;
    _j1min = j1min;
    _j1max = j1max;
  }
	public void setl1(int k1min, int k1max) {
	  setl1(k1min,k1max,k1min,k1max);
	}

  /** 
	 * Sets the constraints on the horizontal time shift lag
   * @param k2min lower bound on shift derivative u'
   * @param k2max upper bound on shift derivative u'
   * @param j2min lower bound on shift second derivative u''
   * @param j2max upper bound on shift second derivative u''
	 */
	public void setl2(int k2min, int k2max, int j2min, int j2max) {
    Check.argument(k2max>=0,"k2max>=0");
    Check.argument(k2min<=0,"k2min<=0");
    Check.argument(j2max>=0,"j2max>=0");
    Check.argument(j2min<=0,"j2min<=0");
    _k2min = k2min;
    _k2max = k2max;
    _j2min = j2min;
    _j2max = j2max;
	}
	public void setl2(int k2min, int k2max) {
	  setl2(k2min,k2max,k2min,k2max);
	}

	/** 
	 * Sets the constraints on the horizontal time shift lag
   * @param k3min lower bound on shift derivative u'
   * @param k3max upper bound on shift derivative u'
   * @param j3min lower bound on shift second derivative u''
   * @param j3max upper bound on shift second derivative u''
	 */
	public void setl3(int k3min, int k3max, int j3min, int j3max) {
    Check.argument(k3max>=0,"k3max>=0");
    Check.argument(k3min<=0,"k3min<=0");
    Check.argument(j3max>=0,"j3max>=0");
    Check.argument(j3min<=0,"j3min<=0");
    _k3min = k3min;
    _k3max = k3max;
    _j3min = j3min;
    _j3max = j3max;
	}
	public void setl3(int k3min, int k3max) {
	  setl3(k3min,k3max,k3min,k3max);
	}

  /**
   * Sets the exponent used to compute alignment errors |f-g|^e.
   * The default exponent is 2.
   * @param e the exponent.
   */
  public void setErrorExponent(double e) {
    _epow = (float)e;
  }

//////////////////////////////////////////////////////////////////////////////
// 1D to 2D 

  /**
   * Returns normalized alignment errors for all samples and fractional lags.
   * @param f array[nf] for the sequence f[if1].
   * @param g array[ng2][ng1] for the image g[ig2][ig1].
   * @return array[n1][nl] of alignment errors.
   */
  public float[][][] computeErrors1d2d(float[] f, float[][] g) {
    int nf = f.length;
    int ods1 = _ods1;
    int ods2 = _ods2;
    int nl1 = _nl1;
    int nl2 = _nl2;
    int ng1 = g[0].length;
    int ng2 = g.length;
    int ng1i = (ng1-1)*ods1+1;
    int ng2i = (ng2-1)*ods2+1;
    float[][] gi = new float[ng2i][ng1i];
		double fg2 = _shifts2.getFirst();
		double dg1 = _shifts1.getDelta();
		double dg2 = _shifts2.getDelta();
		int fg2i = (int)(fg2);
    SincInterpolator si = new SincInterpolator();
		for (int i2=0; i2<ng2i; ++i2) {
			for (int i1=0; i1<ng1i; ++i1) 
			  gi[i2][i1] = si.interpolate(ng1,1f,0f,ng2,1f,fg2,g,dg1*i1,fg2+dg2*i2);
		}
    float[][][] e = new float[nf][nl2][nl1];
    for (int ii=0; ii<nf; ++ii) {
      for (int il2=0; il2<nl2; ++il2) {
        for (int il1=0,ig1=ods1*ii+il1; il1<nl1; ++il1,++ig1) 
          e[ii][il2][il1] = error(f[ii],gi[il2][ig1]);
			}
    }
    normalizeErrors(e);
    return e;
  }

	public float[][][] accumulate(float[][][] e) {
    int nf = e.length;
    int nl2 = e[0].length;
    int nl1 = e[0][0].length;
		_m2 = new int[nf][nl2];
		_m1 = new int[nf][nl1];
		float[][][] d = new float[nf][nl2][nl1];
	  accumulate(_k1min,_k1max,_j1min,_j1max,_k2min,_k2max,_j2min,_j2max,e,d,_m1,_m2);
		return d;
	}

  private static void accumulate(
	  int k1min, int k1max, int j1min, int j1max,
		int k2min, int k2max, int j2min, int j2max, 
		float[][][] e, float[][][] d, int[][] m1, int[][] m2)
  {
    int nf = e.length;
    int nl2 = e[0].length;
    int nl1 = e[0][0].length;
		int nl2m = nl2-1;
		int nl1m = nl1-1;
		// Compute d[0]
		for (int il2=0; il2<nl2; ++il2) {
		  for (int il1=0; il1<nl1; ++il1) 
			  d[0][il2][il1] = e[0][il2][il1];
		}
		// Compute d[1]
		for (int il2=0; il2<nl2; ++il2) {
			for (int il1=0; il1<nl1; ++il1) {
		    float dm = d[0][il2][il1];
				for (int k2=k2min; k2<=k2max; ++k2) {
				  for (int k1=k1min; k1<=k1max; ++k1) {
					  int ik2 = il2+k2;
					  int ik1 = il1+k1;
				  	if (0<=ik2 && ik2<nl2) {
						  if (0<=ik1 && ik1<nl1) {
							  float dv = d[0][ik2][ik1];
							  if (dv<dm) {
							    dm = dv;
							  	m2[1][il2] = k2;
							  	m1[1][il1] = k1;
							  }
							}
						}
					}
				}
				d[1][il2][il1] = dm + e[1][il2][il1];
			}
		}
		// Compute all other d[]
		float em = max(e)*nf;
		for (int ii=2; ii<nf; ++ii) {
		  int im1 = ii-1;
      for (int il2=0; il2<nl2; ++il2) {
		  	for (int il1=0; il1<nl1; ++il1) {
		      float dm = em;
		  		for (int k2=k2min; k2<=k2max; ++k2) {
		  		  for (int k1=k1min; k1<=k1max; ++k1) {
		  			  int ik2 = il2+k2;
		  		  	if (0<=ik2 && ik2<nl2) {
		  			  	int ik1 = il1+k1;
		  				  if (0<=ik1 && ik1<nl1) {
		  					  float dv = d[im1][ik2][ik1];
		  					  if (dv<dm) {
								  	int k2test = k2-m2[im1][ik2];
						      	if (j2min<=k2test && k2test<=j2max) {
								  		int k1test = k1-m1[im1][ik1];
						      		if (j1min<=k1test && k1test<=j1max) {
		  					        dm = dv;
		  					      	m2[ii][il2] = k2;
		  					      	m1[ii][il1] = k1;
		  					      }
										}
									}
		  					}
		  				}
		  			}
		  		}
		  		d[ii][il2][il1] = dm + e[ii][il2][il1];
		  	}
		  }
		}
	}

	public float[][] findShifts(float[][][] d) {
	  int n = d.length;
	  float[] u1  = new float[n]; 
	  float[] u2  = new float[n];  
		backtrack(_shifts1,_shifts2,d,_m1,_m2,u1,u2);
		return new float[][]{u1,u2};
  }
  
	private static void backtrack(
	  Sampling shifts1, Sampling shifts2, float[][][] d, 
		int[][] m1, int[][] m2, float[] u1, float[] u2) 
	{
    int nf = d.length;
    int nl2 = d[0].length;
    int nl1 = d[0][0].length;
		int nl2m = nl2-1;
		int nl1m = nl1-1;
		int mi1 = 0;
		int mi2 = 0;
		int il1 = 0;
		int il2 = 0;
		int ii = nf-1;
		float ds1 = (float)shifts1.getDelta();
		float ds2 = (float)shifts2.getDelta();
		float dl = d[ii][il2][il1];
		for (int jl2=1; jl2<nl2; ++jl2) {
			for (int jl1=1; jl1<nl1; ++jl1) {
			  float dtest = d[ii][jl2][jl1];
			  if (dtest<dl) {
				  dl = dtest;
					il1 = jl1;
					il2 = jl2;
				}
			}
		}
		u1[ii] = (float)shifts1.getValue(il1);
		u2[ii] = (float)shifts2.getValue(il2);
		while (ii>0) {
		  mi1 = m1[ii][il1];
		  mi2 = m2[ii][il2];
			il1 += mi1;
			il2 += mi2;
			--ii;
		  u1[ii] = (float)shifts1.getValue(il1);
		  u2[ii] = (float)shifts2.getValue(il2);
			if (mi1!=0) {
				float du1 = u1[ii]-u1[ii+1];
				u1[ii] = u1[ii+1]+du1;
			}
			if (mi2!=0) {
				float du2 = u2[ii]-u2[ii+1];
				u2[ii] = u2[ii+1]+du2;
			}
		}
	}

	 /**
   * Computes a sequence warped by applying specified shifts.
   * @param u input array of shifts.
   * @param g input array for the sequence to be warped.
   * @param h output array for the warped sequence.
   */
  public float[] applyShifts( 
    float[] u, float[] g) 
  {
    int n1 = g.length;
    float[] r = add(rampfloat(0f,1f,n1),u);
    int n2 = (int)(ceil(r[n1-1])-floor(r[0]));
    float[] s = new float[n2];
    float[] h = new float[n2];
		SincInterpolator si = new SincInterpolator();
    inverseLinearInterpolation(n1,1f,0f,r,n2,1f,u[0],s,r[0],r[n1-1]);
    for (int i2=0; i2<n2; ++i2) 
      h[i2] = si.interpolate(n1,1.0,0.0,g,s[i2]);
    return h;
  }


	//public float[][] applyShifts(float[] u1, float[] u2, float[] f) {
		

	//}

  ///////////////////////////////////////////////////////////////////////////
  // private fields and utils

  private float _epow = 2; // exponent used for alignment errors |f-g|^e
  private int _nl1,_nl2,_nl3; // number of vertical time lags
  private int _ods1,_ods2,_ods3; // fractional shift factor (ds=1.0/_ods is the shift interval)
  private int _k1min,_k2min,_k3min; // upper strain limit of u'
  private int _k1max,_k2max,_k3max; // lower strain limit of u'
  private int _j1min,_j2min,_j3min; // upper strain limit of u''
  private int _j1max,_j2max,_j3max; // lower strain limit of u''
	private int[][] _m1,_m2,_m3; // array of moves for the 1st, 2nd, and 3rd lags
  private Sampling _shifts1,_shifts2,_shifts3; // sampling of shift values

  private float error(float f, float g) {
    return pow(abs(f-g),_epow);
  }

  /**
   * Normalizes alignment errors to be in range [0,1].
   * @param e input/output array of alignment errors.
   */
  private static void normalizeErrors(float[][][] e) {
    int nl1 = e[0][0].length;
    int nl2 = e[0].length;
    int n1 = e.length;
    float emin = e[0][0][0];
    float emax = e[0][0][0];
    for (int i1=0; i1<n1; ++i1) {
      for (int il2=0; il2<nl2; ++il2) {
        for (int il1=0; il1<nl1; ++il1) {
          float ei = e[i1][il2][il1];
          if (ei<emin) emin = ei;
          if (ei>emax) emax = ei;
				}
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
  private static void shiftAndScale(float emin, float emax, float[][][] e) {
    int nl1 = e[0][0].length;
    int nl2 = e[0].length;
    int n1 = e.length;
    float eshift = emin;
    float escale = (emax>emin)?1.0f/(emax-emin):1.0f;
    for (int i1=0; i1<n1; ++i1) {
      for (int il2=0; il2<nl2; ++il2) {
        for (int il1=0; il1<nl1; ++il1) {
          e[i1][il2][il1] = (e[i1][il2][il1]-eshift)*escale;
				}
      }
    }
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

};




