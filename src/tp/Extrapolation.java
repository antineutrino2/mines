package tp;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Extrapolation using linear regression and predictive error methods.
 * @author Andrew Munoz, Colorado School of Mines
 * @version 01.27.14
 */

public class Extrapolation {

  public Extrapolation() {
    this(1);
  }
  public Extrapolation(float sigma) {
    this(sigma,16);
  }
  public Extrapolation(float sigma, int niter) {
    this.sigma = sigma;
    this.niter = niter;
    sm = new ExponentialSmoother(sigma,niter);
  }

///////////////////////////////////////////////////////////////////////////////
// 1D prediction filter

 /**
  *  1D prediction error filter.
  *  dir0: x0 = a*xm : future from past
  *  dir1: xm = a*x0 : past from future
  *  dir2: predict both 
  *  @param x0 input array with zero mean
  *  @param a  array of coefficients
  *  @return e error 
  */
  public float[] predictionError(float[] x0, float[] a) {
    int n = x0.length;
    float[] xm = new float[n];
    xm[0] = x0[0];
    for (int i=1; i<n; ++i)
      xm[i] = x0[i-1];
    float[] ym = mul(a,xm); 
    float[] y0 = mul(a,x0); 
    float[] em = sub(xm,y0); 
    float[] e0 = sub(x0,ym);
    float[] e2 = add(mul(e0,e0),mul(em,em));
    float de0 = sum(mul(em,x0)); 
    float dem = sum(mul(e0,xm)); 
    float de2 = de0+dem;
    if (de2==0) System.out.println("de=0");
    else System.out.println("de!=0, it is "+Float.toString(de2));
    return e2;
  }
  public float[] calculateCoeff(float[] x0) {
    int n = x0.length;
    float[] xm  = new float[n];
    for (int i=1; i<n; ++i)
      xm[i-1] = x0[i];
    xm[n-1] = x0[n-1];
    float[] r1  = mul(2f,mul(x0,xm));
    float[] r11 = add(mul(x0,x0),mul(xm,xm));
    sm.apply(r1,r1);
    sm.apply(r11,r11);
    float[] a = div(r1,r11);
    return a;
  }
  public float[] getTrend(float[] x) {
    int n = x.length;
    float[] y = new float[n];
    sm.apply(x,y);
    return y;
  }


 /**
   * Computes a linear extrapolation for the ends of a function.
   * Also resamples if dz1 and dz2 are different.
   */
  private float[] extrapolateLinear(
    Sampling sz1, Sampling sz2, float[] l, float pnull)
  {
    int nz1 = sz1.getCount(); 
    int nz2 = sz2.getCount(); 
    double dz1 = sz1.getDelta(); 
    double dz2 = sz2.getDelta(); 
    double fz1 = sz1.getFirst(); 
    double fz2 = sz2.getFirst(); 
    double lz2 = sz2.getLast();
    double lz1 = sz1.getLast();
    float fdz2 = (float)dz2;
    float[] le = new float[nz2];
    float[] lr = l;
    if (dz1!=dz2) {
      lr = resampleLinear(sz1,sz2,l);
      sz1 = new Sampling(lr.length,dz2,fz1);
    }
    // extrapolate beginning
    int nr = (int)(lr.length*0.5);
    int nb = infl((fz1-fz2)/dz2);
    int nm = lr.length;
    int ne = nz2-(nb+nm);
    if (pnull<0) {
      for (int i=0; i<nb; ++i)
        le[i] = pnull;
    }
    else if (nb>0) {
      float[] bb = getRegressionTerms(sz1,lr,nr,1);
      float b0 = l[0];
      float b1 = -bb[1];
      for (int i=0; i<nb; ++i)
        le[nb-i-1] = b0+i*b1*fdz2;
    }
    else if (nb>0) {
      // if time-depth curve, extrapolate to t=z=fz2
      float b0 = (float)fz2;
      float b1 = (l[0]-b0)/nb;
      for (int i=0; i<nb; ++i)
        le[i] = b0+i*b1;
    }
    // fill values
    for (int i=0; i<nm; ++i)
      le[i+nb] = lr[i];
    // extrapolate end
    if (pnull<0) {
      int nf = nb+nm;
      for (int i=0; i<ne; ++i)
        le[i+nf] = pnull;
    } 
    else if (ne>0) {
      int nf = nb+nm;
      float[] bb = getRegressionTerms(sz1,lr,nr,-1);
      float b2 = bb[1];
      float fltz = lr[nm-1];
      for (int i=0; i<ne; ++i)
        le[i+nf] = fltz+i*b2*fdz2;
    }
    return le;
  }
  private float[] extrapolateLinear(
    Sampling sz1, Sampling sz2, float[] l) 
  {
    return extrapolateLinear(sz1,sz2,l,1);
  }
  
  /**
   * Resamples a function using linear interpolation.
   */
  private float[] resampleLinear(Sampling sz1, Sampling sz2, float[] x) { 
    int nz1 = sz1.getCount(); 
    double dz2 = sz2.getDelta(); 
    double dz1 = sz1.getDelta(); 
    double fz1 = sz1.getFirst(); 
    double lz1 = sz1.getLast();
    int nzi = inro((lz1-fz1)/dz2);
    float[] xi = new float[nzi];
    LinearInterpolator li = new LinearInterpolator(); 
    li.setUniform(nz1,dz1,fz1,x);
	  li.interpolate(nzi,dz2,fz1,xi);
    return xi;
  }
  
  /**
   * Solves for simple linear regression terms.
   */
  private float[] getRegressionTerms(Sampling sz, float[] l, int nr, int end) {
    int nl = l.length;
    float[] zvl = floats(sz.getValues());
    int nre = nl-nr;
    int n1 = (end>0)?0:nre;
    int n2 = (end>0)?nr:nl;
    float ym = meanW(l,n1,n2);
    float xm = meanW(zvl,n1,n2);
    float num=0f,den=0f;
    for (int j=n1; j<n2; ++j) {
      float xi = zvl[j]-xm;
      num += (l[j]-ym)*xi;
      den += xi*xi;
    }
    float b1 = num/den;
    float b0 = ym - b1*xm;
    return new float[]{b0,b1};
  }
  private static float meanW(float[] x, int w1, int w2) {
	  int n = w2-w1;
	  float[] xw = copy(n,w1,x);
	  return sum(xw)/n;
  }

/////////////////////////////////////////////////////////////////
// private
  
  private float sigma;
  private int dir,niter;
  private ExponentialSmoother sm;

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

};
