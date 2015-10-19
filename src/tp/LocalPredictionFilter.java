package tp;

import static edu.mines.jtk.util.Parallel.*;
import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.dsp.*;

/**
 * Local Prediction Error Filter for 1D or 2D arrays.
 * @author Andrew Munoz, Colorado School of Mines
 * @version 01.27.14
 */

public class LocalPredictionFilter {

  public LocalPredictionFilter(float sigma) {
    this(sigma,0,1);
  }
  public LocalPredictionFilter(float sigma, int dir) {
    this(sigma,dir,1);
  }
  public LocalPredictionFilter(float sigma, int dir, int niter) {
    this.sigma = sigma;
    this.dir = dir;
    this.niter = niter;
    sm = new ExponentialSmoother(sigma,niter);
  }

// For Benchmarking:
  public void apply(float[] x) {
     int n1 = x.length;
     float[] y = new float[n1];
     float[] a = new float[n1];
     float[] e = new float[n1];
     y = getTrend(x);
     x = sub(x,y);
     a = calculateCoeff(x);
     e = predictionError(x,a);
  }
  public void apply(float[][] x) {
     int n1 = x[0].length; 
     int n2 = x.length;
     float[][] y = new float[n2][n1];
     y = getTrend(x);
     x = sub(x,y);
     float[][][] a = calculateCoeff(x);
     float[][] a1 = a[0]; 
     float[][] a2 = a[1];
     float[][][] e = predictionError(x,a1,a2);
  }

////////////////////////////////////////////////////////////////////////////////////////////
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
    float[] e  = new float[n];
    float[] e0 = new float[n];
    float[] em = new float[n];
    float[] e2 = new float[n];
    float[] y0 = new float[n];
    float[] ym = new float[n];
    float[] xm = new float[n];
    float de0,dem,de2,de=1;
    xm[0] = x0[0];
    for (int i=1; i<n; ++i)
      xm[i] = x0[i-1];
    ym = mul(a,xm); 
    y0 = mul(a,x0); 
    em = sub(xm,y0); 
    e0 = sub(x0,ym);
    if (dir==2) e2 = add(mul(e0,e0),mul(em,em));
    de0 = sum(mul(em,x0)); 
    dem = sum(mul(e0,xm)); 
    de2 = de0+dem;
    if (dir==0) {de = dem; copy(e0,e);} 
    if (dir==1) {de = de0; copy(em,e);}
    if (dir==2) {de = de2; copy(e2,e);}
    if (de==0) System.out.println("de=0");
    else System.out.println("de!=0, it is "+Float.toString(de));
    return e;
  }
  public float[] calculateCoeff(float[] x0) {
    int n = x0.length;
    float[] a   = new float[n];
    float[] xm  = new float[n];
    float[] r1  = new float[n];
    float[] r11 = new float[n];
    xm[0] = x0[0];
    for (int i=1; i<n; ++i)
      xm[i-1] = x0[i];
    if (dir==0) {r1 = mul(x0,xm); r11 = mul(xm,xm);}
    if (dir==1) {r1 = mul(xm,x0); r11 = mul(x0,x0);}
    if (dir==2) {r1 = mul(2f,mul(x0,xm)); r11 = add(mul(x0,x0),mul(xm,xm));}
    sm.apply(r1,r1);
    sm.apply(r11,r11);
    a = div(r1,r11);
    return a;
  }
  public float[] getTrend(float[] x) {
    int n = x.length;
    float[] y = new float[n];
    sm.apply(x,y);
    return y;
  }

////////////////////////////////////////////////////////////////////////////////////////////
// 2D prediction filter

 /**
  *  2D linear prediction filter.
  *  
  *  Removal of all predictions:
  *  dir0: e003 = x00 - (a1*x0m + a2*xm0)
  *  dir1: e0m3 = x0m - (a1*x00 + a2*xmm)
  *  dir2: em03 = xm0 - (a1*xmm + a2*x00)
  *  dir3: emm3 = xmm - (a1*xm0 + a2*x0m)
  *  dir4: find error using all directions
  *  Component directions of prediction (a1 and a2) removed individually too 
  *
  *  @param x00 the input array with zero mean
  *  @return e1  the predicted vertical component of the image shown 
  *  @return e2  the predicted horizontal component of the image shown
  *  @return e3  both of the predicted components removed
  */
  public float[][][] predictionError(float[][] x00, float[][] a1, float[][] a2) {
    int n1 = x00[0].length;
    int n2 = x00.length;
    float[][] e1  = new float[n2][n1];
    float[][] e2  = new float[n2][n1];
    float[][] e3  = new float[n2][n1];
    float[][] esum= new float[n2][n1];
    float[][] x0m = new float[n2][n1];
    float[][] xm0 = new float[n2][n1];
    float[][] xmm = new float[n2][n1];
    float[][] y0m1= new float[n2][n1],y0m2= new float[n2][n1];
    float[][] ym01= new float[n2][n1],ym02= new float[n2][n1];
    float[][] ymm1= new float[n2][n1],ymm2= new float[n2][n1];
    float[][] y001= new float[n2][n1],y002= new float[n2][n1];
    float[][] e001= new float[n2][n1],e002= new float[n2][n1],e003= new float[n2][n1];
    float[][] e0m1= new float[n2][n1],e0m2= new float[n2][n1],e0m3= new float[n2][n1];
    float[][] em01= new float[n2][n1],em02= new float[n2][n1],em03= new float[n2][n1];
    float[][] emm1= new float[n2][n1],emm2= new float[n2][n1],emm3= new float[n2][n1];
    float de001,de0m1,dem01,demm1,de002,de0m2,dem02,demm2,desum,de=1;
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        x0m[i2][i1] = x00[i2  ][i1-1];
        xm0[i2][i1] = x00[i2-1][i1  ];
        xmm[i2][i1] = x00[i2-1][i1-1];
      }
    }
    y001 = sub(x00,mul(a1,x0m)); y002 = sub(x00,mul(a2,xm0));
    y0m1 = sub(x0m,mul(a1,x00)); y0m2 = sub(x0m,mul(a2,xmm));
    ym01 = sub(xm0,mul(a1,xmm)); ym02 = sub(xm0,mul(a2,x00));
    ymm1 = sub(xmm,mul(a1,xm0)); ymm2 = sub(xmm,mul(a2,x0m));
    e001 = sub(x00,y001); e002 = sub(x00,y002); e003 = sub(x00,add(y001,y002));
    e0m1 = sub(x0m,y0m1); e0m2 = sub(x0m,y0m2); e0m3 = sub(x0m,add(y0m1,y0m2));
    em01 = sub(xm0,ym01); em02 = sub(xm0,ym02); em03 = sub(xm0,add(ym01,ym02));
    emm1 = sub(xmm,ymm1); emm2 = sub(xmm,ymm2); emm3 = sub(xmm,add(ymm1,ymm2));
    if (dir==4) {
      float[][] est1  = new float[n2][n1];
      float[][] est2  = new float[n2][n1];
      est1 = add(mul(add(e001,e002),add(e001,e002)),mul(add(emm1,emm2),add(e0m1,emm2)));
      est2 = add(mul(add(e0m1,e0m2),add(e0m1,e0m2)),mul(add(em01,em02),add(em01,em02)));
      esum = add(est1,est2);
    }
    de001 = sum(mul(add(e0m1,em02),x00)); de002 = sum(mul(add(e0m1,em02),x00));
    de0m1 = sum(mul(add(e001,emm2),x0m)); de0m2 = sum(mul(add(e001,emm2),x0m));
    dem01 = sum(mul(add(emm1,e002),xm0)); dem02 = sum(mul(add(emm1,e002),xm0));
    demm1 = sum(mul(add(em01,e0m2),xmm)); demm2 = sum(mul(add(em01,e0m2),xmm));
    desum = de001+de0m1+dem01+demm1+de002+de0m2+dem02+demm2;
    if (dir==0) {de = de001+de002; copy(e001,e1); copy(e002,e2); copy(e003,e3);}
    if (dir==1) {de = de0m1+de0m2; copy(e0m1,e1); copy(e0m2,e2); copy(e0m3,e3);} 
    if (dir==2) {de = dem01+dem02; copy(em01,e1); copy(em02,e2); copy(em03,e3);} 
    if (dir==3) {de = demm1+demm2; copy(emm1,e1); copy(emm2,e2); copy(emm3,e3);} 
    if (dir==4) {de = desum; copy(esum,e1);} 
    if (de==0) System.out.println("de=0");
    else System.out.println("de!=0, it is "+Float.toString(de));
    return new float[][][]{e1,e2,e3};
  }

 /**
  *  Coeffiecients for a 2D linear prediction filter.
  *  Solve the positive-definite system:
  *  
  *  | r11  r12 | |a1|   | rc1 |
  *  |          | |  | = |     |
  *  | r21  r22 | |a2|   | rc2 |
  * 
  *  Solve using Cramer's rule:
  *  a1 = (rc1*r22-rc2*r12)/(r11*r22-r12*r21)
  *  a2 = (rc2*r11-rc2*r21)/(r11*r22-r12*r21)
  *
  *  dir0: x00 = a1*x0m + a2*xm0
  *  dir1: x0m = a1*x00 + a2*xmm
  *  dir2: xm0 = a1*xmm + a2*x00
  *  dir3: xmm = a1*xm0 + a2*x0m
  *  dir4: predict all 
  *  @param x00 the input array with zero mean
  *  @return a1  the coefficient array for the vertical direction
  *  @return a2  the coefficient array for the horizontal direction
  */
  public float[][][] calculateCoeff(float[][] x00) {
    int n1 = x00[0].length; 
    int n2 = x00.length; 
    float[][] a1  = new float[n2][n1];
    float[][] a2  = new float[n2][n1];
    float[][] x0m = new float[n2][n1];
    float[][] xm0 = new float[n2][n1];
    float[][] xmm = new float[n2][n1];
    float[][] r11 = new float[n2][n1];
    float[][] r12 = new float[n2][n1];
    float[][] r21 = new float[n2][n1];
    float[][] r22 = new float[n2][n1];
    float[][] rc1 = new float[n2][n1];
    float[][] rc2 = new float[n2][n1];
    float[][] den = new float[n2][n1];
    float[][] nu1 = new float[n2][n1];
    float[][] nu2 = new float[n2][n1];
    for (int i2=1; i2<n2; ++i2) {
      for (int i1=1; i1<n1; ++i1) {
        x0m[i2][i1] = x00[i2  ][i1-1];
        xm0[i2][i1] = x00[i2-1][i1  ];
        xmm[i2][i1] = x00[i2-1][i1-1];
      }
    } if (dir==0) {
        r11 = mul(x0m,x0m); r22 = mul(xm0,xm0);
        r12 = mul(xm0,x0m); r21 = mul(xm0,x0m);
        rc1 = mul(x00,x0m); rc2 = mul(x00,xm0);
    } if (dir==1) {
        r11 = mul(x00,x00); r22 = mul(xmm,xmm);
        r12 = mul(x00,xmm); r21 = mul(x00,xmm);
        rc1 = mul(x00,x0m); rc2 = mul(xmm,x0m);
    } if (dir==2) {
        r11 = mul(xmm,xmm); r22 = mul(x00,x00);
        r12 = mul(x00,xmm); r21 = mul(x00,xmm);
        rc1 = mul(xmm,xm0); rc2 = mul(x00,xm0);
    } if (dir==3) {
        r11 = mul(xm0,xm0); r22 = mul(x0m,x0m);
        r12 = mul(xm0,x0m); r21 = mul(xm0,x0m);
        rc1 = mul(xmm,xm0); rc2 = mul(xmm,x0m);
    } if (dir==4) {
        float[][] r1t1 = new float[n2][n1];
        float[][] r1t2 = new float[n2][n1];
        r1t1 = add(mul(x00,x00),mul(xm0,xm0));
        r1t2 = add(mul(x0m,x0m),mul(xmm,xmm));
        r11 = add(r1t1,r1t2);
        r12 = mul(2f,add(mul(xm0,x0m),mul(x00,xmm)));
        copy(r11,r22); copy(r12,r21);
        rc1 = mul(2f,add(mul(x00,x0m),mul(xmm,xm0)));
        rc2 = mul(2f,add(mul(x00,xm0),mul(xmm,x0m)));
    }
    // Smooth
    sm.apply(r11,r11); sm.apply(r22,r22);
    sm.apply(r12,r12); sm.apply(r21,r21);
    sm.apply(rc1,rc1); sm.apply(rc2,rc2);
    // Using Cramer's rule:
    den = sub(mul(r11,r22),mul(r12,r21));
    nu1 = sub(mul(rc1,r22),mul(rc2,r12));
    nu2 = sub(mul(rc2,r11),mul(rc1,r21));
    a1 = div(nu1,den);
    a2 = div(nu2,den);
    return new float[][][]{a1,a2};
  }
  public float[][] getTrend(float[][] x) {
    int n2 = x.length;
    int n1 = x[0].length;
    float[][] y = new float[n2][n1];
    sm.apply(x,y);
    return y;
  }

/////////////////////////////////////////////////////////////////
// private
  
  private float sigma;
  private int dir,niter;
  private ExponentialSmoother sm;


/////////////////////////////////////////////////////////////////
// utils

  private void copy(final float[] x, final float[] y) {
    final int n1 = x.length;
    loop(n1,new LoopInt() {   
    public void compute(int i1) {
       y[i1] = x[i1];
    }});
  }
  private void copy(final float[][] x, final float[][] y) {
    final int n2 = x.length;
    loop(n2,new LoopInt() {   
    public void compute(int i2) {
      copy(x[i2],y[i2]);
    }});
  }

};
