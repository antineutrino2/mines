package wt;

import static edu.mines.jtk.util.Parallel.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * A Two-sided Exponential Smoothing Filter.
 * This filter is useful because it can approximate a gaussian smoothing filter
 * with about 16 iterations of the filter.
 * @author Andrew Munoz, Colorado School of Mines
 * @version 2011.11.15
 */

public class ExponentialSmoother {

  public ExponentialSmoother(float sigma) {
    this(sigma,1);
  }
  public ExponentialSmoother(float sigma, int niter) {
    this.sigma = sigma;
    this.niter = niter;
  }

  /**
   * Exponential Filters in 1D and 2D.
   */
  public void apply(float[][] x, float[][] y) {
    float eps = sigma/sqrt(niter);
    float a = getA(eps);
    smooth(a,x,y);
    for (int i=1; i<niter; ++i) 
      smooth(a,y,y);
  }
  public void apply(float[] x, float[] y) {
    float eps = sigma/sqrt(niter);
    float a = getA(eps);
    smooth(a,x,y);
    for (int i=1; i<niter; ++i) 
      smooth(a,y,y);
  }

///////////////////////////////////////////////////////////////////////////////// 
// private

  private float sigma;
  private int niter;

  /**
   * Exponential smoother of a 2D array along both 1st and 2nd dimensions
   * In Parallel.
   */
  private static void smooth(float a, float[][] x, float[][] y) {
    smooth1P(a,x,y);
    smooth2P(a,y,y);
  }

  /**
   * Exponential smoother of a 2D array along only the 1st dimension 
   * In Parallel. 
   */
  private static void smooth1P(final float a, final float[][] x, final float[][] y) {
    final int n = x.length;
    loop(n, new LoopInt() {
      public void compute(int i) {
        smooth(a,x[i],y[i]);
      }
    });
  }

  /**
   * Exponential smoother along the 2nd dimension, with outer loops over the
   * the 2nd dimension and inner loops over the 1st dimension
   * In Parallel. 
   */
  private static void smooth2P(final float a, final float[][] x, final float[][] y) {
    final int n1 = x[0].length;
    final int n2 = x.length;
    final float b = 1.0f-a;
    final int ke = n1%nthread;
    final int kb = (int)(n1/nthread);
    loop(1, nthread+1, new LoopInt() {
      public void compute(int i3) {
        int sta = (i3-1)*kb;
	int end = i3*kb;
        if (i3 == nthread) end += ke;
        for (int i1=sta; i1<end; ++i1) {
          y[0][i1] = x[0][i1];
	}
        for (int i2=1; i2<n2; ++i2) {
          for (int i1=sta; i1<end; ++i1) {
            y[i2][i1] = a*y[i2-1][i1] + b*x[i2][i1];
          }
        }
        for (int i1=sta; i1<end; ++i1) {
          y[n2-1][i1] = (x[n2-1][i1] + a*y[n2-1][i1])/(1.0f+a);
        }
        for (int i2=n2-2; i2>=0; --i2) {
          for (int i1=sta; i1<end; ++i1) {
            y[i2][i1] = a*y[i2+1][i1] + b*y[i2][i1];
          }
        }
      }
    });
  } 

  /**
   * Exponential smoother of a 1D array.
   */
  private static void smooth(float a, float[] x, float[] y) {
    // Improved with the "Simon simplification"
    int n = x.length;
    float b = 1.0f-a;
    float yi = y[0] = x[0];
    for (int i=1; i<n; ++i)
      y[i] = yi = a*yi+b*x[i];
    y[n-1] = yi = (a*yi+x[n-1])/(1.0f+a);  // Combined steps algebraically 
    for (int i=n-2; i>=0; --i)
      y[i] = yi = a*yi+b*y[i];
  }

///////////////////////////////////////////////////////////////////////////////// 
// Utils

 public static float mean(float[] x) {
    int n = x.length;
    float sum = 0.0f;
    for (int i=0; i<n; ++i)
      sum += x[i];
    return sum/n;
  }

  public static float mean(float[][] x) {
    int n = x.length;
    float sum = 0.0f;
    for (int i=0; i<n; ++i)
      sum += mean(x[i]);
    return sum/n;
  }

// Thanks to Luming for helping me derive this using Vieta's formula
  private static float getA(float e) {
    float e2 = pow(e,2);
    float a = (e2+1-sqrt(2*e2+1))/e2;
    return a;
  }

  /**
   * Number of threads used for parallel processing.
   */
  private static int nthread = Runtime.getRuntime().availableProcessors();



}
