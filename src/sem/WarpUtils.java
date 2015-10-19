package sem;

import java.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.interp.CubicInterpolator.*;
import edu.mines.jtk.util.*;

import static edu.mines.jtk.util.ArrayMath.*;

public class WarpUtils {

  /**
   * Stacks a 1 parameter semblance spectrum in time.
   * Use to determine a coarse grid sampling.
   * @param s 1 parameter semblance spectrum
   * @return stk a stacked semblance spectrum
   */
  public static float[] stackSemblance(float[][] s) {
    int n2 = s.length;
    int n1 = s[0].length;
    float[] stk = new float[n1];
    for (int i2=0; i2<n2; ++i2)
      stk = add(stk,s[i2]);
    return stk;
  }
  /**
   * Stacks a 2 parameter semblance spectrum in time.
   * Use to determine a coarse grid sampling.
   * @param s 2 parameter semblance spectrum
   * @return stk a stacked semblance spectrum
   */
  public static float[] stackSemblance(float[][][] s) {
    int n3 = s.length;
    float[] stk = new float[s[0][0].length];
    for (int i3=0; i3<n3; ++i3)
      stk = add(stk,stackSemblance(s[i3]));
    return stk;
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
  public static int[] subsample(int n, int kmin, int wb) {
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

  public static int[] subsample(float[] f, int d, int wb) {
     return subsample(f,d,wb,-1);
  }
  public static int[] subsample(float[] f, int d, int wb, int ng) {
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


  /**
   * Get the {@code grid} as sample indices.
   * @param s Sampling that this {@code grid} subsamples.
   * @param grid the subsample locations defined in terms of Sampling {@code s}.
   * @return the {@code grid} as sample indices.
   * @throws IllegalArgumentException if any of the grid locations are not valid
   *   with the given Sampling {@code s}.
   */
  public static int[][][] gridCoordsToSamples(Sampling s, float[][][] grid) {
    int n2 = grid[0].length;
    int n3 = grid.length;
    int[][][] g = new int[n3][n2][];
    for (int i3=0; i3<n3; i3++)
      for (int i2=0; i2<n2; i2++)
        g[i3][i2] = gridCoordsToSamples(s,grid[i3][i2]);
    return g;
  }

  /**
   * Get the {@code grid} as sample indices.
   * @param s Sampling that this {@code grid} subsamples.
   * @param grid the subsample locations defined in terms of Sampling {@code s}.
   * @return the {@code grid} as sample indices.
   * @throws IllegalArgumentException if any of the grid locations are not valid
   *   with the given Sampling {@code s}.
   */
  public static int[][] gridCoordsToSamples(Sampling s, float[][] grid) {
    int n2 = grid.length;
    int[][] g = new int[n2][];
    for (int i2=0; i2<n2; i2++)
      g[i2] = gridCoordsToSamples(s,grid[i2]);
    return g;
  }

  /**
   * Get the {@code grid} as sample indices.
   * @param s Sampling that this {@code grid} subsamples.
   * @param grid the subsample locations defined in terms of Sampling {@code s}.
   * @return the {@code grid} as sample indices.
   * @throws IllegalArgumentException if any of the grid locations are not valid
   *   with the given Sampling {@code s}.
   */
  public static int[] gridCoordsToSamples(Sampling s, float[] grid) {
    Almost a = new Almost();
    float f = (float)s.getFirst();
    float l = (float)s.getLast();
    int ng = grid.length;
    int[] t = new int[ng]; // temp sample indices
    int count = 0;
    int is = -1; // save last index
    for (int ig=0; ig<ng; ig++) {
      float v = grid[ig];
      if (a.ge(v,f) && a.le(v,l)) {
        int i = s.indexOfNearest(v);
        if (i!=is) { // no duplicate entries
          t[count] = i;
          count++;
        }
        is = i;
      } else {
        print("Error: value "+v+" is out of bounds! First="+f+", Last="+l);
      }
    }
    if (count!=ng) {
      print("Grid values:"); dump(grid);
      throw new IllegalArgumentException(
        "Error: Only "+count+" of "+ng+" input grid coordinates are valid "+
        "with the specified sampling "+s.toString());
    }
    return copy(count,t);
  }

 
  //public static float[] interpUS(Sampling sf, float[] u1, float[] ust) {
  //  int n1 = ust.length;
  //  float[] us = new float[n1];
  //  float[] x = getValues(sf);
  //  CubicInterpolator ci = new CubicInterpolator(Method.MONOTONIC,x,ust);
  //  for (int i1=0; i1<n1; i1++)
  //    us[i1] = ci.interpolate(x[i1]+u1[i1]);
  //  return us;
  //}

  public static float[] interpUS(Sampling sf, float[] u1, float[] ust) {
    int n1 = ust.length;
    float[] us = new float[n1];
    float[] xf = getValues(sf);
    float[] xu = new float[n1];
    for (int i1=0; i1<n1; i1++)
      xu[i1] = xf[i1]+u1[i1];
    CubicInterpolator ci = new CubicInterpolator(Method.MONOTONIC,xu,ust);
    for (int i1=0; i1<n1; i1++)
      us[i1] = ci.interpolate(xf[i1]);
    return us;
  }

  public static double[] getSubStrainMin(
      Sampling su, float[] u, double[] rmin)
  {
    float dui = 1.0f/(float)su.getDelta();
    int n = u.length;
    int nm1 = n-1;
    double[] rminNew = new double[n];
    rminNew[ 0 ] = rmin[ 0 ]-(u[ 1 ]-u[  0  ])*dui; // forward diff
    rminNew[nm1] = rmin[nm1]-(u[nm1]-u[nm1-1])*dui; // backward diff
    for (int i1=1; i1<nm1; i1++)
      rminNew[i1] = rmin[i1]-(u[i1+1]-u[i1-1])*0.5*dui;
    return rminNew;
  }

  public static double[] getSubStrainMax(
      Sampling su, float[] u, double[] rmax)
  {
    float dui = 1.0f/(float)su.getDelta();
    int n = u.length;
    int nm1 = n-1;
    double[] rmaxNew = new double[n];
    rmaxNew[ 0 ] = rmax[ 0 ]-(u[ 1 ]-u[  0  ])*dui; // forward diff
    rmaxNew[nm1] = rmax[nm1]-(u[nm1]-u[nm1-1])*dui; // backward diff
    for (int i1=1; i1<nm1; i1++)
      rmaxNew[i1] = rmax[i1]-(u[i1+1]-u[i1-1])*0.5*dui;
    return rmaxNew;
  }

  ////////////////////////////////////////////////////////////////////////////
  // Private

  private static float[] getValues(Sampling s) {
    int n = s.getCount();
    double f = s.getFirst();
    double d = s.getDelta();
    float[] x = new float[n];
    for (int i=0; i<n; i++)
      x[i] = (float)(f+i*d);
    return x;
  }

 private static void print(String s) {
    System.out.println(s);
  }

  private static class MinMax {
    float emin,emax;
    MinMax(float emin, float emax) {
      this.emin = emin;
      this.emax = emax;
    }
  }

}
