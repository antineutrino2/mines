package dtw;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.interp.CubicInterpolator.*;
import edu.mines.jtk.util.*;

import static edu.mines.jtk.util.ArrayMath.*;

public class WarpUtils {

  /**
   * Computes the NRMS value of the warped PS trace {@code h} and the PP trace
   * {@code f}. The NRMS value is the RMS of the difference between {@code f}
   * and {@code h} divided by the average RMS of {@code f} and {@code h}.
   * @param n1Max the length of {@code f} and {@code h} to use for this
   *  computation.
   * @param f the PP trace.
   * @param h the warped PS trace.
   * @return the NRMS value, this measure is between 0.0 and 2.0.
   */
  public static float computeNrms(int n1Max, float[] f, float[] h) {
    int n1f = f.length;
    int n1h = h.length;
    Check.argument(n1Max<=n1f,"n1Max<=n1f");
    Check.argument(n1Max<=n1h,"n1Max<=n1h");
    float[] fs = copy(n1Max,f);
    float[] hs = copy(n1Max,h);
    float scale = 1.0f/n1Max;
    float[] d = sub(hs,fs);
    float frms = sqrt(sum(mul(fs,fs))*scale);
    float hrms = sqrt(sum(mul(hs,hs))*scale);
    float drms = sqrt(sum(mul( d, d))*scale);
    float nrms = (2.0f*drms)/(frms+hrms);
    return nrms;
  }

  /**
   * Computes the NRMS value of the warped PS trace {@code h} and the PP trace
   * {@code f}. The NRMS value is the RMS of the difference between {@code f}
   * and {@code h} divided by the average RMS of {@code f} and {@code h}.
   * @param n1Max the length of {@code f} and {@code h} to use for this
   *  computation.
   * @param f the PP traces.
   * @param h the warped PS traces.
   * @return the NRMS value, this measure is between 0.0 and 2.0.
   */
  public static float computeNrms(
      final int n1Max, final float[][] f, final float[][] h)
  {
    int n1f = f[0].length;
    int n1h = h[0].length;
    Check.argument(n1Max<=n1f,"n1Max<=n1f");
    Check.argument(n1Max<=n1h,"n1Max<=n1h");
    int n2 = f.length;
    float scale = 1.0f/(n1Max*n2);
    float[] rms = Parallel.reduce(n2,new Parallel.ReduceInt<float[]>() {
      @Override
      public float[] compute(int i2) {
        float[] fhdSq = new float[3];
        for (int i1=0; i1<n1Max; i1++) {
          float fv = f[i2][i1];
          float hv = h[i2][i1];
          fhdSq[0] += fv*fv;
          fhdSq[1] += hv*hv;
          fhdSq[2] += (hv-fv)*(hv-fv);
        }
        return fhdSq;
      }

      @Override
      public float[] combine(float[] fhdSq1, float[] fhdSq2) {
        float[] rms = new float[3];
        rms[0] = fhdSq1[0]+fhdSq2[0];
        rms[1] = fhdSq1[1]+fhdSq2[1];
        rms[2] = fhdSq1[2]+fhdSq2[2];
        return rms;
      }
    });
    float frms = sqrt(rms[0]*scale);
    float hrms = sqrt(rms[1]*scale);
    float drms = sqrt(rms[2]*scale);
    float nrms = (2.0f*drms)/(frms+hrms);
    return nrms;
  }

  /**
   * Computes the NRMS value of the warped PS trace {@code h} and the PP trace
   * {@code f}. The NRMS value is the RMS of the difference between {@code f}
   * and {@code h} divided by the average RMS of {@code f} and {@code h}.
   * @param n1Max the length of {@code f} and {@code h} to use for this
   *  computation.
   * @param f the PP traces.
   * @param h the warped PS traces.
   * @return the NRMS value, this measure is between 0.0 and 2.0.
   */
  public static float computeNrms(
      final int n1Max, final float[][][] f, final float[][][] h)
  {
    int n1f = f[0][0].length;
    int n1h = h[0][0].length;
    Check.argument(n1Max<=n1f,"n1Max<=n1f");
    Check.argument(n1Max<=n1h,"n1Max<=n1h");
    final int n3 = f.length;
    final int n2 = f[0].length;
    final int n23 = n2*n3;
    float scale = 1.0f/(n1Max*n23);
    float[] rms = Parallel.reduce(n23,new Parallel.ReduceInt<float[]>() {
      @Override
      public float[] compute(int i23) {
        int i2 = i23%n2;
        int i3 = i23/n2;
        float[] fhdSq = new float[3];
        for (int i1=0; i1<n1Max; i1++) {
          float fv = f[i3][i2][i1];
          float hv = h[i3][i2][i1];
          fhdSq[0] += fv*fv;
          fhdSq[1] += hv*hv;
          fhdSq[2] += (hv-fv)*(hv-fv);
        }
        return fhdSq;
      }

      @Override
      public float[] combine(float[] fhdSq1, float[] fhdSq2) {
        float[] rms = new float[3];
        rms[0] = fhdSq1[0]+fhdSq2[0];
        rms[1] = fhdSq1[1]+fhdSq2[1];
        rms[2] = fhdSq1[2]+fhdSq2[2];
        return rms;
      }
    });
    float frms = sqrt(rms[0]*scale);
    float hrms = sqrt(rms[1]*scale);
    float drms = sqrt(rms[2]*scale);
    float nrms = (2.0f*drms)/(frms+hrms);
    return nrms;
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

  /**
   * Computes the slope from the given {@code vpvs} and compression constant.
   * @param vpvs
   * @param c
   * @return the slope
   */
  public static float getSlope(float vpvs, float c) {
    return (vpvs-2.0f*c+1.0f)/(2.0f*c);
  }

  /**
   * Computes time shifts between PP and PS2 images, here referred to as 'U2'
   * time shifts.
   * @param s1
   * @param u1
   * @param sS
   * @param uS
   * @return
   */
  public static float[] computeU2(
      Sampling s1, float[] u1, Sampling sS, float[] uS)
  {
    int n1 = u1.length;
    float[] x1 = getValues(s1);
    CubicInterpolator ci =
      new CubicInterpolator(Method.LINEAR,getValues(sS),uS);
    float[] u2 = new float[n1];
    for (int i1=0; i1<n1; i1++)
      u2[i1] = u1[i1]+ci.interpolate(x1[i1]+u1[i1]);
    return u2;
  }

  public static float[][] ps1ToPpTime(
      Sampling sf, float[] f, Sampling sg, float[][] g)
  {
    int n2 = g.length;
    float[][] gp = new float[n2][];
    for (int i2=0; i2<n2; i2++)
      gp[i2] = ps1ToPpTime(sf,f,sg,g[i2]);
    return gp;
  }

  public static float[] ps1ToPpTime(
      Sampling sf, float[] u1, Sampling sg, float[] g)
  {
    int n1 = u1.length;
    int ng = g.length;
    float[] xg = new float[ng];
    for (int ig=0; ig<ng; ig++)
      xg[ig] = (float)sg.getValue(ig);
    CubicInterpolator ci = new CubicInterpolator(Method.LINEAR,xg,g);
    float[] gp = new float[n1];
    for (int i1=0; i1<n1; i1++)
      gp[i1] = ci.interpolate((float)(sf.getValue(i1)+u1[i1]));
    return gp;
  }

  /**
   * Computes an array of VpVs ratios an array of shifts u
   * using a backward difference approximation.
   * The relationship is defined as vpvs(t) = 1+2*(du/dt)
   * @param u
   * @return computed vpvs values.
   */
  public static float[] vpvsBd(Sampling su, float[] u) {
    float dui = 1.0f/(float)su.getDelta();
    int n = u.length;
    float[] vpvs = new float[n];
    vpvs[0] = 1.0f + 2.0f*(u[1]-u[0])*dui; // at i1=0, forward difference
    for (int i1=1; i1<n; ++i1)
      vpvs[i1] = 1.0f + 2.0f*(u[i1]-u[i1-1])*dui;
    return vpvs;
  }

  /**
   * Computes an array of VpVs ratios an array of shifts u. This relationship
   * is defined as vpvs(t) = 1+2*(du/dt).
   * @param u the shifts.
   * @return computed interval Vp/Vs values.
   */
  public static float[][][] vpvs(Sampling su, float[][][] u) {
    return vpvs(su,u,1.0f);
  }

  /**
   * Computes an array of VpVs ratios an array of shifts u. This relationship
   * is defined as vpvs(t) = 1+2*(du/dt).
   * @param u the shifts.
   * @param c the initial scale factor applied.
   * @return computed interval Vp/Vs values.
   */
  public static float[][][] vpvs(Sampling su, float[][][] u, float c) {
    int n3 = u.length;
    float[][][] vpvs = new float[n3][][];
    for (int i3=0; i3<n3; i3++)
      vpvs[i3] = vpvs(su,u[i3],c);
    return vpvs;
  }

  /**
   * Computes an array of VpVs ratios an array of shifts u. This relationship
   * is defined as vpvs(t) = 1+2*(du/dt).
   * @param u the shifts.
   * @return computed interval Vp/Vs values.
   */
  public static float[][] vpvs(Sampling su, float[][] u) {
    return vpvs(su,u,1.0f);
  }

  /**
   * Computes an array of VpVs ratios an array of shifts u. This relationship
   * is defined as vpvs(t) = 1+2*(du/dt).
   * @param u the shifts.
   * @param c the initial scale factor applied.
   * @return computed interval Vp/Vs values.
   */
  public static float[][] vpvs(Sampling su, float[][] u, float c) {
    int n2 = u.length;
    float[][] vpvs = new float[n2][];
    for (int i2=0; i2<n2; i2++)
      vpvs[i2] = vpvs(su,u[i2],c);
    return vpvs;
  }

  /**
   * Computes an array of VpVs ratios an array of shifts u. This relationship
   * is defined as vpvs(t) = 1+2*(du/dt).
   * @param u the shifts.
   * @return computed interval Vp/Vs values.
   */
  public static float[] vpvs(Sampling su, float[] u) {
    return vpvs(su,u,1.0f);
  }

  /**
   * Computes an array of VpVs ratios an array of shifts u. This relationship
   * is defined as vpvs(t) = 1+2*(du/dt).
   * @param u the shifts.
   * @param c the initial scale factor applied.
   * @return computed interval Vp/Vs values.
   */
  public static float[] vpvs(Sampling su, float[] u, float c) {
    int n = u.length;
    float twoC = 2.0f*c;
    float twoCm1 = twoC-1.0f;
    float[] vpvs = new float[n];
    float[] du = firstDerivative(su,u);
    for (int i1=0; i1<n; i1++)
      vpvs[i1] = twoCm1 + twoC*du[i1];
    return vpvs;
  }

  /**
   * Computes gammaS, a measure of the time delay between two split shear 
   * waves, from shifts between PP-PS1 {@code u1} and PP-PS2 {@code u2}.
   * @param sf PP time Sampling.
   * @param u1 shifts between PP-PS1 in PP time.
   * @param u2 shifts between PP-PS2 in PP time.
   * @return gammaS
   */
  public static float[] gammaSu12(Sampling sf, float[] u1, float[] u2) {
    int n1 = u1.length;
    Check.argument(n1==sf.getCount(),"u1 consistent with sampling");
    Check.argument(n1==u2.length,"u1.length==u2.length");
    float[] du1 = firstDerivative(sf,u1);
    float[] du2 = firstDerivative(sf,u2);
    float[] gs = new float[n1];
    for (int i1=0; i1<n1; i1++)
      gs[i1] = (2.0f*(du2[i1]-du1[i1])) / (1.0f+2.0f*du1[i1]);
    return gs;
  }

  /**
   * Computes gammaS, a measure of the time delay between two split shear 
   * waves, from shifts between PP-PS1 {@code u1} and PS1-PS2 {@code uS}.
   * @param sf PP time Sampling.
   * @param u1 shifts between PP-PS1 in PP time.
   * @param uS shifts between PS1-PS2 in PP time.
   * @return gammaS
   */
  public static float[] gammaSu1S(Sampling sf, float[] u1, float[] uS) {
    int n1 = u1.length;
    Check.argument(n1==sf.getCount(),"u1 consistent with sampling");
    Check.argument(n1==uS.length,"u1.length==uS.length");
    float[] du1 = firstDerivative(sf,u1);
    float[] duS = firstDerivative(sf,uS);
    float[] gs = new float[n1];
    for (int i1=0; i1<n1; i1++)
      gs[i1] = (2.0f*duS[i1]) / (1.0f+2.0f*du1[i1]);
    return gs;
  }

  /**
   * Computes gammaS, a measure of the time delay between two split shear 
   * waves, from shifts between PP-PS2 {@code u2} and PS1-PS2 {@code uS}.
   * @param sf PP time Sampling.
   * @param u2 shifts between PP-PS2 in PP time.
   * @param uS shifts between PS1-PS2 in PP time.
   * @return gammaS
   */
  public static float[] gammaSu2S(Sampling sf, float[] u2, float[] uS) {
    int n1 = u2.length;
    Check.argument(n1==sf.getCount(),"u2 consistent with sampling");
    Check.argument(n1==uS.length,"u2.length==uS.length");
    float[] du2 = firstDerivative(sf,u2);
    float[] duS = firstDerivative(sf,uS);
    float[] gs = new float[n1];
    for (int i1=0; i1<n1; i1++)
      gs[i1] = (2.0f*duS[i1]) / (1.0f+2.0f*(du2[i1]-duS[i1]));
    return gs;
  }

  /**
   * Applies the shifts {@code u} to the trace {@code g}.
   * @param sf the sampling that {@code g} is being warped to.
   * @param u the shifts to apply to {@code g}.
   * @param sg the sampling of {@code g}.
   * @param g the trace to be warped.
   * @return the warped trace.
   */
  public static float[] applyShifts(
      Sampling sf, final float[] u,
      Sampling sg, final float[] g)
  {
    final int n1 = u.length;
    final int ng = g.length;
    final double dg = sg.getDelta();
    final double df = sf.getDelta();
    final double fg = sg.getDelta();
    final double ff = sf.getDelta();
    final float[] hf = new float[n1];
    final SincInterp si = new SincInterp();
    double v = ff;
    for (int i1=0; i1<n1; i1++, v=ff+i1*df)
      hf[i1] = si.interpolate(ng,dg,fg,g,(float)v+u[i1]);
    return hf;
  }

  /**
   * Applies the shifts {@code u} to the trace {@code g}.
   * @param sf the sampling that {@code g} is being warped to.
   * @param u the shifts to apply to {@code g}.
   * @param sg the sampling of {@code g}.
   * @param g the trace to be warped.
   * @return the warped image.
   */
  public static float[][] applyShifts(
      Sampling sf, final float[][] u,
      Sampling sg, final float[][] g)
  {
    final int n1 = u[0].length;
    final int n2 = u.length;
    final int ng = g[0].length;
    final double dg = sg.getDelta();
    final double df = sf.getDelta();
    final double fg = sg.getDelta();
    final double ff = sf.getDelta();
    final float[][] hf = new float[n2][n1];
    final SincInterp si = new SincInterp();
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      double v = ff;
      for (int i1=0; i1<n1; i1++, v=ff+i1*df) {
        hf[i2][i1] = si.interpolate(ng,dg,fg,g[i2],(float)v+u[i2][i1]);
      }
    }});
    return hf;
  }

  /**
   * Applies the shifts {@code u} to the trace {@code g}.
   * @param sf the sampling that {@code g} is being warped to.
   * @param u the shifts to apply to {@code g}.
   * @param sg the sampling of {@code g}.
   * @param g the trace to be warped.
   * @return the warped image.
   */
  public static float[][][] applyShifts(
      Sampling sf, final float[][][] u,
      Sampling sg, final float[][][] g)
  {
    final int n1 = u[0][0].length;
    final int n2 = u[0].length;
    final int n3 = u.length;
    final int ng = g[0][0].length;
    final double dg = sg.getDelta();
    final double df = sf.getDelta();
    final double fg = sg.getDelta();
    final double ff = sf.getDelta();
    final float[][][] hf = new float[n3][n2][n1];
    final SincInterp si = new SincInterp();
    int n23 = n3*n2;
    Parallel.loop(n23,new Parallel.LoopInt() {
    public void compute(int i23) {
      int i2 = i23%n2;
      int i3 = i23/n2;
      double v = ff;
      for (int i1=0; i1<n1; i1++, v=ff+i1*df) {
        hf[i3][i2][i1] = si.interpolate(
          ng,dg,fg,g[i3][i2],(float)v+u[i3][i2][i1]);
      }
    }});
    return hf;
  }

  public static float[][] compositeShifts(
      final Sampling sf, final float[][] u1, final float[][] u2)
  {
    int n2 = u1.length;
    final float[][] uc = new float[n2][];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      uc[i2] = compositeShifts(sf,u1[i2],u2[i2]);
    }});
    return uc;
  }

  public static float[] compositeShifts(
      Sampling sf, float[] u1, float[] u2)
  {
    int n1 = u1.length;
    float[] x = getValues(sf);
    return compositeShifts(n1,x,u1,u2);
  }

  public static void compress(float c, float[] f) {
    SincInterp si = new SincInterp();
    int n1 = f.length;
    float[] g = new float[n1];
    scale(si,n1,f,c,g);
  }

  public static void compress(float c, float[][] f) {
    SincInterp si = new SincInterp();
    int n1 = f[0].length;
    int n2 = f.length;
    float[] g = new float[n1];
    for (int i2=0; i2<n2; i2++)
      scale(si,n1,f[i2],c,g);
  }

  public static void stretch(float c, float[] f) {
    SincInterp si = new SincInterp();
    int n1 = f.length;
    float[] g = new float[n1];
    scale(si,n1,f,1.0f/c,g);
  }

  public static void scale(
      SincInterp si, int n1, float[] f, float c, float[]g)
  {
    si.interpolate(n1,1.0,0.0,f,n1,c,0.0,g);
    copy(g,f);
  }

  /**
   * Scale shifts by compression factor c.
   * @param sf the trace sampling corresponding to shifts.
   * @param u array of shifts to compress.
   * @param c scale factor. To stretch, c, to compress, 1/c.
   * @return scaled shifts.
   */
  public static float[] getScaledShifts(Sampling sf, float[] u, float c) {
    int n = u.length;
    float[] uc = new float[n];
    for (int i=0; i<n; i++) {
      float fv = (float)sf.getValue(i);
      uc[i] = (fv+u[i])*c-fv;
    }
    return uc;
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

  /**
   * Returns an extrapolated array with length {@code n1} using the last value
   * of the input array {@code f}.
   * @param n1 the length of the output extrapolated array.
   * @param f the input array.
   * @return a constantly extrapolated array.
   */
  public static float[] extrapolate(int n1, float[] f) {
    int n1f = f.length;
    float v = f[n1f-1];
    float[] ef = new float[n1];
    copy(n1f,f,ef);
    for (int i1=n1f; i1<n1; i1++)
      ef[i1] = v;
    return ef;
  }

  /**
   * Returns an extrapolated array with length {@code n1} using the last value
   * of the input array {@code f}.
   * @param n1 the length of the output extrapolated array.
   * @param f the input array.
   * @return a constantly extrapolated array.
   */
  public static float[][] extrapolate(int n1, float[][] f) {
    int n2 = f.length;
    float[][] ef = new float[n2][];
    for (int i2=0; i2<n2; i2++)
      ef[i2] = extrapolate(n1,f[i2]);
    return ef;
  }

  /**
   * Returns an extrapolated array with length {@code n1} using the last value
   * of the input array {@code f}.
   * @param n1 the length of the output extrapolated array.
   * @param f the input array.
   * @return a constantly extrapolated array.
   */
  public static float[][][] extrapolate(int n1, float[][][] f) {
    int n3 = f.length;
    float[][][] ef = new float[n3][][];
    for (int i3=0; i3<n3; i3++)
      ef[i3] = extrapolate(n1,f[i3]);
    return ef;
  }

  /**
   * Normalizes values to be in range [0,1].
   * @param e input/output array.
   */
  public static void normalizeErrors(float[][] e) {
    int nl = e[0].length;
    int n1 = e.length;
    float emin =  Float.MAX_VALUE;
    float emax = -Float.MAX_VALUE;
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
   * Normalizes values to be in range [0,1].
   * @param e input/output array.
   */
  public static void normalizeErrors(
      float[][] e, float ignoreMin, float ignoreMax)
  {
    int nl = e[0].length;
    int n1 = e.length;
    float emin =  Float.MAX_VALUE;
    float emax = -Float.MAX_VALUE;
    for (int i1=0; i1<n1; ++i1) {
      for (int il=0; il<nl; ++il) {
        float ei = e[i1][il];
        if (ei<emin && ei>ignoreMin) emin = ei;
        if (ei>emax && ei<ignoreMax) emax = ei;
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

  public static void normalizeErrors(
      float[][][] e, final float ignoreMin, final float ignoreMax)
  {
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
          if (ei<emin && ei>ignoreMin) emin = ei;
          if (ei>emax && ei<ignoreMax) emax = ei;
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

  public static void normalizeErrors(
      float[][][][] e, final float ignoreMin, final float ignoreMax)
  {
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
              if (ei<emin && ei>ignoreMin) emin = ei;
              if (ei>emax && ei<ignoreMax) emax = ei;
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

  public static void normalize(
      float[][] f, final float nmin, final float nmax)
  {
    final int n1 = f[0].length;
    final int n2 = f.length;
    final float[][] ff = f;
    final float vmin = min(f);
    final float vmax = max(f);
    final float range = vmax-vmin;
    final float nrange = nmax-nmin;
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int i1=0; i1<n1; ++i1) {
        float vi = ff[i2][i1];
        ff[i2][i1] = nrange*(vi-vmin)/range + nmin;
      }
    }});
  }

  public static void normalize(
      float[][][] f, final float nmin, final float nmax)
  {
    final int n1 = f[0][0].length;
    final int n2 = f[0].length;
    final int n3 = f.length;
    final float[][][] ff = f;
    final float vmin = min(f);
    final float vmax = max(f);
    final float range = vmax-vmin;
    final float nrange = nmax-nmin;
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2) {
        for (int i1=0; i1<n1; ++i1) {
          float vi = ff[i3][i2][i1];
          ff[i3][i2][i1] = nrange*(vi-vmin)/range + nmin;
        }
      }
    }});
  }

  /**
   * Returns errors in an array with lag the slowest dimension.
   * Useful only for visualization of errors. Other methods in this
   * class assume that lag is the fastest dimension in arrays of errors.
   * @param e array[n1][nl] of errors.
   * @return transposed array[nl][n1] of errors.
   */
  @Deprecated
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
  @Deprecated
  public static float[][][] transposeLag(final float[][][] e) {
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final float[][][] t = new float[nl][n2][n1];
    Parallel.loop(nl,new Parallel.LoopInt() {
    public void compute(int il) {
      for (int i2=0; i2<n2; ++i2)
        for (int i1=0; i1<n1; ++i1)
          t[il][i2][i1] = e[i2][i1][il];
    }});
    return t;
  }

  @Deprecated
  public static float[][][] transposeLag12(final float[][][] e) {
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final float[][][] t = new float[n2][nl][n1];
    Parallel.loop(n2,new Parallel.LoopInt() {
    public void compute(int i2) {
      for (int il=0; il<nl; ++il)
        for (int i1=0; i1<n1; ++i1)
          t[i2][il][i1] = e[i2][i1][il];
    }});
    return t;
  }

  @Deprecated
  public static float[][][] transposeLag23(final float[][][] e) {
    final int nl = e[0][0].length;
    final int n1 = e[0].length;
    final int n2 = e.length;
    final float[][][] t = new float[n1][n2][nl];
    Parallel.loop(n1,new Parallel.LoopInt() {
    public void compute(int i1) {
      for (int i2=0; i2<n2; ++i2)
        t[i1][i2] = e[i2][i1];
    }});
    return t;
  }

  @Deprecated
  public static float[][][][] transposeLag12(final float[][][][] e) {
    final int nl = e[0][0][0].length;
    final int n1 = e[0][0].length;
    final int n2 = e[0].length;
    final int n3 = e.length;
    final float[][][][] t = new float[n3][n2][nl][n1];
    Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
      for (int i2=0; i2<n2; ++i2)
        for (int il=0; il<nl; ++il)
          for (int i1=0; i1<n1; ++i1)
            t[i3][i2][il][i1] = e[i3][i2][i1][il];
    }});
    return t;
  }

  ////////////////////////////////////////////////////////////////////////////
  // Private

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

  private static float[] compositeShifts(
      int n1, float[] x, float[] u1, float[] u2)
  {
    float[] uc = new float[n1];
    CubicInterpolator ci = new CubicInterpolator(Method.LINEAR,x,u2);
    for (int i1=0; i1<n1; i1++)
      uc[i1] = u1[i1]+ci.interpolate(x[i1]+u1[i1]);
    return uc;
  }

  private static float[] getValues(Sampling s) {
    int n = s.getCount();
    double f = s.getFirst();
    double d = s.getDelta();
    float[] x = new float[n];
    for (int i=0; i<n; i++)
      x[i] = (float)(f+i*d);
    return x;
  }

  private static float[] firstDerivative(Sampling s, float[] f) {
    int n = f.length;
    int nm1 = n-1;
    float[] g = new float[n];
    float di = 1.0f/(float)s.getDelta();
    float di2 = 0.5f*di;
    g[ 0 ] = (f[ 1 ]-f[  0  ])*di; // forward diff
    g[nm1] = (f[nm1]-f[nm1-1])*di; // backward diff
    for (int i1=1; i1<nm1; i1++)
      g[i1] = (f[i1+1]-f[i1-1])*di2;
    return g;
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
