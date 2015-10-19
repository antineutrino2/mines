package sem;

import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.interp.*;

/**
 * Interpolation for warping.
 * @author Stefan Compton
*/

public class ShiftInterp {

  /**
   * Constants for interpolation.
   */
  public enum Method {
    LINEAR,
    MONOTONIC,
    SPLINE
  }

  /**
   * Default constructor.
   */
  public ShiftInterp() {}

  /**
   * Set the interpolation method for this instance.
   */
  public void setMethod(Method method) {
    _method = method;
  }

  /**
   * Interpolates subsampled shifts {@code u[ng]} to uniformly sampled shifts
   * {@code ui[n]}, where {@code n} is the length of the given sampling
   * {@code s}.
   * @param s Sampling for interpolated shifts.
   * @param g sparse grid coordinates.
   * @param u sparse shifts.
   * @return the interpolated shifts.
   */
  public float[] interpolate(Sampling s, float[] g, float[] u) {
    int ng = g.length;
    int n = s.getCount();
    float[] ui = new float[n];
    CubicInterpolator ci = new CubicInterpolator(getMethod(),ng,g,u);
    double d = s.getDelta();
    double f = s.getFirst();
    double v = f;
    for (int i=0; i<n; i++, v=f+i*d)
      ui[i] = ci.interpolate((float)v);
    return ui;
  }

  /**
   * Interpolates subsampled shifts u[ng2][ng1] to uniformly sampled shifts
   * ui[n2][n1], where {@code n2} and {@code n1} are the lengths of the given
   * samplings {@code s2} and {@code s1}. Interpolation is done in two passes.
   * First interpolation is done along the second dimension where subsampling
   * must be regular. The second pass does interpolation along the first
   * dimension where subsampling may be irregular.
   * </p>
   * The locations of the subsampled shifts {@code u} are defined by the input
   * arrays {@code g1} and {@code g2}. The {@code g2} coordinates are assumed to
   * be consistent for all n1 locations. That is, g2 coordinates are the same at
   * every i1, for i1=0,...,n1-1.
   * </p>
   * The {@code g1} coordinate array defines the subsampled locations of the
   * fastest dimension of shifts u, which may vary for each i2, as long as the
   * number of subsampled coordinates is the same for all i2.
   * @param s1 first dimension sampling for interpolated shifts.
   * @param s2 second dimension sampling for interpolated shifts.
   * @param g1 first dimension sparse grid coordinates.
   * @param g2 second dimension sparse grid coordinates.
   * @param u sparse shifts.
   * @return the interpolated shifts ui[n2][n1].
   */
  public float[][] interpolate(
      Sampling s1, Sampling s2, float[][] g1, float[] g2, float[][] u)
  {
    CubicInterpolator.Method m2 = CubicInterpolator.Method.LINEAR;
    CubicInterpolator.Method m1 = getMethod();
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    int ng1 = g1[0].length;
    int ng2 = g2.length;
    double d1 = s1.getDelta();
    double d2 = s2.getDelta();
    double f1 = s1.getFirst();
    double f2 = s2.getFirst();
    float[][] ui = new float[n2][n1];

    // Interpolate along the second dimension.
    float[] u2 = new float[ng2];
    float[][] ui2 = new float[n2][ng1];
    for (int i1=0; i1<ng1; i1++) {
      for (int i2=0; i2<ng2; i2++)
        u2[i2] = u[i2][i1];
      CubicInterpolator ciu = new CubicInterpolator(m2,g2,u2);
      double v2 = f2;
      for (int i2=0; i2<n2; i2++, v2=f2+i2*d2)
        ui2[i2][i1] = ciu.interpolate((float)v2);
    }

    // Interpolate along the first dimension.
    for (int i2=0; i2<n2; i2++) {
      CubicInterpolator ci = new CubicInterpolator(m1,g1[i2],ui2[i2]);
      double v1 = f1;
      for (int i1=0; i1<n1; i1++, v1=f1+i1*d1)
        ui[i2][i1] = ci.interpolate((float)v1);
    }
    return ui;
  }

  /**
   * Interpolates subsampled shifts u[ng3][ng2][ng1] to uniformly sampled
   * shifts ui[_n3][_n2][_ne1]. The locations of the subsampled shifts u are
   * defined by the input arrays g1, g2, and g3. The indices in the g3 array
   * must be the same at every i2, for i2=0,...,n2-1 and i1, for i1=0,...n1-1.
   * The indices in the g2 array must be the same at every i3, for
   * i3=0,...,n3-1, and i1, for i1=0,...,n1-1.
   * </p>
   * The g1 coordinate array defines the subsampled locations of the
   * fastest dimension of shifts u, for all n3 and n2 indices.
   * </p>
   * The interpolation is done in two passes. First interpolation
   * is done in the second and third dimension where subsampling
   * must be regular. The second pass does interpolation in the
   * first dimension where subsampling may be irregular.
   * @param g1 first dimension sparse grid coordinates.
   * @param g2 second dimension sparse grid indices.
   * @param g3 second dimension sparse grid indices.
   * @param u sparse shifts.
   * @param interp1 interpolation method for the i1 (fast) dimension.
   * @param interp23 interpolation method for the i2 (middle) and i3
   *  (slow) dimensions.
   * @return the interpolated shifts ui[_n3][_n2][_ne1].
   */
  public float[][][] interpolate(
      Sampling s1, Sampling s2, Sampling s3,
      float[][][] g1, float[] g2, float[] g3, float[][][] u)
  {
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    int n3 = s3.getCount();
    int ng1 = g1[0][0].length;
    int ng2 = g2.length;
    int ng3 = g3.length;
    double d1 = s1.getDelta();
    double d2 = s2.getDelta();
    double d3 = s3.getDelta();
    double f1 = s1.getFirst();
    double f2 = s2.getFirst();
    double f3 = s3.getFirst();

    CubicInterpolator.Method m1 = getMethod();
    // Interpolate along the second and third dimension.
    float[][] u23 = new float[ng3][ng2];
    float[][][] ui23 = new float[n3][n2][ng1];
    for (int i1=0; i1<ng1; i1++) {
      for (int i3=0; i3<ng3; i3++)
        for (int i2=0; i2<ng2; i2++)
          u23[i3][i2] = u[i3][i2][i1];
      BilinearInterpolator2 bli = new BilinearInterpolator2(g2,g3,u23);
      double v3 = f3;
      for (int i3=0; i3<n3; i3++, v3=f3+i3*d3) {
        float v3f = (float)v3;
        double v2 = f2;
        for (int i2=0; i2<n2; i2++, v2=f2+i2*d2)
          ui23[i3][i2][i1] = bli.interpolate((float)v2,v3f);
      }
    }

    // Interpolate along the first dimension.
    float[][][] ui = new float[n3][n2][n1];
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        CubicInterpolator ci =
            new CubicInterpolator(m1,g1[i3][i2],ui23[i3][i2]);
        double v1 = f1;
        for (int i1=0; i1<n1; i1++, v1=f1+i1*d1)
          ui[i3][i2][i1] = ci.interpolate((float)v1);
      }
    }
    return ui;
  }

  ////////////////////////////////////////////////////////////////////////////
  // Private

  private Method _method;

  private CubicInterpolator.Method getMethod() {
    CubicInterpolator.Method m;
    switch (_method) {
      case LINEAR:    m = CubicInterpolator.Method.LINEAR; break;
      case MONOTONIC: m = CubicInterpolator.Method.MONOTONIC; break;
      case SPLINE:    m = CubicInterpolator.Method.SPLINE; break;
      default: throw new IllegalArgumentException(
        _method.toString()+" is not a recognized interpolation method.");
    }
    return m;
  }

}
