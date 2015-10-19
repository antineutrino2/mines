package sem;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class Subsample {

  /**
   * Subsamples indices of an array of length n, ensuring that
   * the subsampling includes the first and last index. To meet
   * this requirement, the interval of the sparse indices may
   * be greater than or equal to {@code d},
   * @param n the length of the array to subsample
   * @param d the minimum subsampled interval
   * @return an integer array of sparse grid indices.
   */
  public static int[] subsample(int n, int d) {
    if (d>=n)
      d = max(n-1,1);
    int m = 1+(n-1)/d;
    double dd = (double)(n-1)/(double)(m-1);
    int[] g = new int[m];
    for (int ig=0; ig<m; ++ig)
      g[ig] = (int)(ig*dd+0.5);
    return g;
  }

  /**
   * Subsamples indices of an array of floats preferentially
   * selecting the indices of f with the highest absolute value
   * amplitudes. This selection ensures that the first and
   * last sample are always included in the subsampled indices,
   * and that the interval between selected indices is greater
   * than or equal to d.
   * @param f the input array to preferentially subsample.
   * @param d the constraint on the subsample interval.
   * @return an integer array of sparse grid indices. The
   *  length of this array is unknown. It will include as
   *  many indices as possible that satisfy the interval
   *  constraint d.
   */
  public static int[] subsample(float[] f, int d) {
    return subsample(f,d,-1);
  }

  /**
   * Subsamples indices of an array of floats preferentially
   * selecting the indices of f with the highest absolute value
   * amplitudes. This selection ensures that the first and
   * last sample are always included in the subsampled indices,
   * and that the number of subsamples is equal to {@code ng}. The
   * interval between selected indices decreased from {@code d } as
   * needed to meet this requirement.
   * @param f the input array to preferentially subsample.
   * @param d the constraint on the subsample interval.
   * @param ng the number of samples included in the subsampled array.
   * @return an integer array of sparse grid indices. The
   *  length of this array is {@code ng} and the interval {@code d}
   *  may be modified to meet this requirement.
   */
  public static int[] subsample(float[] f, int d, int ng) {
    return subsample(f.length,f,d,ng);
  }

  public static int[] subsample(int nf, float[] f, int d, int ng) {
    int im = nf-1; // max index
    int[] i = rampint(0,1,nf);
    float[] fa = copy(nf,abs(f));
    quickIndexSort(fa,i);
   // StringBuffer sb = new StringBuffer();
   // for (int ii=0; ii<10; ii++)
   //   sb.append(i[im-ii]+",");
   // System.out.println("Top 10: "+sb.toString());

    // Get maximum amplitude indices, adjusting the interval
    // d so that the subsampled array contains ng samples.
    List<Integer> gl = null;
    if (ng!=-1) {
      for (; d>=1; d--) {
        gl = getList(i,nf,im,d,ng);
        if (gl.size()==ng) { // size matches desired sparse array size, quit
          System.out.println("final d="+d);
          break;
        }
      }
    } else {
      gl = getList(i,nf,im,d,ng);
    }
    Collections.sort(gl);
    int nl = gl.size();
    int[] g = new int[nl];
    for (int il=0; il<nl; il++) {
      g[il] = gl.get(il);
    }
    return g;
  }

  private static List<Integer> getList(int[] i, int nf, int im, int d, int ng) {
    List<Integer> gl = new ArrayList<>();
    gl.add(0); gl.add(im);
    for (int ii=0; ii<nf; ii++) {
      int ia = i[im-ii]; // next highest amplitude index
      int s = gl.size();
      if (s==ng) // size matches desired sparse array size, quit
        break;
      boolean add = true;
      for (int ig=0; ig<s; ig++) {
        if (abs(ia-gl.get(ig))<d) {
          add = false;
          break;
        }
      }
      if (add)
        gl.add(ia);
    }
    return gl;
  }

  public static float[] indicesToSampling(Sampling s, int[] ig) {
    int n = ig.length;
    float[] sg = new float[n];
    for (int i=0; i<n; i++)
      sg[i] = (float)s.getValue(ig[i]);
    return sg;
  }

  public static int[] subsampleEnvelope(float[][][] f, int d) {
    return subsampleEnvelope(f,d,-1);
  }

  public static int[] subsampleEnvelope(int nf, float[][][] f, int d) {
    return subsampleEnvelope(nf,f,d,-1);
  }

  public static int[] subsampleEnvelope(float[][][] f, int d, int ng) {
    return subsampleEnvelope(f[0][0].length,f,d,ng);
  }

  public static int[] subsampleEnvelope(int nf, float[][][] f, int d, int ng) {
    int n1 = f[0][0].length;
    int n2 = f[0].length;
    int n3 = f.length;
    float[] es = new float[n1];
    for (int i3=0; i3<n3; i3++) {
      for (int i2=0; i2<n2; i2++) {
        float[] g = hilbertTransform(n1,f[i3][i2]);
        envelopeSum(n1,f[i3][i2],g,es);
      }
    }
    return subsample(nf,es,d,ng);
  }

  public static int[] subsampleEnvelope(float[][] f, int d) {
    return subsampleEnvelope(f,d,-1);
  }

  public static int[] subsampleEnvelope(int nf, float[][] f, int d) {
    return subsampleEnvelope(nf,f,d,-1);
  }

  public static int[] subsampleEnvelope(float[][] f, int d, int ng) {
    return subsampleEnvelope(f[0].length,f,d,ng);
  }

  public static int[] subsampleEnvelope(int nf, float[][] f, int d, int ng) {
    int n1 = f[0].length;
    int n2 = f.length;
    float[] es = new float[n1];
    for (int i2=0; i2<n2; i2++) {
      float[] g = hilbertTransform(n1,f[i2]);
      envelopeSum(n1,f[i2],g,es);
    }
    return subsample(nf,es,d,ng);
  }

  public static int[] subsampleEnvelope(float[] f, int d) {
    return subsampleEnvelope(f,d,-1);
  }

  public static int[] subsampleEnvelope(int nf, float[] f, int d) {
    return subsampleEnvelope(nf,f,d,-1);
  }

  public static int[] subsampleEnvelope(float[] f, int d, int ng) {
    return subsampleEnvelope(f.length,f,d,ng);
  }

  public static int[] subsampleEnvelope(int nf, float[] f, int d, int ng) {
    int n1 = f.length;
    float[] es = new float[n1];
    float[] g = hilbertTransform(n1,f);
    envelopeSum(n1,f,g,es);
    return subsample(nf,es,d,ng);
  }

  private static void envelopeSum(int n, float[] f, float[] g, float[] es) {
    for (int i=0; i<n; i++)
      es[i] += sqrt(f[i]*f[i]+g[i]*g[i]);
  }

  private static float[] hilbertTransform(int n, float[] f) {
    float[] g = copy(f);
    _htf.apply(n,f,g);
    return g;
  }

  public static void main(String[] args) {
    if (args.length != 2) {
      System.out.println("Usage: java Decompose numberOfSamples delta");
      System.exit(0);
    }
    int n = Integer.parseInt(args[0]);
    int d = Integer.parseInt(args[1]);
    int[] g = subsample(n,d);
    dump(g);
  }

  private static HilbertTransformFilter _htf = new HilbertTransformFilter();

}
