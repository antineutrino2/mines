package dtw;

import static edu.mines.jtk.util.ArrayMath.*;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

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
  public static int[] subsample(int n, int d, float[] f) {
	if (f.length>1) {
		return subsample(d,f);
		//return subsample(d,f,subsample(n,d).length);
	}
		return subsample(n,d); 
	}
 
  public static int[] subsample(int d, float[] f) {
    int n = f.length;
    int nm = n-1;
    int[] i = rampint(0,1,n);
    float[] fa = abs(sub(1.0f,f));
    quickIndexSort(fa,i);
    
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
  public static int[] subsample(int d, float[] f, int ng) {
    int n = f.length;
    int nm = n-1;
    int[] i = rampint(0,1,n);
    float[] fa = abs(sub(1.0f,f));
    quickIndexSort(fa,i);
    
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


///////////////////////////////////////////////////////////////////////
// For max amps

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

  public static Map<Integer,float[][]> getMXB(
      int[] g, float rmin, float rmax, boolean up) 
  {
    Map<Integer,float[][]> mxbMap= new HashMap<Integer,float[][]>();
    int ks = up? 1:-1;
    int ms = up?-1: 1;
    int ng = g.length;
    for (int ig=1; ig<ng; ig++) {
      int dg = g[ig]-g[ig-1];
      if (mxbMap.containsKey(dg))
        continue;
      int kmin = (int)ceil( rmin*dg);
      int kmax = (int)floor(rmax*dg);
      int nk = kmax+1;
//      System.out.println("dg="+dg+", kmin="+kmin+", kmax="+kmax+", nk="+nk);
      float[][] mx = new float[nk][dg+1];
      for (int k=kmin; k<=kmax; k++) {
        float m = (float)k/dg*ms;
        for (int x=0; x<=dg; x++) {
          mx[k][x] = m*x+ks*k;
        }  
      }
      mxbMap.put(dg,mx);
    }
    return mxbMap;
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

}
