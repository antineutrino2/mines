package coda;


public class FloatList {
  public int n;
  public float[] a = new float[1024];
  public void add(float f) {
    if (n==a.length) {
      float[] t = new float[2*n];
      System.arraycopy(a,0,t,0,n);
      a = t;
    }
    a[n++] = f;
  }
  public void add(float[] f) {
    int n = f.length;
    for (int i=0; i<n; ++i) 
      add(f[i]);
  }
  public float[] trim() {
    if (n==0)
      return null;
    float[] t = new float[n];
    System.arraycopy(a,0,t,0,n);
    return t;
  }
  public float[] trim(float tau) {
    if (n==0)
      return null;
    float[] t = new float[n];
    int j=0;
    for (int i=0; i<n; ++i) {
      if (a[i]>tau) {
        t[j] = a[i];
        ++j;
      }
    }
    float[] g = new float[j];
    for (int i=0; i<j; ++i) 
      g[i] = t[i];
    return g;
  }

}
