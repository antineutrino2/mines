package tp;

public class LongList {
  public int n;
  public long[] a = new long[1024];
  public void add(long i) {
    if (n==a.length) {
      long[] t = new long[2*n];
      System.arraycopy(a,0,t,0,n);
      a = t;
    }
    a[n++] = i;
  }
  public long[] trim() {
    if (n==0)
      return null;
    long[] t = new long[n];
    System.arraycopy(a,0,t,0,n);
    return t;
  }
  public boolean contains(long x) {
    for (int i=0; i<n; ++i) {
      if (a[i]==x)
        return true;
    }
    return false;
  }
}
