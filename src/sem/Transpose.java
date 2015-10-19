package sem;

import edu.mines.jtk.util.Parallel;

public class Transpose {

  public static final int block1 = 16;
  public static final int block2 = 16;
  public static final int block3 = 16;
  public static final int block4 = 16;
  public static final int block5 = 16;

  public static void transpose4312(final float[][][][] x, final float[][][][] y) {
    final int n1 = x[0][0][0].length; 
    final int n2 = x[0][0].length;    
    final int n3 = x[0].length;
    final int n4 = x.length;
    Parallel.loop(0,n4,block4, new Parallel.LoopInt() {
      public void compute(int i8) {
        for (int i7=0; i7<n3; i7+=block3) {
        for (int i6=0; i6<n2; i6+=block2) {
        for (int i5=0; i5<n1; i5+=block1) {
          for (int i4=i8; i4<i8+block4 && i4<n4; ++i4) {
          for (int i3=i7; i3<i7+block3 && i3<n3; ++i3) {
          for (int i2=i6; i2<i6+block2 && i2<n2; ++i2) {
          for (int i1=i5; i1<i5+block1 && i1<n1; ++i1) {
                  y[i4][i3][i1][i2] = x[i4][i3][i2][i1];
          }
          }
          }
          }
	  		}
        }
        }
    }});
  }

  public static void transpose312(final float[][][] x, final float[][][] y) {
    final int n1 = x[0][0].length; 
    final int n2 = x[0].length;    
    final int n3 = x.length;
    Parallel.loop(0,n3,block3, new Parallel.LoopInt() {
      public void compute(int i6) {
        for (int i5=0; i5<n2; i5+=block2) {
        for (int i4=0; i4<n1; i4+=block1) {
          for (int i3=i6; i3<i6+block3 && i3<n3; ++i3) {
          for (int i2=i5; i2<i5+block2 && i2<n2; ++i2) {
          for (int i1=i4; i1<i4+block1 && i1<n1; ++i1) {
                  y[i3][i1][i2] = x[i3][i2][i1];
          }
          }
          }
	  		}
        }
      }
    });
  }

  public static void transpose132(final float[][][] x, final float[][][] y) {
    final int n1 = x[0][0].length; 
    final int n2 = x[0].length;    
    final int n3 = x.length;
    Parallel.loop(0,n3,block3, new Parallel.LoopInt() {
      public void compute(int i6) {
        for (int i5=0; i5<n2; i5+=block2) {
        for (int i4=0; i4<n1; i4+=block1) {
          for (int i3=i6; i3<i6+block3 && i3<n3; ++i3) {
          for (int i2=i5; i2<i5+block2 && i2<n2; ++i2) {
          for (int i1=i4; i1<i4+block1 && i1<n1; ++i1) {
                  y[i1][i3][i2] = x[i3][i2][i1];
          }
          }
          }
	  		}
        }
      }
    });
  }

  public static void transpose123(final float[][][] x, final float[][][] y) {
    final int n1 = x[0][0].length; 
    final int n2 = x[0].length;    
    final int n3 = x.length;
    Parallel.loop(0,n3,block3, new Parallel.LoopInt() {
      public void compute(int i6) {
        for (int i5=0; i5<n2; i5+=block2) {
        for (int i4=0; i4<n1; i4+=block1) {
          for (int i3=i6; i3<i6+block3 && i3<n3; ++i3) {
          for (int i2=i5; i2<i5+block2 && i2<n2; ++i2) {
          for (int i1=i4; i1<i4+block1 && i1<n1; ++i1) {
                  y[i1][i2][i3] = x[i3][i2][i1];
          }
          }
          }
	  		}
        }
      }
    });
  }

  public static void transpose213(final float[][][] x, final float[][][] y) {
    final int n1 = x[0][0].length; 
    final int n2 = x[0].length;    
    final int n3 = x.length;
    Parallel.loop(0,n3,block3, new Parallel.LoopInt() {
      public void compute(int i6) {
        for (int i5=0; i5<n2; i5+=block2) {
        for (int i4=0; i4<n1; i4+=block1) {
          for (int i3=i6; i3<i6+block3 && i3<n3; ++i3) {
          for (int i2=i5; i2<i5+block2 && i2<n2; ++i2) {
          for (int i1=i4; i1<i4+block1 && i1<n1; ++i1) {
                  y[i2][i1][i3] = x[i3][i2][i1];
          }
          }
          }
	  		}
        }
      }
    });
  }

  public static void transpose231(final float[][][] x, final float[][][] y) {
    final int n1 = x[0][0].length; 
    final int n2 = x[0].length;    
    final int n3 = x.length;
    Parallel.loop(0,n3,block3, new Parallel.LoopInt() {
      public void compute(int i6) {
        for (int i5=0; i5<n2; i5+=block2) {
        for (int i4=0; i4<n1; i4+=block1) {
          for (int i3=i6; i3<i6+block3 && i3<n3; ++i3) {
          for (int i2=i5; i2<i5+block2 && i2<n2; ++i2) {
          for (int i1=i4; i1<i4+block1 && i1<n1; ++i1) {
                  y[i2][i3][i1] = x[i3][i2][i1];
          }
          }
          }
	  		}
        }
      }
    });
  }

	public static void transpose12(final float[][] x, final float[][] y) {
    final int n1 = x[0].length; 
    final int n2 = x.length;    
    Parallel.loop(0,n2,block2, new Parallel.LoopInt() {
      public void compute(int i4) {
      for (int i3=0; i3<n1; i3+=block1) {
        for (int i2=i4; i2<i4+block2 && i2<n2; ++i2) {
        for (int i1=i3; i1<i3+block1 && i1<n1; ++i1) {
                y[i1][i2] = x[i2][i1];
        }
        }
	  	}
      }
    });
  }

}
