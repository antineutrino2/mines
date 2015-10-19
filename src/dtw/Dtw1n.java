/**
 * Dynamic Time Warping in 1D.
 * Subsample Dynamic Time Warping used for the application 
 * of seismic time to depth conversion
 *
 * @author Andrew Munoz, Colorado School of Mines
 * @version 10.20.2011
 */

package dtw;

import java.awt.Color;
import java.lang.Math.*;
import java.util.Random; 
import java.util.Random.*; 
import java.util.List;
import java.util.LinkedList;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.util.Check;
import static edu.mines.jtk.util.Parallel.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class Dtw1n {

  private int dir=1;

  /** 
   * Cost matrix for the two signals.
   * The cost matrix is found by calculating the absolute value of the
   * difference of the two signals (Manhattan Distance).
   * <p>
   * @param x[] Signal 1 (x-axis)
   * @param y[] Signal 2 (y-axis)
   * @return cm[][] The cost matrix
   */
  public float[][] cm(float[] x, float[] y) {
    int n1 = x.length;
    int n2 = y.length;
    float[][] cm = new float[n2][n1];
    for   (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        cm[i2][i1] = pow(abs(x[i1] - y[i2]),2f); // Manhattan Distance 
        //cm[i2][i1] = pow(abs((x[i1+1]-x[i1-1])/2 - (y[i2+1]-y[i2-1])/2),2f); // DDTW 
      }
    }
    // for DDTW
    //for (int i2=0; i2<n2; ++i2) 
    //    cm[i2][0] = cm[i2][1];
    //for (int i1=0; i1<n1; ++i1) 
    //    cm[0][i1] = cm[1][i1];
    return cm;
  }

  /** 
   * Accumulated cost matrix for the two signals of different length.
   * Shows all of the possible solutions for signal correlation.
   * <p>
   * @param cm[][] Cost matrix
   * @return acm[][] The accumulated cost matrix
   */
  public float[][] acm(float[][] cm, int sl) {
    int n1 = cm[0].length;
    int n2 = cm.length;
    float[][] acm = new float[n2][n1];
    //int sl = 2; // min=2
    int bc = 2;
    float[] sums1 = new float[bc];
    // Set matrix edges according to Mueller:
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<bc; ++i1) {
        sums1[i1] += cm[i2][i1];
        acm[i2][i1] = sums1[i1]; 
    }}
    for (int i2=0; i2<bc; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        acm[i2][i1] = cm[i2][i1];
    }}
    // Accumulate Values
    for (int i2=bc; i2<n2; ++i2) {
      for (int i1=bc; i1<n1; ++i1) {
        // slope constraint (P) == 2 extended to 
        //   an arbitrary number of samples sl (Sakoe, Chiba 1978)
        acm[i2][i1] = minarb(sl,acm,cm,i2,i1);
      }
    }
    return acm;
  }
  private float minarb(int sl, float[][] acm, float[][] cm, int i2, int i1) {
    int n1 = acm[0].length; 
    int n2 = acm.length;
    if (i2<sl || i1<sl) sl=(i2<i1)?i2:i1;
    float eo = cm[i2][i1];
    float m3 = acm[i2-1][i1-1]+eo;
    float m = m3;
    float[] es1 = new float[sl];
    float[] es2 = new float[sl];
    float es1o = eo;
    float es2o = eo;
    for (int i=0; i<sl-1; ++i) {
      // Get the error vals for each sample point
      es1[i] = es1o = es1o + cm[i2-i][i1-i-1];
      es2[i] = es2o = es2o + cm[i2-i-1][i1-i];
    }
    //for (int i=1; i<sl; ++i) {
    //  float m1 = acm[i2-i  ][i1-i-1] + es1[i-1]; 
    //  float m2 = acm[i2-i-1][i1-i  ] + es2[i-1]; 
    //  //if (m1>m3 && m2>m3 && m>m3) m=m3;
    //  if (m3>m2 && m1>m2 && m>m2) m=m2;
    //  if (m3>m1 && m2>m1 && m>m1) m=m1;
    //}
    float m1 = acm[i2-sl  ][i1-sl+1] + es1[sl-2]; 
    float m2 = acm[i2-sl+1][i1-sl  ] + es2[sl-2]; 
    if (m1>m3 && m2>m3 && m>m3) m=m3;
    if (m3>m2 && m1>m2 && m>m2) m=m2;
    if (m3>m1 && m2>m1 && m>m1) m=m1;
    return m;
  }

  /** 
   * Optimal warp path across the accumulated cost matrix for signals of different length.
   * <p>
   * @param float[][] cm cost matrix
   * @param bintv interval for local minumum search
   * @param bst starting sample for b* (usually best around 10% of n)
   * @param tli v constraint- deepest time that the well log could possibly start
   * @param tui v constraint- shallowest time that the well log could possibly start
   * @return All of the possible paths in two arrays
   */
  public float[][][] paths(float[][] cm, int sl) {
    int n1 = cm[0].length;
    int n2 = cm.length;
    float[][] acm = acm(cm,sl);

    int[] pi = rampint(0,1,n1); 
    int np = pi.length;
    float[][] tpi = new float[np][];  // time indicies for each path 
    float[][] zpi = new float[np][];  // depth indicies for each path 
    float[][] pe = new float[np][];  // error values on each path 
    for (int ip=0; ip<np; ++ip) { 
      int i1 = pi[ip];
      int i2 = n2-1;
      float[][] allp = path(i1,i2,sl,cm,acm,pe[ip],tpi[ip],zpi[ip]); // path data
      pe[ip]  = new float[allp[0].length];
      tpi[ip] = new float[allp[1].length];
      zpi[ip] = new float[allp[2].length];
      copyp(allp[0],pe[ip]);
      copyp(allp[1],tpi[ip]);
      copyp(allp[2],zpi[ip]);
    }
    return new float[][][]{pe,tpi,zpi};
  }
  /*
  * Explain the velocity constraints applied here in the documentation
  */
  public float[][][] pathRemoval(
    float[][] pe, float[][] tpi, float[][] zpi, int tui, int tli)
    {
    int np = pe.length;
    int[] ps = new int[np];
    int pc = 0;
    for (int ip=0; ip<np; ++ip) {
      int nj = tpi[ip].length;
      if (tpi[ip][nj-1]>=tui && tpi[ip][nj-1]<=tli && tpi[ip][nj-1]>2) {
        ps[ip] = 1;
        pc+=1;
      }
    }
    if (pc==0) System.out.println("No paths followed the constraints");
    if (pc>0) {
      return pathcopy(pe,tpi,zpi,ps,pc);
    } 
    System.out.println("Size of np is: "+np);
    return new float[][][]{null,null,null};
  }
  public float[][][] pathRemoval(
    float[][] pe, float[][] tpi, float[][] zpi, float[] b, int r)
    {
    int bst = (int)tpi[0][0];
    int ben = (int)tpi[tpi.length-1][0];
    int[] mp = btest(bst,ben,2,b,pe,r);
    int np = pe.length;
    int[] ps = new int[np];
    for (int rk=0; rk<r; ++rk) {
      ps[mp[rk]] = 1;
    }
    return pathcopy(pe,tpi,zpi,ps,r);
  }

  private float[][][] pathcopy(float[][] pe, float[][] tpi, float[][] zpi, int[] ps, int pc) {
    float[][] pe2  = new float[pc][];  
    float[][] tpi2 = new float[pc][];  
    float[][] zpi2 = new float[pc][];  
    int np = pe.length;
    int ip2 = 0;
    for (int ip=0; ip<np; ++ip) {
      if (ps[ip]!=0) {
        pe2[ip2]  = new float[ pe[ip].length];  
        tpi2[ip2] = new float[tpi[ip].length];  
        zpi2[ip2] = new float[zpi[ip].length];  
  	    copyp( pe[ip], pe2[ip2]);
  	    copyp(tpi[ip],tpi2[ip2]);
  	    copyp(zpi[ip],zpi2[ip2]);
  	    ++ip2;
      }
    }
    return new float[][][]{pe2,tpi2,zpi2};
  }
  
////////////////////////////////////////////////////////////////////////////////////////////////
//                                        PRIVATE                                             //
////////////////////////////////////////////////////////////////////////////////////////////////

  /**
   * Path that returns the image of all paths.
   * Uses Sakoe, Chiba (1978) slope constraint of P=1.
   */
  private float[][] path(
  int i1, int i2, int sl, float[][] cm, float[][] acm, 
  float[] pe, float[] tpi, float[] zpi) 
  {
    int n1n2 = acm.length+acm[0].length;
    float[]  pe2 = new float[n1n2];
    float[] tpi2 = new float[n1n2];
    float[] zpi2 = new float[n1n2];
    int j=0;
    float v=0f;
    int nj=2;
    if (i1>1) {
    while (i1>1 && i2>1) 
    {
      if (i2<sl || i1<sl) sl=(i2<i1)?i2:i1;
       pe2[j] = acm[i2][i1];
      tpi2[j] = i1;
      zpi2[j] = i2;
      int mi = parb(sl,acm,cm,i2,i1);
      if(mi==3) {--i1;--i2;++j;}  
      else if(mi==2) {
        --i2; ++j;
	for (int k=1; k<sl; ++k) {
          pe2[j] = acm[i2][i1];
          tpi2[j] = i1; 
          zpi2[j] = i2;
          --i1; --i2; ++j;
	}
      }
      else if(mi==1) {
        --i1; ++j;
	for (int k=1; k<sl; ++k) {
          pe2[j] = acm[i2][i1];
          tpi2[j] = i1; 
          zpi2[j] = i2;
          --i1; --i2; ++j;
	}
      }
    }   
    nj = j;
    }
     pe = new float[nj];
    tpi = new float[nj];
    zpi = new float[nj];
    copyp(nj, pe2, pe);
    copyp(nj,tpi2,tpi);
    copyp(nj,zpi2,zpi);
    return new float[][]{pe,tpi,zpi};
  }
  private int parb(int sl, float[][] acm, float[][] cm, int i2, int i1) {
    int n1 = acm[0].length; 
    int n2 = acm.length;
    if (i2<sl || i1<sl) sl=(i2<i1)?i2:i1;
    float eo = cm[i2][i1];
    float m3 = acm[i2-1][i1-1]+eo;
    float m = m3;
    int mi = 3;
    float[] es1 = new float[sl];
    float[] es2 = new float[sl];
    float es1o = eo;
    float es2o = eo;
    for (int i=0; i<sl-1; ++i) {
      // Get the error vals for each sample point
      es1[i] = es1o += cm[i2-i  ][i1-i-1];
      es2[i] = es2o += cm[i2-i-1][i1-i  ];
    } 
    //for (int i=1; i<sl; ++i) {
    //  float m1 = acm[i2-i  ][i1-i-1] + es1[i-1]; 
    //  float m2 = acm[i2-i-1][i1-i  ] + es2[i-1]; 
    //if (m1>m3 && m2>m3 && m>m3) {mi=3; m=m3;}
    //  if (m3>m2 && m1>m2 && m>m2) {mi=2; m=m2;}
    //  if (m3>m1 && m2>m1 && m>m1) {mi=1; m=m1;}
    //}
    float m1 = acm[i2-sl+1][i1-sl  ] + es1[sl-2]; 
    float m2 = acm[i2-sl  ][i1-sl+1] + es2[sl-2]; 
    if (m1>m3 && m2>m3 && m>m3) mi=3;
    if (m3>m2 && m1>m2 && m>m2) mi=2;
    if (m3>m1 && m2>m1 && m>m1) mi=1;
    return mi;
  }



  /**
   * Returns the b* indicies by looking at the array b.
   * Add references and theory to the Mueller Paper
   */
  private int[] btest(int bst, int ben, int bintv, float[] b, float[][] pe, int r) {
    int n1 = b.length;
    int[] bt = new int[n1];
    float[] bv = new float[n1];
    float[] bvn = new float[b.length];
    boolean test = false;
    int bct=0;
    // Normalize b 
    for (int i=bst; i<ben; ++i) 
      bvn[i] = b[i]/(pe[i-bst].length);
    // Find all b values within the local minimum 
    for (int i1=bst; i1<ben; ++i1) {
    //Check to make sure it's a Local Minumum around an interval
      localmin:
        for (int i2=1; i2<bintv; ++i2) {
          if ((i1+i2) < n1 && (i1-i2) > 0) { 
            if (bvn[i1] < bvn[i1+i2] && bvn[i1] < bvn[i1-i2]) {
              test = true;
            }
            else {
              test = false;
      	      break localmin;
            }
          }
          else if ((i1+i2) < n1 && (i1-i2) < 0 ) {
            if (bvn[i1] < bvn[i1+i2] && bvn[i1] < bvn[1]) {
              test = true;
            }
            else {
              test = false;
      	      break localmin;
            }
          } 
          else if ((i1+i2) > n1 && (i1-i2) > 0 ) {
            if (bvn[i1] < bvn[n1-1] && bvn[i1] < bvn[i1-i2]) {
              test = true;
            }
            else {
              test = false;
      	      break localmin;
            }
          }
        }
      if (test == true) {
        bt[i1] = i1; // Index of b*
       bv[i1] = bvn[i1]; // Value of b* normalized
       ++bct;
       test = false;
      }
    }
    if (r>bct) r=bct;
    int[] bt2 = new int[r];
    for (int rk=0; rk<r; ++rk) {
      float bvm=max(bv);
      for (int i=bst,j=0; i<ben; ++i,++j) {
        if (bv[i]>0 && bv[i]<bvm) {
          bt2[rk] = j;
          bvm=bv[i];
          bv[i]=0;
        }
      }
    }
    return bt2;
  }

  /**
   * Path that returns the image of all paths.
   * Uses Sakoe, Chiba (1978) slope constraint of P=0.
   */
  private float[][] patho(
  int i1, int i2, float[][] cm, float[][] acm, 
  float[] pe, float[] tpi, float[] zpi) 
  {
    int n1n2 = acm.length+acm[0].length;
    float[]  pe2 = new float[n1n2];
    float[] tpi2 = new float[n1n2];
    float[] zpi2 = new float[n1n2];
    int j=0;
    float v=0f;
    while (i1>1 && i2>1) 
    {
      v = acm[i2][i1];
       pe2[j] = v;
      tpi2[j] = i1;
      zpi2[j] = i2;
      float stp1 = acm[i2-1][i1-1];
      float stp2 = acm[i2-1][i1  ];
      float stp3 = acm[i2  ][i1-1];
      float stpm = min(stp1,stp2,stp3);
      if(stpm==stp1) {--i1;--i2;}  
      else if(stpm==stp2) --i2;
      else if(stpm==stp3) --i1;
      ++j;
    }   
    int nj = j;
     pe = new float[nj];
    tpi = new float[nj];
    zpi = new float[nj];
    copyp(nj, pe2, pe);
    copyp(nj,tpi2,tpi);
    copyp(nj,zpi2,zpi);
    return new float[][]{pe,tpi,zpi};
  }

  private static void copyp(final float[] x, final float[] y) {
    final int n1 = x.length;
    loop(n1,new LoopInt() {   
    public void compute(int i1) {
       y[i1] = x[i1];
    }});
  }
  private static void copyp(final int n1, final float[] x, final float[] y) {
    loop(n1,new LoopInt() {   
    public void compute(int i1) {
       y[i1] = x[i1];
    }});
  }

//////////////////////////////////////////////////////////////////// 
//// Thanks to Luming Liang for the following methods           //// 
////////////////////////////////////////////////////////////////////

// Methods are public so I can call them from outside python code

  public float[] makeEvents(int n1, long seed) {
    Random r = new Random(seed);
    float[] f = pow(mul(2.0f,sub(randfloat(r,n1),0.5f)),7.0f);
    return f;
  }

  public float[] addRickerWavelet(double fpeak, float[] f) {
    int n1 = f.length;
    int ih = (int)(3.0/fpeak);
    int nh = 1+2*ih;
    float[] h = new float[nh];
    for (int jh=0; jh<nh; ++jh)
      h[jh] = ricker(fpeak,jh-ih);
    float[] g = new float[n1];
    Conv.conv(nh,-ih,h,n1,0,f,n1,0,g);
    return g;
  }

  private float ricker(double fpeak, double time) {
    double x = PI*fpeak*time;
    double xx = x*x;
    return (float)((1.0-2.0*xx)*exp(-xx));
  }

  public float[] addNoise(float nrms, long seed, float[] f) {
    int n = f.length;
    Random r = new Random(seed);
    nrms *= max(abs(f));
    float[] g = mul(2.0f,sub(randfloat(r,n),0.5f));
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.apply1(g,g); // 1st derivative enhances high-frequencies
    float frms = sqrt(sum(mul(f,f))/n);
    float grms = sqrt(sum(mul(g,g))/n);
    g = mul(g,nrms*frms/grms);
    return add(f,g);
  }

};

//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////
// Old Functions not being used anymore

  /*
   * To find the warp path starting from both axes. 
   * Finds the warp path from b1 and b2
   */
  /*
  public float[][] wpath2(int bintv, int bst, float[][] cm) {
    float[][] wpath1 = new float[cm.length][cm[0].length];
    float[][] wpath2 = new float[cm[0].length][cm.length];
    float[][] wpath3 = new float[cm.length][cm[0].length];
    float[][] cm2 = new float[cm[0].length][cm.length];
    wpath1 = wpath1(bintv,bst,cm);
    for (int i2=0; i2<cm.length; ++i2) {
      for (int i1=0; i1<cm[0].length; ++i1) {    // <<<< needs speedup
        cm2[i1][i2] = cm[i2][i1];
      }
    }
    wpath2 = wpath1(bintv,bst,cm2);
    for (int i2=0; i2<cm.length; ++i2) {
      for (int i1=0; i1<cm[0].length; ++i1) {
        wpath3[i2][i1] = wpath2[i1][i2] + wpath1[i2][i1];
      }
    }
    return wpath3;
  }*/


	  /*else {
	    // FIX
	    // slope constraint > 1
            float[][] acmst = new float[sst][sst];
	    for (int j2=0; j2<sst; ++j2) {
	      for (int j1=0; j1<sst; ++j1) {
	        if ((i1-j1+1)>0 && (i2-j2+1)>0 && (i1-j1+2)<n1 && (i2-j2+2)<n2) {
                  acmst[j2][j1] = acm[i2-j2+1][i1-j1+1];
	          acmst[j2][j1] += cm[i2-j2+2][i1-j1+2];
	        }
	      }
	    }
            acm[i2][i1] = min(acmst) + cm[i2][i1];
	  }*/

  ////////////////////////////////////////////////////////
    // Plot b values
    /*SimplePlot plot = new SimplePlot();
    plot.addPoints(b);
    plot.setSize(1000,500);
    plot.setTitle("b function");*/
  ////////////////////////////////////////////////////////
    // Find the global minumum index of b first
    /*float bglob = b[1];
    int bglobi  = 1;
    for (int i1=2; i1<n1-1; ++i1) {
      if (b[i1] <= bglob) {
        bglobi = i1;
	bglob = b[i1];
      }
    }
    System.out.println("The Global b index is: "+bglobi);
    System.out.println("The b interval size is "+bintv+" samples");*/

    // Find the minimum path which has the shortest accumulated distance
    /*float mpath = path[0];
    float[][] minpath = new float[n2][n1];
    for (int i=0; i<bsize; ++i) {
      if ((path[i]) <= mpath) {
        mindex = i;
	mpath = path[i];
      }
      System.out.println("path "+i+": "+path[i]+" of b index: "+bi[i]);
    }
    path(bi[mindex],n2-1,mindex,acm,minpath);
    //System.out.println("min b* index is: "+mindex);
   
   // For background of image to make it more visible
    /*float acmmax = max(acm);
    for (int i2=0; i2<n2; ++i2) {
      for (int i1=0; i1<n1; ++i1) {
        if (wpath[i2][i1] == 0f) {
	  wpath[i2][i1] = acmmax + 100f;
	}
      }
    }*/
    
    //return minpath;

/*
  public float[][][] interpolations(float[][] t, float[][] d, Sampling st, Sampling sd) 
    {
    int m = 2; // slope
    int nb = t.length;
    int n = t[0].length;
    float t1; float d1;
    float t1p1; float d1p1;
    int t1i; int d1i;
    int t1ip1; int d1ip1;
    int[] blens = new int[nb];
    float[][] tl = new float[nb][2*n];
    float[][] dl = new float[nb][2*n];
    for (int b=0; b<nb; ++b) {
      int i1=0; int i2=0; 
      float std = (float)st.getDelta();
      float sdd = (float)sd.getDelta();
      float sdf = (float)sd.getFirst();
      tl[b][0] = ((int)(t[b][i1]))*std; 
      dl[b][0] = ((int)(d[b][i2]))*sdd+sdf;
      do {
        t1i = (int)(t[b][i1]);
        d1i = (int)(d[b][i2]);
        t1ip1 = (int)(t[b][i1+1]);
        d1ip1 = (int)(d[b][i2+1]);
        t1 = t1i*std;
        d1 = sdf+d1i*sdd;
        t1p1 = t1ip1*std;
        d1p1 = sdf+d1ip1*sdd;
        if (t1i!=t1ip1+1) { // slope = 1/2 
          tl[b][i1] = t1-(t1-t1p1)/m; 
          tl[b][i1+1] = t1p1; 
          dl[b][i2] = d1-(d1-d1p1)/m;
	  i1+=m; ++i2;
        }
        else if (d1i!=d1ip1+1) { // slope = 2
          tl[b][i1] = t1-(t1-t1p1)/m; 
          dl[b][i2] = d1-(d1-d1p1)/m;
          dl[b][i2+1] = d1p1;
	  i2+=m; ++i1;
	}
	else {
          tl[b][i1] = t1p1; 
          dl[b][i2] = d1p1;
        ++i1; ++i2;
	}
      } while (t1p1>0.0 || d1p1>sdf);
      blens[b] = max(i1,i2);
    }
    int n2 = max(blens);
    float[][] tl2 = new float[nb][n2];
    float[][] dl2 = new float[nb][n2];
    for (int b=0; b<nb; ++b) {
      copyp(n2,tl[b],tl2[b]);
      copyp(n2,dl[b],dl2[b]);
    }
    return new float[][][]{tl,dl};
  }
*/

/*
    int[] brt = new int[np];
    int bct = 0;
    int t1 = n1+n2;
    for (int b=0; b<np; ++b) {
      // Find where time=1 (2 for 2ms sampling)
      for (int j=0; j<n1+n2; ++j) {
        if (times[b][j]<=2.0) {
	  t1 = j;
	  break;
	}
      }
      iloop:
      for (int i1=0; i1<n1+n2; ++i1) {
        int ti = (int)times[b][i1];
        for (int j=i1+1; j<t1-1; ++j) {
          if (times[b][j]==times[b][j+1]+1) { // Check if at the end of the path using only time indicies
            continue;
          }
	  else continue iloop;
        }
	if (ti<=tli && ti>=tui) { // Constraint for paths below v=1.5 km/s limit and above upper limit
	  for (int i2=2; i2<n1+n2-3; ++i2) {
	    int ti2 = (int)times[b][i2  ];
            int tp1 = (int)times[b][i2+1];
            int tp2 = (int)times[b][i2+2];
	    if (ti2==tp1 && ti2==tp2 && ti2>0) { // remove paths below time 0
              brt[b] = b+1;
              bct += 1;
              break iloop;
            }
	  }
	  break iloop; // If the curve passes, test moves on to the next b 
	}
        else {  // If the curve fails the velocity test then it is put on the remove list below
          brt[b] = b+1;
          bct += 1;
          break iloop;
	}
      }
    }
    
  private float ysolve(float yo, float y1, float x, float xo, float x1) {
    float y = yo + (x-xo)*((y1-yo)/(x1-xo));
    return y;
  }
*/



