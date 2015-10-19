/** 
 * Computes semblance using various NMO equations. 
 * Thanks to Simon Luo for helping me start this code.
 * @author Andrew Munoz, CSM
 * @version 10.26.13
 */

package sem;

import java.util.Random;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class SemblanceOld {

	public SemblanceOld(
    Sampling st, Sampling sx, Sampling sv)
  {
		this(st,sx,sv,null,0.0);
  }

	public SemblanceOld(
    Sampling st, Sampling sx, Sampling sv, double fpeak)
  {
		this(st,sx,sv,null,fpeak);
  }

	public SemblanceOld(
    Sampling st, Sampling sx, Sampling sv, Sampling se)
  {
		this(st,sx,sv,se,0.0);
  }


	public SemblanceOld(
    Sampling st, Sampling sx, Sampling sv, Sampling se, double fpeak)
  {
    _st = st;
    _sx = sx;
    _sv = sv;
    _se = se;
		if (fpeak>0.0) _fpeak = fpeak;
		_sigma = 1.0/(_fpeak*_st.getDelta());
  	setSmoother();
  }

  public float[][] applyHyp(float[][] p) {
    int nv = _sv.getCount();
    float[][] s = new float[nv][];
    for (int iv=0; iv<nv; ++iv) {
      double v = _sv.getValue(iv);
      float[][] q = nmoHyp(v,_sf,_st,_sx,p);
      s[iv] = semblance(q);
    }
    return s;
  }

  public float[][] applyHypDt(float[][] p) {
    int ntv = _sv.getCount();
    float[][] s = new float[ntv][];
    for (int iv=0; iv<ntv; ++iv) {
      double dt = _sv.getValue(iv);
      float[][] q = nmoHypDt(dt,_sf,_st,_sx,p);
      s[iv] = semblance(q);
    }
    return s;
  }

 public float[][][] applyNonHyp(final float[][] p) {
    final int nv = _sv.getCount();
    final int ne = _se.getCount();
    final float[][][] s = new float[nv][ne][];
		Parallel.loop(0,nv, new Parallel.LoopInt() {
			public void compute(int iv) {	
    		for (int ie=0; ie<ne; ++ie) {
    		  double v = _sv.getValue(iv);
					double e = _se.getValue(ie);
    		  float[][] q = nmoNonHypE(v,e,_sf,_st,_sx,p);
    		  s[iv][ie] = semblance(q);
    		}
			}
		});
		float[][][] sr = new float[s[0][0].length][ne][nv];
		transposeP(s,sr);
   	return sr;
  }

	public float[][][] applyNonHyp(final float[][] p, final Sampling sh) {
    final int nv = _sv.getCount();
    final int nh =  sh.getCount();
    final float[][][] s = new float[nv][nh][];
		Parallel.loop(0,nv, new Parallel.LoopInt() {
			public void compute(int iv) {	
    		for (int ih=0; ih<nh; ++ih) {
    		  double v = _sv.getValue(iv);
					double h =  sh.getValue(ih);
    		  float[][] q = nmoNonHypH(v,h,_sf,_st,_sx,p);
    		  s[iv][ih] = semblance(q);
    		}
			}
		});
		float[][][] sr = new float[s[0][0].length][nh][nv];
		transposeP(s,sr);
   	return sr;
  }

	public void setSigma(double sigma) {
		_sigma = sigma;
		setSmoother();
	}

	public void setStretchMute(double sf) {
		_sf = sf;
	}

	public float[][] flattenGatherHyp(float[][] p, float[] vnmo) {
		int nt = _st.getCount();
    double dt = _st.getDelta();
    double ft = _st.getFirst();
    int nx = _sx.getCount();
    float[] t = new float[nt];
    float[][] q = new float[nx][nt];
    SincInterp si = new SincInterp();
    for (int ix=0; ix<nx; ++ix) {
   	  double x = _sx.getValue(ix);
   	  for (int it=0; it<nt; ++it) {
   	    double t0 = _st.getValue(it);
   	    double t1 = sqrt(t0*t0+(x*x)/(vnmo[it]*vnmo[it])); 
        t[it] = ((t1-t0)/t0>_sf)?0f:(float)t1;
   	  }
   	  si.interpolate(nt,dt,ft,p[ix],nt,t,q[ix]);
   	}
		return q;
	}

	public float[][] flattenGatherNonHypE(float[][] p, float[] vnmo, float[] e) {
		int nt = _st.getCount();
    double dt = _st.getDelta();
    double ft = _st.getFirst();
    int nx = _sx.getCount();
    float[] t = new float[nt];
    float[][] q = new float[nx][nt];
    SincInterp si = new SincInterp();
    for (int ix=0; ix<nx; ++ix) {
      double x = (float)_sx.getValue(ix);
			double xx = x*x;
      for (int it=0; it<nt; ++it) {
        double t0 = _st.getValue(it);
        double xxg = xx/(vnmo[it]*vnmo[it]);
      	double t1 = sqrt(t0*t0+xxg-
								2.0*e[it]*xx*xxg/(t0*t0*vnmo[it]*vnmo[it]+(1.0+2.0*e[it])*xx));
        t[it] = ((t1-t0)/t0>_sf)?0f:(float)t1;
      }
      si.interpolate(nt,dt,ft,p[ix],nt,t,q[ix]);
    }
		return q;
	}

	public float[][] flattenGatherNonHypH(float[][] p, float[] vnmo, float[] vhor) {
		int nt = _st.getCount();
    double dt = _st.getDelta();
    double ft = _st.getFirst();
    int nx = _sx.getCount();
    float[] t = new float[nt];
    float[][] q = new float[nx][nt];
    SincInterp si = new SincInterp();
    for (int ix=0; ix<nx; ++ix) {
      double x = _sx.getValue(ix);
			double xx = x*x;
      for (int it=0; it<nt; ++it) {
				double vnmo2 = vnmo[it]*vnmo[it];
				double vhor2 = vhor[it]*vhor[it];
      	double xxg = xx/vnmo2;
        double t0 = _st.getValue(it);
      	double t1 = sqrt(t0*t0+xxg-
								(vhor2-vnmo2)*xx*xxg/(t0*t0*vnmo2*vnmo2+vhor2*xx));
        t[it] = ((t1-t0)/t0>_sf)?0f:(float)t1;
      }
      si.interpolate(nt,dt,ft,p[ix],nt,t,q[ix]);
    }
		return q;
	}

	public float[] pickMaxSemblance(float[][] sr) {
		int nv = sr.length;
		int nt = sr[0].length;
		double fv = _sv.getFirst();
		double dv = _sv.getDelta();
		float[] v = new float[nt];
		float[][] srt = new float[nt][nv];
		transposeP(sr,srt);
		int[] ind = new int[1];
		for (int it=0; it<nt; ++it) {
			float ms = max(srt[it],ind);
			v[it] = (float)(ind[0]*dv+fv);
		}
		return v;
	}

	public float[][] pickMaxSemblance(float[][][] sr, Sampling sa) {
		int nt = sr.length;
		int na = sr[0].length;
		int nv = sr[0][0].length;
		double fv = _sv.getFirst();
		double dv = _sv.getDelta();
		double fa =  sa.getFirst();
		double da =  sa.getDelta();
		float[] v = new float[nt];
		float[] a = new float[nt];
		int[] ind = new int[2];
		for (int it=0; it<nt; ++it) {
			float ms = max(sr[it],ind);
			v[it] = (float)(ind[0]*dv+fv);
			a[it] = (float)(ind[1]*da+fa);
		}
		return new float[][]{v,a};
	}


  ///////////////////////////////////////////////////////////////////////////

  private double _sigma = 2.0; // smoother half-width
  private double _fpeak = 30.0; // peak frequency of data
  private double _sf    = 1000.0; // stretch factor mute
  private Sampling _st,_sx,_sv,_se;
	private RecursiveExponentialFilter _ref;

  private static float[][] nmoHyp(
    double vnmo, double sf, Sampling st, Sampling sx, float[][] p)
  {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    int nx = sx.getCount();
    float[] t = new float[nt];
    float[][] q = new float[nx][nt];
    SincInterp si = new SincInterp();
    for (int ix=0; ix<nx; ++ix) {
      double x = sx.getValue(ix);
      double xxg = (x*x)/(vnmo*vnmo);
      for (int it=0; it<nt; ++it) {
        double t0 = st.getValue(it);
        double t1 = sqrt(t0*t0+xxg); 
        t[it] = ((t1-t0)/t0>sf)?0f:(float)t1;
      }
      si.interpolate(nt,dt,ft,p[ix],nt,t,q[ix]);
    }
    return q;
  }

  private static float[][] nmoHypDt(
    double dtv, double sf, Sampling st, Sampling sx, float[][] p)
  {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    int nx = sx.getCount();
    float[] t = new float[nt];
    float[][] q = new float[nx][nt];
    SincInterp si = new SincInterp();
    for (int ix=0; ix<nx; ++ix) {
      double x = sx.getValue(ix);
      for (int it=0; it<nt; ++it) {
        double t0 = st.getValue(it);
        double vnmo = sqrt((x*x)/(t0*t0-(dtv+t0)*(dtv+t0)));
        double xxg = (x*x)/(vnmo*vnmo);
        double t1 = sqrt(t0*t0+xxg); 
        t[it] = ((t1-t0)/t0>sf)?0f:(float)t1;
      }
      si.interpolate(nt,dt,ft,p[ix],nt,t,q[ix]);
    }
    return q;
  }

  private static float[][] nmoNonHypE(
    double vnmo, double e, double sf, Sampling st, Sampling sx, float[][] p)
  {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    int nx = sx.getCount();
    float[] t = new float[nt];
    float[][] q = new float[nx][nt];
    SincInterp si = new SincInterp();
    for (int ix=0; ix<nx; ++ix) {
      double x = sx.getValue(ix);
			double xx = x*x;
      double xxg = xx/(vnmo*vnmo);
			double a = 2.0*e*xx*xxg;
			double b = (1.0+2.0*e)*xx;
      for (int it=0; it<nt; ++it) {
        double t0 = st.getValue(it);
      	double t1 = sqrt(t0*t0+xxg);//-a/(t0*t0*vnmo*vnmo+b));
        t[it] = ((t1-t0)/t0>sf)?0f:(float)t1;
      }
      si.interpolate(nt,dt,ft,p[ix],nt,t,q[ix]);
    }
    return q;
  }

  private static float[][] nmoNonHypH(
    double vnmo, double vhor, double sf, Sampling st, Sampling sx, float[][] p)
  {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    int nx = sx.getCount();
    float[] t = new float[nt];
    float[][] q = new float[nx][nt];
    SincInterp si = new SincInterp();
    for (int ix=0; ix<nx; ++ix) {
      double x = sx.getValue(ix);
			double xx = x*x;
			double vnmo2 = vnmo*vnmo;
			double vhor2 = vhor*vhor;
      double xxg = xx/vnmo2;
			double a = (vhor2-vnmo2)*xx*xxg;
      for (int it=0; it<nt; ++it) {
        double t0 = st.getValue(it);
      	double t1 = sqrt(t0*t0+xxg-a/(t0*t0*vnmo2*vnmo2+vhor2*xx));
        t[it] = ((t1-t0)/t0>sf)?0f:(float)t1;
      }
      si.interpolate(nt,dt,ft,p[ix],nt,t,q[ix]);
    }
    return q;
  }

	private void setSmoother() {
	  _ref = new RecursiveExponentialFilter(_sigma);
    _ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_SLOPE);
		// For testing only
    //_ref.setEdges(RecursiveExponentialFilter.Edges.INPUT_ZERO_VALUE);
    //_ref.setEdges(RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE);	
	}

	/**
	 * Thanks to Simon Luo for this method
	 */
  private float[] semblance(float[][] q) {
    int nx = q.length;
    int nt = q[0].length;
    float[] sn = new float[nt];
    float[] sd = new float[nt];
    for (int ix=0; ix<nx; ++ix) {
      for (int it=0; it<nt; ++it) {
        float qi = q[ix][it];
        if (qi==0.0f) continue;
        sn[it] += qi;
        sd[it] += qi*qi;
      }
    }
    mul(sn,sn,sn);
    mul(nx,sd,sd);
    _ref.apply(sn,sn);
    _ref.apply(sd,sd);
    float[] s = sn;
    for (int it=0; it<nt; ++it)
      s[it] = (sd[it]>0.0)?sn[it]/sd[it]:0.0f;
    return s;
  }

  private static void transposeP(final float[][][] x, final float[][][] y) {
    final int n1 = x[0][0].length; 
    final int n2 = x[0].length;    
    final int n3 = x.length;
    final int block1 = 16;//(int)(0.08*n1); // About the same speed, maybe slightly faster
    final int block2 = 16;//(int)(0.08*n2);
    final int block3 = 16;//(int)(0.08*n3);
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
	private static void transposeP(final float[][] x, final float[][] y) {
    final int n1 = x[0].length; 
    final int n2 = x.length;    
    final int block1 = 16;//(int)(0.08*n1); // About the same speed, maybe slightly faster
    final int block2 = 16;//(int)(0.08*n2);
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


};
