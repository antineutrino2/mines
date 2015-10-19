/** 
 * Computes semblance using hyperbolic and non-hyperbolic NMO approximations.
 * Thanks to Simon Luo for helping me start this code.
 * @author Andrew Munoz, CSM
 * @version 01.30.14
 */

package sem;

import java.util.Random;

import edu.mines.jtk.mosaic.SimplePlot;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class Semblance {

  public Semblance() {
  }

  /**
   * Constructs a Semblance.
   * Semblance smoother in time is assumed to be sigma=2.0
   * @param st0 time sampling of cmp gathers
   * @param sx  offset sampling of cmp gathers
   */
	public Semblance(Sampling st0, Sampling sx) {
    this(st0,sx,0.0);
  }

  /**
   * Constructs a Semblance.
   * Use fpeak to adjust size of time-smoothing window for 
   * semblance calculations.
   * @param st0 time sampling of cmp gathers
   * @param sx  offset sampling of cmp gathers
   * @param fpeak peak frequency of seismic 
   */
	public Semblance(Sampling st0, Sampling sx, double fpeak) {
    _st0 = st0;
    _sx = sx;
    _fpeak = fpeak;
    if (fpeak>0.0) _sigma = 1.0/(_fpeak*st0.getDelta());
  	setSmoother();
  }

  // loop over dtv
  public float[][] applyHypDt(float[][] p, Sampling stv) {
    float xlim = (float)max(abs(_sx.getFirst()),abs(_sx.getLast()));
    return applyHypDt(p,stv,xlim);
  }
  public float[][] applyHypDt(float[][] p, Sampling stv, float xlim) {
    int n1 = _st0.getCount();
    int ntv = stv.getCount();
    double dt = _st0.getDelta();
    double xmax = _sx.getFirst();
    double x2 = xmax*xmax;
    float[][] s = new float[ntv][];
    float[] v = new float[n1];
    for (int it=0; it<ntv; ++it) {
      //for (int io=1; io<n1; ++io)
      //  v[io] = (float)sqrt(x2/(pow((io*dt-stv.getValue(it)),2)-io*dt));
      //v[0] = v[1];
      double vv = sqrt(x2/(pow((-stv.getValue(it)),2)));
      float[][] q = nmoV(vv,_mute,_st0,_sx,p,xlim);
      s[it] = semblance(q);
    }
    return s;
  }
  // loop over q 
  public float[][] applyHypQ(float[][] p, float vmin, float vmax, float dq) {
    double qmin = 1.0/(vmax*vmax);
    double qmax = 1.0/(vmin*vmin);
    int nq = inro((qmax-qmin)/dq);
    Sampling sq = new Sampling(nq,dq,qmin);
    float[][] s = new float[nq][];
    float[][] v = new float[nq][];
    for (int iq=0; iq<nq; ++iq) {
      double qv = sq.getValue(iq);
      float[][] q = nmoQ1(qv,_mute,_st0,_sx,p);
      s[iq] = semblance(q);
    }
    return s;
  }
  // loop over v
  public float[][] applyHypV(float[][] p, float vmin, float vmax, float dv) {
    float xlim = (float)(abs(_sx.getFirst())>abs(_sx.getLast())?_sx.getFirst():_sx.getLast());
    return applyHypV(p,vmin,vmax,dv,xlim);
  }
  public float[][] applyHypV(float[][] p, float vmin, float vmax, float dv, float xlim) {
    int nv = (int)((vmax-vmin)/dv);
    float[][] s = new float[nv][];
    float[][] v = new float[nv][];
    for (int iv=0; iv<nv; ++iv) {
      float vnmo = vmin+iv*dv;
      float[][] q = nmoV(vnmo,_mute,_st0,_sx,p,xlim);
      s[iv] = semblance(q);
      if (iv==100) SimplePlot.asPoints(s[iv]);
    }
    return s;
  }


// Non-hyperbolic methods 
  public float[][][] applyNonHypV(
    final float[][] p, final Sampling sv, final Sampling se)
  {
    final int nv = sv.getCount();
    final int ne = se.getCount();
    final float[][][] s = new float[nv][ne][];
		Parallel.loop(0,nv, new Parallel.LoopInt() {
			public void compute(int iv) {	
    		for (int ie=0; ie<ne; ++ie) {
    		  double v = sv.getValue(iv);
					double e = se.getValue(ie);
    		  float[][] q = nmoNonHypE(v,e,_mute,_st0,_sx,p);
    		  s[iv][ie] = semblance(q);
    		}
			}
		});
		float[][][] sr = new float[s[0][0].length][ne][nv];
		transposeP(s,sr);
   	return sr;
  }

  // Non-hyperbolic methods 
  public float[][][] applyNonHypV(
    final float[][] p, final Sampling sv, final Sampling se,
    final float[] vl, final float[] vu, 
    final float[] el, final float[] eu)
  {
    final int nv = sv.getCount();
    final int ne = se.getCount();
    final float[][][] s = new float[nv][ne][];
		Parallel.loop(0,nv, new Parallel.LoopInt() {
			public void compute(int iv) {	
    		for (int ie=0; ie<ne; ++ie) {
    		  double v = sv.getValue(iv);
					double e = se.getValue(ie);
    		  float[][][] qm = nmoNonHypE(v,e,_mute,_st0,_sx,p,vl,vu,el,eu);
    		  s[iv][ie] = semblance(qm);
    		}
			}
		});
		float[][][] sr = new float[s[0][0].length][ne][nv];
		transposeP(s,sr);
   	return sr;
  }

	public float[][][] applyNonHypQ(
    final float[][] p, final Sampling q1, final Sampling q2)
  {
    final int nq1 = q1.getCount();
    final int nq2 = q2.getCount();
    final float[][][] s = new float[nq1][nq2][];
		Parallel.loop(0,nq1, new Parallel.LoopInt() {
			public void compute(int i1) {	
    		double v = 1.0/sqrt(q1.getValue(i1));
    		for (int i2=0; i2<nq2; ++i2) {
					double e = 0.5*q2.getValue(i2)/q1.getValue(i1);
    		  float[][] q = nmoNonHypE(v,e,_mute,_st0,_sx,p);
    		  s[i1][i2] = semblance(q);
    		}
			}
		});
		float[][][] sr = new float[s[0][0].length][nq2][nq1];
		transposeP(s,sr);
   	return sr;
  }

//  public float[][][] applyNonHypDt(
//    final float[][] p, final Sampling dtv, final Sampling dte)
//  {
//    final int ndtv = dtv.getCount();
//    final int ndte = dte.getCount();
//    final float[][][] s = new float[ndtv][ndte][];
//    Parallel.loop(0,ndtv, new Parallel.LoopInt() {
//      public void compute(int i1) {	
//        for (int i2=0; i2<ndte; ++i2) {
//          double v = 1.0/sqrt(dtv.getValue(i1));
//        	double e = 0.5*dte.getValue(i2)/q1.getValue(i1);
//          float[][] q = nmoNonHypE(v,e,_mute,_st0,_sx,p);
//          s[iv][ih] = semblance(q);
//        }
//      }
//    });
//    float[][][] sr = new float[s[0][0].length][nq2][nq1];
//    transposeP(s,sr);
//    return sr;
//  }


	public void setSigma(double sigma) {
		_sigma = sigma;
		setSmoother();
	}

	public void setStretchMute(double mute) {
		_mute = mute;
	}


	public float[][] flattenGatherHyp(float[][] p, float[] vnmo) {
		int nt = _st0.getCount();
    double dt = _st0.getDelta();
    double ft = _st0.getFirst();
    int nx = _sx.getCount();
    float[] t = new float[nt];
    float[][] q = new float[nx][nt];
    SincInterp si = new SincInterp();
    for (int ix=0; ix<nx; ++ix) {
   	  double x = _sx.getValue(ix);
   	  for (int it=0; it<nt; ++it) {
   	    double t0 = _st0.getValue(it);
   	    double t1 = sqrt(t0*t0+(x*x)/(vnmo[it]*vnmo[it])); 
        t[it] = ((t1-t0)/t0>_mute)?0f:(float)t1;
   	  }
   	  si.interpolate(nt,dt,ft,p[ix],nt,t,q[ix]);
   	}
		return q;
	}

	public float[][] flattenGatherNonHypE(float[][] p, float[] vnmo, float[] e) {
		int nt = _st0.getCount();
    double dt = _st0.getDelta();
    double ft = _st0.getFirst();
    int nx = _sx.getCount();
    float[] t = new float[nt];
    float[][] q = new float[nx][nt];
    SincInterp si = new SincInterp();
    for (int ix=0; ix<nx; ++ix) {
      double x = (float)_sx.getValue(ix);
			double xx = x*x;
      for (int it=0; it<nt; ++it) {
        double t0 = _st0.getValue(it);
        double xxg = xx/(vnmo[it]*vnmo[it]);
      	double t1 = sqrt(t0*t0+xxg-
								2.0*e[it]*xx*xxg/(t0*t0*vnmo[it]*vnmo[it]+(1.0+2.0*e[it])*xx));
        t[it] = ((t1-t0)/t0>_mute)?0f:(float)t1;
      }
      si.interpolate(nt,dt,ft,p[ix],nt,t,q[ix]);
    }
		return q;
	}

	public float[][] flattenGatherNonHypH(float[][] p, float[] vnmo, float[] vhor) {
		int nt = _st0.getCount();
    double dt = _st0.getDelta();
    double ft = _st0.getFirst();
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
        double t0 = _st0.getValue(it);
      	double t1 = sqrt(t0*t0+xxg-
								(vhor2-vnmo2)*xx*xxg/(t0*t0*vnmo2*vnmo2+vhor2*xx));
        t[it] = ((t1-t0)/t0>_mute)?0f:(float)t1;
      }
      si.interpolate(nt,dt,ft,p[ix],nt,t,q[ix]);
    }
		return q;
	}

	public float[] pickMaxSemblance(float[][] sr, Sampling sv) {
		int nv = sr.length;
		int nt = sr[0].length;
		double fv = sv.getFirst();
		double dv = sv.getDelta();
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

	public float[][] pickMaxSemblance(float[][][] sr, Sampling sa, Sampling sv) {
		int nt = sr.length;
		int na = sr[0].length;
		int nv = sr[0][0].length;
		double fv =  sv.getFirst();
		double dv =  sv.getDelta();
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

  public static float[][] localRMSnorm(float[][] x, float s) {
	 	int n2 = x.length;
	 	int n1 = x[0].length;
		float[][] y = new float[n2][n1];
		for (int i2=0; i2<n2; ++i2)
				y[i2] = localRMSnorm(x[i2],s);
		return y;
	}
  public static float[] localRMSnorm(
     float[] x, float s)
   {
     int n = x.length;
     float[] y = new float[n];
     float[] yy = new float[n];
     float[] xx = mul(x,x);
		 RecursiveExponentialFilter es = new RecursiveExponentialFilter(s);
		 //es.setEdges(RecursiveExponentialFilter.Edges.OUTPUT_ZERO_VALUE);
     es.apply(xx,yy);
     for (int i=0; i<n; ++i)
       y[i] = x[i]/sqrt(yy[i]);
     return y;
   }




  ///////////////////////////////////////////////////////////////////////////

  private double _sigma = 2.0; // smoother half-width
  private double _fpeak = 30.0; // peak frequency of data
  private double _mute    = 1000.0; // stretch factor mute
  private Sampling _st,_sx,_sv,_st0,_se;
	private RecursiveExponentialFilter _ref;

  private static float[][] nmoQ1(
    double q1, double mute, Sampling st, Sampling sx, float[][] p)
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
      double xxg = x*x*q1;
      for (int it=0; it<nt; ++it) {
        double t0 = st.getValue(it);
        double t1 = sqrt(t0*t0+xxg); 
        t[it] = ((t1-t0)/t0>mute)?0f:(float)t1;
      }
      si.interpolate(nt,dt,ft,p[ix],nt,t,q[ix]);
    }
   return q;
  }
  private static float[][] nmoV(
    double vnmo, double mute, Sampling st, Sampling sx, float[][] p, float xmax)
  {
    return nmoV(fillfloat((float)vnmo,st.getCount()),mute,st,sx,p,xmax);
  }
  private static float[][] nmoV(
    float[] vnmo, double mute, Sampling st, Sampling sx, float[][] p, float xmax)
  {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    int nx = sx.getCount();
    int dir = (sx.getFirst()>0)?1:-1;
    int iex = (dir<0)?nx:nx-sx.indexOfNearest(xmax);
    int ibx = (dir<0)?sx.indexOfNearest(xmax):0;
    float[] t = new float[nt];
    float[][] q = new float[nx][nt];
    SincInterp si = new SincInterp();
    for (int ix=ibx; ix<iex; ++ix) {
      double x = sx.getValue(ix);
      for (int it=0; it<nt; ++it) {
        double xxg = (x*x)/(vnmo[it]*vnmo[it]);
        double t0 = st.getValue(it);
        double t1 = sqrt(t0*t0+xxg); 
        t[it] = ((t1-t0)/t0>mute)?0f:(float)t1;
      }
      si.interpolate(nt,dt,ft,p[ix],nt,t,q[ix]);
    }
    return q;
  }
 
  private static float[][] nmoNonHypE(
    double vnmo, double e, double mute, Sampling st, Sampling sx, float[][] p)
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
      	double t1 = sqrt(t0*t0+xxg-a/(t0*t0*vnmo*vnmo+b));
        t[it] = ((t1-t0)/t0>mute)?0f:(float)t1;
      }
      si.interpolate(nt,dt,ft,p[ix],nt,t,q[ix]);
    }
    return q;
  }

  private static float[][][] nmoNonHypE(
    double vnmo, double e, double mute, Sampling st, Sampling sx, float[][] p,
    float[] vl, float[] vu, float[] el, float[] eu)
  {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    int nx = sx.getCount();
    float[] t = new float[nt];
    float[][] m = new float[1][nt]; // mask
    float[][] q = new float[nx][nt];
    SincInterp si = new SincInterp();
    for (int ix=0; ix<nx; ++ix) {
      double x = sx.getValue(ix);
			double xx = x*x;
      double xxg = xx/(vnmo*vnmo);
			double a = 2.0*e*xx*xxg;
			double b = (1.0+2.0*e)*xx;
      for (int it=0; it<nt; ++it) {
        if (vl[it]<=vnmo && vnmo<=vu[it] && el[it]<=e && e<=eu[it]) 
          m[0][it] = 1f;
        double t0 = st.getValue(it);
      	double t1 = sqrt(t0*t0+xxg-a/(t0*t0*vnmo*vnmo+b));
        t[it] = ((t1-t0)/t0>mute)?0f:(float)t1;
      }
      si.interpolate(nt,dt,ft,p[ix],nt,t,q[ix]);
    }
    return new float[][][]{q,m};
  }


  private static float[][] nmoNonHypH(
    double vnmo, double vhor, double mute, Sampling st, Sampling sx, float[][] p)
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
        t[it] = ((t1-t0)/t0>mute)?0f:(float)t1;
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
  private float[] semblance(float[][][] qm) {
    float[][] q = qm[0];
    float[][] m = qm[1];
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
      s[it] = (m[0][it]>0.0)?((sd[it]>0.0)?sn[it]/sd[it]:0.0f):0.0f;
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

	private int inro(float x) {
		return round(x);
	}
	private int inro(double x) {
		return (int)round(x);
	}
	private int infl(float x) {
		return (int)floor(x);
	}
	private int infl(double x) {
		return (int)floor(x);
	}
	private int ince(float x) {
		return (int)ceil(x);
	}
	private int ince(double x) {
		return (int)ceil(x);
	}


};
