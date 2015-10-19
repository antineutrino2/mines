/** 
 * Models a gather using various NMO equations. 
 * Thanks to Simon Luo for helping me start this code.
 * @author Andrew Munoz, CSM
 * @version 10.26.13
 */

package sem;

import java.util.Random;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class Gather {
	
	/**
	 * @param st time sampling
	 * @param sx time sampling
	 * @param sv nmo velocity sampling
	 */
	public Gather(Sampling st, Sampling sx, Sampling sv) {
		this(st,sx,sv,null);
	}
	
	/**
	 * @param st time sampling
	 * @param sx time sampling
	 * @param sv nmo velocity sampling
	 * @param se eta sampling
	 */
	public Gather(Sampling st, Sampling sx, Sampling sv, Sampling se) {
		if (se!=null) _hyp = false;
			_st = st;
			_sx = sx;
			_sv = sv;
			_se = se;
	}

	/**
	 * Models a gather with single layer hyperbolic moveout equation.
	 * Uses a ricker wavelet with a certian peak frequency.
	 * @param fpeak peak frequency
	 * @param snr signal to noise ration
	 * @param vnmo Vnmo velocities
	 * @author Simon Luo
	 */
	public float[][] modelGatherHyp(
		double fpeak, double snr, float[] vnmo) 
	{
	  int nt = _st.getCount();
    double dt = _st.getDelta();
    double ft = _st.getFirst();
    double lt = _st.getLast();
    int nx = _sx.getCount();
    double dx = _sx.getDelta();
    double fx = _sx.getFirst();
    float[][] p = new float[nx][nt];
    double thalf =  1.0/fpeak;
    Random random = new Random(314);
    for (int jt=0; jt<nt; ++jt) {
      if (random.nextDouble()<0.1) // if reflection, ...
        continue;
      double t0 = _st.getValue(jt);
      double a0 = 2.0*random.nextDouble()-1.0;
      double v0 = vnmo[_st.indexOfNearest(t0)];
      double gamma = 1.0/(v0*v0);
      for (int ix=0; ix<nx; ++ix) {
        double x = (float)_sx.getValue(ix);
        double t = sqrt(t0*t0+x*x*gamma);
        int itlo = max(0,(int)((t-thalf-ft)/dt));
        int ithi = min(nt-1,(int)((t+thalf-ft)/dt));
        for (int it=itlo; it<=ithi; ++it) {
          double twave = _st.getValue(it)-t;
          p[ix][it] += (float)(a0*ricker(fpeak,twave));
        }
      }
    }
    return addRandomNoise((float)snr,p); // noise
	}

	/**
	 * Models a gather with single-layer non-hyperbolic moveout equation.
	 * Specifically, the equation defined by Alkhaliha and Tsvankin (1995) which 
	 * uses the quartic coefficient from the Taylor series expansion of travel
	 * times A4. This coefficient assumes to have a negligable effect from shear
	 * wave vertical velocities and is controlled primarily by eta and Vnmo. 
	 * This method also uses a ricker wavelet with a certian peak frequency.
	 * @param fpeak peak frequency
	 * @param snr signal to noise ration
	 * @param vnmo Vnmo velocities
	 * @param eta eta = (epsilon - delta)/(1 + 2*delta)
	 */
	public float[][] modelGatherNonHyp(
		double fpeak, double snr, float[] vnmo, float[] eta) 
	{
	  int nt = _st.getCount();
    double dt = _st.getDelta();
    double ft = _st.getFirst();
    double lt = _st.getLast();
    int nx = _sx.getCount();
    double dx = _sx.getDelta();
    double fx = _sx.getFirst();
    float[][] p = new float[nx][nt];
    double thalf =  1.0/fpeak;
    Random random = new Random(314);
    for (int jt=0; jt<nt; ++jt) {
      if (random.nextDouble()<0.1) // if reflection, ...
        continue;
      double t0 = _st.getValue(jt);
      double a0 = 2.0*random.nextDouble()-1.0;
      double v0 = vnmo[_st.indexOfNearest(t0)];
      double e0 = 2.0*eta[_st.indexOfNearest(t0)];
      double gamma1 = 1.0/(v0*v0);
      double gamma2 = gamma1*gamma1;
      for (int ix=0; ix<nx; ++ix) {
        double x = (float)_sx.getValue(ix);
				// less accurate for long offsets
        //double t = sqrt(t0*t0+x*x*gamma1 -
				//					 e0*gamma2/(t0*t0); 
				// more accurate for long offsets
				double b = t0*t0*v0*v0+(1.0+e0)*x*x;
				if (b==0.0) b = 1.0;
				double a =  e0*x*x*x*x*gamma1/b;
        double t = sqrt(t0*t0+x*x*gamma1 - a);
        int itlo = max(0,(int)((t-thalf-ft)/dt));
        int ithi = min(nt-1,(int)((t+thalf-ft)/dt));
        for (int it=itlo; it<=ithi; ++it) {
          double twave = _st.getValue(it)-t;
          p[ix][it] += (float)(a0*ricker(fpeak,twave));
        }
      }
    }
    return addRandomNoise((float)snr,p); // noise
	}

	public float[] makeLinearParameter(
		double vmin, double vmax) 
  {
    int nt = _st.getCount();
    double dt = _st.getDelta();
    double ft = _st.getFirst();
    double lt = _st.getLast();
    float[] v = new float[nt];
    for (int it=0; it<nt; ++it) 
      v[it] = (float)((it*vmax+(nt-1-it)*vmin)/(nt-1));
    return v;
  }

	public float[] makeRandomParameter(
		double emin, double emax, double sigma) 
  {
    int nt = _st.getCount();
    double dt = _st.getDelta();
    double ft = _st.getFirst();
    double lt = _st.getLast();
    float[] e = new float[nt];
    Random random = new Random(333);
		double range = emax-emin;
    for (int it=0; it<nt; ++it)
      e[it] = (float)(emin + random.nextDouble()*0.5*range);
    float[] es = smooth(e,sigma);
    return es;
  }

	//////////////////////////////////////////////////////////////////////////////

	private Sampling _st,_sx,_sv,_se;
	private boolean _hyp = true;

	private float[] smooth(float[] x, double sigma) {
    float[] xs = new float[x.length];
		RecursiveExponentialFilter rgf = new RecursiveExponentialFilter(sigma);
		rgf.setEdges(RecursiveExponentialFilter.Edges.OUTPUT_ZERO_SLOPE);
    rgf.apply(x,xs);
		for (int i=1; i<16; ++i) 
    	rgf.apply(xs,xs);
		return xs;
	}

	/////////////////////////////////////////////////////////////////////////////
	/**
	 * Utilities. 
	 * @author Simon Luo
	 */
	private float[] makeLinearVelocity(
    double vmin, double vmax, Sampling st) 
  {
    int nt = st.getCount();
    double dt = st.getDelta();
    double ft = st.getFirst();
    double lt = st.getLast();
    float[] v = new float[nt];
    for (int it=0; it<nt; ++it) {
      v[it] = (float)((it*vmax+(nt-1-it)*vmin)/(nt-1));
    }
    return v;
  }
  private float[][] addRandomNoise(float r, float[][] p) {
    int nt = p[0].length;
    int nx = p.length;
    //float pmax = max(abs(p)); // peak signal
    float prms = sqrt(sum(mul(p,p))/nt/nx); // rms of signal
    Random random = new Random(3);
    float[][] s = sub(randfloat(random,nt,nx),0.5f);
    RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1.0);
    rgf.apply0X(s,s); // noise, bandlimited in time only
    float srms = sqrt(sum(mul(s,s))/nt/nx); // rms of noise
    return add(mul(prms/(srms*r),s),p); // r = rms-signal / rms-noise
  }
  private float ricker(double fpeak, double time) {
    double x = PI*fpeak*time;
    double xx = x*x;
    return (float)((1.0-2.0*xx)*exp(-xx));
  }
};
