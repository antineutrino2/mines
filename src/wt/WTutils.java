package wt;

/**
 * A set of methods to support automatically tying well logs to seismic data.
 * The primary purpose of this is for speedup over python.
 * @author Andrew Munoz, Colorado School of Mines
 * @version 09.11.12
 */

import java.util.Random;
import java.util.Arrays;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;


public class WTutils {

  /**
   * Generates a simple synthetic seismogram.
   * The method computes a sum of shifted wavelets (ricker with peak frequency fp)
   * scaled by the reflection coefficients r found from the well logs. 
   * @param fp The peak frequency of the seismic trace
   * @param r The reflectivity log
   * @param toz The time-depth function t(z) found by integrating the velocity log 
   * @param st The seismic sampling
   * @param sz The log sampling
   * @return sy2 The resulting synthetic seismogram
   */
  public float[] simpleSeismogram(
     Sampling sz, Sampling st, float fp, float[] r, float[] toz) 
  {
	  int nz = sz.getCount();
	  int nt = st.getCount();
    float dz = (float)sz.getDelta();
    float dt = (float)st.getDelta();
    float fz = (float)sz.getFirst();
    float ft = (float)st.getFirst();
    float[] ti = rampfloat(toz[0],dt,nt);
    float[] sy = new float[nt];
    double sig = 10.0/(fp*PI);
    int lt=0;
    for (int iz=0; iz<nz; ++iz) {
      for (int it=0; it<nt; ++it) {
        if ((toz[iz]-sig)<=ti[it] && ti[it]<=(toz[iz]+sig)) {
          sy[it] += r[iz]*ricker(fp,toz[iz]-ti[it]);
          if (abs(sy[it])>0.0001f) lt=it;
        }
      }
    }
    float[] sy2 = new float[lt];
    copy(lt,sy,sy2);
    return sy2;
  }
  private float ricker(double fpeak, double time) {
    double x = PI*fpeak*time;
    double xx = x*x;
    return (float)((1.0-2.0*xx)*exp(-xx));
  }
  private float morlet(double fpeak, double x) {
    double c = 1/sqrt(2*PI);
    return (float)(c*exp(-x*x/2)*cos(2*PI*fpeak*x));
  }

	public float[] convolutionSeismogram(
		float[] r, float[] toz, float[] w, Sampling st, Sampling sz) 
	{
    int nz = toz.length;
		int nw = w.length;
    float dt = (float)st.getDelta();
    float dz = (float)sz.getDelta();
    float fz = (float)sz.getFirst();
    float ft = (float)st.getFirst();
    int ntu = (int)((max(toz)-ft)/dt);
    float[] ti = rampfloat(toz[0],dt,ntu);
    float[] sy = new float[ntu];
		float[] rt = new float[ntu];
		// Convert r to time
		LinearInterpolator li = new LinearInterpolator();
		li.setUniform(nz,dz,0.0,r);
		li.interpolate(ntu,ti,rt);
		// Convolve rt with wavelet
    Conv.conv(nw,0,w,ntu,0,rt,ntu,0,sy);
		return sy;
	}

	public float[] convolutionSeismogramRicker(
		float[] r, float[] toz, float fp, Sampling st, Sampling sz, int ntu) 
	{
    int nz = toz.length;
    float dt = (float)st.getDelta();
    float ft = (float)st.getFirst();
		float dz = (float)sz.getDelta();
    //int ntu = (int)((max(toz)-ft)/dt);
    float[] ti = rampfloat(toz[0],dt,ntu);
    float[] sy = new float[ntu];
		float[] rt = new float[ntu];
		// Convert r to time
		LinearInterpolator li = new LinearInterpolator();
		li.setUniform(nz,dz,0.0,r);
		li.interpolate(ntu,ti,rt);
		// Convolve rt with wavelet
		sy = addRickerWavelet(3.0/fp,rt);
		return sy;
	}
	
	/**
	 * A cosine taper for seismograms. 
	 * Applies a cosine taper to a specified time window using a Hanning Window.
	 * @param st time sampling
	 * @param tb first time in taper window
	 * @param te last time in taper window
	 * @param x signal to be tapered
	 * @return full signal with tapered window
	 */
	public float[] applyTaper(Sampling st, float tb, float te, float[] x) {
		int nt = x.length;
		float dt = (float)st.getDelta();
		float ft = (float)st.getFirst();
		int tbi = (int)floor((tb-ft)/dt);
		int tei = (int)ceil( (te-ft)/dt);
		tbi = (tbi>0 )?tbi:0 ;
		tei = (tei<nt)?tei:nt;
		int ntw = tei-tbi+1-1;
		float[] xt = copy(x);
		for (int it=tbi,j=ntw; it<tei; ++it,--j)
			xt[it] = x[it]*0.5f*(1.0f-cos(2.0f*FLT_PI*j/(ntw)));	
		return xt;	
	}

	/**
	 * A exponential dampner.
	 * Applies a exponential dampner to the ends of a signal
	 * @param st time sampling
	 * @param tb end of first window
	 * @param te start of last window
	 * @param x signal to be dampened
	 * @return full signal with dampened window
	 */
	public float[] applyDampner(Sampling st, float tb, float te, float[] x) {
		int nt = x.length;
		float dt = (float)st.getDelta();
		float ft = (float)st.getFirst();
		int tbi = (int)floor((tb-ft)/dt);
		int tei = (int)ceil( (te-ft)/dt);
		float[] xt = copy(x);
		for (int it=tbi,j=0; it>=0; --it,++j)
			xt[it] = x[it]*exp(-0.1f*j);	
		for (int it=tei,j=0; it<nt; ++it,++j)
			xt[it] = x[it]*exp(-0.1f*j);	
		return xt;	
	}


  /**
   * Integrates a velocity log to get a time-depth function t(z).
   * @param v The velocity log
   * @param sz The log sampling
   * @return toz The time-depth function t(z)
   */
  public float[] tzVlog(
    float[] v, float fz, float dz) 
  {
    int nz = v.length;
		int nm = 50;
		MedianFinder mf = new MedianFinder(nm);
		float vm = mf.findMedian(copy(nm,v));
    float to = fz*(2f/vm);//+1f/2.7432f);// 2.74-m/s = 9000-ft/s, static corr vel
    //float to = fz*(2f/v[0]-1f/35f;
    float[] toz = new float[nz];
    float dz2 = dz*2f;
    float tc = 0.0f;
		toz[0] = tc = to;
    for (int i=1; i<nz; ++i) 
      toz[i] = tc = tc+dz2/v[i];
    return toz;
  }

	/**
   * Gets a depth-time function z(t) from t(z).
	 * Uses inverse linear interpolation.
   * @return zt The depth-time function z(t)
   */
  public float[] getzt(
    float[] tz, Sampling st, Sampling sz) 
  {
		int nt = st.getCount();
		float dt = (float)st.getDelta();
		float ft = (float)st.getFirst();
		int nz = sz.getCount();
		float dz = (float)sz.getDelta();
		float fz = (float)sz.getFirst();
		float[] zt = new float[nt];	
		inverseLinearInterpolation(nz,dz,fz,tz,nt,dt,tz[0],zt,tz[0],tz[nz-1]);
    return zt;
  }


	public float[][] getTimeDepths(Sampling st, float[][] vi) {	
		int nx = vi.length;
		int nt = vi[0].length;
		float dt = (float)st.getDelta();
		float[][] zt = new float[nx][nt];
		float ci = 0.0f;
		float c  = 0.0f;
		for (int ix=0; ix<nx; ++ix) {
			zt[ix][0] = c = ci;
			for (int it=1; it<nt; ++it) 
				zt[ix][it] = c = c + 0.5f*(dt*it)*vi[ix][it];
		}
		return zt;
	}

  /**
   * Computes seismic interval velocity from a time-depth function.
   * Computes by applying a central difference system to
   * the time-depth function t(z). 
   * @param toz t(z) values used for time to depth conversion
   * @param sz  sampling in depth
   * @param st  sampling in time
   * @return v  seismic interval velocity
   */
   public float[] intervalVelocity(float[] toz, Sampling sz) {
     int nz = toz.length;
     int nzm = nz-1;
     float dz = (float)sz.getDelta();
     float[] v = new float[nz];
     for (int iz=1; iz<nzm; ++iz)
       v[iz] = (4*dz)/(toz[iz+1]-toz[iz-1]);
     return v;
   }
  
  /**
   * Computes centered difference to approximate derivative.
   * @param x  sequence
   * @param dx sequence's sampling
   * @return y approximated derivative of sequence
   */
  public float[] centerDiff(float[] x, float dx) {
    int n = x.length;
		int nm = n-1;
    float[] y = new float[n];
    for (int i=1; i<nm; ++i)
      y[i] = (x[i+1]-x[i-1])/(2*dx);
    y[nm] = y[nm-1];
		y[0] = y[1];
    return y;
  }

  /**
   * Computes forward difference to approximate derivative.
   * @param x  sequence
   * @param dx sequence's sampling
   * @return y approximated derivative of sequence
   */
  public float[] forwardDiff(float[] x, float dx) {
    int n = x.length;
		int nm = n-1;
    float[] y = new float[n];
    for (int i=0; i<nm; ++i)
      y[i] = (x[i+1]-x[i])/(dx);
    y[nm] = y[nm-1];
    return y;
  }

	/**
   * Computes backwards difference to approximate derivative.
   * @param x  sequence
   * @param dx sequence's sampling
   * @return y approximated derivative of sequence
   */
  public float[] backwardsDiff(float[] x, float dx) {
    int n = x.length;
		int nm = n-1;
    float[] y = new float[n];
    for (int i=1; i<n; ++i)
      y[i] = (x[i]-x[i-1])/(dx);
    y[0] = y[1];
    return y;
  }

  /**
   * Get 2D seismic horizons from three 1D arrays within a specified slice in x.
   * Output used to plot on seismic 2D vertical slices. 
   * @param x1 Horizon values in z direction
   * @param x2 Horizon values in y direction
   * @param x3 Horizon values in x direction
   * @param s3 Horizon sampling in x direction
   * @param s1 Horizon sampling in z direction
   * @param x3c Well coordinate index in x
   */
  public float[][] horizon2d(float[] x1, float[] x2, float[] x3, 
   Sampling s3, Sampling s1, int x3c)
  {
    int n = x1.length;
    double d3 = s3.getDelta();
    double d1 = s1.getDelta();
    double f1 = s1.getFirst();
    double od3 = 1.0/d3;
    float[] x1d = new float[n];
    float[] x2d = new float[n];
    int x3i,xi=0,j=0;
    for (int i=0; i<n; ++i) {
      x3i = (int)(x3[i]*od3);
      if (x3i==x3c) {
        x1d[j] = x1[i];
        x2d[j] = x2[i];
	++j;
      }
    }
    float[] x1t = new float[j];
    float[] x2t = new float[j];
    copy(j,x1d,x1t);
    copy(j,x2d,x2t);
    return new float[][]{x1t,x2t};
  }
	

	/**
   * Get 2D seismic horizons from three 1D arrays within a specified x2 and x3 curves.
   * Output used to plot on seismic 2D vertical slices. 
   * @param x1 Horizon values in z direction
   * @param x2 Horizon values in y direction
   * @param x3 Horizon values in x direction
   * @param s2 Horizon sampling in y direction
   * @param s1 Horizon sampling in z direction
   * @param x2s Well coordinate curve in y
   * @param x3s Well coordinate curve in x
   */
  public float[][] horizon2d(float[] x1, float[] x2, float[] x3, 
   Sampling s2, Sampling s1, float[] x2s, float[] x3s)
  {
		int np = x1.length;
		int ns = x2s.length;
		int npm = np-1;
		int ii=1;
		float[] x3n = new float[ns];
		float[] x2n = new float[ns];
		float[] x1n = new float[ns];
		float[] xsc = rampfloat(x2s[0],(float)s2.getDelta(),ns);
		for (int is=0; is<ns; ++is) {
			float x2j = x2s[is];
			float x3j = x3s[is];
			for (int ip=0; ip<npm; ++ip) {
				 float x2hm = x2[ip];
				 float x2hp = x2[ip+1];
				 if (x2hm <= x2j && x2j < x2hp) {
				 	float x3hm = x3[ip];
				 	float x3hp = x3[ip+1];
				 	if (x3hm <= x3j && x3j < x3hp) {
						if (x2hm>x2n[ii-1]) {
							//x3n[ii] = x3hm;
							x2n[ii] = x2hm;
							x1n[ii] = x1[ip];
							++ii;
						}
					}
				}
			}
    }	
		--ii;
		//x3n = copy(ii,1,x3n);
		x2n = copy(ii,1,x2n);
		x1n = copy(ii,1,x1n);
		int[] i2 = rampint(0,1,ii);
		quickIndexSort(x2n,i2);
		float[] x1nn = new float[ii];
		float[] x2nn = new float[ii];
		for (int ij=0; ij<ii; ++ij) {
			x1nn[ij] = x1n[i2[ij]];
			x2nn[ij] = x2n[i2[ij]];
		}
		CubicInterpolator ci = new CubicInterpolator(CubicInterpolator.Method.LINEAR,x2nn,x1nn);
		float[] x1s = ci.interpolate(xsc);
		return new float[][]{x1s,xsc};
		//return new float[][]{x1n,x2n};
  }


  /**
   * Gets the intersecting depth or time value of a horizon with a well.
   * @param x1 Horizon values in z direction
   * @param x2 Horizon values in y direction
   * @param x3 Horizon values in x direction
   * @param s3 Horizon sampling in x direction
   * @param s2 Horizon sampling in y direction
   * @param x3c Well coordinate index in x
   * @param x2c Well coordinate index in y
   */
   public float getds(float[] x1, float[] x2, float[] x3, 
    Sampling s3, Sampling s2, int x3c, int x2c)
   {
     int n = x1.length;
     double d3 = s3.getDelta();
     double d2 = s2.getDelta();
     double f3 = s3.getFirst();
     double f2 = s2.getFirst();
     double od3 = 1.0/d3;
     double od2 = 1.0/d2;
     int x3i,x2i,xi=0;
     for (int i=0; i<n; ++i) {
       x3i = (int)((x3[i]-f3)*od3);
       x2i = (int)((x2[i]-f2)*od2);
       if (x3i==x3c) {
        if (x2i==x2c || x2i==x2c+1) xi = i;}
     }
     if (xi==0) System.out.println("Horizon Fail");
     return x1[xi];
   }

  /**
   * Performs a local RMS normalization on a sequence.
   * This method is extremely useful when avoiding a global 
   * normalization. It applies an exponential smoother on a 
   * window of samples, s, around a sample in the sequence
   * and normalizes that value by the smoothed samples.
   * @param x The input sequence
   * @param s The window size, sigma, to smooth samples
   * @return y The normalized sequence
   */
   public float[] localRMSnorm(
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
	public float[][][] localRMSnorm(final float[][][] x, final float s) {
	 	final int n3 = x.length;
	 	final int n2 = x[0].length;
	 	final int n1 = x[0][0].length;
		final float[][][] y = new float[n3][n2][n1];
		Parallel.loop(n3,new Parallel.LoopInt() {
    	public void compute(int i3) {
				for (int i2=0; i2<n2; ++i2)
					y[i3][i2] = localRMSnorm(x[i3][i2],s);
		}});
		return y;
	}


  /**
   * Returns the maximum horizontal deviation distance from the surface.
   * @param x The array of x coordinates
   * @param y The array of y coordinates
   * @return The maximum horizontal deviation distance
   */
  public float maxDeviation(
    float[] x, float[] y)
  {
    int nx = x.length;
    int ny = y.length;
    float xmax=0,ymax=0;
    for (int ix=0; ix<nx; ++ix)
      if (x[ix]-x[0] > x[ix-1]-x[0]) xmax=x[ix]-x[0];
    for (int iy=0; iy<ny; ++iy)
      if (y[iy]-y[0] > y[iy-1]-y[0]) xmax=y[iy]-y[0];
    return sqrt(xmax*xmax+ymax*ymax);
  }

  /**
   * A despike filter.
   * This filter defines a spike as a point that has an amplitude
   * greater than l percent of the adjacent values, and removes it 
   * by setting spike amplitudes to the median amplitude of the nearest 
   * 2m+1 samples. For non-spikes, the system retains the values of the 
   * original sequence.
   * @param m filter half-width.
   * @param l filter clip percentage in terms of adjacent samples.
   * @param sx input sequence x[n].
   * @return output sequence y[n].
   */
  public float[] despike(int m, float l, float[] x) {
    int nt = x.length;
    int ntm = nt-1;
    int m2 = 2*m+1;
    int mp = (m>0)?m+1:0;
    float xvp,xi1m,xi1p;
    int i1m=0,i1p=0;
    float[] y = new float[nt];
    float[] r = new float[m2];
		float mean = sum(x)/nt;
		float std = sqrt(sum(pow(sub(x,mean),2.0f))/nt);
    for (int i1=0; i1<nt; ++i1) {
      if (i1>0 && i1<ntm) {
        i1m = i1-1; i1p = i1+1;
      } 
      else if (i1==0) i1m = i1+1;
      else            i1p = i1-1;
      //xvp = abs(x[i1]*l);
      //xi1m = abs(x[i1m]);
      //xi1p = abs(x[i1p]);
      xvp = abs(x[i1]);
      xi1m = mean-std*l;
      xi1p = mean+std*l;
      //if (xi1m<xvp && xi1p<xvp) {
      if (xi1m>xvp || xi1p<xvp) {
        for (int i2=i1-m,i3=0; i2<=i1+m; ++i2,++i3) 
          r[i3] = ((i2>0)?((i2<nt)?x[i2]:0):0); // Extrapolate 0 off the ends
        Arrays.sort(r);
        y[i1] = r[mp];
      }
      else y[i1] = x[i1];
    }
		return y;
  }
   /**
    * Converts ints to floats
    */
   public float[][] floats(int[][] x) {
     int n1 = x.length;
     int n2 = x[0].length;
     float[][] y = new float[n1][n2];
     for (int i1=0; i1<n1; ++i1) {
       for (int i2=0; i2<n2; ++i2)
         y[i1][i2] = (float)x[i1][i2];
     }
     return y;
   }
   public float[] floats(int[] x) {
     int n1 = x.length;
     float[] y = new float[n1];
     for (int i1=0; i1<n1; ++i1) {
       y[i1] = (float)x[i1];
     }
     return y;
   }

  /**
   * Sorts an array a1 by indexes x and outputs and array b1
   */
  public void asort(float[] a1, int[] x, float[] b1) {
    int n = a1.length;
    for (int i=0; i<n; ++i) 
      b1[i] = a1[x[i]];
    }

	/**
   * Interpolates missing data using values for nearby samples.
   * @param missing value that represents missing data.
   * @param x input sequence x[n] with missing values.
   * @return output sequence y[n] with missing values interpolated.
   */
  public static float[] interpolateMissingValues(float missing, float[] x) {
    int nt = x.length;
    float[] y = new float[nt];
    int i=1;
    for (int n=0; n<nt; ++n) {
      float x0 = x[n];
      if (x0!=missing) {
        y[n] = x0;
      } else {
			while (x[n+i]<missing) ++i; 
			float x1 = x[n+i];
			if (n>0) {
			  float xm = y[n-1];
			  y[n] = xm + (x1-xm)/((i+n)-(n-1));
			} else { 
			// extrapolate the slope of the first two values to get n=0
  		  	float x2 = x[i+1]; 
			  y[n] = x1 - i*(x2-x1); 
			}
			i=1;
      }
    }
    return y;
  }

  public float[] addWavelet(float[] w, float[] f) {
    int n1 = f.length;
		int nw = w.length;
    float[] g = new float[n1];
    Conv.conv(nw,0,w,n1,0,f,n1,0,g);
    return g;
  }

  public float[] addMorletWavelet(double fpeak, float[] f) {
    int n1 = f.length;
    int ih = (int)(3.0/fpeak);
    int nh = 1+2*ih;
    float[] h = new float[nh];
    for (int jh=0; jh<nh; ++jh)
      h[jh] = morlet(fpeak,jh-ih);
    float[] g = new float[n1];
    Conv.conv(nh,-ih,h,n1,0,f,n1,0,g);
    return g;
  }

	public float[][] getSlice(final int iz, final float[][][] x) {
		final int n3 = x.length;
		final int n2 = x[0].length;
		final float[][] y = new float[n3][n2];
		Parallel.loop(n3,new Parallel.LoopInt() {
    public void compute(int i3) {
			for (int i2=0; i2<n2; ++i2)
				y[i3][i2] = x[i3][i2][iz];
    }});
		return y;
	}

 /**
  * Computes the NRMS value 
  * @return the NRMS value, this measure is between 0.0 and 2.0.
  */
 public float computeNrms(float[] fs, float[] gs) {
 	 int n1 = fs.length;
   float[] d = sub(gs,fs);
   float frms = sqrt(sum(mul(fs,fs))/(n1));  
   float grms = sqrt(sum(mul(gs,gs))/(n1));
   float drms = sqrt(sum(mul( d, d))/(n1));
   float nrms = (2.0f*drms)/(frms+grms);
   return nrms;
 }
  public float computeNrms(float[][] fs, float[][] gs) {
 	 int n1 = fs[0].length;
   float[][] d = sub(gs,fs);
   float frms = sqrt(sum(mul(fs,fs))/(n1));  
   float grms = sqrt(sum(mul(gs,gs))/(n1));
   float drms = sqrt(sum(mul( d, d))/(n1));
   float nrms = (2.0f*drms)/(frms+grms);
   return nrms;
 }
 public float computeNrms(float[][][] fs, float[][][] gs) {
 	 int n1 = fs[0][0].length;
   float[][][] d = sub(gs,fs);
   float frms = sqrt(sum(mul(fs,fs))/(n1));  
   float grms = sqrt(sum(mul(gs,gs))/(n1));
   float drms = sqrt(sum(mul( d, d))/(n1));
   float nrms = (2.0f*drms)/(frms+grms);
   return nrms;
 }

	public float[][] syntheticImageFill(
		int xi, int ns, int f, float[] s, float[][] g) 
	{
		for (int i=0; i<ns; ++i) 
			g[xi][f+i] = s[i];
		return g;
	}
	public float[][][] syntheticImageFill(
		int x2i, int x3i, int ns, int f, float[] s, float[][][] g) 
	{
		for (int i=0; i<ns; ++i) 
			g[x3i][x2i][f+i] = s[i];
		return g;
	}
	public float[][][] syntheticImageFill(
		float[] x2, float[] x3, int ns, int f, float[] s, float[][][] g) 
	{
		for (int i=0; i<ns; ++i) 
			g[inro(x3[i])][inro(x2[i])][f+i] = s[i];
		return g;
	}
	public float[] logImageFill(
		int nt, int nz, float pnull,
		float[] l, float[] tz, float[] ta, float[] g, int f) 
	{
		int ntm = nt-1;
		float zc = 0.0f;
		float ssum = 0.0f;
		float vel = pnull;
		float[] sl = div(1.0f,l);
		for (int it=0; it<ntm; ++it) {
			for (int iz=0; iz<nz; ++iz) {
				float tzv = tz[iz];
				if (tzv>ta[it] && tzv<=ta[it+1]) {
					ssum += sl[iz];
					zc += 1.0f;
				}
				if (tzv>ta[it+1]) 
					break;
			}
			if (ssum>0.0f) 
				vel = 1.0f/(ssum/zc);
			g[it] = vel;
			zc=0f;ssum=0f;vel=pnull;
		}
		return g;
	}
	public float[] logImageFill(
		int nt, int nz, float pnull,
		float[] l, float[] tz, float[] ta, float[] g) 
	{
	return logImageFill(nt,nz,pnull,l,tz,ta,g,0);
	}

	public float[][] logImageFill(
		int xi, int nt, int nz, float pnull,
		float[] l, float[] tz, float[] ta, float[][] g) 
	{
		int ntm = nt-1;
		float zc = 0.0f;
		float ssum = 0.0f;
		float vel = pnull;
		float[] sl = div(1.0f,l);
		for (int it=0; it<ntm; ++it) {
			for (int iz=0; iz<nz; ++iz) {
				float tzv = tz[iz];
				if (tzv>ta[it] && tzv<=ta[it+1]) {
					ssum += sl[iz];
					zc += 1.0f;
				}
				if (tzv>ta[it+1]) 
					break;
			}
			if (ssum>0.0f) 
				vel = 1.0f/(ssum/zc);
			g[xi][it] = vel;
			zc=0f;ssum=0f;vel=pnull;
		}
		return g;
	}
	public float[][][] logImageFill(
		int x2i, int x3i, int nt, int nz, float pnull,
		float[] l, float[] tz, float[] ta, float[][][] g) 
	{
		int ntm = nt-1;
		float zc = 0.0f;
		float ssum = 0.0f;
		float vel = pnull;
		float[] sl = div(1.0f,l);
		for (int it=0; it<ntm; ++it) {
			for (int iz=0; iz<nz; ++iz) {
				float tzv = tz[iz];
				if (tzv>ta[it] && tzv<=ta[it+1]) {
					ssum += sl[iz];
					zc += 1.0f;
				}
				if (tzv>ta[it+1]) 
					break;
			}
			if (ssum>0.0f) 
				vel = 1.0f/(ssum/zc);
			g[x3i][x2i][it] = vel;
			zc=0f;ssum=0f;vel=pnull;
		}
		return g;
	}

  public float[][][] imageTimeToDepth3(
    Sampling s1, Sampling s2, Sampling s3, 
    float[][][] timg, float[][][] tzs) 
  {
	  int nz = tzs[0][0].length;
	  int n2 = s2.getCount();
	  int n3 = s3.getCount();
	  float[][][] dimg = new float[n3][n2][nz];
	  SincInterp si = new SincInterp();
	  double[] s2v = s2.getValues();
	  double[] s3v = s3.getValues();
	  for (int i3=0; i3<n3; ++i3)
	    for (int i2=0; i2<n2; ++i2)
	      for (int i1=0; i1<nz; ++i1)
	  			dimg[i3][i2][i1] = si.interpolate(s1,s2,s3,timg,tzs[i3][i2][i1],s2v[i2],s3v[i3]);
	  return dimg;
  }


	public int inro(float x) {
		return (int)round(x);
	}
	public int inro(double x) {
		return (int)round(x);
	}
	public int infl(float x) {
		return (int)floor(x);
	}
	public int infl(double x) {
		return (int)floor(x);
	}
	public int ince(float x) {
		return (int)ceil(x);
	}
	public int ince(double x) {
		return (int)ceil(x);
	}

//////////////////////////////////////////////////////////////////// 
//// Thanks to Luming Liang for the following methods           //// 
////////////////////////////////////////////////////////////////////

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

 /**
   * Inverse linear interpolation.
   * @author Dave Hale, Colorado School of Mines
   * @version 1989.06.02
   * @author Translated to java by Andrew Munoz, Colorado School of Mines
   * @version 2012.09.15
   *
   * Compute a regularly-sampled, monotonically increasing function x(y) from a 
   * regularly-sampled, monotonically increasing function y(x) by inverse linear 
   * interpolation.
   * @param nx	number of samples of y(x)
   * @param dx	x sampling interval; dx &gt 0.0 is required
   * @param fx	first x
   * @param y	array[nx] of y(x) values; y[0] &lt y[1] &lt ... &lt y[nx-1] required
   * @param ny	number of samples of x(y)
   * @param dy	y sampling interval; dy &gt 0.0 is required
   * @param fy	first y
   * @param xylo  x value assigned to x(y) when y is less than smallest y(x)
   * @param xyhi  x value assigned to x(y) when y is greater than largest y(x)
   * @return x  array[ny] of x(y) values
   */
  private void inverseLinearInterpolation(
    int nx, float dx, float fx, float[] y, 
    int ny, float dy, float fy, float[] x, float xylo, float xyhi) 
  { 
    int nxi,nyo,jxi1,jxi2,jyo;
    float dxi,fxi,dyo,fyo,fyi,yo,xi1,yi1,yi2,yid,q; 
    nxi = nx; dxi = dx; fxi = fx;
    nyo = ny; dyo = dy; fyo = fy;
    fyi = y[0];
    // loop over output y less than smallest input y
    for (jyo=0,yo=fyo; jyo<nyo; jyo++,yo+=dyo) {
      if (yo>=fyi) break;
      x[jyo] = xylo;
    }
    // loop over output y between smallest and largest input y
    if (jyo==nyo-1 && yo==fyi) {
      x[jyo++] = fxi;
      yo += dyo;
    }
    jxi1 = 0;
    jxi2 = 1;
    xi1 = fxi;
    while (jxi2<nxi && jyo<nyo) {
      yi1 = y[jxi1];
      yi2 = y[jxi2];
      if (yi1<=yo && yo<=yi2) {
        yid = abs(yi2-yi1);
	x[jyo++] = (yid>0.0f)?xi1+dxi*(yo-yi1)/(yi2-yi1):xi1;
        yo += dyo;
    } else {
        jxi1++;
        jxi2++;
        xi1 += dxi;
      }
    }
    // loop over output y greater than largest input y
    while (jyo<nyo) x[jyo++] = xyhi;
  }

};
