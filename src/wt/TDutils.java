package wt;

/**
 * A set of methods for TD conversion.
 * The primary purpose of this is for speedup 
 * over python.
 * @author Andrew Munoz, Colorado School of Mines
 * @version 2.03.12
 */

import java.util.Random;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;


public class TDutils {

  /**
   * Stretches well to trace or trace to well. 
   * Linearly interpolate a function with a certian non-uniform sampling in one dim to a 
   * uniformed sampling in another dim. (z to t)
   * 
   * @param zot z(t) values used for z to t conversion
   * @param f   values in z
   * @param st  sampling in t
   * @return fs,toz the stretched values and the unstretching values
   */
   public float[] stretch(float[] zot, float[] f, Sampling st) {
     int nf = f.length;
     float dt = (float)st.getDelta(); // The desired sampling interval
     int n2 = (int)((zot[nf-1]/dt)-(zot[0]/dt)); // The desired new length of uniformly sampled func
     float[] fs = new float[n2];
     float[] ut = rampfloat(zot[0],dt,n2); // x values in the interpolation
     int j=0;
     for (int i=1; i<n2-1; ++i) {
       // Linear decimation? (squeezing)
       while (ut[i]>zot[j+1] || ut[i]<zot[j]) ++j;
       // Linear Interpolation (stretching)
       fs[i] = ysolve(f[j],f[j+1],ut[i],zot[j],zot[j+1]);
     }
     fs[0]=f[0];
     fs[n2-1]=f[nf-1];
     return fs;
   }
   private float ysolve(float yo, float y1, float x, float xo, float x1) {
     float y = yo + (x-xo)*((y1-yo)/(x1-xo));
     return y;
   }

  /**
   * Finds vertical two-way time to depth.
   * Uses the time to tau relationship tm(t) and fixed point iteration
   * to solve for tn(z) using the initial t(z) used to find tm. 
   * <p>
   * @param sz Sampling in depth
   * @param st Sampling in vertical two-way time
   * @param stau Sampling in tau (of the synthetic)
   * @param taut tau to time relationship found using dtw tau(t)
   * @param ttau t(tau)
   * @param tz original tau to depth relationship using the sonic log
   * @return tn time to depth relationship 
   *
   */
  public float[][] tdcurve(Sampling sz, Sampling st, Sampling stau, float[] tpi, float[] taupi, float[] tauz) {
    int nz = sz.getCount();
    int nt = st.getCount();
    int ntau = stau.getCount();
    int nj = tpi.length;
    float dz = (float)sz.getDelta();
    float dt = (float)st.getDelta(); 
    float dtau = (float)stau.getDelta(); 
    float fz = (float)sz.getFirst();
    float ft = (float)st.getFirst();
    float ftau = (float)stau.getFirst();
    float ft2 = tpi[nj-1]*dt+ft;
    float[] tm = new float[nj];
    float[] tpv = new float[nj];
    float[] zpv = new float[nj];
    int k=0;
    //  Remove the horizontal and vertical repeating values so function is always increasing 
    for (int j=0; j<nj-1; ++j,++k) {
      if (tpi[j]!=tpi[j+1] && taupi[j]!=taupi[j+1]) { // slope 1
        tpv[k] = tpi[j]*dt;
        zpv[k] = taupi[j]*dtau+ftau;
        continue;
      }
      if (tpi[j]!=tpi[j+1] && taupi[j]==taupi[j+1]) { // slope h
        tpv[k] = ((tpi[j]*dt)+(tpi[j+1]*dt))/2;
        zpv[k] = taupi[j]*dtau+ftau;
        ++j;
        continue;
      }
      if (tpi[j]==tpi[j+1] && taupi[j]!=taupi[j+1]) { // slope v
        tpv[k] = tpi[j]*dt;
        zpv[k] = ((taupi[j]*dtau+ftau)+(taupi[j+1]*dtau+ftau))/2;
        ++j;
        continue;  
      }
      if (j==nj-2) {
        tpv[k] = tpi[nj-1]*dt;
        zpv[k] = taupi[nj-1]*dtau;
      }
    }
    // Reverse TPV
    float[] tpvr = new float[k];
    float[] zpvr = new float[k];
    for (int i=k-1,j=0; i>=0; --i,++j) {
      tpvr[j] = tpv[i];
      zpvr[j] = zpv[i];
    }
    // Compute t(tau)
    float[] ttau = new float[ntau];
    CubicInterpolator cttau = new CubicInterpolator(CubicInterpolator.Method.LINEAR,k,zpvr,tpvr);
    for (int itau=0; itau<ntau; ++itau) 
      ttau[itau] = cttau.interpolate(ftau+itau*dtau);

    // Compute tau(t) (for warping sy to tr)
    float ftn = tpi[tpi.length-1];
    float ltn = tpi[0];
    int ntn = (int)(ltn-ftn);
    float[] taut = new float[ntn];
    Sampling ss = new Sampling(ntn,dt,ftn*dt);
    CubicInterpolator ctaut = new CubicInterpolator(CubicInterpolator.Method.LINEAR,k,tpvr,zpvr);
    for (int it=0; it<ntn; ++it) 
      taut[it] = ctaut.interpolate(ftn*dt+it*dt);
    
    // Compute t(z) from t(tau=tau(z))
    float[] tz = new float[nz];
    int nz2=0;
    int iz2=0;
    LinearInterpolator li = new LinearInterpolator();
    li.setUniform(ntau,dtau,ftau,ttau);
    for (int iz=0; iz<nz; ++iz) 
      tz[iz] = li.interpolate(tauz[iz]);
    

    // Plots for testing

    //SimplePlot si = new SimplePlot();
    //si.addPoints(sz,tz);
    //si.addTitle("tz");
    //si.setHLabel("Depth (km)");
    //si.setVLabel("Time (s)");
    //SimplePlot si3 = new SimplePlot();
    //si3.addPoints(stau,ttau);
    //si3.addTitle("ttau");
    //si3.setHLabel("tau (s)");
    //si3.setVLabel("t (s)");
    ////SimplePlot si4 = new SimplePlot();
    ////si4.addPoints(ss,taut);
    ////si4.addTitle("taut");
    ////si4.setHLabel("t (s)");
    ////si4.setVLabel("tau (s)");
    //SimplePlot si5 = new SimplePlot();
    //si5.addPoints(sz,tauz);
    //si5.addTitle("tauz");
    //si5.setHLabel("Depth (km)");
    //si5.setVLabel("tau (s)");
    //SimplePlot si6 = new SimplePlot();
    //si6.addPoints(sz,tauz).setLineColor(java.awt.Color.RED);
    //si6.addPoints(sz,tz);
    //si6.addTitle("tz");
    //si6.setHLabel("Depth (km)");
    //si6.setVLabel("t (s)");

    return new float[][]{ttau,taut,tz};
  }
 

  public float[] applyTWavelet(float[] r, Sampling st, Sampling sz, float fp, float[] toz) {
    int nz = toz.length;
    float dt = (float)st.getDelta();
    float dz = (float)sz.getDelta();
    float fz = (float)sz.getFirst();
    int ntu = (int)(1.5*nz);
    float[] ti = rampfloat(toz[0],dt,ntu);
    float[] zi = rampfloat(fz,dz,nz);
    double sig = 1/(fp*PI);
    double s4 = sig*4;
    float[] sy = new float[ntu];
    int lt=0;
    for (int iz=0; iz<nz-1; ++iz) {
      for (int it=0; it<ntu; ++it) {
        if ((toz[iz]-s4)<=ti[it] && ti[it]<=(toz[iz]+s4)) {
          sy[it] += r[iz]*ricker(fp,toz[iz]-ti[it]);
          lt=it;
        }
      }
    }
    float[] sy2 = new float[lt];
    copy(lt,sy,sy2);
    return sy2;
  }
  public float ricker(double fpeak, double time) {
    double x = PI*fpeak*time;
    double xx = x*x;
    return (float)((1.0-2.0*xx)*exp(-xx));
  }

  public float[] maketauz(float[] v, Sampling sz) {
    int nz = v.length;
    float fz = (float)sz.getFirst();
    float dz = (float)sz.getDelta();
    // Find the median of the top values of v
    int nfm = 30;
    float[] vfs = new float[nfm];
    copy(nfm,v,vfs);
    MedianFinder mf = new MedianFinder(nfm);
    float md = mf.findMedian(vfs);
    float to = 2f*(fz/md);
    float[] toz = new float[nz];
    float[] ov = div(1f,v);
    toz[0] = to;
    for (int i=1; i<nz; ++i) {
      toz[i] = toz[i-1]+(ov[i])*dz*2f;
    }
    return toz;
  }

  public float[] addRickerWavelet2(Sampling st, double fpeak, float[] zot, float[] f) {
    int n1 = f.length;
    int nzt = zot.length;
    double dt = st.getDelta();
    double sig = 1/(fpeak*PI*dt);
    int dr = (int)(4.0*sig);
    float[] g = new float[nzt];
    for (int zt=0; zt<nzt; ++zt) {
      float tz = zot[zt];
      int tzms = (int)((tz-dr));
      int tzps = (int)((tz+dr));
      for (int ti=tzms; ti<tzps; ++ti) {
        //g[zt+ti+dr] += ricker((tz-ti)/sig);
      }
      //SimplePlot si = new SimplePlot();
      //si.addPoints(g);
    }
    return g;
  }



  /**
   * Gets T and D values from indicies, reverses path order, and throws out samples.
   * Creates an array of time and depth values from the index arrays.
   * The path order of values is reversed so that they are increasing.
   * Path values are also removed to take out any bias or monotonicity. 
   * @param tpi array of indicies that correspond to time
   * @param zpi array of indicies that correspond to depth 
   * @param sz  sampling in depth
   * @param st  sampling in time
   * @return 
   */
   public float[][] reverseInputThrow(float[] tpi, float[] zpi, Sampling st, Sampling sz) {
     int np = tpi.length; //=zpi.length
     int nt = st.getCount();
     double dt = st.getDelta();
     double ft = st.getFirst();
     int nz = sz.getCount();
     double dz = sz.getDelta();
     double fz = sz.getFirst();
     float[] tp = new float[np];
     float[] zp = new float[np];
     int i1=1; int n1=1;
     tp[0] = (float)(ft+dt*tpi[np-1]);
     zp[0] = (float)(fz+dz*zpi[np-1]);
     // Reverse the order, replace with values, and throw out vert/horiz bias 
     for (int i=1; i<np-1; ++i) {
       int ii = np-i-1;
       if ((tpi[ii]==tpi[ii+1] || tpi[ii]==tpi[ii-1]) 
        || (zpi[ii]==zpi[ii+1] || zpi[ii]==zpi[ii-1])) 
         continue;
       else {
         tp[i1] = (float)(ft+dt*tpi[ii]);
         zp[i1] = (float)(fz+dz*zpi[ii]);
         i1+=1; n1+=1;
       }
     }
     tp[i1] = (float)(ft+dt*tpi[0]);
     zp[i1] = (float)(fz+dz*zpi[0]);
     float[] tp2 = new float[n1];
     float[] zp2 = new float[n1];
     copy(n1,tp,tp2);
     copy(n1,zp,zp2);
     return new float[][]{tp2,zp2};
   }

  /**
   * Interpolation of velocity curves.
   * Use the CubicInterpolator to interpolate the
   * current velocity. 
   * @param zot z(t) values used for depth to time conversion
   * @param toz t(z) values used for time to depth conversion
   * @param sz  sampling in depth
   * @param st  sampling in time
   * @return zot,toz interpoloated z(t) and t(z)
   */
   public float[][] interpolation(float[] tpv, float[] zpv, Sampling st, Sampling sz) {
     int nt = st.getCount();
     double dt = st.getDelta();
     double ft = st.getFirst();
     int nz = sz.getCount();
     double dz = sz.getDelta();
     double fz = sz.getFirst();
     monotonicityTest(tpv);
     monotonicityTest(zpv);
     int n1 = tpv.length; //=zpv.length
     float[] toz = new float[nt];
     float[] zot = new float[nz];
     CubicInterpolator ctoz = new CubicInterpolator(CubicInterpolator.Method.LINEAR,n1,tpv,zpv);
     CubicInterpolator czot = new CubicInterpolator(CubicInterpolator.Method.LINEAR,n1,zpv,tpv);
     for (int it=0; it<nt; ++it) 
       toz[it] = ctoz.interpolate((float)(ft+it*dt));
     for (int iz=0; iz<nz; ++iz)  
       zot[iz] = czot.interpolate((float)(fz+iz*dz));
     return new float[][]{toz,zot};
   }
   private void monotonicityTest(float[] x) {
     // Prints location of non-monotonic points
     int n = x.length;
     for (int i=0; i<n-1; ++i) {
       if (x[i]>=x[i+1]) 
         System.out.format("not montonic between %d and %d with a val of %f and %f\n",i,i+1,x[i],x[i+1]);
     }
   }
  
  /**
   * Velocity.
   * Calculates the velocity using the time and 
   * depth values for depth to time and time to
   * depth conversion. Uses centered derivative.
   * @param zot z(t) values used for depth to time conversion
   * @param toz t(z) values used for time to depth conversion
   * @param sz  sampling in depth
   * @param st  sampling in time
   * @return vz,vt velocity for T2D and D2T conversion
   */
   public float[][] velocities(float[] zot, float[] toz, Sampling sz, Sampling st) {
     int nt = toz.length;
     int nz = zot.length;
     double dz = sz.getDelta();
     double dt = st.getDelta();
     float[] vz = new float[nz];
     float[] vt = new float[nt];
     for (int iz=1; iz<nz-1; ++iz) 
       vz[iz] = (float)((zot[iz+1]-zot[iz-1])/(dt));
     for (int it=1; it<nt-1; ++it)
       vt[it] = (float)((4*dz)/(toz[it+1]-toz[it-1]));
     return new float[][]{vz,vt};
   }

  /**
   * RMS error.
   * Find the RMS error for each path in an array of paths
   * @param tpe  array of path errors 
   * @return rms array of rms error for each path
   */
  public float[] RMS(float[][] x) {
    int np = x.length;
    float[] rms = new float[np];
    for (int ip=0; ip<np; ++ip) {
      int nj = x[ip].length;
      float rs=0;
      for (int j=0; j<nj; ++j) 
        rs += x[ip][j]*x[ip][j];
      rms[ip] = sqrt(rs/nj);
    }
    return rms;
  }

  /**
   * Get seismic horizons from 3 1D arrays in specified x slice
   */
  public float[][] horizon2d(float[] x1, float[] x2, float[] x3, 
   Sampling s3, Sampling s1, int x3c, float[] tz, String z)
  {
    int n = x1.length;
    float dc = 0f;
    if (z=="z") dc = .102f; // km- datum correction
    double d3 = s3.getDelta();
    double d1 = s1.getDelta();
    double f1 = s1.getFirst()+dc;
    float[] x1d = new float[n];
    float[] x2d = new float[n];
    int x3i;
    int xi=0, j=0;
    for (int i=0; i<n; ++i) {
      x3i = (int)(x3[i]/d3);
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
    if (z=="z") {
      float[] x1z = new float[j];
      for (int i1=0; i1<j; ++i1) {
        for (int i2=1; i2<tz.length; ++i2) {
          if ((f1+(i2-1)*d1)<=x1d[i1] && x1d[i1]<=(f1+i2*d1)) {
	    x1z[i1] = tz[i2]; 
      	  }
        }
      }
      return new float[][]{x1z,x2t};
    }
    return new float[][]{x1t,x2t};
  }
  /**
   */
   public void printds(float[] x1, float[] x2, float[] x3, 
    Sampling s3, Sampling s2, int x3c, int x2c, String horz)
   {
     int n = x1.length;
     double d3 = s3.getDelta();
     double d2 = s2.getDelta();
     double f3 = s3.getFirst();
     double f2 = s2.getFirst();
     int x3i,x2i,xi=0;
     for (int i=0; i<n; ++i) {
       x3i = (int)((x3[i]-f3)/d3);
       x2i = (int)((x2[i]-f2)/d2);
       if (x3i==x3c) {
        if (x2i==x2c || x2i==x2c+1) xi = i;}
     }
     if (xi==0) System.out.println("Horizon Fail");
     System.out.println(horz+"="+x1[xi]);
   }



     

  /**
   * Correlation.
   * Compute the cross correlation of two signals and 
   * return the rms avg of values of the cross correlation
   * @param f signal 1 
   * @param g signal 2
   * @param sf signal 1 sampling
   * @param sg signal 2 sampling
   * @return rms the rms avg value of the cross correlation. 
   */
   public float correlation(float[] f, float[] g, Sampling sf, Sampling sg) {
     int nf = sf.getCount();
     int ng = sg.getCount();
     int ff = (int)(sg.getFirst()/sf.getDelta());
     int fg = 0;
     int nz = max(nf,ng);
     int fz = 0;
     float[] z = new float[nz];
     Conv c = new Conv();
     c.xcor(nf,ff,f,ng,fg,g,nz,fz,z);
     float rms=0f;
     for (int i=0; i<nz; ++i)
       rms += z[i]*z[i];
     rms /= nz;
     rms = sqrt(rms);
     return rms;
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


