package wt;

import java.util.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.Cfloat.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Synthetic seismogram of a vertically propogating plane waves. 
 * Computes a 1-D seismogram in using plane waves, and can model  
 * kinematic effects including muliples, attenuation, and dispersion. 
 * Applicable for any vertically heterogeneous media and can contain 
 * any number of zones with different velocity, density, and Q. 
 * Any number of sources and recievers can be modeled at any depth, and
 * isotropic (vibroseis, airgun) or anisotropic (dynamite) sources can 
 * be specified as well as velocity (geophone) or pressure (hydrophone) 
 * recievers. 
 * The attenuation is assumed constant. For Q \leq 8 the error (from 
 * the small angle approximation) becomes greater than 1%, but this 
 * quality factor is not encountered often for most geophysical
 * applications. A complex frequency IFFT must be used to compute the
 * seismogram. The wavefields are computed using the propogator matrix 
 * method of upgoing and downgoing wavefields.
 *
 * The theory comes from a derivation from Dr. Dave Hale's notes (CSM) and 
 * two papers:
 * 
 * Kjartansson, E. (1979), "Constant Q-wave Propagation and Attenuation", 
 * Journal of Geophysical Research, Vol. 84 No. B9, pp. 4737-4747
 *
 * Gilbert, F., and Backus, G. (1966), "Propagator Matricies in Elastic 
 * Wave and Vibration Problems", Geophysics, Vol. XXXI No. 2, pp. 326-332
 *
 * @author Andrew Munoz, Colorado School of Mines
 * @version 28.01.2013
 */

public class SyntheticSeismogram {

	// Used to import arrays of velocity, density, and q
	public void setZones(Sampling sz, float[] v, float[] d, float[] q) { 
		int nz = sz.getCount();
		double dz = sz.getDelta();
		double fz = sz.getFirst();
  	checklengths(nz,v,d,q); 
		addZone(0.0,v[0],d[0],q[0]);
		for (int iz=0; iz<nz; ++iz) 
			addZone(iz*dz+fz,v[iz],d[iz],q[iz]);
	}
	public void setZones(Sampling sz, float[] v, float[] d, float q) { 
	  setZones(sz,v,d,fillfloat(q,v.length));
	}

	public void addReciever(double z) {
		_recievers.add((float)z);
	}

  // Sets the reflection coefficient for the surface
	public void setSRC(double src) {
		_r1 = (float)src;
	}

	// Sets the oversampling factor to compute nfft for the FFT
	public void setOversampling(double c) {
		_oversampling = (float)c;
	}
	
	// Sets the decay rate for the complex frequency
	public void setDecay(double decay) {
		_decay = (float)decay;
	}

	// Turn off interbeded multiples
	public void turnOffMultiples(boolean multiples) {
		if (multiples) _mult = 0.001f;
	}

	// Adds a ricker wavelet to the seismogram with a certian peak frequency fp [Hz]
	public void addRicker(double fp) {
		_ricker = (float)fp;
	}


	// Generates a layer with a constant depth, velocity, density and q.
	public void addZone(double z, double v, double d, double q) {
    addZone((float)z,(float)v,(float)d,(float)q);
  }
  public void addZone(float z, float v, float d, float q) {
		Zone zone = _zones.get(z);
		if (zone!=null) {
			zone.v = v;
			zone.d = d;
			zone.q = q;
		} 
		else {
			zone = new Zone(z,v,d,q,0.0f,false);
			_zones.put(z,zone);
		}
	  for (Float hk=_zones.higherKey(z); hk!=null; hk=_zones.higherKey(hk)) {
      Zone zk = _zones.get(hk);
      if (zk.source) {
        zk.v = v;
        zk.d = d;
        zk.q = q;
      } else {
        break;
      }
    }
	}

  public void addSource(double z, double s) {
		addSource((float)z,(float)s);
	}
  public void addSource(float z, float s) {
		Zone zone = _zones.get(z);
		if (zone!=null) {
			zone.s = s;
		} 
		else {
		  Float zk = _zones.floorKey(z);
		  //if (zk==null) zk = _zones.ceilingKey(z);
			Zone zz = _zones.get(zk);
		  zone = new Zone(z,zz.v,zz.d,zz.q,s,true);
		  _zones.put(z,zone);
		}
	}

  /**
   * 1D Synthetic Seismogram. 
	 * @param st time sampling of seismic
   * @return seismograms 
   */
  public float[][] getSeismogram(Sampling st) {
    int nt = st.getCount();
    float dt = (float)st.getDelta();
		int ntm = nt-1;

		// Frequency Sampling
    int nfft = FftReal.nfftFast((int)(_oversampling*nt));
    int nw = nfft/2+1;
    float dw = 2.0f*FLT_PI/(nfft*dt); // rad/sample
		float df = dw/(2f*FLT_PI); // cycles/s
		float f0 = 0.5f/dt; // nyquist cycles/s
		float w0 = 2.0f*FLT_PI*f0; // nyquist rad/sample
		Sampling sw = new Sampling(nw,df,0.0);
 
		// complex seismogram with real and im parts in frequency
		int nr = _recievers.size();
		float[][] cseis  = new float[nr][2*nw];
		float[][] cseisr = new float[nr][  nw];
		float[][] cseisi = new float[nr][  nw];

		// Surface
    Float zs = _zones.firstKey();

		// i
		Cfloat i0 = new Cfloat(0f,1f);
		Cfloat c1 = new Cfloat(1f,0f);

    // imaginary part of the complex frequency
    float eps = -log(_decay)/(dt*nt);//dt*nfft

		for (int iw=0,iwr=0,iwi=1; iw<nw; ++iw,iwr+=2,iwi+=2) {
			Cfloat cw = new Cfloat(iw*dw,eps); // complex frequency
			Cfloat ciw = i0.times(cw); // i times frequency
			Cfloat nciw = ciw.neg(); 
			
			Zone zu,zl;
			Float uk,lk;

			// upgoing wavefield in lower halfspace is zero
    	Cfloat cd1 = new Cfloat(1f,0f);
    	Cfloat cu1 = new Cfloat(0f,0f);
    	
			// No source in lower halfspace
    	Cfloat cds = new Cfloat(0f,0f);
    	Cfloat cus = new Cfloat(0f,0f);

			Cfloat b11l = cus;
			Cfloat b21l = cus;

			// Compute wavefield of last layer
			lk = _zones.lastKey();
			zl = _zones.get(lk);

			float gamma = atan(1f/zl.q)/FLT_PI;
			 // complex velocity
			zl.cv = pow(nciw.over(w0),gamma).times(zl.v*cos(gamma*0.5f*FLT_PI));
			 // complex wavenumber
   		zl.ck = cw.over(zl.cv);
			 // complex impedance
		  zl.cim = zl.cv.times(zl.d); 

			// loop over layers last to first to compute dn and u1
			for (uk=_zones.lowerKey(lk); uk!=null; lk=uk,uk=_zones.lowerKey(lk)) {
				zl = _zones.get(lk);
				zu = _zones.get(uk);

				// Compute wavefield of each layer
			  gamma = atan(1f/zl.q)/FLT_PI;
			   // complex velocity
			  zu.cv = pow(nciw.over(w0),gamma).times(zu.v*cos(gamma*0.5f*FLT_PI));
			   // complex wavenumber
        zu.ck = cw.over(zu.cv);
			   // complex impedance upper layer
		    zu.cim = zu.cv.times(zu.d);
				 // complex reflectivity
				zu.cr = zu.cim.minus(zl.cim).over(zu.cim.plus(zl.cim));
				zu.cr = zu.cr.times(_mult);
        zu.cr.negEquals();
				
				
				// Compute propagation matrix
				Cfloat crs = c1.over(zu.cr.plus(c1));
				Cfloat cex = ciw.over(zu.cv).times(zl.z-zu.z);
				Cfloat cep = crs.times(exp(cex));
				Cfloat cem = crs.times(exp(cex.neg()));
				Cfloat b11 = cem;
				Cfloat b22 = cep;
				Cfloat b12 = zu.cr.times(cem);
				Cfloat b21 = zu.cr.times(cep);
		
				// Recursively compute wavefields without source to get to z=1
				Cfloat cdt = cd1;
        cd1 = b11.times(cdt).plus(b12.times(cu1));
        cu1 = b21.times(cdt).plus(b22.times(cu1));
			  
				// Compute upgoing and downgoing parts (with isotropic source)
				cdt = cds;
      	cds = b11.times(cdt).plus(b12.times(cus));
      	cus = b21.times(cdt).plus(b22.times(cus));
				if (zl.s>0.0f) {
					cds.minusEquals(cem.times(zl.s));
					cus.plusEquals(cep.times(zl.s));
				}

				// Save B11 and B21
				b11l = b11;
				b21l = b21;
			}	
		
		  // Amplitude of source at the surface
			zl = _zones.get(lk);
		  Cfloat css = new Cfloat(zl.s,0f);

		  // Compute dn and un
			lk = _zones.lastKey();
			zl = _zones.get(lk);
		  //zl.cd = (css.minus(cus.times(_r1)).minus(cds)).over(b11l.plus(b21l.times(_r1)));
    	zl.cd = cus.times(_r1).minus(cds).plus(css).over(cd1.minus(cu1.times(_r1)));
		  zl.cu = new Cfloat(0f,0f);
			
			gamma = atan(1f/zl.q)/FLT_PI;
			zl.cv = pow(nciw.over(w0),gamma).times(zl.v);
      zl.ck = cw.over(zl.cv);
		  zl.cim = zl.cv.times(zl.d);
			
			// loop over layers last to first to compute all dj and uj
			for (uk=_zones.lowerKey(lk); uk!=null; lk=uk,uk=_zones.lowerKey(lk)) {
				zl = _zones.get(lk);
				zu = _zones.get(uk);

				// Compute propagation matrix
				Cfloat crs = c1.over(zu.cr.plus(c1));
				Cfloat cex = ciw.over(zu.cv).times(zl.z-zu.z);
				Cfloat cep = crs.times(exp(cex));
				Cfloat cem = crs.times(exp(cex.neg()));
				Cfloat b11 = cem;
				Cfloat b22 = cep;
				Cfloat b12 = zu.cr.times(cem);
				Cfloat b21 = zu.cr.times(cep);
		
       	zu.cd = b11.times(zl.cd).plus(b12.times(zl.cu));
       	zu.cu = b21.times(zl.cd).plus(b22.times(zl.cu));

				// Compute upgoing and downgoing parts (with isotropic source)
				if (zl.s!=0.0f) {
        	zu.cd.minusEquals(cem.times(zl.s));
        	zu.cu.plusEquals(cep.times(zl.s));
				}
			}
			int ir = 0;
			for (Float zr: _recievers) {
				uk = _zones.floorKey(zr);
				// If reciever is at interface exactly
        if (!zr.equals(zs) && _zones.containsKey(zr)) {
					zl = _zones.get(zr);
					zu = _zones.get(_zones.lowerKey(zr));
					Cfloat cex = i0.times(zu.ck).times(zl.z-zu.z);
      		Cfloat ced = exp(cex);
      		Cfloat ceu = exp(cex.neg());
					// Average upper and lower wavefields
      		Cfloat csj = zu.cd.times(ced).plus(zu.cu.times(ceu)
											 .plus(zl.cd).plus(zl.cu)).times(0.5f);
      		cseis[ir][iwr] = cseisr[ir][iw] = csj.r;
      		cseis[ir][iwi] = cseisi[ir][iw] = csj.i;
				}
				// Find wavefield at the reciever located in a zone or at surface
				else if (zr.equals(zs) || uk!=null) {
					zu = _zones.get(uk);
					Cfloat cex = i0.times(zu.ck).times(zr-zu.z);
      		Cfloat ced = exp(cex);
      		Cfloat ceu = exp(cex.neg());
      		Cfloat csj = zu.cd.times(ced).plus(zu.cu.times(ceu));
      		cseis[ir][iwr] = cseisr[ir][iw] = csj.r;
      		cseis[ir][iwi] = cseisi[ir][iw] = csj.i;
				}
			++ir;
			}
		}
		/*
		// Anti-alias filter from Dave's code
		  int npole = 12;
      float pio2 = FLT_PI/2.0f;
      float w3db = pio2/dt;
      for (int ipole=0; ipole<npole; ++ipole) {
        Cfloat cpower = new Cfloat(0.0f,-(2*ipole+1)*pio2/npole);
        Cfloat cpole = exp(cpower).times(w3db);
        for (int iw=0,iwr=0,iwi=1; iw<nw; ++iw,iwr+=2,iwi+=2) {
          Cfloat cw = new Cfloat(iw*dw,eps);
          Cfloat ch = cpole.over(cpole.minus(cw));
          for (int ir=0; ir<nr; ++ir) {
            float csr = cseis[ir][iwr];
            float csi = cseis[ir][iwi];
            cseis[ir][iwr] = csr*ch.r-csi*ch.i;
            cseis[ir][iwi] = csr*ch.i+csi*ch.r;
          }
        }
      }*/

		if (_ricker>0.0f) {
			float sigma = sqrt(2.0f)/(2.0f*_ricker*FLT_PI);
			float ascal = sqrt(2.0f*FLT_PI)*sigma;
			float sigex = sigma*sigma/2.0f;
			System.out.format("sigma = %f \n",sigma);
			System.out.format("peak f = %f \n",2.0f/sigma);
    	for (int iw=0,iwr=0,iwi=1; iw<nw; ++iw,iwr+=2,iwi+=2) {
				Cfloat cw = new Cfloat(iw*dw,eps);
				Cfloat cws = cw.times(cw);
				Cfloat exs = neg(cws.times(sigex));
				Cfloat rik = cws.times(ascal).times(exp(exs));
				for (int ir=0; ir<nr; ++ir) {
					cseis[ir][iwr] *= rik.r;
					cseis[ir][iwi] *= rik.i;
				}
			}
			/*
			double fp = 3.0/_ricker;
    	int ih = (int)(3.0/fp);
    	int nh = 1+2*ih;
    	float[] h = new float[nh];
    	for (int jh=0; jh<nh; ++jh) {
				double x = PI*fp*(jh-ih);
    		double xx = x*x;
    		h[jh] = (float)((1.0-2.0*xx)*exp(-xx));
			}
			SimplePlot si = new SimplePlot();
			si.addPoints(h);
 			int nwfft = FftReal.nfftFast((int)(_oversampling*nh));
			int nww = nwfft/2+1;
			int nwpad = nwfft+2; 
    	float[] g = new float[nwpad];
			FftReal fft = new FftReal(nwfft);
			constantExtrap(h,g);
			copy(nh,h,g);
			float[] gr = new float[nww];
			float[] gi = new float[nww];
			fft.realToComplex(1,g,g);
    	for (int iw=0,iwr=0,iwi=1; iw<nww; ++iw,iwr+=2,iwi+=2) {
				for (int ir=0; ir<nr; ++ir) {
					cseis[ir][iwr] *= g[iwr];
					cseis[ir][iwi] *= g[iwi];
					gr[iw] = g[iwr];
					gi[iw] = g[iwi];
					
				}
			}
			Sampling sww = new Sampling(nw,dw,0.0);
			plot(sww,gr,"frequency cycles/s","real part");
			plot(sww,gi,"frequency cycles/s","real part");
			*/
		}

		System.out.format("nw = %d \n",nw);
		System.out.format("nt = %d \n",nt);
		System.out.format("nfft = %d \n",nfft);
		System.out.format("dw = %f \n",dw);
		System.out.format("df = %f \n",df);
		System.out.format("w0 = %f \n",w0);
		System.out.println("f0 = "+f0);
		System.out.println("mf = "+(nw*df));

		
		//plot(sw,cseisr[0],"frequency cycles/s","real part");
		//plot(sw,cseisi[0],"frequency cycles/s","imaginary part");

		// IFFT
		float[][] seis = new float[nr][nt];
		FftReal fft = new FftReal(nfft);
		float[] sfft = new float[nfft];
		float omult = 1.0f/(_mult*nfft);
    for (int ir=0; ir<nr; ++ir) {
      fft.complexToReal(-1,cseis[ir],sfft);
      for (int it=0; it<nt; ++it)
        seis[ir][it] = sfft[it]*exp(eps*it*dt)*omult;
    }
		return seis;
  }

///////////////////////////////////////////////////////////////////////////////////////////
// private

	private TreeMap<Float,Zone> _zones = new TreeMap<Float,Zone>(); // zones
	private TreeSet<Float> _recievers = new TreeSet<Float>(); // recievers
	private float _r1 = 1.0f; // surface reflectivity (set to 0 to remove surface reflection)
	private float _oversampling = 2.0f;
	private float _decay = 0.01f; // decay factor to reduce wraparound after IFFT
	private float _mult = 1.0f;
	private float _ricker=0.0f; // peak frequency for ricker wavelet

	/**
	 * A Zone is the area below and above an interface at depth z.
	 */
	private class Zone {
		public Zone(float z, float v, float d, float q, float s, boolean source) {
			this.z = z;
			this.v = v;
			this.d = d;
			this.q = q;
			this.s = s;
			this.source=source;
		}
		float z; // Depth
		float v; // Velocity
		float d; // Density
		float q; // Attenuation factor
		float s; // Source amplitude
		Cfloat cd; // Downgoing wave
		Cfloat cu; // Upgoing wave
		Cfloat ck; // Complex wavenumber
		Cfloat cr; // Complex reflectivity
		Cfloat cv; // Complex velocity
		Cfloat cim; // Complex impedance
		boolean source; // source layer
	}

  private void checklengths(
    int n, float[] v, float[] d, float[] q) 
  {
    Check.argument(v.length==n,"v length == nz");
    Check.argument(d.length==n,"d length == nz");
    Check.argument(q.length==n,"q length == nz");
  }

	private static void constantExtrap(float[] x, float[] z) {
	  int nx = x.length;
		int nz = z.length;
		int nm = nz-nx;
		int nmh = nx+(int)(nm/2.0);
		float x0 = x[0];
		float xn = x[nx-1];
	  for (int i=nx; i<nmh; ++i)
	    z[i] = xn;
	  for (int i=nmh; i<nz; ++i)
	    z[i] = x0;
  }

  
///////////////////////////////////////////////////////////////////////////////////////////
// testing


  public static void main(String[] args) {
    test1();
	}

  public static void test1() {
	  // Length of Logs
		int nz = 2;
		float dz = 5.0f;//0.1524f;
		float fz = 0.0f; 
		Sampling sz = new Sampling(nz,dz,fz);
		// Length of Synthetic Seismogram
		int nt = 101;
		double dt = 1.0;
		double ft = 0.0;
		Sampling st = new Sampling(nt,dt,ft);
		// Reflection coefficient at surface
		float r1 = 1.0f; 
		// Create logs
		SyntheticSeismogram ssm = new SyntheticSeismogram();
		ssm.addZone(0.0f,1.0f,1.0f,1.0e6f);
		ssm.addZone(20.0f,2.0f,2.0f,10);
		ssm.addSource(10.0f,1.0f);
		ssm.addReciever(10.0f);
		float[][] seis = ssm.getSeismogram(st);
		// Plot
		plot(st,seis[0],"time (s)","synthetic");

		// Spectrum
		//Spectrum spc = new Spectrum();
		//spc.apply(st,seis,"Seismic");
	}
	
	private static void plot(Sampling s, float[] x, String vaxis, String haxis) {
	  SimplePlot si = new SimplePlot(SimplePlot.Origin.UPPER_LEFT);
		si.addPoints(s,x);
		si.setSize(500,1000);
		si.setHLabel(haxis);
		si.setVLabel(vaxis);
	}


//////////////////////////////////////////////////////////////////////////////////////////////////////
// Spectrum Plotting Class

	private static class Spectrum {
	
	   public static void apply(Sampling xs, float[] x, String title) {
	
	    // Time sampling.
	    int nt = xs.getCount();
	    double dt = xs.getDelta();
	    double ft = xs.getFirst();
	 
	    // Frequency Sampling
	    int nfft = FftReal.nfftSmall(2*nt);
	    int nf = nfft/2 + 1;
	    double df = 1.0/(nfft*dt);
	    double ff = 0.0;
	    Sampling fs = new Sampling(nf,df,ff);
	
	    // Fourier-Transform
	    FftReal fft = new FftReal(nfft);
	    int nf2 = nf*2;
	    float[] y = new float[nf2];
	    copy(nt,x,y);
	    fft.realToComplex(-1,y,y);
	
	    // Adjust phase for possibly non-zero time of first sample.
	    float[] wft = rampfloat(0.0f,-2.0f*FLT_PI*(float)(df*ft),nf);
	    y = cmul(y,cmplx(cos(wft),sin(wft)));
	
	    // Amplitude Spectrum, normalized
	    float[] af = cabs(y);
	    float amax = max(max(af),FLT_EPSILON);
	    af = mul(1.0f/amax,af);
	    
	    double fmax = 0;
	    for (int k=0; k<nf; ++k) {
	      if (af[k]>=0.15) {
		  fmax = ff+k*df;
	      }
	    }
	    fmax *= 1.2;
	    
	    // Peak frequency
	    double pfr = 0;
	    int kmx = 0;
	    float[] af2 = new float[nf];
	    copy(af,af2);
	    //ExponentialSmoother es = new ExponentialSmoother(0.5f);
	    //es.apply(af,af2);
	    for (int k=0; k<nf; ++k) {
	      if (af2[kmx]<af2[k]) {
	        pfr = ff+k*df;
		kmx=k;
	      }
	    }
	    int rt = (int)(nt/(2*fmax));
	
	    //plot(fs,af,fmax,"f",title);
	    //plot(fs,af2,fmax,"f",title);
	    //System.out.format("max f of "+title+" is: %f\n",fmax);
	    //System.out.format("peak f of "+title+" is: %f\n",pfr);
	
	    // Phase spectrum, in cycles
	    float[] pf = carg(y);
	    pf = mul(0.5f/FLT_PI,pf);
	    plot(fs,pf,fmax,"p",title);
	
	    plotPanels(fs,af,pf,fmax,title);
	    //plotPanels(fs,af2,pf,fmax,title+" smoothed");
	  }
	
	//////////////////////////////////////////////////////
	// private for Spectrum
	
	  private static void plotPanels(Sampling fs, float[] f, float[] p, double fmax, String title) {
	    PlotPanel pp = new PlotPanel(2,1);
	    PointsView fl = pp.addPoints(0,0,fs,f);
	    PointsView pl = pp.addPoints(1,0,fs,p);
	    pp.setVLabel(0,"Amplitude");
	    pp.setVLabel(1,"Phase (cycles)");
	    pp.setHLimits(0.0,fmax);
	    pp.setVLimits(0,0.0,1.0);
	    pp.setVLimits(1,-0.5,0.5);
	    pp.setHLabel("Frequency (cycles/sample)");
	    pp.setTitle(title);
	    PlotFrame frame = new PlotFrame(pp);
	    frame.setSize(900,1200);
	    frame.setDefaultCloseOperation(PlotFrame.EXIT_ON_CLOSE);
	    frame.setVisible(true);
	  }
	    
	
	  private static void plot(Sampling xs, float[] x, double fmax, String axis, String title) {
	    SimplePlot sp = new SimplePlot();
	    sp.addPoints(xs,x);
	    if (axis=="t") {
	      sp.setHLabel("Time (ms)");
	      sp.setVLabel("Amplitide");
	    }
	    else if (axis=="f") {
	      sp.setHLabel("Frequency (cycles/sample)");
	      sp.setVLabel("Amplitide");
	      sp.setHLimits(0.0,fmax);
	    }
	    else if (axis=="d") {
	      sp.setHLabel("Depth (m)");
	      sp.setVLabel("Amplitide");
	    }
	    else if (axis=="p") {
	      sp.setHLabel("Frequency (cycles/sample)");
	      sp.setVLabel("Phase (cycles)");
	      sp.setVLimits(-0.5,0.5);
			}
		}
	}
};

