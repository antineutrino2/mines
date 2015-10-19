package tp;

public class LogSamples {

  public float[] vl,dl,gl;
  public float xl,yl,zl;

  public LogSamples() {
    } 

  /**
   * Gets all samples f(x1,x2,x3) for the specified well log curve.
	 * For null samples, the values are a constant exrapolation before and after
	 * the first and last recorded samples. If there are null samples between the
	 * first and last sample, we just take an average of nearby samples.
   * @return array of arrays {v,d,g,x1,x2,x3} of all samples 
   */
	public void fixLogsForSynthetics(WellLog wlog) {
		int n = wlog.n;
		float[] v = wlog.v;
		float[] d = wlog.d;
		float[] g = wlog.g;
		float[] x1 = wlog.x1;
		float[] x2 = wlog.x2;
		float[] x3 = wlog.x3;
		float NULL_VALUE = wlog.NULL_VALUE;
		float zf = 0.0f; float yf = 0.0f; float xf = 0.0f; 
		int kf=0; int kl=0;
		boolean first = false;
		boolean last = false;
		// find first and last samples of depth and velocity
		for (int i=0,j=n-1; i<n; ++i,--j) {
      float vi = v[i];
      float di = d[i];
      if (vi!=NULL_VALUE && di!=NULL_VALUE && first==false) {
				zf = x1[i];
				yf = x2[i];
				xf = x3[i];
				kf = i;
				first = true;
			}
      float vj = v[j];
      float dj = d[j];
      if (vj!=NULL_VALUE && dj!=NULL_VALUE && last==false) {
				kl = j+1;
				last = true;
			}
		}
    kl -= 20; // remove last 20 samples
		int nn = kl-kf;
		float[] ve = new float[nn];
		float[] de = new float[nn];
		float[] ge = new float[nn];
		for (int i=kf,j=0; i<kl; ++i,++j) {
			ve[j] = v[i];
			de[j] = d[i];
			ge[j] = g[i];
		}
		vl = denull(ve,NULL_VALUE);
		dl = denull(de,NULL_VALUE);
		gl = denull(ge,NULL_VALUE);
		xl = xf;
		yl = yf;
		zl = zf;
	}
	private static float[] denull(float[] x, float missing) {
    // interpolates missing values
		int nt = x.length;
    float[] y = new float[nt];
		int i=1; boolean kill=false;
		for (int n=0; n<nt; ++n) {
      float x0 = x[n];
      if (x0!=missing) {
        y[n] = x0;
      } 
			else {
				while (x[n+i]==missing) {
					++i; 
					if (n+i>=nt) {
						kill = true;
						break;
					}
				}
				if (kill) {
					y[n] = y[n-1];
					break;
				}
				float x1 = x[n+i];
				if (n>0) {
	 			 	float xm = y[n-1];
	  			y[n] = xm + (x1-xm)/((i+n)-(n-1));
				}
				else { 
        	float x2 = x[i+1]; 
	  			y[n] = x1 - i*(x2-x1); 
				}
			i=1;
     	}
		}
		return y;
	}
};
