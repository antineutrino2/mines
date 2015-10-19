package wt;
/**
 * Decimates a 1D signal.
 * Finds the nyquist frequency and cuts the high frequency samples
 * in the Fourier Domain.
 * @author Andrew Munoz, Colorado School of Mines
 * @version 13.01.2012
 */
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.mosaic.SimplePlot;
import static edu.mines.jtk.util.Parallel.*;
import static edu.mines.jtk.util.ArrayMath.*;


public class Decimator {

  public static int apply(Sampling xs, float[] x) {
    int n = x.length;

    double dt = xs.getDelta();
    
    // Frequency Sampling
    int nfft = FftReal.nfftSmall(2*n);
    int nf = nfft/2 + 1;
    double df = 1.0/(nfft*dt);
    double ff = 0.0;
    Sampling fs = new Sampling(nf,df,ff);

    // Fourier-Transform
    FftReal fft = new FftReal(nfft);
    int nf2 = nf*2;
    float[] y = new float[nf2];
    copy(x,y);
    fft.realToComplex(-1,y,y);

    // Amplitude Spectrum, normalized
    float[] af = cabs(y);
    float amax = max(max(af),FLT_EPSILON);
    af = mul(1.0f/amax,af);
    
    int ny = 0;
    for (int k=0; k<nf; ++k) {
      if (af[k]>=0.15) {
	  ny = k;
      }
    }
    float fmax = ny*1.2f;
    //plot(af,"f");
    System.out.println("max k is: "+ny);
    
    // Decimate x
    int rt = (int)(n/(2*fmax));

    /*
    int n2 = (int)(n/rt);
    float[] x2 = new float[n2];
    for (int i=0,j=0; i<n2; ++i,j+=rt)
      x2[i] = x[j];
    double srn = rt*xs.getDelta();
    double vup = 4*srn/(0.001);
    double vlw = srn/(0.001);
    System.out.println("m is: "+rt);
    System.out.format("new n is: %d %n%n",n2);
    System.out.format("new dz is: %f km %n",srn);
    System.out.format("Velocity bounds are %f to %f km/s %n%n",vlw,vup);
    //plot(x2,"d"); 
    //plot(x,"d");
    //return x2;
    */
    return rt;
  }

//////////////////////////////////////////////////////
// private

  private static void copy(final float[] x, final float[] y) {
    final int n1 = x.length;
    loop(n1,new LoopInt() {   
    public void compute(int i1) {
       y[i1] = x[i1];
    }});
  }
  private static void copy2(final float[] x, final float[] y) {
    final int n1 = y.length;
    loop(n1,new LoopInt() {   
    public void compute(int i1) {
       y[i1] = x[i1];
    }});
  }


  private static void plot(float[] x, String axis) {
    SimplePlot sp = new SimplePlot();
    sp.addPoints(x);
    if (axis=="t") {
      sp.setHLabel("Time (ms)");
      sp.setTitle("Trace Amplitudes");
    }
    else if (axis=="f") {
      sp.setHLabel("Frequency (cycles/sample)");
      sp.setTitle("Frequency-Amplitude spectrum");
    }
    else if (axis=="d") {
      sp.setHLabel("Depth (m)");
      sp.setTitle("Well Log in depth");
    }
    else {
      sp.setHLabel("Samples");
    }
    sp.setVLabel("Amplitide");
    sp.setSize(700,500);
  }

//////////////////////////////////////////////////////
// testing
/*
  public static void main(String[] args) {
    // Make the function
    int n = 100;
    float[] x = new float[n];
    float pi = FLT_PI;
    for (int i=0; i<n; ++i) 
      //x[i] = sin(i*pi/10f)/(i*pi/10f);
      //x[i] = exp(-(i*i)/500f);
      x[i] = sin(i*pi/10f);

    // Interpolate
    SincInterpolator si = new SincInterpolator();
    si.setExtrapolation(SincInterpolator.Extrapolation.CONSTANT);
    si.setUniform(n,1.0,0.0,x);
    int c = 100;
    float oc = 1f/c;
    float[] y = new float[n*c];
    for (int i=0; i<n*c; ++i)
      y[i] = si.interpolate(i*oc);
   
    // Apply the fft
    //Sampling xs = new Sampling(n);
    //apply(xs,x,1.0f);
    Sampling ys = new Sampling(n*c);
    apply(ys,y,1.0f);
    
    
    // Decimation Test
    int ny = 26;
    float bu = 2f;
    int rf = (int)(n*c/(2*ny*bu));
    int n2 = (int)(n*c/rf);
    float[] y2 = new float[n2];
    for (int i=0,j=0; i<n2; ++i,j+=rf)
      y2[i] = y[j];
    Sampling ys2 = new Sampling(n2);
    apply(ys2,y2,1.0f);
    
    }
    */

};
