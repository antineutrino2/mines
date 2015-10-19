package wt;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.Parallel.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Computes and plots amplitude spectrum.
 * @author Andrew Munoz, Colorado School of Mines
 * @version 2011.11.15
 */


public class Spectrum {

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
    ExponentialSmoother es = new ExponentialSmoother(0.5f);
    es.apply(af,af2);
    for (int k=0; k<nf; ++k) {
      if (af2[kmx]<af2[k]) {
        pfr = ff+k*df;
	kmx=k;
      }
    }
    int rt = (int)(nt/(2*fmax));

    //plot(fs,af,fmax,"f",title);
    //plot(fs,af2,fmax,"f",title);
    System.out.format("max f of "+title+" is: %f\n",fmax);
    System.out.format("peak f of "+title+" is: %f\n",pfr);
    System.out.format("decimation of "+title+" is: %d\n",rt);

    // Phase spectrum, in cycles
    float[] pf = carg(y);
    pf = mul(0.5f/FLT_PI,pf);
    //plot(fs,pf,fmax,"p",title);

    plotPanels(fs,af,pf,fmax,title);
    plotPanels(fs,af2,pf,fmax,title+" smoothed");
  }

//////////////////////////////////////////////////////
// private

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
      sp.setHLimits(0.0,fmax);
    }
    else {
      sp.setHLabel("Samples");
      sp.setVLabel("Amplitide");
    }
    sp.setTitle(title);
    sp.setSize(700,500);
  }
  
}; 
