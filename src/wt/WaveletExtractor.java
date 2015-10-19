/**
 * A statistical wavelet extractor.
 * Extracts a wavelet statistically from anywhere in the seismic 
 * using autocorrelation method and writes the wavelet to file for future access. 
 * Only returns a zero-phased wavelet.
 *
 * @author Andrew Munoz, Colorado School of Mines
 * @version 2011.10.27
 */

package wt;

import edu.mines.jtk.dsp.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class WaveletExtractor {

  /**
   * Get a statistical Wavelet from the given traces and in a certian sample window
   * @param wmin window sample minimum 
   * @param wmax window sample maximum
   * @param wvl wavelet sample length
   * @param trace trace in terms of trace number and amplitude with depth
   * @return wavelet 
   */

  public float[] getWavelet(
		int wvl, int wmin, int wmax, float[][] trace) 
	{ 
    int ns = trace[0].length; // Number of z samples
    int nt = trace.length; // Number of Traces
    int tpl = 10; // Taper sample length
    float[][] window = new float[nt][ns];
    float[][] winwav = new float[nt][wvl];
    float[]  wavelet = new float[wvl];
		Fft fft = new Fft(wvl);
    for (int it=0; it<nt; ++it) { // Increment the traces
      for (int is=wmin; is<wmax; ++is) {
        window[it][is] = trace[it][is];
				// Taper the window
				if (10>((wmax-wmin)/4)) tpl = (wmax-wmin)/4;
				if (is<tpl && (tpl-is)>0) window[it][is] *= 1/(tpl-is);
				if (is>(wmax-tpl) && (is-(wmax-tpl))>0) window[it][is] *= 1/(is-(wmax-tpl));
      }
      // Compute the autocorrelation of the window- length of ac is 1/2 wavelet length
      Conv.xcor(ns,0,window[it],ns,0,window[it],(wvl/2),0,winwav[it]);
      // Get the amplitude spectrum
      winwav[it] = fft.applyForward(winwav[it]);
			winwav[it] = csqrt(winwav[it]);
			winwav[it] = fft.applyInverse(winwav[it]); 
    }
    // Average the wavelets from each trace
    for (int it=0; it<nt; ++it) {
    	for (int is=0; is<wvl; ++is) 
        wavelet[is] += winwav[it][is]; 
    }
		RecursiveGaussianFilter rgf = new RecursiveGaussianFilter(1);
		rgf.apply0(wavelet,wavelet);
		float wx = max(wavelet);
		wavelet = mul(wavelet,1.0f/wx);
    // Make the wavelet centered and symmetrical
    for (int is=0; is<(wvl/2); ++is) 
      wavelet[is+(wvl/2)] = wavelet[is];
    for (int is=0; is<(wvl/2); ++is) 
      wavelet[is] = wavelet[wvl-is-1];
    return wavelet;
  }

};

