/** 
 * Computes nmo equations for grid testing.
 * For further details, refer to:
 *  Brahim Abbad, B. U., and D. Rappin, 2009, Automatic
 *    nonhyperbolic velocity analysis: Geophysics, 74, U1-U12.
 * @author Andrew Munoz, CSM
 * @version 02.17.14
 */

package sem;

import java.util.Random;

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import static edu.mines.jtk.util.ArrayMath.*;

public class NmoSamplingTest {

  public NmoSamplingTest(
    Sampling s1, Sampling sx, 
    float vmin, float vmax, 
    float emin, float emax)
  {
    this.s1 = s1;
    this.sx = sx;
    n1 = s1.getCount();
    nx = sx.getCount();
    d1 = s1.getDelta();
    dx = sx.getDelta();
    f1 = s1.getFirst();
    fx = sx.getFirst();
    this.vmin = vmin;
    this.vmax = vmax;
    this.emin = vmin;
    this.emax = emax;
    xmax = (float)(abs(sx.getFirst())>abs(sx.getLast())?
                              sx.getFirst():sx.getLast());
  }

  public void nmoVnVh(float dv, float dh) {
    float vhmin = vmin*sqrt(1f+2f*emin);
    float vhmax = vmax*sqrt(1f+2f*emax);
    int nv = round((vmax-vmin)/dv);
    int nh = round((vhmax-vhmin)/dh);
    int na = nv*nh;
    float[] gvn = new float[na];
    float[] gvh = new float[na];
    float[] gdtv = new float[na];
    float[] gdth = new float[na];
    int i1=0,i2=0,i3=0,i4=0;
    for (int ih=0; ih<nh; ++ih) {
      for (int iv=0; iv<nv; ++iv) {
        gvn[++i1]  = vmin+iv*dv;
        gvh[++i2]  = vhmin+ih*dh;
        gdtv[++i3] = 1;
        gdth[++i4] = 1;


      }
    }
    int n1v = n1*nv;
    int n1h = n1*nh;
    float[] gtvn = new float[n1v];
    float[] gtvh = new float[n1h];
    float[] gtdtvn = new float[n1v];
    float[] gtdtvh = new float[n1h];

  }

  public void nmoVnE(float dv, float de) {

  }

  public void nmoDtvDth(float dtv, float dth) {
    float vhmin = vmin*sqrt(1f+2f*emin);
    float vhmax = vmax*sqrt(1f+2f*emax);

  }

  public void nmoDtvDte(float dtv, float dte) {

  }

  public void nmoDtvDthS(float dtvs, float dths) {
    float vhmin = vmin*sqrt(1f+2f*emin);
    float vhmax = vmax*sqrt(1f+2f*emax);

  }

  public void nmoDtvDteS(float dtvs, float dths) {

  }


  public void nmoQ1Q2(float dq1, float dq2) {

  }


////////////////////////////////////////////////////////////////////////////////
// private

  private Sampling s1,sx;
  private int n1,nx;
  private float vmin,vmax,emin,emax,xmax;
  private double d1,dx,f1,fx;



};

