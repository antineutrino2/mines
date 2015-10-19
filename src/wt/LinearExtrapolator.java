package wt;

import edu.mines.jtk.util.*;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.SincInterp;
import edu.mines.jtk.dsp.SincInterp.Extrapolation;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Computes a linear extrapolation of a 1D curve in a least squares sense. 
 * 
 * @author Andrew Munoz, Colorado School of Mines
 * @version 11.23.13
 */

public class LinearExtrapolator {

  public LinearExtrapolator(int m) {
    _m = m;
  }

  public float[] extrapolateFirst(
    double fx, double fy, 
    Sampling sx, double[] y) 
  {
    return extrapolateFirst(sx.getDelta(),fx,fy,sx.getValues(),y);
  }
  
  public float[] extrapolateLast(
    double lx, double ly, 
    Sampling sx, double[] y) 
  {
    return extrapolateLast(sx.getDelta(),lx,ly,sx.getValues(),y);
  }

  public float[] extrapolateFirst(
    double dx,
    double fx, double fy, 
    double[] x, double[] y) 
  { 
    double lx = x[x.length-1];
    double ly = y[y.length-1];
    return extrapolate(dx,fx,fy,lx,ly,x,y);
  }

  public float[] extrapolateLast(
    double dx,
    double lx, double ly, 
    double[] x, double[] y) 
  { 
    double fx = x[0];
    double fy = y[0];
    return extrapolate(dx,fx,fy,lx,ly,x,y);
  }

  public float[] extrapolate(
    double dx,
    double fx, double fy, 
    double lx, double ly, 
    double[] x, double[] y) 
  {
    int nb = x.length;
    Check.argument(nb==y.length,"x and y lengths must be equal");
    int nf =    (int)floor(fx/dx);
    int nl = nb-(int)ceil(lx/dx);
    int ne = nf+nb+nl;
    double[] xn = new double[ne]; 
    
    return new float[nb];
  }



  private int _m;

};
