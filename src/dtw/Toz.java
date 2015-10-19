package dtw;

/**
 * @author Andrew Munoz, Colorado School of Mines
 * @version 2012.09.06
 */

import edu.mines.jtk.dsp.*;
import edu.mines.jtk.interp.CubicInterpolator;
import edu.mines.jtk.mosaic.SimplePlot;
import static edu.mines.jtk.util.ArrayMath.*;

public class Toz {

  /**
   * Computes t(z) for functions t1 and t2 with samplings. 
   *
   * @param u The shift field from lags computed from DTW
   * @param itoz The initial t(z) function typically from logs
   * @param s1 Sampling of the longer sequence (g for well ties)
   * @param s2 Sampling of the shorter sequence (f for well ties)
   * @return t(z) updated from DTW shifts
   */
  public float[] computeToz(
    float[] u, float[] itoz, Sampling s1, Sampling s2) 
  {
    int n1 = s1.getCount();
    int n2 = s2.getCount();
    int nz = itoz.length;
    float d1 = (float)s1.getDelta();
    float d2 = (float)s2.getDelta();
    float f1 = (float)s1.getDelta();
    float f2 = (float)s2.getDelta();
    
    float[] tt = add(rampfloat(f2,d2,n2),u);
    float[] tt2 = new float[n2];
    tt2[0]=tt[0];
    int j=0;
    for (int i=1; i<n2; ++i) {
      if (tt[i]<=tt[i-1]) tt2[i]=-1;
      else {tt2[i]=tt[i]; ++j;}}
    float[] t2 = new float[j];
    for (int i=1,i2=0; i<n2; ++i)
      if (tt2[i]>0) {t2[i2]=tt2[i]; ++i2;}
    t2[0]=tt[0];
    float[] t1 = rampfloat(f2,d2,j);
    
    for (int i=0; i<j; ++i)
      System.out.format("t1 is %f t2 is %f \n",t1[i],t2[i]);
    
    float[] t1ot2 = new float[n2];
    CubicInterpolator ct1ot2 = new CubicInterpolator(CubicInterpolator.Method.LINEAR,n2,t2,t1);
    for (int it2=0; it2<n2; ++it2) 
      t1ot2[it2] = ct1ot2.interpolate(f2+it2*d2);
    
    float[] toz = new float[nz];
    LinearInterpolator li = new LinearInterpolator();
    li.setUniform(n2,d2,f2,t1ot2);
    for (int iz=0; iz<nz; ++iz) 
      toz[iz] = li.interpolate(itoz[iz]);
   
   //For testing:
    SimplePlot si = new SimplePlot();
    si.addPoints(toz);
    si.setTitle("t(z)");
    si.setSize(300,1000);
    
    SimplePlot si2 = new SimplePlot();
    si2.addPoints(t1ot2);
    si2.setTitle("t1(t2)");
    si2.setSize(300,1000);
    
    return toz;
  }

};
