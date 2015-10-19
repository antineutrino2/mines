package anm;

import edu.mines.jtk.sgl.*;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import static edu.mines.jtk.util.ArrayMath.*;

// THANKS TO STEFAN COMPTON FOR PROVIDING THIS

public class TriangleGroupArray {

  public static TriangleGroup[] makeArray(
      float[][][] f, Sampling sx, Sampling sy, double[] clips) 
  {
    int n3 = f.length;
    //System.out.println("n3="+n3);
    int n2 = f[0].length;
    int n1 = f[0][0].length;
    TriangleGroup[] tga = new TriangleGroup[n3];
    ColorMap cmap = new ColorMap(ColorMap.JET);
    cmap.setValueRange(clips[0],clips[1]);
    for (int i3=0; i3<n3; i3++) {
      float[] rgb = cmap.getRgbFloats(flatten(f[i3]));
      int p = 0;
      float[][] r = zerofloat(n1,n2);
      float[][] g = zerofloat(n1,n2);
      float[][] b = zerofloat(n1,n2);
      for (int i2=0; i2<n2; i2++) {
        for (int i1=0; i1<n1; i1++) {
          r[i2][i1] = rgb[ p ];
          g[i2][i1] = rgb[p+1];
          b[i2][i1] = rgb[p+2];
          p += 3;
        }
      }
      TriangleGroup tg = new TriangleGroup(true,sx,sy,f[i3],r,g,b);
      tg.setColorMap(cmap,clips[0],clips[1]);
      tga[i3] = tg;
    }
    return tga;
  }
  
}
