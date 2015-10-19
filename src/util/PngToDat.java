/** @author Luming Liang
 */

package util;

import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.*;
import javax.imageio.ImageIO;

import edu.mines.jtk.mosaic.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.io.ArrayOutputStream;

import static edu.mines.jtk.util.ArrayMath.*;

public class PngToDat {
  public PngToDat() {
  }
  public float[][] getData(String fName) {
    File pngRef = new File(fName);
    BufferedImage imageRef; 
    try {
      imageRef = ImageIO.read(pngRef); 
    } catch (IOException e) {
      throw new RuntimeException(e);
    }
    int n2 = imageRef.getWidth();
    int n1 = imageRef.getHeight();
    float[][] image = new float[n2][n1];
    int r, g, b, rgb;
    for (int i2=0; i2<n2; ++i2)
      for (int i1=0; i1<n1; ++i1) {
        rgb = imageRef.getRGB(i2,i1);
        r = (rgb>>16)&255;
        g = (rgb>>8)&255;
        b = (rgb)&255;
        image[i2][i1] = (r+g+b)/255.0f/3;
      }
    SimplePlot sp = SimplePlot.asPixels(image);
    sp.setSize(new java.awt.Dimension(n2+100,n1+100));
    sp.addColorBar();
    System.out.println(n2+"  "+n1);
    float factor = 500.0f/min(n2,n1);
    if (factor<=1) return image;
    int nn2 = (int)(n2*factor);
    int nn1 = (int)(n1*factor);
    SincInterpolator si = new SincInterpolator();
    float[][] y = new float[nn2][nn1];
    for (int i2=0; i2<nn2; ++i2)
      for (int i1=0; i1<nn1; ++i1) 
        y[i2][i1] = si.interpolate(n1,1.0f,0.0f,n2,1.0f,0.0f,image,i1*n1/nn1,i2*n2/nn2);
    // regularize
    //float regF = 1.0f/max(abs(y));
    //y = mul(y,-regF);
    SimplePlot spy = SimplePlot.asPixels(y);
    spy.addColorBar();
    spy.setSize(new java.awt.Dimension(nn2+100,nn1+100));
    return y;
  }

  // Write image
  public static void WriteImg(float[][] image, String name) {
    try {
      ArrayOutputStream aos = new ArrayOutputStream(name+".dat");
      aos.writeFloats(image);
      aos.close();
      }	catch(Exception e){
      e.printStackTrace();
    }
  }

  public static void main(String[] args) {
    PngToDat ptd = new PngToDat();
    float[][] image = ptd.getData("andrew2.jpg");
    WriteImg(image, "andrew"+image.length+"_"+image[0].length);
  }
}

