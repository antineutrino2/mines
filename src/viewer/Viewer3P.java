package viewer;

import java.awt.Color;
import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.KeyEvent;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.event.MouseMotionAdapter;
import java.awt.event.MouseMotionListener;
import java.awt.image.IndexColorModel;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;

import javax.swing.JMenuItem;
import javax.swing.KeyStroke;

import edu.mines.jtk.awt.Mode;
import edu.mines.jtk.awt.ModeManager;
import edu.mines.jtk.awt.ModeMenuItem;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.io.ArrayInputStream;
import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.ContoursView;
import edu.mines.jtk.mosaic.PlotPanelPixels3.Orientation;
import edu.mines.jtk.mosaic.PlotPanel;
import edu.mines.jtk.mosaic.Tile;
import edu.mines.jtk.mosaic.PlotPanelPixels3;
import edu.mines.jtk.mosaic.PlotPanelPixels3.AxesPlacement;
import edu.mines.jtk.mosaic.PointsView;
import static edu.mines.jtk.util.ArrayMath.*;
import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.SimpleFloat3;

/**
 * Wraps PlotPanelPixels3 with some convenient options. 
 */
public class Viewer3P {

  /**
   * Construct a 3 panel view of f, with default orientation,
   * {@link Orientation#X1DOWN_X2RIGHT} and default Samplings.
   * @param f
   */
  public Viewer3P(float[][][] f) {
    this(f,null);
  }
  
  /**
   * Construct a 3 panel view of f, with specified orientation
   * and default Samplings.
   * @param f
   */
  public Viewer3P(float[][][] f, Orientation o) {
    this(null,null,null,f,o);
  }
  
  /**
   * Construct a 3 panel view of f, with default orientation,
   * {@link Orientation#X1DOWN_X2RIGHT} and specified Samplings.
   * @param f
   */
  public Viewer3P(Sampling s1, Sampling s2, Sampling s3, float[][][] f) {
    this(s1,s2,s3,f,null);
  }
  
  /**
   * Construct a 3 panel view of f, with specified orientation
   * and Samplings.
   * @param f
   */
  public Viewer3P(
      Sampling s1, Sampling s2, Sampling s3, float[][][] f, 
      Orientation orientation)
  {
    _s1 = (s1==null)?new Sampling(f[0][0].length):s1;
    _s2 = (s2==null)?new Sampling(f[0].length   ):s2;
    _s3 = (s3==null)?new Sampling(f.length      ):s3;
    _orientation = (orientation==null )?Orientation.X1DOWN_X2RIGHT:orientation;
    if (orientation==Orientation.X1DOWN_X2RIGHT) {
      _transpose23 = false;
    } else if (orientation==Orientation.X1DOWN_X3RIGHT) {
      _transpose23 = true;
    } else if (orientation==Orientation.X1RIGHT_X2UP) {
      _transpose23 = true;
    } else if (orientation==Orientation.X1RIGHT_X3UP) {
      _transpose23 = false;
    }
    _n1 = _s1.getCount();
    _n2 = _s2.getCount();
    _n3 = _s3.getCount();
    _k1 = _n1/2;
    _k2 = _n2/2;
    _k3 = _n3/2;
    _l1Min = _s1.getFirst();
    _l2Min = _s2.getFirst();
    _l3Min = _s3.getFirst();
    _l1Max = _s1.getLast();
    _l2Max = _s2.getLast();
    _l3Max = _s3.getLast();
    // Make initial panel, displaying the middle frame.
    _pp = new PlotPanelPixels3(_orientation,AxesPlacement.LEFT_BOTTOM,
        _s1,_s2,_s3,f);
    _pp.setSlices(_k1,_k2,_k3);
    _pp.getMosaic().setHeightElastic(0,100);
    _pp.getMosaic().setHeightElastic(1,200);
    _pv1 = new PixelsView[]{
        _pp.getPixelsView12(),
        _pp.getPixelsView13(),
        _pp.getPixelsView23()
    };
    _vf = new ViewerFrame(_pp,_pv1);
    SliceMode sm = new SliceMode(_vf.getModeManager());
    _vf.addToMenu(new ModeMenuItem(sm));
    
    // Add SliceFrame dialog to options menu.
    JMenuItem changeSlices = new JMenuItem("Change Slices");
    final SliceFrame sf = new SliceFrame(this);
    changeSlices.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        sf.getFrame().setVisible(true);
      }
    });
    _vf.addToMenu(changeSlices);
    
    // Add LimitsFrame3P dialog to the options menu.
    JMenuItem changeLimits = new JMenuItem("Change Limits");
    final LimitsFrame3P lf3p = new LimitsFrame3P(this);
    changeLimits.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        lf3p.getFrame().setVisible(true);
      }
    });
    _vf.addToMenu(changeLimits);
  }
  
  public Sampling getSampling1() {
    return _s1;
  }
  public Sampling getSampling2() {
    return _s2;
  }
  public Sampling getSampling3() {
    return _s3;
  }
  public double getL1Min() {
    return _l1Min;
  }
  public double getL2Min() {
    return _l2Min;
  }
  public double getL3Min() {
    return _l3Min;
  }
  public double getL1Max() {
    return _l1Max;
  }
  public double getL2Max() {
    return _l2Max;
  }
  public double getL3Max() {
    return _l3Max;
  }
  public int getN1() {
    return _n1;
  }
  public int getN2() {
    return _n2;
  }
  public int getN3() {
    return _n3;
  }
  public void setSlices(int k1, int k2, int k3) {
    _pp.setSlices(k1,k2,k3);
    _k1 = k1;
    _k2 = k2;
    _k3 = k3;
  }
  
  public int[] getSlices() {
    return new int[]{_k1,_k2,_k3};
  }
  
  public void updateSlices(int k1, int k2, int k3) {
    _pp.setSlices(k1,k2,k3);
    setSlicesPixels(k1, k2, k3,_sf3);
    setSlicesPoints(k1,k2,k3);
    _k1 = k1;
    _k2 = k2;
    _k3 = k3;
  }
  
  /**
   * Sets the height elastic for the specified row. If extra 
   * height is available in this mosaic, it is allocated to 
   * the specified row of tiles in proportion to the specified 
   * height elastic. For fixed-height rows, the height elastic 
   * should be zero. The default height elastic is 100.
   * @param irow the row index.
   * @param heightElastic the height elastic.
   */
  public void setHeightElastic(int irow, int heightElastic) {
    _pp.getMosaic().setHeightElastic(irow, heightElastic);
  }
  
  /**
   * Sets the height minimum for the specified row. All tiles in the 
   * specified row will have height not less than the specified minimum. 
   * Height minimums are used to compute the preferred height of this mosaic.
   * The default height minimum is 100.
   * @param irow the row index.
   * @param heightMinimum the height minimum.
   */
  public void setHeightMinimum(int irow, int heightMinimum) {
    _pp.getMosaic().setHeightMinimum(irow, heightMinimum);
  }
  
  /**
   * Sets the width elastic for the specified column. If extra width is 
   * available in this mosaic, it is allocated to the specified column 
   * of tiles in proportion to the specified width elastic. 
   * For fixed-width columns, the width elastic should be zero.
   * The default width elastic is 100.
   * @param icol the column index.
   * @param widthElastic the width elastic.
   */
  public void setWidthElastic(int icol, int widthElastic) {
    _pp.getMosaic().setWidthElastic(icol, widthElastic);
  }
  
  /**
   * Sets the width minimum for the specified column. All tiles in the 
   * specified column will have width not less than the specified minimum. 
   * Width minimums are used to compute the preferred width of this mosaic.
   * The default width minimum is 100.
   * @param icol the column index.
   * @param widthMinimum the width minimum.
   */
  public void setWidthMinimum(int icol, int widthMinimum) {
    _pp.getMosaic().setWidthMinimum(icol, widthMinimum);
  }
  
  public void addPixels(float[][][] f) {
    Check.argument(f.length==_s3.getCount(),
        "f.length is not consistent with sampling");
    Check.argument(f[0].length==_s2.getCount(),
        "f[0].length is not consistent with sampling");
    Check.argument(f[0][0].length==_s1.getCount(),
        "f[0][0].length is not consistent with sampling");
    _sf3 = new SimpleFloat3(f);
    PixelsView p212 = _pp.addPixels(1,0,_s1,_s2,slice12(_sf3));
    PixelsView p213 = _pp.addPixels(1,1,_s1,_s3,slice13(_sf3));
    PixelsView p223 = _pp.addPixels(0,0,_s2,_s3,slice23(_sf3));
    p212.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
    p213.setOrientation(PixelsView.Orientation.X1DOWN_X2RIGHT);
    p223.setOrientation(PixelsView.Orientation.X1RIGHT_X2UP);
    _pv2 = new PixelsView[]{p212,p213,p223};
    _vf.addOptions(_pv2,"overlay");
  }
  
  public void addPoints(Map<Integer,float[][][]>[] maps) {
    addPoints(maps,"rO",8.0f);
  }

  public void addPoints(
      Map<Integer,float[][][]>[] maps, String style, float size) 
  {
    _havePoints = true;
    PointsView pt12 = null; 
    PointsView pt13 = null;
    PointsView pt23 = null;
    _i1Map = maps[0];
    _i2Map = maps[1];
    _i3Map = maps[2];
    float[][][] i1Floats = _i1Map.get(_k1);
    float[][][] i2Floats = _i2Map.get(_k2);
    float[][][] i3Floats = _i3Map.get(_k3);
    if (i1Floats!=null) {
//      pt23 = _pp.addPoints(0,0,i1Floats[0],i1Floats[1]);
      pt23 = _pp.addPoints(0,0,i1Floats[1],i1Floats[0]);
      pt23.setStyle(style);
      pt23.setMarkSize(size);
    }
    if (i2Floats!=null) {
      pt13 = _pp.addPoints(1,1,i2Floats[0],i2Floats[1]);
      pt13.setStyle(style);
      pt13.setMarkSize(size);
    }
    if (i3Floats!=null) {
      pt12 = _pp.addPoints(1,0,i3Floats[0],i3Floats[1]);
      pt12.setStyle(style);
      pt12.setMarkSize(size);
    }
    _pts = new PointsView[]{pt12,pt13,pt23};
//    _vf.addPointsOverlay(_pts[0],"points","2");
  }

  public static Map<Integer,float[][][]>[] getSparseCoordsMap(
      int[][][] g1, int[] g2, int[] g3, float d1, float d2, float d3)
  {
    int ng1 = g1[0][0].length;
    int ng2 = g2.length;
    int ng3 = g3.length;
    Map<Integer,float[][][]>[] maps = new HashMap[3];
    Map<Integer,float[][][]> i3Map = new HashMap<Integer, float[][][]>();
    Map<Integer,float[][][]> i2Map = new HashMap<Integer, float[][][]>();
    Map<Integer,float[][][]> i1Map = new HashMap<Integer, float[][][]>();

    Map<Integer,List<Float>> i12Map = new HashMap<Integer,List<Float>>();
    Map<Integer,List<Float>> i13Map = new HashMap<Integer,List<Float>>();
    for (int i3=0; i3<ng3; i3++) {
      float[][] x211 = new float[ng2][ng1];
      float[][] x212 = new float[ng2][ng1];
      for (int i2=0; i2<ng2; i2++) {
        for (int i1=0; i1<ng1; i1++) {
          int i1g = g1[g3[i3]][g2[i2]][i1];
          if (i12Map.containsKey(i1g)) {
            i12Map.get(i1g).add(g2[i2]*d2);
            i13Map.get(i1g).add(g3[i3]*d3);
          } else {
            List<Float> l12 = new ArrayList<Float>();
            l12.add(g2[i2]*d2);
            i12Map.put(i1g,l12);
            List<Float> l13 = new ArrayList<Float>();
            l13.add(g3[i3]*d3);
            i13Map.put(i1g,l13);
          }
          x211[i2][i1] = (float)i1g*d1;
          x212[i2][i1] = (float)g2[i2]*d2;
        }
      }
      i3Map.put(g3[i3],new float[][][]{x211,x212});
    }
    
    for (int i2=0; i2<ng2; i2++) {
      float[][] x131 = new float[ng3][ng1];
      float[][] x132 = new float[ng3][ng1];
      for (int i3=0; i3<ng3; i3++) {
        for (int i1=0; i1<ng1; i1++) {
          x131[i3][i1] = (float)g1[g3[i3]][g2[i2]][i1]*d1;
          x132[i3][i1] = (float)g3[i3]*d3;
        }
      }
      i2Map.put(g2[i2],new float[][][]{x131,x132});
    }
    
    Iterator<Integer> it = i12Map.keySet().iterator();
    while (it.hasNext()) {
      int i1 = it.next();
      List<Float> l12 = i12Map.get(i1);
      List<Float> l13 = i13Map.get(i1);
      int nl2 = l12.size();
      int nl3 = l13.size();
      float[][] x232 = new float[nl3][nl2];
      float[][] x233 = new float[nl3][nl2];
      for (int i=0; i<nl3; i++) {
        x232[i][i] = l12.get(i);
        x233[i][i] = l13.get(i);
      }
      i1Map.put(i1,new float[][][]{x232,x233});
    }
    maps[0] = i1Map;
    maps[1] = i2Map;
    maps[2] = i3Map;
    return maps;
  }
  
  /**
   * Sets the plot title.
   * @param title the title; null for no title.
   */
  public void setTitle(String title) {
    _title = title;
    _pp.setTitle(_title);
  }

  /**
   * Sets the label for axis 1.
   * @param label the label.
   */
  public void setLabel1(String label) {
    _pp.setLabel1(label);
  }
  
  /**
   * Sets the label for axis 2.
   * @param label the label.
   */
  public void setLabel2(String label) {
    _pp.setLabel2(label);
  }
  
  /**
   * Sets the label for axis 3.
   * @param label the label.
   */
  public void setLabel3(String label) {
    _pp.setLabel3(label);
  }
  
  /**
   * Sets the limits for axis 1.
   * @param min the minimum value.
   * @param max the maximum value.
   */
  public void setLimits1(double min, double max) {
    _l1Min = min;
    _l1Max = max;
    if (_orientation==Orientation.X1DOWN_X2RIGHT) {
      _pp.setVLimits(1,min,max);
    } else if (_orientation==Orientation.X1DOWN_X3RIGHT) {
      _pp.setVLimits(1,min,max);
    } else if (_orientation==Orientation.X1RIGHT_X2UP) {
      _pp.setHLimits(0,min,max);
    } else if (_orientation==Orientation.X1RIGHT_X3UP) {
      _pp.setHLimits(0,min,max);
    }
  }
  
  /**
   * Sets the limits for axis 2.
   * @param min the minimum value.
   * @param max the maximum value.
   */
  public void setLimits2(double min, double max) {
    _l2Min = min;
    _l2Max = max;
    if (_orientation==Orientation.X1DOWN_X2RIGHT) {
      _pp.setHLimits(0,min,max);
    } else if (_orientation==Orientation.X1DOWN_X3RIGHT) {
      _pp.setHLimits(1,min,max);
      _pp.setVLimits(0,min,max);
    } else if (_orientation==Orientation.X1RIGHT_X2UP) {
      _pp.setVLimits(1,min,max);
    } else if (_orientation==Orientation.X1RIGHT_X3UP) {
      _pp.setHLimits(1,min,max);
      _pp.setVLimits(0,min,max);
    }
  }
  
  /**
   * Sets the limits for axis 3.
   * @param min the minimum value.
   * @param max the maximum value.
   */
  public void setLimits3(double min, double max) {
    _l3Min = min;
    _l3Max = max;
    if (_orientation==Orientation.X1DOWN_X2RIGHT) {
      _pp.setHLimits(1,min,max);
      //_pp.setVLimits(0,min,max,PlotPanel.Orientation.X1RIGHT_X2UP);
    } else if (_orientation==Orientation.X1DOWN_X3RIGHT) {
      _pp.setHLimits(0,min,max);
    } else if (_orientation==Orientation.X1RIGHT_X2UP) {
      _pp.setHLimits(1,min,max);
      _pp.setVLimits(0,min,max);
    } else if (_orientation==Orientation.X1RIGHT_X3UP) {
      _pp.setVLimits(1,min,max);
    }
  }
  
  public void setHInterval(double interval) {
    _pp.setHInterval(interval);
  }

  public void setVInterval(double interval) {
    _pp.setVInterval(interval);
  }
  
  public void setClips1(float clipMin, float clipMax) {
    _pp.setClips(clipMin,clipMax);
  }
  
  public void setClips2(float clipMin, float clipMax) {
    if (_pv2==null)
      throw new IllegalStateException(
          "Second PixelsView has not been added to plot panel.");
    for (PixelsView pv: _pv2)
      pv.setClips(clipMin,clipMax);
  }

  public void setLineColor(Color color) {
    _pp.setLineColor(color);
  }
  
  public void setColorModel1(IndexColorModel colorModel) {
    _pp.setColorModel(colorModel);
  }
  
  public void setColorModel2(IndexColorModel colorModel) {
    for (PixelsView pv: _pv2)
      pv.setColorModel(colorModel);
  }
  
  public void addColorBar(String label) {
    _pp.addColorBar(label);
//    _pp.setColorBarWidthMinimum(100);
  }
  
  public void setVFormat(int irow, String format) {
    _pp.setVFormat(irow,format);
  }
  
  public void setVInterval(int irow, float interval) {
    _pp.setVInterval(irow,interval);
  }
  
  public void setColorBarWidthMinimum(int widthMinimum) {
    _pp.setColorBarWidthMinimum(widthMinimum);
  }
  
  public void setSize(int width, int height) {
    _vf.setSize(width,height);
  }
  
  public void setFontSizeForPrint(double fontSize, double plotWidth) {
    _vf.setFontSizeForPrint(fontSize,plotWidth);
  }
  
  public void setFontSizeForSlide(double fracWidth, double fracHeight) {
    _vf.setFontSizeForSlide(fracWidth, fracHeight);
  }
  
  public void paintToPng(double dpi, double win, String fileName) {
    _vf.paintToPng(dpi,win,fileName);
  }
  
  public void show() {
    _vf.setVisible(true);
  }
  
  public static void main(String[] args) throws IOException {
    if (args.length != 4) {
      System.out.println("usage: java Viewer datasetPath n1 n2 n3");
      System.exit(0);
    }
    ArrayInputStream ais = new ArrayInputStream(args[0]);
    int n1 = Integer.parseInt(args[1]);
    int n2 = Integer.parseInt(args[2]);
    int n3 = Integer.parseInt(args[3]);
    float[][][] f = new float[n3][n2][n1]; 
    ais.readFloats(f);
    ais.close();
    Viewer3P v = new Viewer3P(f);
    v.show();
  }

	public void addPts(PointsView pv) {
		_pp.getPixelsView23().getTile().addTiledView(pv);
	}

	public void addContours
		(Sampling s1, Sampling s2, Sampling s3, float[][][] x)
	{
		ContoursView cv = _pp.addContours(1,0,s1,s2,x[_k3]);
		cv.setLineColor(Color.WHITE);
	}
		

  ////////////////////////////////////////////////////////////////////////////
  // Private
  
  int _n3;
  int _n2;
  int _n1;
  int _k3;
  int _k2;
  int _k1;
  private double _l1Min;
  private double _l2Min;
  private double _l3Min;
  private double _l1Max;
  private double _l2Max;
  private double _l3Max;
  private Sampling _s1;
  private Sampling _s2;
  private Sampling _s3;
  private Orientation _orientation;
  private String _title = "";
  private ViewerFrame _vf;
  private PlotPanelPixels3 _pp;
  private PixelsView[] _pv1;
  private PixelsView[] _pv2;
  private PointsView[] _pts;
  private SimpleFloat3 _sf3;
  private Map<Integer,float[][][]> _i1Map;
  private Map<Integer,float[][][]> _i2Map;
  private Map<Integer,float[][][]> _i3Map;
  private boolean _havePoints = false;
  private boolean _transpose23;
  
  private void setSlicesPixels(int k1, int k2, int k3, SimpleFloat3 sf3) {
    if (_pv2==null)
      return;
    if (_k1!=k1) {
      if (_transpose23)
        _pv2[2].set(_s3,_s2,slice23(sf3));
      else
        _pv2[2].set(_s2,_s3,slice23(sf3));
    }
    if (_k2!=k2)
      _pv2[1].set(_s1,_s3,slice13(sf3));
    if (_k3!=k3)
      _pv2[0].set(_s1,_s2,slice12(sf3));
  }
  
  private void setSlicesPoints(int k1, int k2, int k3) {
    if (!_havePoints) {
      return;
    }
    float[][][] i1Floats = _i1Map.get(_k1);
    float[][][] i2Floats = _i2Map.get(_k2);
    float[][][] i3Floats = _i3Map.get(_k3);
//    System.out.printf("k1=%d, k2%d, k3=%d%n,",_k1,_k2,_k3);
    if (i1Floats!=null) {
//      System.out.println("found i1floats at k1="+_k1);
      if (_pts[2]==null) {
        PointsView pt23 = _pp.addPoints(0,0,i1Floats[1],i1Floats[0]);
//        pt23 = _pp.addPoints(0,0,i1Floats[0],i1Floats[1]);
        pt23.setStyle("rO");
        pt23.setMarkSize(8.0f);
        _pts[2] = pt23;
      } else {
        _pts[2].set(i1Floats[1],i1Floats[0]);
      }
      _pp.addTiledView(0,0,_pts[2]);
    } else {
      if (_pts[2]!=null) {
        _pp.remove(_pts[2]);
      }
    }
    
    if (i2Floats!=null) {
//      System.out.println("found i2floats at k2="+_k2);
      if (_pts[1]==null) {
        PointsView pt13 = _pp.addPoints(1,1,i2Floats[0],i2Floats[1]);
        pt13.setStyle("rO");
        pt13.setMarkSize(8.0f);
        _pts[1] = pt13;
      } else {
        _pts[1].set(i2Floats[0],i2Floats[1]);       
      }
      _pp.addTiledView(1,1,_pts[1]);
    } else {
      if (_pts[1]!=null) {
        _pp.remove(_pts[1]);
      }
    }
    
    if (i3Floats!=null) {
//      System.out.println("found i3floats at k3="+_k3);
      if (_pts[0]==null) {
        PointsView pt12 = _pp.addPoints(1,0,i3Floats[0],i3Floats[1]);
        pt12.setStyle("rO");
        pt12.setMarkSize(8.0f);
        _pts[0] = pt12;
      } else {
        _pts[0].set(i3Floats[0],i3Floats[1]);
      }
      _pp.addTiledView(1,0,_pts[0]);
    } else {
      if (_pts[0]!=null) {
        _pp.remove(_pts[0]);
      }
    }
  }
  
  private float[][] slice12(SimpleFloat3 sf3) {
    float[][] f12 = new float[_n2][_n1];
    sf3.get12(_n1,_n2,0,0,_k3,f12);
    return f12;
  }
  
  private float[][] slice13(SimpleFloat3 sf3) {
    float[][] f13 = new float[_n3][_n1];
    sf3.get13(_n1,_n3,0,_k2,0,f13);
    return f13;
  }
  
  private float[][] slice23(SimpleFloat3 sf3) {
    float[][] f23 = new float[_n3][_n2];
    sf3.get23(_n2,_n3,_k1,0,0,f23);
    return (_transpose23)?transpose(f23):f23;
  }
  
  private class SliceMode extends Mode {

    protected SliceMode(ModeManager manager) {
      super(manager);
      setName("Slice Mode");
      setMnemonicKey(KeyEvent.VK_S);
      setAcceleratorKey(KeyStroke.getKeyStroke(KeyEvent.VK_S,0));
      setShortDescription("Click/drag in any tile to change slices");
    }

    @Override
    protected void setActive(Component component, boolean active) {
      if (component instanceof Tile) {
        if (active) {
          component.addMouseListener(_ml);
          component.addMouseMotionListener(_mml);
        } else {
          component.removeMouseListener(_ml);
          component.removeMouseMotionListener(_mml);
        }
      }
    }
    
    private static final long serialVersionUID = 1L;
    
    // Handles mouse pressed and released events.
    private MouseListener _ml = new MouseAdapter() {
      public void mousePressed(MouseEvent e) {
        boolean print = false;
        if (e.isControlDown())
          print = true;
        int x = e.getX();
        int y = e.getY();
        Object source = e.getSource();
        if (source instanceof Tile) {
          Tile tile = (Tile)source;
          setSlices(tile,x,y,print);
        }
      }
      public void mouseReleased(MouseEvent e) {} // Do nothing.
    };

    // Handles mouse dragged events.
    private MouseMotionListener _mml = new MouseMotionAdapter() {
      public void mouseDragged(MouseEvent e) {
        boolean print = false;
        if (e.isControlDown())
          print = true;
        int x = e.getX();
        int y = e.getY();
        Object source = e.getSource();
        if (source instanceof Tile) {
          Tile tile = (Tile)source;
          setSlices(tile,x,y,print);
        }
      }
    };

    private void setSlices(Tile tile, int x, int y, boolean print) {
      int row = tile.getRowIndex();
      int col = tile.getColumnIndex();
//      System.out.format("row=%d, col=%d\n",row,col);
      double worldX = tile.pixelToWorldHorizontal(x);
      double worldY = tile.pixelToWorldVertical(y);
      int k1=_k1,k2=_k2,k3=_k3;
      if (_orientation==Orientation.X1DOWN_X2RIGHT) {
        if (row==0) {
          k2 = _s2.indexOfNearest(worldX);
          k3 = _s3.indexOfNearest(worldY);
        } else if (col==0) {
          k1 = _s1.indexOfNearest(worldY);
          k2 = _s2.indexOfNearest(worldX);
        } else {
          k1 = _s1.indexOfNearest(worldY);
          k3 = _s3.indexOfNearest(worldX);
        }
      }
      if (print)
        System.out.format(
            "worldX=%4.2f, worldY=%4.2f, k1=%3d, k2=%3d, k3=%3d\n",
            worldX,worldY,k1,k2,k3);
      updateSlices(k1,k2,k3);
    }
  }
    
}
