package util;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.image.IndexColorModel;
import java.io.IOException;

import javax.swing.DefaultBoundedRangeModel;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.io.ArrayInputStream;
import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.PixelsView.Interpolation;
import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanel;
import edu.mines.jtk.mosaic.PlotPanel.Orientation;
import edu.mines.jtk.mosaic.PointsView;
import edu.mines.jtk.mosaic.TiledView;
import edu.mines.jtk.util.Check;
import edu.mines.jtk.util.Clips;

public class Viewer {

  public enum CMaps {
    GRAY("gray",ColorMap.GRAY),
    JET("jet",ColorMap.JET),
    HUE("hue",ColorMap.HUE),
    BWR("blue white red",ColorMap.BLUE_WHITE_RED),
    PRISM("prism",ColorMap.PRISM),
    AT("alpha test",ColorMap.JET);
    
    private CMaps(String s, IndexColorModel cim) {
      _s = s;
      _cim = cim;
    }
    
    @Override
    public String toString() {
      return _s;
    }
    
    public static CMaps fromString(String map) {
      for (CMaps cmaps : values()) {
        if (map.equals(cmaps._s))
          return cmaps;
      }
      return null;
    }
    
    private String _s;
    private IndexColorModel _cim;
  }
  
  public Viewer(float[][] f) {
    this(f,null);
  }
  
  public Viewer(float[][] f, Orientation o) {
    this(null,null,f,o);
  }
  
  public Viewer(Sampling s1, Sampling s2, float[][] f) {
    this(s1,s2,f,null);
  }
  
  public Viewer(Sampling s1, Sampling s2, float[][] f, Orientation o) {
    _s1 = (s1==null)?new Sampling(f[0].length):s1;
    _s2 = (s2==null)?new Sampling(f.length   ):s2;
    o  = (o==null )?Orientation.X1DOWN_X2RIGHT:o;
    
    // Make initial panel.
    _pp = new PlotPanel(o);
    _pv1 = _pp.addPixels(_s1,_s2,f);
    
    JMenuBar menuBar = new JMenuBar();
    _options = new JMenu("Options");
    addInterpolationOption(_options,_pv1);
    addClipOptions(_options,_pv1,null);
    addColorOptions(_options,_pv1,null);
    menuBar.add(_options);
    
    // Add everything to the PlotFrame, and display.
    _pf = new PlotFrame(_pp);
    _pf.setJMenuBar(menuBar);
    _pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
  }
  
  public Viewer(float[][][] f) {
    this(f,null);
  }
  
  public Viewer(float[][][] f, Orientation o) {
    this(null,null,null,f,o);
  }
  
  public Viewer(Sampling s1, Sampling s2, Sampling s3, float[][][] f) {
    this(s1,s2,s3,f,null);
  }
  
  public Viewer(
      Sampling s1, Sampling s2, Sampling s3, float[][][] f, Orientation o)
  {
    _f = f;
    _s1 = (s1==null)?new Sampling(f[0][0].length):s1;
    _s2 = (s2==null)?new Sampling(f[0].length   ):s2;
    _s3 = (s3==null)?new Sampling(f.length      ):s3;
    o  = (o==null )?Orientation.X1DOWN_X2RIGHT:o;
    
    // Make initial panel, displaying the middle frame.
    int n3 = _s3.getCount();
    _i3 = n3/2;
    int r3 = (int)_s3.getValue(_i3);
    _pp = new PlotPanel(o);
    _title = String.valueOf(_i3);
    _pp.setTitle(_title);
    _pv1 = _pp.addPixels(_s1,_s2,f[_i3]);
    Clips clips = new Clips(f);
    _pv1.setClips(clips.getClipMin(),clips.getClipMax());
    
    JMenuBar menuBar = new JMenuBar();
    _options = new JMenu("Options");
    addInterpolationOption(_options,_pv1);
    addClipOptions(_options,_pv1,null);
    addColorOptions(_options,_pv1,null);
    menuBar.add(_options);
    
    SliderListener sl = new SliderListener();
    DefaultBoundedRangeModel brm = new DefaultBoundedRangeModel(
        r3,0,(int)_s3.getFirst(),(int)_s3.getLast());
    JSlider slider = new JSlider(brm);
    slider.setMajorTickSpacing(n3/10);
    slider.setMinorTickSpacing(n3/150);
    slider.setPaintLabels(true);
    slider.setPaintTicks(true);
    slider.addChangeListener(sl);
      
    _pf = new PlotFrame(_pp);
    _pf.add(slider,BorderLayout.SOUTH);
    _pf.setJMenuBar(menuBar);
    _pf.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
  }
  
  public void addPixels(float[][] f) {
    Check.argument(f.length==_s2.getCount(),
        "f.length is not consistent with sampling");
    Check.argument(f[0].length==_s1.getCount(),
        "f[0].length is not consistent with sampling");
    _pv2 = _pp.addPixels(_s1,_s2,f);
    updateOptions(_pv2,"2");
  }
  
  public void addPixels(float[][][] f) {
    Check.argument(f.length==_s3.getCount(),
        "f.length is not consistent with sampling");
    Check.argument(f[0].length==_s2.getCount(),
        "f[0].length is not consistent with sampling");
    Check.argument(f[0][0].length==_s1.getCount(),
        "f[0][0].length is not consistent with sampling");
    _g = f;
    _pv2 = _pp.addPixels(_s1,_s2,f[_i3]);
    updateOptions(_pv2,"2");
  }
  
  public void addPoints(float[] x2) {
    Check.argument(_s1.getCount()==x2.length,
        "x2.length is not consistend with sampling");
    _pt1 = _pp.addPoints(_s1,x2);
    _pt1.setLineColor(Color.WHITE);
    addRemoveOptions(_options,_pt1,"pt1");
  }
  
  public void addPoints(float[][] x2) {
    Check.argument(_s1.getCount()==x2[0].length,
        "x2.length is not consistend with sampling");
    _p = x2;
    _pt1 = _pp.addPoints(_s1,x2[_i3]);
    _pt1.setLineColor(Color.WHITE);
    addRemoveOptions(_options,_pt1,"pt1");
  }
  
  public void addPoints(float[] x1, float[] x2) {
    _pt2 = _pp.addPoints(x1,x2);
    _pt2.setStyle("rO");
    addRemoveOptions(_options,_pt2,"pt2");
  }
  
  public void addPoints2(float[][] x1, float[][] x2) {
    _pt2 = _pp.addPoints(x1,x2);
    _pt2.setStyle("rO");
    addRemoveOptions(_options,_pt2,"pt2");
  }
  
  public void addPoints3(float[][][] x1, float[][][] x2) {
    _x13 = x1;
    _x23 = x2;
    _pt2 = _pp.addPoints(x1[_i3],x2[_i3]);
    _pt2.setStyle("rO");
    addRemoveOptions(_options,_pt2,"pt2");
  }
  
  public void addPoints(float[][] x1, float[][] x2) {
    Check.argument(x1.length==_f.length,"x1.length==f.length");
    Check.argument(x2.length==_f.length,"x2.length==f.length");
    _x1 = x1;
    _x2 = x2;
    _pt2 = _pp.addPoints(x1[_i3],x2[_i3]);
    _pt2.setStyle("rO");
    addRemoveOptions(_options,_pt2,"pt2");
  }

  public void setTitle(String title) {
    _title = title;
    if (_i3!=Integer.MIN_VALUE)
      _pp.setTitle(_title+" "+_i3);
    else
      _pp.setTitle(_title);
  }

  public void setHLabel(String label) {
    _pp.setHLabel(label);
  }
  
  public void setVLabel(String label) {
    _pp.setVLabel(label);
  }
  
  public void setHLimits(double hmin, double hmax) {
    _pp.setHLimits(hmin,hmax);
  }
  
  public void setVLimits(double vmin, double vmax) {
    _pp.setVLimits(vmin,vmax);
  }
  
  public void setClips1(float clipMin, float clipMax) {
    _pv1.setClips(clipMin,clipMax);
  }
  
  public void setClips2(float clipMin, float clipMax) {
    if (_pv2==null)
      throw new IllegalStateException(
          "Second PixelsView has not been added to plot panel.");
    _pv2.setClips(clipMin,clipMax);
  }

  public void setColorModel1(IndexColorModel colorModel) {
    _pv1.setColorModel(colorModel);
  }
  
  public void setColorModel2(IndexColorModel colorModel) {
    _pv2.setColorModel(colorModel);
  }
  
  public void addColorBar(String label) {
    _pp.addColorBar(label);
    _pp.setColorBarWidthMinimum(100);
  }
  
  public void setSize(int width, int height) {
    _pf.setSize(width,height);
  }
  
  public void show() {
    _pf.setVisible(true);
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
    Viewer v = new Viewer(f);
    v.show();
  }

  ////////////////////////////////////////////////////////////////////////////
  // Private
  
  private PlotFrame _pf;
  private PlotPanel _pp;
  private PixelsView _pv1;
  private PixelsView _pv2;
  private PointsView _pt1;
  private PointsView _pt2;
  private Sampling _s1;
  private Sampling _s2;
  private Sampling _s3;
  private JMenu _options;
  private String _title = "";
  private float[][][] _f;
  private float[][][] _g;
  private float[][] _p;
  private float[][] _x1;
  private float[][] _x2;
  private float[][][] _x13;
  private float[][][] _x23;
  private int _i3 = Integer.MIN_VALUE;

  private void updateOptions(PixelsView pv, String label) {
    addClipOptions(_options,pv,label);
    addColorOptions(_options,pv,label);
    addAlphaOptions(_options,pv,label);
  }
  
  private void addClipOptions(
      JMenu options, final PixelsView pv, String label)
  {
    String name = (label==null)?"Change Clips":"Change Clips ("+label+")";
    JMenuItem changeClips = new JMenuItem(name);
    changeClips.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        float clipMin = pv.getClipMin();
        float clipMax = pv.getClipMax();
        new ClipFrame(clipMin,clipMax,pv);
      }
    });
    options.add(changeClips);
  }
  
  private void addAlphaOptions(
      JMenu options, final PixelsView pv, String label)
  {
    String name = (label==null)?"Change Alpha":"Change Alpha("+label+")";
    JMenuItem changeAlpha = new JMenuItem(name);
    changeAlpha.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        float a = pv.getColorModel().getAlpha(0)/255.0f;
        new AlphaFrame(a,pv);
      }
    });
    options.add(changeAlpha);
  }
  
  private void addInterpolationOption(
      JMenu options, final PixelsView pv)
  {
    String name = "Change Interpolation";
    JMenu changeInterp = new JMenu(name);
    JMenuItem nearest = new JMenuItem("Nearest Neighbor");
    JMenuItem linear  = new JMenuItem("Linear");
    ChangeInterpolationListener cil = new ChangeInterpolationListener(pv);
    nearest.addActionListener(cil);
    linear.addActionListener(cil);
    changeInterp.add(nearest);
    changeInterp.add(linear);
    options.add(changeInterp);
  }
  
  private void addColorOptions(
      JMenu options, final PixelsView pv, String label)
  {
    String name = (label==null)?"Change Colormap":"Change Colormap ("+label+")";
    JMenu changeCmap = new JMenu(name);
    JMenuItem gray = new JMenuItem(CMaps.GRAY.toString());
    JMenuItem jet = new JMenuItem(CMaps.JET.toString());
    JMenuItem bwr = new JMenuItem(CMaps.BWR.toString());
    JMenuItem hue = new JMenuItem(CMaps.HUE.toString());
    JMenuItem prism = new JMenuItem(CMaps.PRISM.toString());
    ChangeColorMapListener ccml = new ChangeColorMapListener(pv);
    gray.addActionListener(ccml);
    jet.addActionListener(ccml);
    bwr.addActionListener(ccml);
    hue.addActionListener(ccml);
    prism.addActionListener(ccml);
    changeCmap.add(gray);
    changeCmap.add(jet);
    changeCmap.add(bwr);
    changeCmap.add(hue);
    changeCmap.add(prism);
    options.add(changeCmap);
  }
  
  private void addRemoveOptions(
      JMenu options, final TiledView tv, String label)
  {
    String name = "Add/Remove "+label;
    JMenuItem addRemove = new JMenuItem(name);
    addRemove.addActionListener(new AddRemoveListener(tv));
    options.add(addRemove);
  }
  
  private class AddRemoveListener implements ActionListener {
    public AddRemoveListener(TiledView tv) {
      _tv = tv;
      _tvLive = true;
    }
    @Override
    public void actionPerformed(ActionEvent e) {
      if (_tvLive) {
        _pp.remove(_tv);
        _tvLive = false;
      } else {
        _pp.addTiledView(_tv);
        _tvLive = true;
      }
    }
    private TiledView _tv;
    private boolean _tvLive;
  }
  
  private class SliderListener implements ChangeListener {
    @Override
    public void stateChanged(ChangeEvent e) {
      JSlider source = (JSlider)e.getSource();
      int i3 = (int)source.getValue();
      _pv1.set(_f[i3]);
      if (_g!=null)
        _pv2.set(_g[i3]);
      if (_p!=null)
        _pt1.set(_s1,_p[i3]);
      if (_x1!=null && _x2!=null)
        _pt2.set(_x1[i3],_x2[i3]);
      if (_x13!=null && _x23!=null)
        _pt2.set(_x13[i3],_x23[i3]);
      _i3 = i3;
      _pp.removeTitle();
      _pp.setTitle(_title+" "+_i3);  
    }
  }
  
  private static class ChangeColorMapListener implements ActionListener {
    public ChangeColorMapListener(PixelsView pv) {_pv = pv;}
    @Override
    public void actionPerformed(ActionEvent e) {
      CMaps source = CMaps.fromString(e.getActionCommand());
      float a = _pv.getColorModel().getAlpha(0)/255.0f;
      IndexColorModel icm;
      switch (source) {
        case GRAY:  icm = ColorMap.setAlpha(CMaps.GRAY._cim, a); break;
        case JET:   icm = ColorMap.setAlpha(CMaps.JET._cim,  a); break;
        case BWR:   icm = ColorMap.setAlpha(CMaps.BWR._cim,  a); break;
        case HUE:   icm = ColorMap.setAlpha(CMaps.HUE._cim,  a); break;
        case PRISM: icm = ColorMap.setAlpha(CMaps.PRISM._cim,a); break;
        default: throw new IllegalArgumentException(
            source+" is not a valid color map");
      }
      _pv.setColorModel(icm);
    }
    private PixelsView _pv;
  }
  
  private static class ChangeInterpolationListener implements ActionListener {
    public ChangeInterpolationListener(PixelsView pv) {_pv = pv;}
    @Override
    public void actionPerformed(ActionEvent e) {
      String interp = e.getActionCommand();
      if (interp.equals("Nearest Neighbor"))
        _pv.setInterpolation(Interpolation.NEAREST);
      else
        _pv.setInterpolation(Interpolation.LINEAR);
    }
    private PixelsView _pv;
  }
}
