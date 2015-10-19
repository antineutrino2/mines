package viewer;

import java.awt.Component;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.MouseListener;
import java.awt.image.IndexColorModel;
import java.io.File;

import javax.swing.AbstractAction;
import javax.swing.Action;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.KeyStroke;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.mosaic.ColorBar;
import edu.mines.jtk.mosaic.PixelsView;
import edu.mines.jtk.mosaic.PixelsView.Interpolation;
import edu.mines.jtk.mosaic.PlotFrame;
import edu.mines.jtk.mosaic.PlotPanel;
import edu.mines.jtk.mosaic.PlotPanelPixels3;
import edu.mines.jtk.mosaic.PointsView;
import edu.mines.jtk.mosaic.TiledView;

/**
 * A PlotFrame that includes interactive features such
 * as options to change the color model, pixel interpolation,
 * and clips. These options can also be changed for a pixel 
 * overlay.
 */
public class ViewerFrame extends PlotFrame {

  /**
   * Constructs a ViewerFrame for the {@code panel} with options
   * to modify settings for the {@link PixelsView} array. Changing
   * settings will effect all PixelsViews in this array. For
   * instance, if constructing a ViewerFrame for a 
   * {@link PlotPanelPixels3}, all three PixelsViews can be passed
   * into the {@code pv} array.
   * </p>
   * Other pixels can be added with the 
   * {@link #addOptions(PixelsView[], String)} method.
   * @param panel
   * @param pv
   */
  public ViewerFrame(PlotPanel panel, PixelsView[] pv) {
    super(panel);
    _panel = panel;
    JMenuBar menuBar = new JMenuBar();
    _options = new JMenu("Options");
    addInterpolationOption(pv);
    addClipOptions(pv,null);
    addColorOptions(pv,null);
    addSaveOption(this);
    menuBar.add(_options);
    
    setJMenuBar(menuBar);
    addMouseListener(_ml);
    setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
  }
  
  /**
   * Adds a Component {@code c} to the options menu
   * of this frame.
   * @param c the component to add.
   */
  public void addToMenu(Component c) {
    _options.add(c);
  }
  
  /**
   * Adds setting for the {@link PixelsView} array to the
   * options menu with the {@code label}.
   * @param pv 
   * @param label the label for the added settings in
   *  the options menu or {@code null} for no label.
   */
  public void addOptions(PixelsView[] pv, String label) {
    addClipOptions(pv,label);
    addColorOptions(pv,label);
    addAlphaOptions(pv,label);
    addRemoveOptions(pv[0],label,"3"); //TODO what about mulitple PixelsViews?
  }
  
  /**
   * Adds settings for the {@link PointsView} to the 
   * options menu with the {@code label}.
   * @param ptv
   * @param label label the label for the added settings in
   *  the options menu or {@code null} for no label.
   */
  public void addOptions(PointsView ptv, String label, String key) {
//    addRemoveOptions(ptv,label,"2");
    addRemoveOptions(ptv,label,key);
  }
  
  ////////////////////////////////////////////////////////////////////////
  // Private
  
  private JMenu _options;
  private PlotPanel _panel;
  private static final long serialVersionUID = 1L;
  
  private void addInterpolationOption(final PixelsView[] pv) {
    String name = "Change Interpolation";
    JMenu changeInterp = new JMenu(name);
    JMenuItem nearest = new JMenuItem("Nearest Neighbor");
    JMenuItem linear  = new JMenuItem("Linear");
    ChangeInterpolationListener cil = new ChangeInterpolationListener(pv);
    nearest.addActionListener(cil);
    linear.addActionListener(cil);
    changeInterp.add(nearest);
    changeInterp.add(linear);
    _options.add(changeInterp);
  }
  
  private void addAlphaOptions(final PixelsView[] pv, String label) {
    String name = (label==null)?"Change Alpha":"Change Alpha("+label+")";
    JMenuItem changeAlpha = new JMenuItem(name);
    changeAlpha.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        float a = pv[0].getColorModel().getAlpha(0)/255.0f;
        new AlphaFrame(a,pv);
      }
    });
    _options.add(changeAlpha);
  }
  
  private void addClipOptions(final PixelsView[] pv, String label) {
    String name = (label==null)?"Change Clips":"Change Clips ("+label+")";
    JMenuItem changeClips = new JMenuItem(name);
    changeClips.addActionListener(new ActionListener() {
      @Override
      public void actionPerformed(ActionEvent e) {
        float clipMin = pv[0].getClipMin();
        float clipMax = pv[0].getClipMax();
        new ClipFrame(clipMin,clipMax,pv);
      }
    });
    _options.add(changeClips);
  }
  
  private void addColorOptions(final PixelsView[] pv, String label) {
    String name = (label==null)?"Change Colormap":"Change Colormap ("+label+")";
    JMenu changeCmap = new JMenu(name);
    JMenuItem gray = new JMenuItem(ColorModel.GRAY.toString());
    JMenuItem jet = new JMenuItem(ColorModel.JET.toString());
    JMenuItem bwr = new JMenuItem(ColorModel.BWR.toString());
    JMenuItem hue = new JMenuItem(ColorModel.HUE.toString());
    JMenuItem prism = new JMenuItem(ColorModel.PRISM.toString());
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
    _options.add(changeCmap);
  }

  private void addSaveOption(final PlotFrame pf) {
    Action saveToPngAction = new AbstractAction("Save to PNG") {
      private static final long serialVersionUID = 1L;
      public void actionPerformed(ActionEvent event) {
        JFileChooser fc = new JFileChooser(System.getProperty("user.dir"));
        fc.showSaveDialog(pf);
        File file = fc.getSelectedFile();
        if (file!=null) {
          String filename = file.getAbsolutePath();
          if (!filename.toLowerCase().endsWith(".png"))
            filename = filename+".png";
          paintToPng(300,6,filename);
        }
      }
    };
    JMenuItem saveToPng = new JMenuItem(saveToPngAction);
    _options.add(saveToPng);
  }
  
  private void addRemoveOptions(
      final TiledView tv, String label, String key)
  {
    String name = "Add/Remove "+label;
    JMenuItem addRemove = new JMenuItem(name);
    addRemove.addActionListener(new AddRemoveListener(tv));
    addRemove.setAccelerator(KeyStroke.getKeyStroke(key));
    _options.add(addRemove);
  }
  
  private class AddRemoveListener implements ActionListener {
    public AddRemoveListener(TiledView tv) {
      _tv = tv;
      _tvLive = true;
    }
    @Override
    public void actionPerformed(ActionEvent e) {
      if (_tvLive) {
        _panel.remove(_tv);
        _tvLive = false;
      } else {
        _panel.addTiledView(_tv);
        _tvLive = true;
      }
    }
    private TiledView _tv;
    private boolean _tvLive;
  }
  
  private static class ChangeColorMapListener implements ActionListener {
    public ChangeColorMapListener(PixelsView[] pv) {_pv = pv;}
    @Override
    public void actionPerformed(ActionEvent e) {
      ColorModel source = ColorModel.fromString(e.getActionCommand());
      float a = _pv[0].getColorModel().getAlpha(0)/255.0f;
      IndexColorModel model;
      switch (source) {
        case GRAY:  model = ColorModel.GRAY.getIndexColorModel(); break;
        case JET:   model = ColorModel.JET.getIndexColorModel(); break;
        case BWR:   model = ColorModel.BWR.getIndexColorModel(); break;
        case HUE:   model = ColorModel.HUE.getIndexColorModel(); break;
        case PRISM: model = ColorModel.PRISM.getIndexColorModel(); break;
        default: throw new IllegalArgumentException(
            source+" is not a valid ColorModel");
      }
      IndexColorModel icm = ColorMap.setAlpha(model,a);
      for (PixelsView pv : _pv)
        pv.setColorModel(icm);
    }
    private PixelsView[] _pv;
  }
  
  private static class ChangeInterpolationListener implements ActionListener {
    public ChangeInterpolationListener(PixelsView[] pv) {_pv = pv;}
    @Override
    public void actionPerformed(ActionEvent e) {
      String interp = e.getActionCommand();
      for (PixelsView pv : _pv) {
        if (interp.equals("Nearest Neighbor"))
          pv.setInterpolation(Interpolation.NEAREST);
        else
          pv.setInterpolation(Interpolation.LINEAR);  
      }
    }
    private PixelsView[] _pv;
  }

  private MouseListener _ml = new MouseAdapter() {
    public void mousePressed(MouseEvent e) {
      if (e.isControlDown() && e.isAltDown()) {
        ColorBar cbar = null;
        for (Component c : _panel.getComponents()) {
          if (c instanceof ColorBar)
            cbar = (ColorBar)c;
        }
        if (cbar!=null) {
          int cbw = cbar.getWidth();
          System.out.println("ColorBar width = "+cbw);
        }
      }
    }
    public void mouseReleased(MouseEvent e) {} // Do nothing.
  };
}
