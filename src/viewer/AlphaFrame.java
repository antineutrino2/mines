package viewer;

import java.awt.image.IndexColorModel;

import javax.swing.DefaultBoundedRangeModel;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.mosaic.PixelsView;

public class AlphaFrame extends JPanel implements ChangeListener {
  
  public AlphaFrame(float alpha, PixelsView[] pv) {
    _pv = pv;
    DefaultBoundedRangeModel brm = 
        new DefaultBoundedRangeModel((int)(alpha*100),0,0,100);
    JSlider slider = new JSlider(brm);
    slider.setMajorTickSpacing(25);
    slider.setMinorTickSpacing(5);
    slider.setPaintLabels(true);
    slider.setPaintTicks(true);
    slider.addChangeListener(this);
    add(slider);
    
    JFrame frame = new JFrame("Transparency");
    frame.add(this);
    frame.pack();
    frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
    frame.setVisible(true);
  }

  @Override
  public void stateChanged(ChangeEvent e) {
    JSlider source = (JSlider)e.getSource();
    int perc = (int)source.getValue();
    IndexColorModel icm = _pv[0].getColorModel();
    for (PixelsView pv : _pv)
      pv.setColorModel(ColorMap.setAlpha(icm,perc/100.0));
  }
  
  private PixelsView[] _pv;
  
}
