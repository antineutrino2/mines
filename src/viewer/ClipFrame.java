package viewer;

import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.text.NumberFormat;

import javax.swing.JFormattedTextField;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

import edu.mines.jtk.mosaic.PixelsView;

public class ClipFrame extends JPanel implements PropertyChangeListener {
  
  public ClipFrame(float minValue, float maxValue, PixelsView[] pv) {
    _minValue = minValue;
    _maxValue = maxValue;
    _pv = pv;
    NumberFormat clipFormat = NumberFormat.getNumberInstance();
    _min = new JFormattedTextField(clipFormat);
    _min.setColumns(10);
    _min.setValue(_minValue);
    _min.addPropertyChangeListener("value",this);
    _max = new JFormattedTextField(clipFormat);
    _max.setColumns(10);
    _max.setValue(_maxValue);
    _max.addPropertyChangeListener("value",this);
    JLabel minLabel = new JLabel("min:");
    JLabel maxLabel = new JLabel("max:");
    add(minLabel);
    add(_min);
    add(maxLabel);
    add(_max);
    
    JFrame frame = new JFrame("Enter Clips");
    frame.add(this);
    frame.pack();
    frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
    frame.setVisible(true);
  }

  @Override
  public void propertyChange(PropertyChangeEvent e) {
    Object source = e.getSource();
    if (source == _min) {
      float clipMin = ((Number)_min.getValue()).floatValue();
      for (PixelsView pv : _pv)
        pv.setClips(clipMin,_maxValue);
      _minValue = clipMin;
    }
    if (source == _max) {
      float clipMax = ((Number)_max.getValue()).floatValue();
      for (PixelsView pv : _pv)
        pv.setClips(_minValue,clipMax);
      _maxValue = clipMax;
    }
  }
  
  public static void main(String[] args) {
    new ClipFrame(-1.0f,1.0f,null);
  }

  private float _minValue;
  private float _maxValue;
  private PixelsView[] _pv;
  private JFormattedTextField _min;
  private JFormattedTextField _max;
}
