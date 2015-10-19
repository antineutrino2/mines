package viewer;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.text.NumberFormat;

import javax.swing.DefaultBoundedRangeModel;
import javax.swing.JFormattedTextField;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

import edu.mines.jtk.util.ArrayMath;

public class SliceFrame extends JPanel implements PropertyChangeListener {
  
  /**
   * Constructs a JFrame with options to change the slices for
   * a {@link PlotPanelPixels3}. Slices can be changed via a text 
   * field or JSlider. The frame is constructed here, but is not 
   * made visible by default. Use the {@link #getFrame()} method 
   * to make the frame visible when required. 
   * @param v the Viewer3P instance that contains the 
   *  {@link PlotPanelPixels3} instance.
   */
  public SliceFrame(Viewer3P v) {
    super(new GridBagLayout());
    SliderListener sl = new SliderListener();
    GridBagConstraints c = new GridBagConstraints();
    c.insets = new Insets(0,10,0,0);
    
    _v = v;
    int[] slices = v.getSlices();
    _k1Val = slices[0];
    _k2Val = slices[1];
    _k3Val = slices[2];
    _n1 = v.getN1();
    _n2 = v.getN2();
    _n3 = v.getN3();
    
    NumberFormat sliceFormat = NumberFormat.getIntegerInstance();
    sliceFormat.setParseIntegerOnly(true);
    sliceFormat.setGroupingUsed(false);
    
    JLabel k1Label = new JLabel("k1:");
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;
    c.gridx = 0;
    c.gridy = 0;
    add(k1Label,c);
    _k1 = new JFormattedTextField(sliceFormat);
    _k1.setColumns(5);
    _k1.setValue(_k1Val);
    _k1.addPropertyChangeListener("value",this);
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;
    c.gridx = 1;
    c.gridy = 0;
    add(_k1,c);
    DefaultBoundedRangeModel brm1 = new DefaultBoundedRangeModel(
        _k1Val,0,0,_n1-1);
    _sliderK1 = new JSlider(brm1);
    int ts1 = 1;
    while ((_n1/ts1)>N_TICKS) ts1++; 
    _sliderK1.setMajorTickSpacing(ts1);
    _sliderK1.setPaintLabels(true);
    _sliderK1.setPaintTicks(true);
    _sliderK1.addChangeListener(sl);
    c.weightx = 1.0;
    c.gridwidth = 4;
    c.gridx = 2;
    c.gridy = 0;
    add(_sliderK1,c);
    
    JLabel k2Label = new JLabel("k2:");
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;
    c.gridx = 0;
    c.gridy = 1;
    add(k2Label,c);
    _k2 = new JFormattedTextField(sliceFormat);
    _k2.setColumns(5);
    _k2.setValue(_k2Val);
    _k2.addPropertyChangeListener("value",this);
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;    
    c.gridx = 1;
    c.gridy = 1;
    add(_k2,c);
    DefaultBoundedRangeModel brm2 = new DefaultBoundedRangeModel(
        _k2Val,0,0,_n2-1);
    _sliderK2 = new JSlider(brm2);
    int ts2 = 1;
    while ((_n2/ts2)>N_TICKS) ts2++;
    _sliderK2.setMajorTickSpacing(ts2);
    _sliderK2.setPaintLabels(true);
    _sliderK2.setPaintTicks(true);
    _sliderK2.addChangeListener(sl);
    c.weightx = 1.0;
    c.gridwidth = 4;    
    c.gridx = 2;
    c.gridy = 1;
    add(_sliderK2,c);
    
    JLabel k3Label = new JLabel("k3:");
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;
    c.gridx = 0;
    c.gridy = 2;
    add(k3Label,c);
    _k3 = new JFormattedTextField(sliceFormat);
    _k3.setColumns(5);
    _k3.setValue(_k3Val);
    _k3.addPropertyChangeListener("value",this);
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;
    c.gridx = 1;
    c.gridy = 2;
    add(_k3,c);
    DefaultBoundedRangeModel brm3 = new DefaultBoundedRangeModel(
        _k3Val,0,0,_n3-1);
    _sliderK3 = new JSlider(brm3);
    int ts3 = 1;
    while ((_n3/ts3)>N_TICKS) ts3++;
    _sliderK3.setMajorTickSpacing(ts3);
    _sliderK3.setPaintLabels(true);
    _sliderK3.setPaintTicks(true);
    _sliderK3.addChangeListener(sl);
    c.weightx = 1.0;
    c.gridwidth = 4;
    c.gridx = 2;
    c.gridy = 2;
    add(_sliderK3,c);
    
    _frame = new JFrame("Enter Slices:");
    _frame.add(this);
    _frame.setSize(500,180);
    _frame.setDefaultCloseOperation(JFrame.DISPOSE_ON_CLOSE);
  }
  
  /**
   * Returns the frame.
   * @return the frame.
   */
  public JFrame getFrame() {
    return _frame;
  }

  @Override
  public void propertyChange(PropertyChangeEvent e) {
    Object source = e.getSource();
    if (source==_k1) {
      int k1 = ((Number)_k1.getValue()).intValue();
      k1 = (k1<0)?0:(k1>=_n1)?_n1-1:k1;
      _v.updateSlices(k1,_k2Val,_k3Val);
      _sliderK1.setValue(k1);
      _k1Val = k1;
    }
    if (source==_k2) {
      int k2 = ((Number)_k2.getValue()).intValue();
      k2 = (k2<0)?0:(k2>=_n2)?_n2-1:k2;
      _v.updateSlices(_k1Val,k2,_k3Val);
      _sliderK2.setValue(k2);
      _k2Val = k2;
    }
    if (source==_k3) {
      int k3 = ((Number)_k3.getValue()).intValue();
      k3 = (k3<0)?0:(k3>=_n3)?_n2-1:k3;
      _v.updateSlices(_k1Val,_k2Val,k3);
      _sliderK3.setValue(k3);
      _k3Val = k3;
    }
  }
  
  public static void main(String[] args) {
    new SliceFrame(new Viewer3P(ArrayMath.zerofloat(1999,500,14)));
  }

  private JFrame _frame;
  private Viewer3P _v;
  private JFormattedTextField _k1;
  private JFormattedTextField _k2;
  private JFormattedTextField _k3;
  private int _k1Val;
  private int _k2Val;
  private int _k3Val;
  private int _n1;
  private int _n2;
  private int _n3;
  private JSlider _sliderK1;
  private JSlider _sliderK2;
  private JSlider _sliderK3;
  private static final long serialVersionUID = 1L;
  private static final int N_TICKS = 10;

  private class SliderListener implements ChangeListener {
    @Override
    public void stateChanged(ChangeEvent e) {
      Object source = e.getSource();
      if (source==_sliderK1) {
        int k1 = _sliderK1.getValue();
        _v.updateSlices(k1,_k2Val,_k3Val);
        _k1.setValue(k1);
        _k1Val = k1;
      }
      if (source==_sliderK2) {
        int k2 = _sliderK2.getValue();
        _v.updateSlices(_k1Val,k2,_k3Val);
        _k2.setValue(k2);
        _k2Val = k2;
      }
      if (source==_sliderK3) {
        int k3 = _sliderK3.getValue();
        _v.updateSlices(_k1Val,_k2Val,k3);
        _k3.setValue(k3);
        _k3Val = k3;
      }
    }
  }
  
}
