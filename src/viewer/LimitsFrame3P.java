package viewer;

import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.beans.PropertyChangeEvent;
import java.beans.PropertyChangeListener;
import java.text.NumberFormat;

import javax.swing.JFormattedTextField;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JPanel;

import edu.mines.jtk.mosaic.PlotPanelPixels3;
import edu.mines.jtk.util.ArrayMath;

public class LimitsFrame3P extends JPanel implements PropertyChangeListener {
  
  /**
   * Constructs a JFrame with options to set the limits for each
   * axis of a {@link PlotPanelPixels3}. The frame is constructed here,
   * but is not made visible by default. Use the {@link #getFrame()}
   * method to make the frame visible when required. 
   * @param v the Viewer3P instance that contains the 
   *  {@link PlotPanelPixels3} instance.
   */
  public LimitsFrame3P(Viewer3P v) {
    super(new GridBagLayout());
    GridBagConstraints c = new GridBagConstraints();
    c.insets = new Insets(0,10,0,0);
    
    _v = v;
    _l1MinVal = v.getL1Min();
    _l2MinVal = v.getL2Min();
    _l3MinVal = v.getL3Min();
    _l1MaxVal = v.getL1Max();
    _l2MaxVal = v.getL2Max();
    _l3MaxVal = v.getL3Max();
    
    NumberFormat limitsFormat = NumberFormat.getNumberInstance();
    limitsFormat.setGroupingUsed(false);
    
    JLabel k1Label = new JLabel("Limits 1:");
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;
    c.gridx = 0;
    c.gridy = 0;
    add(k1Label,c);
    _l1Min = new JFormattedTextField(limitsFormat);
    _l1Min.setColumns(5);
    _l1Min.setValue(_l1MinVal);
    _l1Min.addPropertyChangeListener("value",this);
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;
    c.gridx = 1;
    c.gridy = 0;
    add(_l1Min,c);
    _l1Max = new JFormattedTextField(limitsFormat);
    _l1Max.setColumns(5);
    _l1Max.setValue(_l1MaxVal);
    _l1Max.addPropertyChangeListener("value",this);
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;
    c.gridx = 2;
    c.gridy = 0;
    add(_l1Max,c);
        
    JLabel k2Label = new JLabel("Limits 2:");
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;
    c.gridx = 0;
    c.gridy = 1;
    add(k2Label,c);
    _l2Min = new JFormattedTextField(limitsFormat);
    _l2Min.setColumns(5);
    _l2Min.setValue(_l2MinVal);
    _l2Min.addPropertyChangeListener("value",this);
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;    
    c.gridx = 1;
    c.gridy = 1;
    add(_l2Min,c);
    _l2Max = new JFormattedTextField(limitsFormat);
    _l2Max.setColumns(5);
    _l2Max.setValue(_l2MaxVal);
    _l2Max.addPropertyChangeListener("value",this);
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;
    c.gridx = 2;
    c.gridy = 1;
    add(_l2Max,c);
    
    JLabel k3Label = new JLabel("Limits 3:");
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;
    c.gridx = 0;
    c.gridy = 2;
    add(k3Label,c);
    _l3Min = new JFormattedTextField(limitsFormat);
    _l3Min.setColumns(5);
    _l3Min.setValue(_l3MinVal);
    _l3Min.addPropertyChangeListener("value",this);
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;    
    c.gridx = 1;
    c.gridy = 2;
    add(_l3Min,c);
    _l3Max = new JFormattedTextField(limitsFormat);
    _l3Max.setColumns(5);
    _l3Max.setValue(_l3MaxVal);
    _l3Max.addPropertyChangeListener("value",this);
    c.fill = GridBagConstraints.HORIZONTAL;
    c.weightx = 0.0;
    c.gridwidth = 1;
    c.gridx = 2;
    c.gridy = 2;
    add(_l3Max,c);
    
    _frame = new JFrame("Enter Limits");
    _frame.add(this);
    _frame.setSize(280,150);
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
    if (source==_l1Min) {
      double l1MinVal = ((Number)_l1Min.getValue()).doubleValue();
      _v.setLimits1(l1MinVal,_l1MaxVal);
      _l1MinVal = l1MinVal;
    }
    if (source==_l1Max) {
      double l1MaxVal = ((Number)_l1Max.getValue()).doubleValue();
      _v.setLimits1(_l1MinVal,l1MaxVal);
      _l1MaxVal = l1MaxVal;
    }
    if (source==_l2Min) {
      double l2MinVal = ((Number)_l2Min.getValue()).doubleValue();
      _v.setLimits2(l2MinVal,_l2MaxVal);
      _l2MinVal = l2MinVal;
    }
    if (source==_l2Max) {
      double l2MaxVal = ((Number)_l2Max.getValue()).doubleValue();
      _v.setLimits2(_l2MinVal,l2MaxVal);
      _l2MaxVal = l2MaxVal;
    }
    if (source==_l3Min) {
      double l3MinVal = ((Number)_l3Min.getValue()).doubleValue();
      _v.setLimits3(l3MinVal,_l3MaxVal);
      _l3MinVal = l3MinVal;
    }
    if (source==_l3Max) {
      double l3MaxVal = ((Number)_l3Max.getValue()).doubleValue();
      _v.setLimits3(_l3MinVal,l3MaxVal);
      _l3MaxVal = l3MaxVal;
    }
  }
  
  public static void main(String[] args) {
    new LimitsFrame3P(new Viewer3P(ArrayMath.zerofloat(1999,500,14)));
  }

  private final JFrame _frame;
  private final Viewer3P _v;
  private JFormattedTextField _l1Min; // input for limits 1 min.
  private JFormattedTextField _l1Max; // input for limits 1 max.
  private JFormattedTextField _l2Min; // input for limits 2 min.
  private JFormattedTextField _l2Max; // input for limits 2 max.
  private JFormattedTextField _l3Min; // input for limits 3 min.
  private JFormattedTextField _l3Max; // input for limits 3 max.
  private double _l1MinVal; // current limits 1 minimum value.
  private double _l1MaxVal; // current limits 1 maximum value.
  private double _l2MinVal; // current limits 2 minimum value.
  private double _l2MaxVal; // current limits 2 maximum value.
  private double _l3MinVal; // current limits 3 minimum value.
  private double _l3MaxVal; // current limits 3 maximum value.
  private static final long serialVersionUID = 1L;
  
}
