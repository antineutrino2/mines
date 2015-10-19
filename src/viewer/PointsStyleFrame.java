package viewer;

import java.awt.Color;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;

import javax.swing.JComboBox;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

import edu.mines.jtk.mosaic.PointsView;

/*package*/ class PointsStyleFrame extends JFrame {
  
  /*package*/ enum Colors {
    BLACK(Color.BLACK),
    WHITE(Color.WHITE),
    GRAY(Color.GRAY),
    RED(Color.RED),
    GREEN(Color.GREEN),
    BLUE(Color.BLUE),
    YELLOW(Color.YELLOW),
    ORANGE(Color.ORANGE),
    MAGENTA(Color.MAGENTA),
    PINK(Color.PINK); 
    
    private Colors(Color color) {
      _color = color;
    }
    
    private Color _color;
  }

  /*package*/ PointsStyleFrame(PointsView[] ptv, boolean isLine) {
    _ptv = ptv;
    _panel = new JPanel(new GridBagLayout());
    GridBagConstraints c = new GridBagConstraints();
    
    JRadioButton line = new JRadioButton("Line");
    line.setSelected(isLine);
    line.addActionListener(new ActionListener() {
    public void actionPerformed(ActionEvent e) {
        
    }});
    c.fill = GridBagConstraints.HORIZONTAL;
    c.gridx = 0;
    c.gridy = 0;
    _panel.add(line,c);
    
    JRadioButton marks = new JRadioButton("Marks");
    marks.setSelected(!isLine);
    c.fill = GridBagConstraints.HORIZONTAL;
    c.gridx = 1;
    c.gridy = 0;
    _panel.add(marks,c);
    
    JComboBox lineStyles = new JComboBox(PointsView.Line.values());
    c.fill = GridBagConstraints.HORIZONTAL;
    c.gridx = 0;
    c.gridy = 1;
    _panel.add(lineStyles,c);
    
    JComboBox markStyles = new JComboBox(PointsView.Mark.values());
    c.fill = GridBagConstraints.HORIZONTAL;
    c.gridx = 1;
    c.gridy = 1;
    _panel.add(markStyles,c);
    
//    JLabel sizeLabel = new JLabel("Size:");
//    c.fill = GridBagConstraints.HORIZONTAL;
//    c.anchor = GridBagConstraints.CENTER;
//    c.weightx = 10;
//    c.gridx = 0;
//    c.gridy = 2;
//    _panel.add(sizeLabel,c);
    
    JComboBox size = new JComboBox(_sizes);
    size.setEditable(true);
    c.fill = GridBagConstraints.HORIZONTAL;
    c.gridx = 1;
    c.gridy = 2;
    _panel.add(size,c);
    
    JComboBox colors = new JComboBox(Colors.values());
    c.fill = GridBagConstraints.HORIZONTAL;
    c.gridx = 1;
    c.gridy = 3;
    _panel.add(colors,c);
    
    add(_panel);
    pack();
    setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
    setVisible(true);
  }
  
  public static void main(String[] args) {
    new PointsStyleFrame(null,true);
  }
  
  private PointsView[] _ptv;
  private JPanel _panel;
  private Float[] _sizes = new Float[]{1.0f,2.0f,4.0f,6.0f,8.0f,10.0f,12.0f};
  private static final long serialVersionUID = 1L;

}
