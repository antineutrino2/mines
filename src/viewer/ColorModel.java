package viewer;

import java.awt.image.IndexColorModel;

import edu.mines.jtk.awt.ColorMap;

/**
 * IndexColorModels for plotting.
 */
public enum ColorModel {
  GRAY("gray",ColorMap.GRAY),
  JET("jet",ColorMap.JET),
  HUE("hue",ColorMap.HUE),
  BWR("blue white red",ColorMap.BLUE_WHITE_RED),
  PRISM("prism",ColorMap.PRISM),
  AT("alpha test",ColorMap.JET);

  /**
   * Constructor.
   * @param name the name of IndexColorModel.
   * @param icm the IndexColorModel.
   */
  private ColorModel(String name, IndexColorModel icm) {
    _name = name;
    _icm = icm;
  }
  
  public IndexColorModel getIndexColorModel() {
    return _icm;
  }

  @Override
  public String toString() {
    return _name;
  }

  /**
   * Gets the enum constant from the name.
   * @param name
   * @return the enum, or null.
   */
  public static ColorModel fromString(String name) {
    for (ColorModel colorModel : values()) {
      if (name.equals(colorModel._name))
        return colorModel;
    }
    return null;
  }

  private String _name;
  private IndexColorModel _icm;
}