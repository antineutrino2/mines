package anm;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import edu.mines.jtk.sgl.*;

import javax.swing.Timer;

// THANKS TO STEFAN COMPTON FOR PROVIDING THIS

public class AnimateSimpleFrame implements ActionListener {

  public AnimateSimpleFrame(SimpleFrame sf, TriangleGroup[] tga, int delay) {
    _sf = sf;
    _tga = tga;
    _ntg = _tga.length;
    _view = _sf.getViewCanvas().getView();
    _world = new World();
    _timer = new Timer(delay,this);
    _timer.start();
  }
  
  @Override
  public void actionPerformed(ActionEvent arg0) {
    if (_tgc!=null)
      _world.removeChild(_tgc);
    _tgc = _tga[_index];
    _world.addChild(_tga[_index]);
    _view.setWorld(_world);
    _index++;
    if (_index>=_ntg)
      _index = 0;
  }

  private Timer _timer;
  private SimpleFrame _sf;
  private View _view;
  private World _world;
  private TriangleGroup[] _tga;
  private TriangleGroup _tgc = null;
  private int _ntg;
  private int _index = 0;

}
