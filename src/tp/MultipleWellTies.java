package tp;

import java.io.*;
import java.util.*;
import wt.SeismicModel1D;
import dtw.DynamicWarpingWT;
import edu.mines.jtk.io.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.awt.ColorMap;
import edu.mines.jtk.dsp.Sampling;
import edu.mines.jtk.dsp.SincInterp;
import edu.mines.jtk.mosaic.SimplePlot;
import edu.mines.jtk.dsp.SincInterp.Extrapolation;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Compute a multiple well ties for Teapot Dome wells.
 * 
 * TODO
 * 
 * @author Andrew Munoz, Colorado School of Mines
 * @version 01.11.14
 */

public class MultipleWellTies {

  // This value represents missing values in all well curve data.
  public static final float NULL_VALUE = -999.2500f;

  // Static correction velocity for the seismic survey (km/s)
  public static final float STATIC_VEL = 2.7432f; // 9000 ft/s


 /**
   * Multiple single-well ties.
   * @param ids log API numbers
   * @param dz  log sampling interval
   * @param g   seismic image
   * @param s1  x1 time sampling of the seismic 
   * @param wset file path to binary well log set
   * @param csmDir file path to data
   */
  public MultipleWellTies(
    long[] ids, double dz, Sampling s1, String wset, String csmDir) 
  {
    this(ids,dz,null,s1,null,wset,csmDir);
  }

  /**
   * Multiple well ties in 2D.
   * @param ids log API numbers
   * @param dz  log sampling interval
   * @param g   seismic image
   * @param s1  x1 time sampling of the seismic 
   * @param s2  x2 spatial sampling of the seismic 
   * @param wset file path to binary well log set
   * @param csmDir file path to data
   */
  public MultipleWellTies(
    long[] ids, double dz, float[][] g, 
    Sampling s1, Sampling s2, String wset, String csmDir) 
  {
    int n = ids.length;
    for (int i=0; i<n; ++i) {
      WellTie well = new WellTie(ids[i],dz,s1,wset); 
      addWell(well);
    }
    _csmDir = csmDir;
    _s1 = s1; _n1 = s1.getCount(); _d1 = s1.getDelta(); _f1 = s1.getFirst();
    if (s2!=null) {
      _s2 = s2; _n2 = s2.getCount(); _d2 = s2.getDelta(); _f2 = s2.getFirst();
    }
    _dz = dz;
    _ids = ids;
    _g2 = g;
    _static = new LongList();
    computeMaxWellDistance();
  }

  /**
   * Multiple well ties in 3D.
   * @param ids log API numbers
   * @param dz  log sampling interval
   * @param g   seismic image
   * @param s1  x1 time sampling of the seismic 
   * @param s2  x2 spatial sampling of the seismic 
   * @param s3  x3 spatial sampling of the seismic 
   * @param wset file path to binary well log set
   * @param csmDir file path to data
   */
  public MultipleWellTies(
    long[] ids, double dz, float[][][] g, 
    Sampling s1, Sampling s2, Sampling s3, String wset, String csmDir) 
  {
    int n = ids.length;
    for (int i=0; i<n; ++i) {
      WellTie well = new WellTie(ids[i],dz,s1,wset); 
      addWell(well);
    }
    _csmDir = csmDir;
    _s1 = s1; _n1 = s1.getCount(); _d1 = s1.getDelta(); _f1 = s1.getFirst();
    _s2 = s2; _n2 = s2.getCount(); _d2 = s2.getDelta(); _f2 = s2.getFirst();
    _s3 = s3; _n3 = s3.getCount(); _d3 = s3.getDelta(); _f3 = s3.getFirst();
    _dz = dz;
    _ids = ids;
    _g3 = g;
    _static = new LongList();
    computeMaxWellDistance();
  }


  /**
   * Computes synthetic seismograms.
   * This method uses vertically propagating plane waves to compute
   * synthetic seismograms in a attenuative medium (constant q). 
   * This method computes internal and free surface multiples.
   * @param dataDir the directory for storing computed synthetic seismograms
   * @param fp peak frequency of synthetic seismogram
   * @param q  constant attenuation factor
   * @param ph constant phase rotation
   */
  public void makeSyntheticSeismograms(
    String dataDir, double fp, double q, double ph) 
  {
    for (WellTie well : wellties) {
      if (_nsyn) well.makeNewSeismograms();
      well.makePropagatorSeismogram(dataDir,fp,q,ph);
    }
  }


  /**
   * Computes an initial synthetic image in 2D.
   * Uses previously created synthetic seismograms to compute a 
   * synthetic image by interpolating the synthetic seismograms 
   * between the wells guided by seismic image structure.
   * @param bgs smoothness of blended neighbor interpolation
   * @param p0 scalar emphasizing the amplitude of tensors
   * @param p1 scalar emphasizing the linearity of tensors
   * @return a synthetic image
   */
  public float[][] makeInitialSyntheticImage2( 
    float bgs, float p0, float p1) 
  {
    float pnull = -999.99f;
    float[][] f2 = new float[_n2][_n1];
    float[][] fv = fillfloat(pnull,_n1,_n2);
    for (WellTie well : wellties) {
      float x2 = well.x2f;
      int ix2 = inro((x2-_f2)/_d2);
      float[] f = (_swttm)?well.h:well.f;
      int ns = f.length;
      Sampling sf = (_swttm)?well.sh:well.sf;
      int ssyi = inro((sf.getFirst()-_f1)/_d1);
      syntheticImageFill(ix2,ns,ssyi,f,fv);
    }
    interpolate2(pnull,bgs,p0,p1,_s1,fv,f2);
    return f2;
  }

  /**
   * Computes an initial synthetic image in 3D.
   * Uses previously created synthetic seismograms to compute a 
   * synthetic image by interpolating the synthetic seismograms 
   * between the wells guided by seismic image structure.
   * @param bgs smoothness of blended neighbor interpolation
   * @param p0 scalar emphasizing the amplitude of tensors
   * @param p1 scalar emphasizing the linearity of tensors
   * @param p1 scalar emphasizing the planarity of tensors
   * @return a synthetic image
   */
  public float[][][] makeInitialSyntheticImage3(
    float bgs, float p0, float p1, float p2) 
  {
    float pnull = -999.99f;
    float[][][] f3 = new float[_n3][_n2][_n1];
    float[][][] fv = fillfloat(pnull,_n1,_n2,_n3);
    String fname = getfnameisi("synimg",bgs,p0,p1,p2);
    if (exists(fname) && !_misyni && _csmDir!=null) {
      System.out.println("Reading initial synthetic image...");
      readImage(fname,f3);
    } else {
      System.out.println("Making initial synthetic image...");
      for (WellTie well : wellties) {
        float x2 = well.x2f;
        float x3 = well.x3f;
        int ix2 = inro((x2-_f2)/_d2);
        int ix3 = inro((x3-_f3)/_d3);
        float[] f = (_swttm)?well.h:well.f;
        int ns = f.length;
        Sampling sf = (_swttm)?well.sh:well.sf;
        int ssyi = inro((sf.getFirst()-_f1)/_d1);
        syntheticImageFill(ix2,ix3,ns,ssyi,f,fv);
      }
      interpolate3(pnull,bgs,p0,p1,p2,_s1,fv,f3);
      writeImage(fname,f3);
    }
    return f3;
  }


  /**
   * Computes single-well ties. 
   * Use smooth dynamic time warping to compute a well tie. Refer 
   * to DynamicWarpingWT class for extra details.
   * @param r1min lower bound on strain in time.
   * @param r1max upper bound on strain in time.
   * @param smin minimum time shift between synthetic seismogram and seismic 
   * @param smax maximum time shift between synthetic seismogram and seismic
	 * @param dr1 time smoothing parameter
   */
  public void computeSingleWellTies(
    float r1min, float r1max, 
    float smin, float smax, float dr1)
  {
    for (WellTie well : wellties) {
      float x2 = well.x2f;
      float x3 = well.x3f;
      int ix2 = inro((x2-_f2)/_d2);
      int ix3 = inro((x3-_f3)/_d3);
      float[][] gt = getNearbyTraces(_g3,ix2,ix3);
      if (_fphase) well.findPhase();
      well.computeSingleWellTie(gt[0],r1min,r1max,smin,smax,dr1);
    }
  }

  /**
   * Computes simultaneous multiple-well ties in 2D. 
   * Use smooth dynamic image warping to compute simultaneous
   * multiple-well ties. For more details, refer 
   * to DynamicWarpingWT class. 
   * @param r1min lower bound on strain in time.
   * @param r1max upper bound on strain in time.
   * @param r2min lower bound on strain in x2.
   * @param r2max upper bound on strain in x2.
	 * @param dr1 time smoothing parameter
	 * @param dr2 distance smoothing parameter
   * @param smin minimum time shift between synthetic image and seismic 
   * @param smax maximum time shift between synthetic image and seismic
   * @param f    the initial synthetic image
   * @return the warped synthetic image
   */
  public float[][] computeMultipleWellTies2(
    float r1min, float r1max, 
    float r2min, float r2max, 
    float dr1,   float dr2,
    float smin, float smax, float[][] f)
  {
    DynamicWarpingWT dw = new DynamicWarpingWT(smin,smax,_s1,_s2);
    dw.setStrainLimits(r1min,r1max,r2min,r2max);
    dw.setSmoothing(dr1,dr2);
    float[][] u = dw.findShifts(_s1,f,_s1,_g2);
    float[][] h = dw.applyShifts(f,u);
    //float[][] u = dw.findShifts(_s1,_g2,_s1,f);
    //float[][] h = dw.applyShiftsX(f,u);
    _u2 = u;
    for (WellTie well : wellties) {
      Sampling sf = well.sf;
      if (_swttm) 
        sf = well.sh;
      double fs = sf.getFirst();
      float x2 = well.x2f;
      int ix2 = inro((x2-_f2)/_d2);
      int nx1 = sf.getCount();
      int jx1 = inro((fs-_f1)/_d1);
      float[] f1 = copy(nx1,jx1,f[ix2]);
      float[] u1 = copy(nx1,jx1,u[ix2]);
      float[] h1 = dw.applyShifts(f1,u1,sf);
      Sampling sh = new Sampling(h1.length,_d1,fs+u1[0]);
      if (_swttm) 
        sh = new Sampling(h1.length,_d1,fs+u1[0]);
      //int jx1 = infl((sf.getFirst()-_f1)/_d1);
      //float[] u1 = copy(nx1,jx1,u[ix2]);
      //float[] ua = add(floats(_s1.getValues()),u[ix2]);
      //float fh = ua[jx1];
      //int jx2 = infl((fh-_f1)/_d1);
      //float[] h1 = copy(nx1,jx2,h[ix2]);
      //Sampling sh = new Sampling(h1.length,_d1,fh);
      if (_swttm) 
        well.tz0 = copy(well.tz1);
      well.h = h1;
      well.sh = sh;
      well.update(u1,(float)fs);
    }
    return h;
  }

  /**
   * Computes simultaneous multiple-well ties in 3D. 
   * Use smooth dynamic image warping to compute simultaneous
   * multiple-well ties. For more details, refer 
   * to DynamicWarpingWT class. 
   * @param r1min lower bound on strain in time.
   * @param r1max upper bound on strain in time.
   * @param r2min lower bound on strain in x2.
   * @param r2max upper bound on strain in x2.
   * @param r3min lower bound on strain in x3.
   * @param r3max upper bound on strain in x3.
	 * @param dr1 time smoothing parameter
	 * @param dr2 distance smoothing parameter
	 * @param dr3 distance smoothing parameter
   * @param smin minimum time shift between synthetic image and seismic 
   * @param smax maximum time shift between synthetic image and seismic
   * @param f    the initial synthetic image
   * @return the warped synthetic image
   */
  public float[][][] computeMultipleWellTies3(
    float r1min, float r1max, 
    float r2min, float r2max, 
    float r3min, float r3max, 
    float dr1,   float dr2, float dr3,
    float smin, float smax, float[][][] f)
  {
    DynamicWarpingWT dw = new DynamicWarpingWT(smin,smax,_s1,_s2,_s3);
    dw.setStrainLimits(r1min,r1max,r2min,r2max,r3min,r3max);
    dw.setSmoothing(dr1,dr2,dr3);
    float[][][] u = new float[_n3][_n2][_n1];
    String fname = getfnamew(r1min,r1max,r2min,r2max,r3min,r3max,
                              dr1,dr2,dr3,smin,smax);
    if (exists(fname) && !_mwarps && _csmDir!=null) {
      System.out.println("Reading shifts...");
      readImage(fname,u);
    } else {
      System.out.println("Making shifts...");
      u = dw.findShifts(_s1,f,_s1,_g3);
      writeImage(fname,u);
    }
    float[][][] h = dw.applyShifts(f,u);
    _u3 = u;
    for (WellTie well : wellties) {
      Sampling sf = well.sf;
      if (_swttm) 
        sf = well.sh;
      double fs = sf.getFirst();
      float x2 = well.x2f;
      float x3 = well.x3f;
      int ix2 = inro((x2-_f2)/_d2);
      int ix3 = inro((x3-_f3)/_d3);
      int nx1 = sf.getCount();
      int jx1 = inro((fs-_f1)/_d1);
      float[] f1 = copy(nx1,jx1,f[ix3][ix2]);
      float[] u1 = copy(nx1,jx1,u[ix3][ix2]);
      float[] h1 = dw.applyShifts(f1,u1,sf);
      Sampling sh = new Sampling(h1.length,_d1,u1[0]+fs);
      if (_swttm) {
        sh = new Sampling(h1.length,_d1,fs+u1[0]);
      }
      well.h = h1;
      well.sh = sh;
      well.update(u1,(float)fs);
    }
    return h;
  }

  /**
   * Interpolates log velocities in 2D.
   * Interpolates between the wells with the 
   * initial time-depth function.
   */
  public float[][] interpolateInitialVelocity2(
    float bgs, float p0, float p1) 
  {
    return interpolateLogs2(bgs,p0,p1,"vi");
  }

  /**
   * Interpolates updated velocities in 2D.
   * Interpolates between the wells with the 
   * updated time-depth function.
   */
  public float[][] interpolateUpdatedVelocity2(
    float bgs, float p0, float p1) 
  {
    return interpolateLogs2(bgs,p0,p1,"vn");
  }

  /**
   * Interpolates log velocities in 2D. 
   * Interpolates between the wells with the 
   * updated time-depth function.
   */
  public float[][] interpolateInitialVelocityTied2(
    float bgs, float p0, float p1) 
  {
    return interpolateLogs2(bgs,p0,p1,"vs");
  }

  /**
   * Interpolates log velocities in 3D.
   * Interpolates between the wells with the 
   * initial time-depth function.
   */
  public float[][][] interpolateInitialVelocity3(
    float bgs, float p0, float p1, float p2) 
  {
    return interpolateLogs3(bgs,p0,p1,p2,"vi");
  }

  /**
   * Interpolates updated velocities in 3D.
   * Interpolates between the wells with the 
   * updated time-depth function.
   */
  public float[][][] interpolateUpdatedVelocity3(
    float bgs, float p0, float p1, float p2) 
  {
    return interpolateLogs3(bgs,p0,p1,p2,"vn");
  }

  /**
   * Interpolates log velocities in 3D.
   * Interpolates between the wells with the 
   * updated time-depth function.
   */
  public float[][][] interpolateInitialVelocityTied3(
    float bgs, float p0, float p1, float p2) 
  {
    return interpolateLogs3(bgs,p0,p1,p2,"vs");
  }

  /**
   * Interpolates average velocities in 2D.
   * Interpolates in depth.
   */
  public float[][] interpolateAverageVelocity2(
    float bgs, float p0, float p1, Sampling sd)
  {
    return interpolateTimeDepths2(bgs,p0,p1,sd,"vavg");
  }
  /**
   * Interpolates updated time-depth functions in 2D.
   * Interpolates in depth.
   */
  public float[][] interpolateUpdatedTimeDepths2(
    float bgs, float p0, float p1, Sampling sd)
  {
    return interpolateTimeDepths2(bgs,p0,p1,sd,"tz1");
  }
  /**
   * Interpolates initial time-depth functions in 2D.
   * Interpolates in depth.
   */
  public float[][] interpolateInitialTimeDepths2(
    float bgs, float p0, float p1, Sampling sd)
  {
    return interpolateTimeDepths2(bgs,p0,p1,sd,"tz0");
  }

  /**
   * Interpolates average velocities in 3D.
   * Interpolates in depth.
   */
  public float[][][] interpolateAverageVelocity3(
    float bgs, float p0, float p1, float p2, Sampling sd)
  {
    return interpolateTimeDepths3(bgs,p0,p1,p2,sd,"vavg");
  }
  /**
   * Interpolates updated time-depth functions in 3D.
   * Interpolates in depth.
   */
  public float[][][] interpolateUpdatedTimeDepths3(
    float bgs, float p0, float p1, float p2, Sampling sd)
  {
    return interpolateTimeDepths3(bgs,p0,p1,p2,sd,"tz1");
  }
  /**
   * Interpolates initial time-depth functions in 3D.
   * Interpolates in depth.
   */
  public float[][][] interpolateInitialTimeDepths3(
    float bgs, float p0, float p1, float p2, Sampling sd)
  {
    return interpolateTimeDepths3(bgs,p0,p1,p2,sd,"tz0");
  }

  
  public float[][] depthConvertSeismic2(float[][] tz) {
    return convertTimeToDepth2(_g2,tz);
  } 
  public float[][][] depthConvertSeismic3(float[][][] tz) {
    return convertTimeToDepth3(_g3,tz);
  }

  /**
   * Converts 2D image to depth.
   */
  public float[][] convertTimeToDepth2(float[][] ft, float[][] tz) {
	  int nz = tz[0].length;
	  float[][] fz = new float[_n2][nz];
	  SincInterp si = new SincInterp();
	  double[] x2 = _s2.getValues();
	  for (int i2=0; i2<_n2; ++i2)
	    for (int i1=0; i1<nz; ++i1)
	  		fz[i2][i1] = si.interpolate(_s1,_s2,ft,tz[i2][i1],x2[i2]);
	  return fz;
  }

  /**
   * Converts 3D image to depth.
   */
  public float[][][] convertTimeToDepth3(
    final float[][][] ft, final float[][][] tz) 
  {
    final Sampling s3 = _s3;
    final Sampling s2 = _s2;
    final Sampling s1 = _s1;
    final int n3 = _n3;
    final int n2 = _n2;
	  final int nz = tz[0][0].length;
	  final float[][][] fz = new float[n3][n2][nz];
	  final SincInterp si = new SincInterp();
	  final double[] x2 = _s2.getValues();
	  final double[] x3 = _s3.getValues();
    Parallel.loop(n3,new Parallel.LoopInt() {
      public void compute(int i3) {
	      for (int i2=0; i2<n2; ++i2)
	        for (int i1=0; i1<nz; ++i1)
	  	  		fz[i3][i2][i1] = si.interpolate(s1,s2,s3,ft,
                               tz[i3][i2][i1],x2[i2],x3[i3]);
     }});
	  return fz;
  }

 public float[] maketz1(Sampling sd, float[] va) {
    int n1 = va.length;
    float[] tz = new float[n1];
    for (int i1=0; i1<n1; ++i1) 
      tz[i1] = 2.0f*(float)sd.getValue(i1)/va[i1];
    return tz;
  }

  public float[][] maketz2(Sampling sd, float[][] va) {
    int n2 = va.length;
    float[][] tz = new float[n2][];
    for (int i2=0; i2<n2; ++i2) 
      tz[i2] = maketz1(sd,va[i2]);
    return tz;
  }

  public float[][][] maketz3(Sampling sd, float[][][] va) {
    int n3 = va.length;
    float[][][] tz = new float[n3][][];
    for (int i3=0; i3<n3; ++i3) 
      tz[i3] = maketz2(sd,va[i3]);
    return tz;
  }

  /////////////////////////////////////////////////////////////////////////////
  // research 
  
  public float[] extrapolateLinearX(Sampling sz,Sampling sd,float[] f1) {
    return extrapolateLinear(sz,sd,f1);
  }
  public void tzeXOn() {
    _tzext = true;
  }
  public void tzeXOff() {
    _tzext = false;
  }

  /////////////////////////////////////////////////////////////////////////////
  // parameter setting and modifications

  public void addWell(WellTie well) {
    wellties.add(well);
  }

  public void removeWell(long id) {
    for (WellTie well : wellties) 
      if(well.id==id) wellties.remove(well);
  }

  public void replaceWell(WellTie nwell) {
    for (WellTie well : wellties) 
      if(well.id==nwell.id) 
        wellties.set(wellties.indexOf(well),nwell);
  }

  public void replaceWells(ArrayList<WellTie> wl) {
    for (WellTie well : wl) 
      replaceWell(well);
  }

  public WellTie getWell(long id) {
    WellTie rwell = new WellTie();
    for (WellTie well : wellties) 
      if(well.id==id) rwell = well;
    return rwell;
  }

  public float[][] getWellCoordinates(long[] ids) {
    FloatList x2 = new FloatList();
    FloatList x3 = new FloatList();
    for (WellTie well : wellties) {
      if(well.id==Arrays.binarySearch(ids,well.id)) {
        x2.add(well.x2f);
        x3.add(well.x3f);
      }
    }
    return new float[][]{x2.trim(),x3.trim()} ;
  }

  public float[][] getWellCoordinates() {
    FloatList x2 = new FloatList();
    FloatList x3 = new FloatList();
    for (WellTie well : wellties) {
      x2.add(well.x2f);
      x3.add(well.x3f);
    }
    return new float[][]{x2.trim(),x3.trim()} ;
  }
 
  public void fixX2(float[] x2n) {
    int i=0;
    for (WellTie well : wellties) {
      well.x2f = x2n[i];
      i+=1;
    }
  }

  public double getMaxz() {
    double mz=0.0;
    for (WellTie well : wellties) {
      double z = well.sz.getLast();
      if (z>mz) 
        mz = z;
    }
    return mz;
  }

  public int getNumberOfWells() {
    return wellties.size();
  }

  public float[][] getShifts2() {
    return _u2;
  }

  public float[][][] getShifts3() {
    return _u3;
  }

  public float[][] getSeismic2() {
    return _g2;
  }

  public float[][][] getSeismic3() {
    return _g3;
  }

  public void makeNewSyntheticSeismograms() {
    _nsyn = true;
  }

  public void findPhases() {
    _fphase = true;
  }

  public void makeNewSyntheticImage() {
    _misyni = true;
  }

 public void makeNewTimeDepthFunction() {
    _mtz = true;
  }

  public void makeNewMultipleWellTie() {
    _mwarps = true;
  }

  public void makeNewInterpolatedLogs() {
    _milogs = true;
  }

  // only use this for research purposes
  public void computeSingleThenMultiple() {
    _swttm = true;
  }

  public void setStatic(long id) {
    _static.add(id);
    WellTie well = getWell(id);
    well.setStatic(true);
    well.maketz0(well.v);
    replaceWell(well);
  }

  public List<WellTie> wellties = new ArrayList<WellTie>();


////////////////////////////////////////////////////////////////////////////////
// private

  private Sampling _s1,_s2,_s3; // seismic samplings
  private int  _n1,_n2,_n3; // seismic samplings
  private double  _d1,_d2,_d3,_f1,_f2,_f3; // seismic samplings
  private double _dz; // log sample rate
  private long[] _ids; // well ids
  private float[][] _g2; // 2D seismic image
  private float[][][] _g3; // 3D seismic image 
  private float[][] _u2; // 2D computed time-shifts
  private float[][][] _u3; // 3D computed time-shifts
  private boolean _tens = true; // use structure oriented tensors
  private boolean _nsyn = false; // make new synthetic seismograms
  private boolean _swttm = false; // compute single well ties then multiple
  private boolean _misyni = false; // compute new synthetic image
  private boolean _fphase = false; // finds the optimal phase of synthetics
  private boolean _mwarps = false; // compute new multiple well tie
  private boolean _milogs = false; // compute new interpolated logs
  private boolean _tzext = false; // t(z) functions being interp and extrap
  private boolean _mtz = false; // compute new interpolated t(z) functions
  private double _slof = 4.0; // sigma used in lof for tensors
  private double _gtmax = 100.0; // maximum distance for interpolation
  private String _csmDir; // root directory to store data
  private String _isifname; // file name for the initial synthetic image
  private String _wfname; // file name for the initial synthetic image
  private LongList _static;


  private float[][] interpolateLogs2(
    float bgs, float p0, float p1, String which)
  {
    // interpolate time-converted logs
    float pnull = -999.99f;
    float[] ta = floats(_s1.getValues());
    float[][] f2 = new float[_n2][_n1];
    float[][] fv = fillfloat(pnull,_n1,_n2);
    for (WellTie well : wellties) {
      float x2 = well.x2f;
      float x3 = well.x3f;
      int ix2 = inro((x2-_f2)/_d2);
      int ix3 = inro((x3-_f3)/_d3);
      float[] l = (which=="vn")?well.v1:well.v;
      float[] tz = (which=="vi")?well.tz0:well.tz1;
      int nz = l.length;
      logImageFill(ix2,_n1,nz,pnull,l,tz,ta,fv);
    }
    interpolate2(pnull,bgs,p0,p1,_s1,fv,f2);
    return f2;
  }

  private float[][][] interpolateLogs3(
    float bgs, float p0, float p1, float p2, String which)
  {
    // interpolate time-converted logs
    float pnull = -999.99f;
    float[] ta = floats(_s1.getValues());
    float[][][] f3 = new float[_n3][_n2][_n1];
    float[][][] fv = fillfloat(pnull,_n1,_n2,_n3);
    String fname = getfnamel(which,bgs,p0,p1,p2);
    if (exists(fname) && !_milogs && _csmDir!=null) {
      System.out.println("Reading "+which+"...");
      readImage(fname,f3);
    } else {
      System.out.println("Making "+which+"...");
      for (WellTie well : wellties) {
        float x2 = well.x2f;
        float x3 = well.x3f;
        int ix2 = inro((x2-_f2)/_d2);
        int ix3 = inro((x3-_f3)/_d3);
        float[] l = (which=="vn")?well.v1:well.v;
        float[] tz = (which=="vi")?well.tz0:well.tz1;
        int nz = l.length;
        logImageFill(ix2,ix3,_n1,nz,pnull,l,tz,ta,fv);
      }
      interpolate3(pnull,bgs,p0,p1,p2,_s1,fv,f3);
      writeImage(fname,f3);
    }
    return f3;
  }

  private float[][] interpolateTimeDepths2(
    float bgs, float p0, float p1, Sampling sd, String which)
  {
    int nz = sd.getCount();
    //float pnull = -999.99f;
    float pnull = 0f;
    float[] tsdv = mul(2f,floats(sd.getValues()));
    float[][] f2 = new float[_n2][nz];
    float[][] fv = fillfloat(pnull,nz,_n2);
    String wh = new String(which);
    if (which=="tvavg") return interpolateVavg(bgs,p0,p1,sd);
    if (wh.startsWith("t")) 
      _tzext = true;
    else
      _tzext = false;
    int nw = wellties.size();
    float[] x2a = new float[nw];int i=0;
    float[] x1a = floats(sd.getValues());
    float[][] fa = fillfloat(pnull,nz,nw);
    for (WellTie well : wellties) {
      Sampling sz = well.sz;
      float x2 = well.x2f;
      int ix2 = inro((x2-_f2)/_d2);
      x2a[i] = x2;
      float[] f1 = (which=="vavg")?computeVavg(sz,well.tz1)
                  :((which=="tz1")?well.tz1:well.tz0);

      float[] fe = extrapolateLinear(sz,sd,f1,pnull);
      //fv[ix2] = (which=="vavg")?maketz1(sd,fe):copy(fe);
      fa[i] = copy(fe);
      ++i;
    }
    _tzext = true;
    //interpolate2(pnull,bgs,p0,p1,sd,fv,f2);
    BicubicInterpolator2 ci2 = new 
      BicubicInterpolator2(BicubicInterpolator2.Method.MONOTONIC,
													 BicubicInterpolator2.Method.MONOTONIC,
													 nz,nw,x1a,x2a,fa); 

  	ci2.interpolate(sd,_s2,f2);
    if (which=="vavg") return maketz2(sd,f2);
    return f2;
  }

  private float[][][] interpolateTimeDepths3(
    float bgs, float p0, float p1, float p2, Sampling sd, String which)
  {
    int nz = sd.getCount();
    float pnull = -999.99f;
    float[][][] f3 = new float[_n3][_n2][nz];
    float[][][] fv = fillfloat(pnull,nz,_n2,_n3);
    String fname = getfnamed(which,bgs,p0,p1,p2);
    if (exists(fname) && !_mtz && _csmDir!=null) {
      System.out.println("Reading "+which+"...");
      readImage(fname,f3);
    } else {
      System.out.println("Interpolating "+which+"...");
      String wh = new String(which);
      if (which=="tvavg") return interpolateVavg(bgs,p0,p1,p2,sd,fname);
      if (wh.startsWith("t")) 
        _tzext = true;
      else
        _tzext = false;
      for (WellTie well : wellties) {
        Sampling sz = well.sz;
        float x2 = well.x2f;
        float x3 = well.x3f;
        int ix2 = inro((x2-_f2)/_d2);
        int ix3 = inro((x3-_f3)/_d3);
        float[] f1 = (which=="vavg")?computeVavg(sz,well.tz1)
                    :((which=="tz1")?well.tz1:well.tz0);
        float[] fe = extrapolateLinear(sz,sd,f1);
        fv[ix3][ix2] = (which=="vavg")?maketz1(sd,fe):copy(fe);
      }
      _tzext = true;
      interpolate3(pnull,bgs,p0,p1,p2,sd,fv,f3);
      writeImage(fname,f3);
    }
    return f3;
  }

  private float[][] interpolateVavg(
    float bgs, float p0, float p1, Sampling sd)
  {
    int nz = sd.getCount();
    float pnull = -999.99f;
    float[] ta = floats(_s1.getValues());
    float[] tsdv = mul(2f,floats(sd.getValues()));
    float[][] f2t = new float[_n2][_n1];
    float[][] fvt = fillfloat(pnull,_n1,_n2);
    float[][] f2z = new float[_n2][nz];
    for (WellTie well : wellties) {
      Sampling sz = well.sz;
      float x2 = well.x2f;
      int ix2 = inro((x2-_f2)/_d2);
      float[] f1 = computeVavg(sz,well.tz1);
      float[] fe = extrapolateLinear(sz,sd,f1);
      float[] tze = maketz1(sd,fe);
      float[] zt = getzt(tze,_s1,sz);
      f2t[ix2] = toTime(sd,fe,zt);
    }
    interpolate2(pnull,bgs,p0,p1,_s1,fvt,f2t);
    //float[][] zt2 = mul(0.5f,mul(floats(_s1.getValues()),f2t));
    //float[][] f2z = toDepth(_s1,sd,_s2,f2t,gettz(zt2,sd,_s1));
    return f2z;
  }

  public float[][][] interpolateVavg(
    float bgs, float p0, float p1, float p2, Sampling sd, String fname)
  {
    int nz = sd.getCount();
    float pnull = -999.99f;
    float[][][] f3 = new float[_n3][_n2][nz];
    float[][][] fv = fillfloat(pnull,nz,_n2,_n3);
    int nw = wellties.size(); int i=0;
    float[][][] fa = new float[nw][nw][nz];
    float[] x2a = new float[nw];
    float[] x3a = new float[nw];
    float[] x1a = floats(sd.getValues());
    for (WellTie well : wellties) {
      Sampling sz = well.sz;
      float x2 = well.x2f;
      float x3 = well.x3f;
      int ix2 = inro((x2-_f2)/_d2);
      int ix3 = inro((x3-_f3)/_d3);
      float[] f1 = computeVavg(sz,well.tz1);
      float[] fe = extrapolateLinear(sz,sd,f1);
      fv[ix3][ix2] = maketz1(sd,fe);
      x2a[i] = x2;
      x3a[i] = x3;
      fa[i][i] = copy(fe);
      ++i;
    }
    _tzext = true;
    //interpolate3(pnull,bgs,p0,p1,p2,sd,fv,f3);
    //writeImage(fname,f3);
    System.out.println("Tri-cubic Interpolation...");
    //float[][] fxy = SimpleGridder3.getGriddedSamples(pnull,sd,_s2,_s3,fv); 
    //TricubicInterpolator3 ci3 = new 
    //  TricubicInterpolator3(TricubicInterpolator3.Method.MONOTONIC,
		//											 TricubicInterpolator3.Method.MONOTONIC,
		//											 TricubicInterpolator3.Method.MONOTONIC,
		//											 //fxy[1],fxy[2],fxy[3],fv); 
		//											 x1a,x2a,x3a,fv); 
    //sort(x2a,x3a);
    TrilinearInterpolator3 ci3 = new TrilinearInterpolator3(x1a,x2a,x3a,fv);
													 //fxy[1],fxy[2],fxy[3],fv); 
  	ci3.interpolate(sd,_s2,_s3,f3);
    return f3;
  }
  // implemented sort
  private static void sort(float[] x1, float[] x2) {
    dump(x1);
    dump(x2);
    int n = x1.length;
    float[] x1t = new float[n];
    float[] x2t = new float[n];
    float x1u = x1[0];
    float x2u = x2[0];
    float x1d = 0f;
    float x2d = 0f;
    x1t[0] = x1d;
    x2t[0] = x2d;
    for (int i=1; i<n; ++i) {
      float x11 = x1[i];
      float x21 = x2[i];
      for (int j=0; j<n; ++j) {
        float x12 = x1[j];
        float x22 = x2[j];
        if (x12>=x1d && x12<x11 && x22>=x1d && x22<x21) {
          x1u = x11;
          x2u = x21;
          x1t[i] = x12;
          x2t[i] = x22;
        }
      }
    }
    copy(x1t,x1);
    copy(x2t,x2);
    dump(x1);
    dump(x2);
  }
  private static float[][] toDepth(
    Sampling st, Sampling sd, Sampling s2, float[][] v, float[][] tz)
  {
    int n2 = v.length;
	  int nz = sd.getCount();
	  int nt = st.getCount();
    float[] x2 = floats(s2.getValues());
	  float[][] fz = new float[n2][nz];
	  SincInterp si = new SincInterp();
	  for (int i2=0; i2<n2; ++i2)
	    for (int i1=0; i1<nz; ++i1)
	 		  fz[i2][i1] = si.interpolate(st,s2,v,tz[i2][i1],x2[i1]);
	  return fz;
  }
  private static float[] toTime(Sampling sz, float[] fz, float[] zt) {
   int nt = zt.length;
	  float[] ft = new float[nt];
	  SincInterp si = new SincInterp();
	  for (int i1=0; i1<nt; ++i1)
	 		ft[i1] = si.interpolate(sz,fz,zt[i1]);
	  return ft;
  }
  private static float[] getzt(float[] tz, Sampling st, Sampling sz) {
		int nt = st.getCount();
		float dt = (float)st.getDelta();
		float ft = (float)st.getFirst();
		int nz = sz.getCount();
		float dz = (float)sz.getDelta();
		float fz = (float)sz.getFirst();
		float[] zt = new float[nt];	
		inverseLinearInterpolation(nz,dz,fz,tz,nt,dt,tz[0],zt,tz[0],tz[nz-1]);
    return zt;
  }
  private static float[] gettz(float[] zt, Sampling sz, Sampling st) {
    return getzt(zt,st,sz);
  }


  /**
   * Computes a linear extrapolation for the ends of a function.
   * Also resamples if dz1 and dz2 are different.
   */
  private float[] extrapolateLinear(
    Sampling sz1, Sampling sz2, float[] l, float pnull)
  {
    int nz1 = sz1.getCount(); 
    int nz2 = sz2.getCount(); 
    double dz1 = sz1.getDelta(); 
    double dz2 = sz2.getDelta(); 
    double fz1 = sz1.getFirst(); 
    double fz2 = sz2.getFirst(); 
    double lz2 = sz2.getLast();
    double lz1 = sz1.getLast();
    float fdz2 = (float)dz2;
    float[] le = new float[nz2];
    float[] lr = l;
    if (dz1!=dz2) {
      lr = resampleLinear(sz1,sz2,l);
      sz1 = new Sampling(lr.length,dz2,fz1);
    }
    // extrapolate beginning
    int nr = (int)(lr.length*0.5);
    int nb = infl((fz1-fz2)/dz2);
    int nm = lr.length;
    int ne = nz2-(nb+nm);
    if (pnull<0) {
      for (int i=0; i<nb; ++i)
        le[i] = pnull;
    }
    else if (nb>0 && !_tzext) {
      float[] bb = getRegressionTerms(sz1,lr,nr,1);
      float b0 = l[0];
      float b1 = -bb[1];
      for (int i=0; i<nb; ++i)
        le[nb-i-1] = b0+i*b1*fdz2;
    }
    else if (nb>0 && _tzext) {
      // if time-depth curve, extrapolate to t=z=0
      float b0 = (float)0;
      float b1 = (l[0]-b0)/nb;
      for (int i=0; i<nb; ++i)
        le[i] = b0+i*b1;
    }
    // fill values
    for (int i=0; i<nm; ++i)
      le[i+nb] = lr[i];
    // extrapolate end
    if (pnull<0) {
      int nf = nb+nm;
      for (int i=0; i<ne; ++i)
        le[i+nf] = pnull;
    } 
    else if (ne>0) {
      int nf = nb+nm;
      float[] bb = getRegressionTerms(sz1,lr,nr,-1);
      float b2 = bb[1];
      float fltz = lr[nm-1];
      for (int i=0; i<ne; ++i)
        le[i+nf] = fltz+i*b2*fdz2;
    }
    return le;
  }
  private float[] extrapolateLinear(
    Sampling sz1, Sampling sz2, float[] l) 
  {
    return extrapolateLinear(sz1,sz2,l,1);
  }
  
  /**
   * Resamples a function using linear interpolation.
   */
  public float[] resampleLinear(Sampling sz1, Sampling sz2, float[] x) { 
    int nz1 = sz1.getCount(); 
    double dz2 = sz2.getDelta(); 
    double dz1 = sz1.getDelta(); 
    double fz1 = sz1.getFirst(); 
    double lz1 = sz1.getLast();
    int nzi = inro((lz1-fz1)/dz2);
    float[] xi = new float[nzi];
    LinearInterpolator li = new LinearInterpolator(); 
    li.setUniform(nz1,dz1,fz1,x);
	  li.interpolate(nzi,dz2,fz1,xi);
    return xi;
  }
  
  /**
   * Solves for simple linear regression terms.
   */
  private float[] getRegressionTerms(Sampling sz, float[] l, int nr, int end) {
    int nl = l.length;
    float[] zvl = floats(sz.getValues());
    int nre = nl-nr;
    int n1 = (end>0)?0:nre;
    int n2 = (end>0)?nr:nl;
    float ym = meanW(l,n1,n2);
    float xm = meanW(zvl,n1,n2);
    float num=0f,den=0f;
    for (int j=n1; j<n2; ++j) {
      float xi = zvl[j]-xm;
      num += (l[j]-ym)*xi;
      den += xi*xi;
    }
    float b1 = num/den;
    float b0 = ym - b1*xm;
    return new float[]{b0,b1};
  }
  private static float meanW(float[] x, int w1, int w2) {
	  int n = w2-w1;
	  float[] xw = copy(n,w1,x);
	  return sum(xw)/n;
  }
  private float[] getTrend(float sigma, int niter, float[] x) {
    int n = x.length;
    float[] y = new float[n];
    ExponentialSmoother es = new ExponentialSmoother(sigma,niter);
    es.apply(x,y);
    return y;
  }

	public static float[] logImageFill(
		int nt, int nz, float pnull,
		float[] l, float[] tz, float[] ta, float[] g) 
	{
		return logImageFill(0,nt,nz,pnull,l,tz,ta,new float[][]{g})[0];
	}
	private static float[][] logImageFill(
		int xi, int nt, int nz, float pnull,
		float[] l, float[] tz, float[] ta, float[][] g) 
	{
		return logImageFill(xi,0,nt,nz,pnull,l,tz,ta,new float[][][]{g})[0];
	}
	private static float[][][] logImageFill(
		int x2i, int x3i, int nt, int nz, float pnull,
		float[] l, float[] tz, float[] ta, float[][][] g) 
	{
		int ntm = nt-1;
		float zc = 0.0f;
		float ssum = 0.0f;
		float vel = pnull;
		float[] sl = div(1.0f,l);
		for (int it=0; it<ntm; ++it) {
			for (int iz=0; iz<nz; ++iz) {
				float tzv = tz[iz];
				if (tzv>ta[it] && tzv<=ta[it+1]) {
					ssum += sl[iz];
					zc += 1.0f;
				}
				if (tzv>ta[it+1]) 
					break;
			}
			if (ssum>0.0f) 
				vel = 1.0f/(ssum/zc);
			g[x3i][x2i][it] = vel;
			zc=0f;ssum=0f;vel=pnull;
		}
		return g;
	}

  public void interpolate2(
    float pnull, float bgs, float p0, float p1,
    Sampling s1, float[][] fv, float[][] f2) 
  {
    float[][] fxy = SimpleGridder2.getGriddedSamples(pnull,s1,_s2,fv); 
    BlendedGridder2 bg;
    if (_tens && !_tzext) {
      Tensors2 ts2 = makeTensors2(_g2,p0,p1);
      bg = new BlendedGridder2(ts2,fxy[0],fxy[1],fxy[2]);
    } else {
      bg = new BlendedGridder2(fxy[0],fxy[1],fxy[2]);
      _tzext = false;
    }
	  bg.setTimeMax(_gtmax);
	  bg.setSmoothness(bgs);
	  float[][] ft = bg.gridNearest(pnull,fv);
	  bg.gridBlended(ft,fv,f2);
  }

  public void interpolate3(
    float pnull, float bgs, float p0, float p1, float p2,
    Sampling s1, float[][][] fv, float[][][] f3) 
  {
    float[][] fxy = SimpleGridder3.getGriddedSamples(pnull,s1,_s2,_s3,fv); 
    BlendedGridder3 bg;
    if (_tens && !_tzext) {
      Tensors3 ts3 = makeTensors3(_g3,p0,p1,p2);
      bg = new BlendedGridder3(ts3,fxy[0],fxy[1],fxy[2],fxy[3]);
    } else {
      bg = new BlendedGridder3(fxy[0],fxy[1],fxy[2],fxy[3]);
      _tzext = false;
    }
	  bg.setTimeMax(_gtmax);
	  bg.setSmoothness(bgs);
	  float[][][] ft = bg.gridNearest(pnull,fv);
	  bg.gridBlended(ft,fv,f3);
  }
    
  private static void syntheticImageFill(
		int xi, int ns, int f, float[] s, float[][] g) 
	{
		for (int i=0; i<ns; ++i) 
			g[xi][f+i] = s[i];
	}
	private static void syntheticImageFill(
		int x2i, int x3i, int ns, int f, float[] s, float[][][] g) 
	{
		for (int i=0; i<ns; ++i) 
			g[x3i][x2i][f+i] = s[i];
	}

  private static float[] computeVavg(Sampling sz, float[] tz) {
    float[] szv = floats(sz.getValues());
    return div(mul(2,szv),tz);
  }

  private static float[][] getNearbyTraces(
    float[][][] g, int x2, int x3, int nr) 
  {
    // returns adjacent traces to well of radius nr
    int n3 = g.length;
    int n2 = g[0].length;
    int n1 = g[0][0].length;
	  int nt = 9+12*(nr-1)+(nr-1)*(nr-1)*4;
    float[][] gt = new float[nt][n1]; 
    int it=1;
    for (int i3=-nr; i3<=nr; ++i3) {
      for (int i2=-nr; i2<=nr; ++i2) {
        if (i3!=0 || i2!=0) {
          if (x2+i2>=0 && x3+i3>=0 && x2+i2<n2 && x3+i3<n3)
            gt[it++] = g[x3+i3][x2+i2];
        }
      }
    }
    gt[0] = g[x3][x2]; // puts closest trace at 0
    return gt;  
  }
  private static float[][] getNearbyTraces(float[][][] g, int x2, int x3) {
    return getNearbyTraces(g,x2,x3,1);
  }

  private void computeMaxWellDistance() {
    int n = wellties.size();
    float[] x2 = new float[n];
    float[] x3 = new float[n];
    int k=0;
    for (WellTie well : wellties) {
      x2[k] = well.x2f;
      x3[k] = well.x3f;
      k+=1;
    }
    float dmax = 0;
    for (int i2=0; i2<n; ++i2) {
      for (int i1=0; i1<n; ++i1) {
        float d = hypot(x2[i1]-x2[i2],x3[i1]-x3[i2]);
        if (d>dmax) 
          dmax = d;
      }
    }
    _gtmax = dmax;
  }
  
  private EigenTensors2 makeTensors2(
    float[][] g, float p0, float p1) 
  {
	  LocalOrientFilter lof = new LocalOrientFilter(_slof);
    EigenTensors2 d = lof.applyForTensors(g);
	  d.invertStructure(p0,p1);
    return d;
  }

  private EigenTensors3 makeTensors3(
    float[][][] g, float p0, float p1, float p2) 
  {
    String fname = _csmDir+"tensors/tpts_"
      +Float.toString(p0)+"_"
      +Float.toString(p1)+"_"
      +Float.toString(p2)+"_"
      +Float.toString((float)_n3)+"_"
      +Float.toString((float)_n2)+"_"
      +Float.toString((float)_n1)+".dat";
    EigenTensors3 d;
    if (exists(fname) && _csmDir!=null) {
      System.out.println("Reading tensors...");
      d = readTensors(fname);
    } else {
      System.out.println("Making tensors...");
	    LocalOrientFilter lof = new LocalOrientFilter(_slof);
      d = lof.applyForTensors(g);
	    d.invertStructure(p0,p1,p2);
      writeTensors(fname,d);
    }
    return d;
  }


	private static int inro(float x) {
		return round(x);
	}
	private static int inro(double x) {
		return (int)round(x);
	}
	private static int infl(float x) {
		return (int)floor(x);
	}
	private static int infl(double x) {
		return (int)floor(x);
	}
	private static int ince(float x) {
		return (int)ceil(x);
	}
	private static int ince(double x) {
		return (int)ceil(x);
	}
  private static float[][] floats(double[][] x) {
     int n2 = x.length;
     int n1 = x[0].length;
     float[][] y = new float[n2][n1];
     for (int i2=0; i2<n2; ++i2) 
       for (int i1=0; i1<n1; ++i1)
         y[i2][i1] = (float)x[i2][i1];
     return y;
   }
   private static float[] floats(double[] x) {
     int n1 = x.length;
     float[] y = new float[n1];
     for (int i1=0; i1<n1; ++i1) 
       y[i1] = (float)x[i1];
     return y;
   }

  private static void inverseLinearInterpolation(
    int nx, float dx, float fx, float[] y, 
    int ny, float dy, float fy, float[] x, float xylo, float xyhi) 
  { 
    int nxi,nyo,jxi1,jxi2,jyo;
    float dxi,fxi,dyo,fyo,fyi,yo,xi1,yi1,yi2,yid,q; 
    nxi = nx; dxi = dx; fxi = fx;
    nyo = ny; dyo = dy; fyo = fy;
    fyi = y[0];
    // loop over output y less than smallest input y
    for (jyo=0,yo=fyo; jyo<nyo; jyo++,yo+=dyo) {
      if (yo>=fyi) break;
      x[jyo] = xylo;
    }
    // loop over output y between smallest and largest input y
    if (jyo==nyo-1 && yo==fyi) {
      x[jyo++] = fxi;
      yo += dyo;
    }
    jxi1 = 0;
    jxi2 = 1;
    xi1 = fxi;
    while (jxi2<nxi && jyo<nyo) {
      yi1 = y[jxi1];
      yi2 = y[jxi2];
      if (yi1<=yo && yo<=yi2) {
        yid = abs(yi2-yi1);
	x[jyo++] = (yid>0.0f)?xi1+dxi*(yo-yi1)/(yi2-yi1):xi1;
        yo += dyo;
    } else {
        jxi1++;
        jxi2++;
        xi1 += dxi;
      }
    }
    // loop over output y greater than largest input y
    while (jyo<nyo) x[jyo++] = xyhi;
  }


  /////////////////////////////////////////////////////////////////////////////
  // I/O 

  private float[][][] readImage(String fname, float[][][] f) {
    try {
      ArrayInputStream ais = new ArrayInputStream(fname);
      ais.readFloats(f);
      ais.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
    return f;
  }

  private void writeImage(String fname, float[][][] f) {
    try {
      ArrayOutputStream aos = new ArrayOutputStream(fname);
      aos.writeFloats(f);
      aos.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  private EigenTensors3 readTensors(String fname) {
    EigenTensors3 d;
    try {
      FileInputStream fis = new FileInputStream(fname);
      ObjectInputStream ois = new ObjectInputStream(fis);
      d = (EigenTensors3)ois.readObject();
      ois.close();
      fis.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
    return d;
  }

  private void writeTensors(String fname, EigenTensors3 d) {
    try {
      FileOutputStream fos = new FileOutputStream(fname);
      ObjectOutputStream oos = new ObjectOutputStream(fos);
      oos.writeObject(d);
      oos.close();
      fos.close();
    } catch (Exception e) {
      throw new RuntimeException(e);
    }
  }

  private static boolean exists(String fname) {
    return new File(fname).exists();
  }

  private String getfnameisi(
    String kind, float bgs, float p0, float p1, float p2) 
  {
    String s1 = Float.toString(bgs);
    String s2 = Float.toString(p0);
    String s3 = Float.toString(p1);
    String s4 = Float.toString(p2);
    String s5 = Float.toString((float)_n1);
    String s6 = Float.toString((float)_n2);
    String s7 = Float.toString((float)_n3);
    String s8 = Double.toString(_f1);
    String s9 = Double.toString(_f2);
    String s10 = Double.toString(_f3);
    String s11 = Boolean.toString(_swttm);
    String s12 = Boolean.toString(_tens);
    String s13 = Float.toString((float)_ids.length);
    String s14 = Float.toString((float)sum(_ids));
    String d = "_";
    String fname = s1+d+s2+d+s3+d+s4+d+s5+d
                   +s6+d+s7+d+s8+d+s9+d+s10+d
                   +s11+d+s12+d+s13+d+s14;
    _isifname = fname;
    return _csmDir+"interps/"+fname+".dat";
  }
  private String getfnamew(
    float r1min, float r1max, float r2min, float r2max, float r3min, 
    float r3max, float dr1, float dr2, float dr3, float smin, float smax)
  {
    String s1 = Float.toString(r1min);
    String s2 = Float.toString(r1max);
    String s3 = Float.toString(r2min);
    String s4 = Float.toString(r2max);
    String s5 = Float.toString(r3min);
    String s6 = Float.toString(r3max);
    String s7 = Float.toString(dr1);
    String s8 = Float.toString(dr2);
    String s9 = Float.toString(dr3);
    String s10 = Float.toString(smin);
    String s11 = Float.toString(smax);
    String d = "_";
    String fname = "u3"+d+_isifname+d
                   +s1+d+s2+d+s3+d+s4+d+s5+d
                   +s6+d+s7+d+s8+d+s9+d+s10+d+s11;
    _wfname = fname;
    return _csmDir+"warps/"+fname+".dat";
  }
  private String getfnamel(
    String kind, float bgs, float p0, float p1, float p2) 
  {
    String s1 = Float.toString(bgs);
    String s2 = Float.toString(p0);
    String s3 = Float.toString(p1);
    String s4 = Float.toString(p2);
    String d = "_";
    String fname = _csmDir+"interps/"+kind+d+_wfname+d
                          +s1+d+s2+d+s3+d+s4;
    return fname+".dat";
  }
  private String getfnamed(
    String kind, float bgs, float p0, float p1, float p2) 
  {
    String s1 = Float.toString(bgs);
    String s2 = Float.toString(p0);
    String s3 = Float.toString(p1);
    String s4 = Float.toString(p2);
    String d = "_";
    String fname = _csmDir+"interps/"+kind+d+_wfname+d
                          +s1+d+s2+d+s3+d+s4;
    return fname+".dat";
  }

};
