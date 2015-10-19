package tp;

import java.util.*;
import java.io.*;
import edu.mines.jtk.io.*;
import edu.mines.jtk.dsp.*;
import edu.mines.jtk.util.*;
import edu.mines.jtk.interp.*;
import edu.mines.jtk.mosaic.*;
import static edu.mines.jtk.util.ArrayMath.*;

/**
 * Finds and returns TP well tops for specified TP well API numbers.
 * The only well tops that are found for each well API are the ones
 * corresponding to the seismic time horzions provided in the teapot
 * dome data set.
 * @author Andrew Munoz, CSM
 * @version 2.8.2013
 */

public class WellTops {

	public WellTops(long[] ids, String[] fmcodes, String filepath) {
		_ids = ids;
		_fms = fmcodes;
		//read the well log file and fill _wells 
		int ns = ids.length;
		int nf = fmcodes.length;
		String fileName = filepath+"TeapotDomeFormationLogTops.txt";
		String datumName = filepath+"TeapotDomeWellDatums.csv";
		// Get tops from file
		try {
      FileInputStream fis = new FileInputStream(fileName);
      Scanner s = new Scanner(fis); 
      while (s.hasNextLine()) {
				String line = s.nextLine();
				String[] fields = line.split(",");	
        if (fields.length>4)
          continue;
				long id = Long.parseLong(fields[0],10);
				String fmtop = fields[2];
				double tdepth = Double.parseDouble(fields[3]);
				for (int ii=0; ii<ns; ++ii) {
					if (id==ids[ii]) {
						Tops tops = _wells.get(id);
						if (tops==null)
							tops = new Tops(0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0);
							_wells.put(id,tops);
						for (int ic=0; ic<nf; ++ic) {
							if (fmtop.equals(fmcodes[ic])) {
								tops = addTop(tops,fmtop,tdepth*0.0003048); // convert from ft to km
							}
						}
					}
				}
			}
		}	catch (IOException e) {
        throw new RuntimeException(e);
		}
		// Get the well datums from file
		try {
      FileInputStream fis = new FileInputStream(datumName);
      Scanner s = new Scanner(fis); 
      while (s.hasNextLine()) {
				String line = s.nextLine();
				String[] fields = line.split(",");	
        if (fields.length>4 || fields.length<2)
          continue;
				if (fields[0].isEmpty() || fields[1].isEmpty())
					continue;
				long id = Long.parseLong(fields[0],10);
				double datum = Double.parseDouble(fields[1]);
				for (int ii=0; ii<ns; ++ii) {
					if (id==ids[ii]) {
						Tops tops = _wells.get(id); 
						tops.datum = datum*0.0003048; // convert from ft to km
					}
				}
			}
		}	catch (IOException e) {
        throw new RuntimeException(e);
		}
	}

	public void writeToFile(String fpath) {
		//String fpath = "/Users/amunoz/Home/data/tp/csm/welllogs/welltops/";
		long[] idw = _ids;
		String[] fms = _fms;
		int nw = idw.length;
		int nf = fms.length;
		for (int iw=0; iw<nw; ++iw) {
			long id = idw[iw];
			Tops tops = _wells.get(id);
			String outDir = fpath+String.valueOf(id);
			File fl = new File(outDir);
			fl.mkdir();
			for (int ic=0; ic<nf; ++ic) {
				String name = fms[ic];
				double fmd = getTop(tops,name);
				if (fmd>0.0) {
					try {
						String outFile = outDir+"/"+name+".dat";
						ArrayOutputStream aos = new ArrayOutputStream(outFile);
						aos.writeDouble(fmd);
						//aos.close();
				}	catch (IOException e) {
        		throw new RuntimeException(e);
					}
				}
			}
			double datum = tops.datum; 
			try {
				ArrayOutputStream aos = new ArrayOutputStream(outDir+"/datum.dat");
				aos.writeDouble(datum);
				//aos.close();
		}	catch (IOException e) {
        		throw new RuntimeException(e);
			}
		}
	}
	
	public double getTop(long id, String name) {
		Tops tops = _wells.get(id);
		if (name.equals("F2WC")) 
			return tops.f2wc;
		else if (name.equals("DKOT"))
			return tops.dkot;
		else if (name.equals("LKOT"))
			return tops.lkot;
		else if (name.equals("CRMT"))
			return tops.crmt;
		else if (name.equals("RDPK"))
			return tops.rdpk;
		else if (name.equals("A Sand"))
			return tops.asnd;
		else if (name.equals("C1 Dolo"))
			return tops.cdol;
		else if (name.equals("PC"))
			return tops.pcbg;
		else if (name.equals("datum"))
			return tops.datum;
		else System.out.println("Not a valid formation code");
		return 0.0f;
	}
	
	public void addTops(long id, double f2wc, double dkot, double lkot, 
			double crmt, double rdpk, double asnd, double cdol, double pcbg)
	{
		Tops tops = _wells.get(id); 
		if (tops!=null) {
			tops.f2wc = f2wc;
      tops.dkot = dkot;
      tops.lkot = lkot;
      tops.crmt = crmt;
      tops.rdpk = rdpk;
      tops.asnd = asnd;
      tops.cdol = cdol;
      tops.pcbg = pcbg;
		} else {
			tops = new Tops(f2wc,dkot,lkot,crmt,rdpk,asnd,cdol,pcbg);
			_wells.put(id,tops);
		}
	}

  private Map<Long,Tops> _wells = new HashMap<Long,Tops>();
	private long[] _ids;
	private String[] _fms;

	private class Tops {
		public Tops(
			double f2wc, double dkot, double lkot, double crmt, 
			double rdpk, double asnd, double cdol, double pcbg) 
		{
			this.f2wc = f2wc;
			this.dkot = dkot;
			this.lkot = lkot;
			this.crmt = crmt;
			this.rdpk = rdpk;
			this.asnd = asnd;
			this.cdol = cdol;
			this.pcbg = pcbg;
		}
		double f2wc; // Frontier 2nd Wall Creek (KF2) "F2WC"
		double dkot; // Dakota Sandstone (FallRiver) "DKOT"
		double lkot; // Lakota Sandstone (Lakota/Morrison) "LKOT"
		double crmt; // Crow Mountian "CRMT"
		double rdpk; // Red Peak "RDPK"
		double asnd; // Tensleep A Sandstone (Tensleep) "A Sand"
		double cdol; // Tensleep C1 Dolomite (TensleepBbase) "C1 Dolo"
		double pcbg; // Pre-Cambrian Basement Granite (Basement) "PC"
		double datum; // KB datum
	}

	private Tops addTop(Tops tops, String name, double v) {
		if (name.equals("F2WC")) {
			tops.f2wc = v;
			return tops;
		}
		else if (name.equals("DKOT"))  {
			tops.dkot = v;
			return tops;
		}
		else if (name.equals("LKOT")) {
			tops.lkot = v;
			return tops;
		}
		else if (name.equals("CRMT")) {
			tops.crmt = v;
			return tops;
		}
		else if (name.equals("RDPK")) {
			tops.rdpk = v;
			return tops;
		}
		else if (name.equals("A Sand")) {
			tops.asnd = v;
			return tops;
		}
		else if (name.equals("C1 Dolo")) {
			tops.cdol = v;
			return tops;
		}
		else if (name.equals("PC")) {
			tops.pcbg = v;
			return tops;
		}
		return tops;
	}

	private double getTop(Tops tops, String name) {
		if (name.equals("F2WC")) 
			return tops.f2wc;
		else if (name.equals("DKOT"))
			return tops.dkot;
		else if (name.equals("LKOT"))
			return tops.lkot;
		else if (name.equals("CRMT"))
			return tops.crmt;
		else if (name.equals("RDPK"))
			return tops.rdpk;
		else if (name.equals("A Sand"))
			return tops.asnd;
		else if (name.equals("C1 Dolo"))
			return tops.cdol;
		else if (name.equals("PC"))
			return tops.pcbg;
		else return 0.0f;
	}

	private void dump(double x) {
		System.out.println(x);
	}
	private void dump(String x) {
		System.out.println(x);
	}

//////////////////////////////////////////////////////////////////
// Tests and write out


  public static void main(String[] args) {
    writeDeepWellTops();
	}

	public static void writeDeepWellTops() {
		// Deep well set with velocity and density logs
		String fpath = "/Users/amunoz/data/tp/doe/WellLogs/";
		long[] dws = {
				 490251095000l, 490251091800l, 490251105400l, 490252305400l, 490251061000l, 
	 			 490251090200l, 490251087700l, 490251104600l, 490251113400l, 490251116100l, 
	 			 490251094400l, 490251096100l, 490251094600l, 490251104700l, 490251099200l, 
	 			 490251091600l, 490251097300l, 490251106400l};
		// Well tops that have seismic time horizons
		String[] fms = {"F2WC","DKOT","LKOT","CRMT","RDPK","A Sand","C1 Dolo","PC"};
		WellTops wtops = new WellTops(dws,fms,fpath);
		wtops.writeToFile("/Users/amunoz/data/tp/csm/welllogs/welltops/");
		System.out.println("Tops and datums written to file");
	}



};
