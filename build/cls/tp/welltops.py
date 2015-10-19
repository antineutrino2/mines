from imports import *

"""
Runs WellTops.java for testing and writing out tops
"""
def main(args):
		# Deep well set with velocity and density logs
		dws = [490251095000, 490251091800, 490251105400, 490252305400, 490251061000, 
	 			 	 490251090200, 490251087700, 490251104600, 490251113400, 490251116100, 
	 			 	 490251094400, 490251096100, 490251094600, 490251104700, 490251099200, 
	 			 	 490251091600, 490251097300, 490251106400]
		# Well tops that have seismic time horizons
		fms = ["F2WC","DKOT","LKOT","CRMT","RDPK","A Sand","C1 Dolo","PC"]
		wtops = WellTops(dws,fms)
		wtops.writeToFile()
		print "Tops and datums written to file"

#############################################################################
# Run the function main on the Swing thread
import sys
class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())

