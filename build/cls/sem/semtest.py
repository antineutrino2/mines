from bputils import *
from imports import *

pngDir = "../../../../../pics/sem/"
#pngDir = "/Users/amunoz/Home/pics/sem/"


def main(args):
  #goHypTest()
  #goNonHypTest()
  data = getData(2790)
  goHyp(data)
  goNonHyp(data)































#------------------------------------------------------------------------------#

class RunMain(Runnable):
  def run(self):
    main(sys.argv)
SwingUtilities.invokeLater(RunMain())
