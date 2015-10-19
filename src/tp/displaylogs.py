"""
Displays log curves of specified well
"""

from tputils import *

setupForSubset("subt_251_4_500")
#setupForSubset("subz_401_4_600")


def main(args):
  displaySubset()

def displaySubset():
  world = World()
  addLogsToWorld(world,"s","v",2.5,2.6,None,3)
  make2DFrame(world)
  

run(main)
