#==============================================================================
# Because v05 of DumPy is not normally in the path on this system (take this
# out to test on other systems)
import sys
sys.path.append("/Users/bkrueger/research/CCSNe/heating_layer/results/sample/dumses/visu")
#==============================================================================

import numpy as np

from dumpy_v05.data.rd_dumses import DumsesData
import SimulationData as SD

if __name__ == "__main__":

   # Create a SimulationInput
   si = SD.SimulationInput("/Users/bkrueger/research/CCSNe/heating_layer/results/sample/data/input")
   print si.__dict__

   # Some quick testing
   a = si.layer_width * (np.arange(31, dtype=float) - 15) / 10
   x = a.reshape((len(a),1,1))
   y = np.zeros((1,1,1))
   z = np.zeros((1,1,1))
   d = np.ones_like(x*y*z)
   g = si.gravity(x, y, z)
   h = si.heating(x, y, z, d)
   for i in xrange(len(a)):
      print "{x:10.3e}  {h:10.3e}  {gx:10.3e}  {gy:10.3e}  {gz:10.3e}".format(
            x=x[i,0,0], h=h[i,0,0],
            gx=g[i,0,0,0], gy=g[i,0,0,1], gz=g[i,0,0,2])

   # Get a DumsesData
   dd = DumsesData(0, filedir="/Users/bkrueger/research/CCSNe/heating_layer/results/sample/data/")

   # Create a SimulationState
   ss = SD.SimulationState(dd, si)

   density = ss.extract("density", "base state")
   pressure = ss.extract("pressure", "base state")
   velocity_x = ss.extract("x velocity", "base state")
   heating = si.heating(ss.x, ss.y, ss.z, density)
   H_Pv = heating / (pressure * velocity_x)
   H_Pv = H_Pv[1:-1,0,0]

   s = ss.extract("specific entropy", "base state")
   dsdx = (s[2:,0,0] - s[:-2,0,0]) / (ss.x[2:,0,0] - ss.x[:-2,0,0])

   mean_error = np.average(np.abs(H_Pv - dsdx))

   print ""
   print "mean error of entropy derivative = {0:10.3e}".format(mean_error)

   print "Let's run some basic tests:"
   print " - loop through all variables to ensure they compute"
   print " - loop through all modes to ensure they compute"
   for variable in SD.SimulationState.known_variables:
      for mode in SD.SimulationState.known_modes:
         print "   - trying " + variable + " in " + mode + " mode"
         var = ss.extract(variable, mode)
   print "success!"

