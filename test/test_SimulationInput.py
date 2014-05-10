import unittest
import numpy as np
import sys
sys.path.append("..")
import site_setup
from SimulationData import SimulationInput, SimulationError

class SimulationInputTest(unittest.TestCase):

   def test_uninitialized(self):
      si = SimulationInput()

      a = (np.arange(31, dtype=float) - 15) * 0.1
      x = a[:,None,None]
      y = a[None,:,None]
      z = a[None,None,:]
      rho = x * y * z

      self.assertRaises(SimulationError, si.heating, x, y, z, rho)

      self.assertRaises(SimulationError, si.gravity, x, y, z)

   def test_run_source_functions(self):
      for input_file_name in ["input_allvalues", "input_novalues"]:
         si = SimulationInput(input_file_name)

         a = (np.arange(31, dtype=float) - 15) * 0.1
         x = a[:,None,None]
         y = a[None,:,None]
         z = a[None,None,:]
         rho = x * y * z

         layer = si._shape_function(x, y, z)
         gravity = si.gravity(x, y, z)
         heating = si.heating(x, y, z, rho)

if __name__ == "__main__":
   unittest.main()
