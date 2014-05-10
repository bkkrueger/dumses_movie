import unittest
import numpy as np
import sys
sys.path.append("..")
import site_setup
from SimulationData import SimulationInput, SimulationError

class SimulationInputTest(unittest.TestCase):

   def test_relations(self):
      """
      Test the SimulationInput class
      """
      for input_file_name in ["input_allvalues", "input_novalues"]:
         si = SimulationInput(input_file_name)

         pres_up = si.dens_up * si._csnd_up**2 / si.gamma
         self.assertAlmostEqual(si.pres_up, pres_up,
               msg="Mismatched pressure scale.")

         grav_coef = - si._Kgrav * si._csnd_up**2 / si.layer_width
         self.assertAlmostEqual(si._grav_coef, grav_coef,
               msg="Mismatched gravity coefficient.")

         heat_coef = si._Kheat * si._mach_up * si._csnd_up**3
         heat_coef /= (si.layer_width * si.gamma)
         self.assertAlmostEqual(si._heat_coef, heat_coef,
               msg="Mismatched heating coefficient.")

   def test_uninitialized(self):
      si = SimulationInput()

      a = (np.arange(31, dtype=float) - 15) * 0.1
      x = a[:,None,None]
      y = a[None,:,None]
      z = a[None,None,:]
      rho = x * y * z

      self.assertRaises(SimulationError, si.heating, x, y, z, rho)

      self.assertRaises(SimulationError, si.gravity, x, y, z)

   def test_source_functions(self):
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

         for i in xrange(layer.shape[0]):
            for j in xrange(layer.shape[1]):
               for k in xrange(layer.shape[2]):
                  self.assertAlmostEqual(gravity[i,j,k,0],
                        si._grav_coef * layer[i,j,k],
                        msg="Gravity scaling is off.")
                  self.assertAlmostEqual(gravity[i,j,k,1], 0.0,
                        msg="Gravity scaling is off.")
                  self.assertAlmostEqual(heating[i,j,k],
                        si._heat_coef * rho[i,j,k] * layer[i,j,k],
                        msg="Heating scaling is off.")

if __name__ == "__main__":
   unittest.main()
