"""
This module performs a series of tests on the SimulationInput object.
"""

import unittest
import numpy as np
import sys
sys.path.append("..")
import site_setup
from SimulationData import SimulationInput, SimulationError

class SimulationInputTest(unittest.TestCase):
   """
   Test the SimulationInput object.
   """

   #===========================================================================

   def test_uninitialized(self):
      """
      An uninitialized SimulationInput should raise an error when the
      source functions are accessed.
      """

      # An uninitialized SimulationInput
      si = SimulationInput()

      # Create some dummy data; values aren't important.
      a = (np.arange(31, dtype=float) - 15) * 0.1
      x = a[:,None,None]
      y = a[None,:,None]
      z = a[None,None,:]
      rho = x * y * z

      # Heating shouldn't work because the scaling parameters aren't there
      self.assertRaises(SimulationError, si.heating, x, y, z, rho)

      # Gravity shouldn't work because the scaling parameters aren't there
      self.assertRaises(SimulationError, si.gravity, x, y, z)

   #===========================================================================

   def test_source_functions(self):
      """
      Initialize the SimulationInput (with and without defaults) and make sure
      the heating and gravity work.
      """
      for input_file_name in ["input_allvalues", "input_novalues"]:
         si = SimulationInput(input_file_name)

         a = (np.arange(31, dtype=float) - 15) * 0.1
         x = a[:,None,None]
         y = a[None,:,None]
         z = a[None,None,:]
         rho = x * y * z

         # Step 1: Do the functions actually work?
         gravity = si.gravity(x, y, z)
         heating = si.heating(x, y, z, rho)

         # Step 2: Do they give results with the right shape?
         self.assertEquals(len(gravity.shape), 4)
         self.assertEquals(len(heating.shape), 3)
         for i in xrange(3):
            self.assertEquals(gravity.shape[i], rho.shape[i])
            self.assertEquals(heating.shape[i], rho.shape[i])
         self.assertEquals(gravity.shape[3], 3)

         # Step 3: Do they scale correctly with the parameters?
         layer = si.shape_function(x, y, z)
         g = - si.K_grav * (si.csnd_up**2 / si.layer_width) * layer
         h = (si.K_heat * ((si.csnd_up**3 * si.mach_up) /
            (si.gamma * si.layer_width)) * rho * layer)
         for i in xrange(g.shape[0]):
            for j in xrange(g.shape[1]):
               for k in xrange(g.shape[2]):
                  self.assertAlmostEquals(h[i,j,k], heating[i,j,k])
                  self.assertAlmostEquals(gravity[i,j,k,0], g[i,j,k])
                  self.assertAlmostEquals(gravity[i,j,k,1], 0)
                  self.assertAlmostEquals(gravity[i,j,k,2], 0)

#==============================================================================

if __name__ == "__main__":
   unittest.main()
