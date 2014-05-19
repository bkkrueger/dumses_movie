"""
This module performs a series of tests on the MovieLimits object.
It will not be exhaustive, but it should cover basic usage and catch
the most glaring errors.
"""

import unittest
import numpy as np
import sys
sys.path.append("..")
import site_setup
from Descriptors import MovieLimits

class MovieLimitsTest(unittest.TestCase):
   """
   Test the SimulationInput object.
   """

   #===========================================================================

   def test_limits(self):
      """
      Make sure the basic functionality of the MovieLimits works.
      """

      # Can we make one?
      ml = MovieLimits()

      # Are the default values (and string conversion) as expected?
      self.assertEqual(str(ml), "[*,*]", msg="Incorrect string.")

      # Test applying the limits
      x = np.arange(7, dtype=float) * 0.5

      N = len(x[ml.test(x)])
      self.assertEqual(N, 7, msg="No limits failed.")

      ml.lo = 1.0
      N = len(x[ml.test(x)])
      self.assertEqual(N, 5, msg="Low limit failed.")

      ml.hi = 2.0
      N = len(x[ml.test(x)])
      self.assertEqual(N, 3, msg="Double limit failed.")

      ml.lo = None
      N = len(x[ml.test(x)])
      self.assertEqual(N, 5, msg="High limit failed.")

#==============================================================================

if __name__ == "__main__":
   unittest.main()
