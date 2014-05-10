import unittest
import numpy as np
import sys
sys.path.append("..")
import site_setup
from SimulationData import SimulationInput, SimulationError, SimulationState

class SimulationStateTest(unittest.TestCase):

   def test_construction(self):
      si = SimulationInput("input_allvalues")

      # Note : SimulationState doesn't actually require a DumsesData, but only
      #        something that meets the following criteria:
      #        -- Attributes: time, x, y, z, rho, rhou, E, rho0, rhou0, E0
      #        -- time should be a scalar
      #        -- coordinate arrays x, y, z should be 1D arrays
      #        -- rho, E should be 3D arrays that match the lengths of the
      #           coordinate arrays
      #        -- rhou should be a 4D array with the first three dimensions
      #           matching the coordinate arrays and the last having length 3
      #        -- rho0, rhou0, E0 should be 1D arrays that match the length of
      #           the x-coordinate array
      class A(object):
         pass
      data = A()
      Nx = 3
      Ny = 4
      Nz = 5
      data.time = 15.3
      data.x = np.arange(Nx, dtype=float)
      data.y = np.arange(Ny, dtype=float)
      data.z = np.arange(Nz, dtype=float)
      data.rho = np.zeros((Nx, Ny, Nz))
      data.rhou = np.zeros((Nx, Ny, Nz, 3))
      data.E = np.zeros((Nx, Ny, Nz))
      data.rho0 = np.zeros(Nx)
      data.rhou0 = np.zeros(Nx)
      data.E0 = np.zeros(Nx)

      # Make sure the construction works
      ss = SimulationState(data, si)

   def test_uninitialized(self):
      si = SimulationInput()
      ss = SimulationState()

      self.assertRaises(SimulationError, ss.extract, "variable")

   def test_variables(self):
      si = SimulationInput("input_allvalues")

      # Note : SimulationState doesn't actually require a DumsesData, but only
      #        something that meets the following criteria:
      #        -- Attributes: time, x, y, z, rho, rhou, E, rho0, rhou0, E0
      #        -- time should be a scalar
      #        -- coordinate arrays x, y, z should be 1D arrays
      #        -- rho, E should be 3D arrays that match the lengths of the
      #           coordinate arrays
      #        -- rhou should be a 4D array with the first three dimensions
      #           matching the coordinate arrays and the last having length 3
      #        -- rho0, rhou0, E0 should be 1D arrays that match the length of
      #           the x-coordinate array
      class A(object):
         pass
      data = A()
      Nx = 3
      Ny = 4
      Nz = 5
      data.time = 15.3
      data.x = np.arange(Nx, dtype=float)
      data.y = np.arange(Ny, dtype=float)
      data.z = np.arange(Nz, dtype=float)
      data.rho = np.ones((Nx, Ny, Nz))
      data.rhou = np.ones((Nx, Ny, Nz, 3))
      data.E = 3*np.ones((Nx, Ny, Nz))
      data.rho0 = np.ones(Nx)
      data.rhou0 = np.ones(Nx)
      data.E0 = np.ones(Nx)

      # Make sure the construction works
      ss = SimulationState(data, si)

      for v in ss.known_variables:
         self.assertSequenceEqual(ss.extract(v).shape, (Nx, Ny, Nz))

   def test_calculation(self):
      # TODO : Create a small data set, hand-verify all the computations for
      #        that small data set, then run the calculations from
      #        SimulationState and verify the results.
      si = SimulationInput("input_allvalues")

      class A(object):
         pass
      data = A()
      Nx = 5
      Ny = 3
      Nz = 1
      data.time = 15.3
      data.x = np.arange(Nx, dtype=float)
      data.y = np.arange(Ny, dtype=float)
      data.z = np.arange(Nz, dtype=float)
      data.rho = np.ones((Nx, Ny, Nz))
      data.rhou = np.ones((Nx, Ny, Nz, 3))
      data.E = 3*np.ones((Nx, Ny, Nz))
      data.rho0 = np.ones(Nx)
      data.rhou0 = np.ones(Nx)
      data.E0 = np.ones(Nx)

      # Make sure the construction works
      ss = SimulationState(data, si)

      for v in ss.known_variables:
         self.assertSequenceEqual(ss.extract(v).shape, (Nx, Ny, Nz))

if __name__ == "__main__":
   unittest.main()
