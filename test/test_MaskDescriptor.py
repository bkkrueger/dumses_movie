"""
This module performs a series of tests on the MaskDescriptor object.
It will not be exhaustive, but it should cover basic usage and catch
the most glaring errors.
"""

import unittest
import copy
import numpy as np
import sys
sys.path.append("..")
import site_setup
from Descriptors import MaskDescriptor, DescriptorError
from SimulationData import SimulationInput, SimulationState

class MaskDescriptorTest(unittest.TestCase):
   """
   Test the MaskDescriptor object.
   """

   #===========================================================================

   def test_construct(self):
      """
      Make sure we can build a mask.
      """

      # Valid values
      good = {"variable"  : "density",
              "mode"      : "pseudocolor: perturbation",
              "threshold" : "19.3",
              "operator"  : "<="}
      # Invalid values
      bad  = {"mode"      : "profile: mean state",
              "threshold" : "nineteen point three",
              "operator"  : "less than or equal to"}

      # Make sure the good values work
      md = MaskDescriptor(good)

      # Make sure the bad values don't work
      for k in bad.keys():
         pdict = copy.copy(good)
         pdict[k] = bad[k]
         self.assertRaises(DescriptorError, MaskDescriptor, pdict)

   #===========================================================================

   def test_apply(self):
      """
      Make sure we can apply a mask.
      """

      si = SimulationInput("input_allvalues")

      # Build a simple DumsesData-like object
      class A(object):
         pass
      data = A()
      Nx = 5
      Ny = 1
      Nz = 1
      data.time = 15.3
      data.x = (np.arange(Nx, dtype=float) - Nx//2)
      data.y = (np.arange(Ny, dtype=float) - Ny//2)
      data.z = (np.arange(Nz, dtype=float) - Nz//2)
      dens = np.arange(Nx*Ny*Nz, dtype=float).reshape((Nx, Ny, Nz)) + 10.0
      eint = np.arange(Nx*Ny*Nz, dtype=float).reshape((Nx, Ny, Nz)) + 10.0
      velx = np.arange(Nx*Ny*Nz, dtype=float).reshape((Nx, Ny, Nz))
      vely = np.zeros((Nx, Ny, Nz))
      velz = np.zeros((Nx, Ny, Nz))
      data.rho = dens
      data.rhou = np.empty((Nx, Ny, Nz, 3))
      data.rhou[...,0] = dens * velx
      data.rhou[...,1] = dens * vely
      data.rhou[...,2] = dens * velz
      data.E = dens * (eint + 0.5 * (velx**2 + vely**2 + velz**2))
      dens0 = dens[:,0,0]
      velx0 = velx[:,0,0]
      eint0 = eint[:,0,0]
      zero0 = np.zeros_like(velx0)
      data.rho0 = dens0
      data.rhou0 = dens0 * velx0
      data.E0 = dens0 * (eint0 + 0.5 * velx0**2)
      dens0 = dens0.reshape((Nx,1,1)).repeat(Ny, axis=1)
      velx0 = velx0.reshape((Nx,1,1)).repeat(Ny, axis=1)
      eint0 = eint0.reshape((Nx,1,1)).repeat(Ny, axis=1)
      zero0 = np.zeros((Nx, Ny, Nz))

      # Build the simulation state
      ss = SimulationState(data, si)

      md = MaskDescriptor({"variable" : "x velocity", "threshold" : 2,
         "operator" : "<="})

      # Apply a mask
      mask = md.apply(ss)
      allow = np.logical_not(mask)
      self.assertEquals(len(ss.extract("entropy")[np.logical_not(mask)]), 2)

#==============================================================================

if __name__ == "__main__":
   unittest.main()
