from __future__ import division
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
      si = SimulationInput("input_allvalues")

      class A(object):
         pass
      data = A()
      Nx = 5
      Ny = 3
      Nz = 1
      data.time = 15.3
      data.x = (np.arange(Nx, dtype=float) - Nx//2) * 0.75
      data.y = (np.arange(Ny, dtype=float) - Ny//2) * 0.75
      data.z = (np.arange(Nz, dtype=float) - Nz//2) * 0.75
      dens = np.empty((Nx, Ny, Nz))
      dens[:,0,0] = [1.0, 1.1, 1.2, 1.3, 1.4]
      dens[:,1,0] = [1.0, 1.1, 1.2, 1.3, 1.4]
      dens[:,2,0] = [1.0, 1.1, 1.2, 1.3, 1.4]
      eint = np.empty((Nx, Ny, Nz))
      eint[:,0,0] = [2.4, 2.3, 2.2, 2.1, 2.0]
      eint[:,1,0] = [2.4, 2.3, 2.2, 2.1, 2.0]
      eint[:,2,0] = [2.4, 2.3, 2.2, 2.1, 2.0]
      velx = np.empty((Nx, Ny, Nz))
      velx[:,0,0] = [0.1, 0.2, 0.3, 0.2, 0.1]
      velx[:,1,0] = [0.3, 0.4, 0.5, 0.4, 0.3]
      velx[:,2,0] = [0.5, 0.6, 0.7, 0.8, 0.9]
      vely = np.empty((Nx, Ny, Nz))
      vely[:,0,0] = [0.1, 0.2, 0.3, 0.2, 0.1]
      vely[:,1,0] = [0.3, 0.4, 0.5, 0.4, 0.3]
      vely[:,2,0] = [0.5, 0.6, 0.7, 0.8, 0.9]
      velz = np.zeros((Nx, Ny, Nz))
      data.rho = dens
      data.rhou = np.empty((Nx, Ny, Nz, 3))
      data.rhou[...,0] = dens * velx
      data.rhou[...,1] = dens * vely
      data.rhou[...,2] = dens * velz
      data.E = dens * (eint + 0.5 * (velx**2 + vely**2 + velz**2))
      dens0 = dens[:,1,0]
      velx0 = velx[:,1,0]
      eint0 = eint[:,1,0]
      zero0 = np.zeros_like(velx0)
      data.rho0 = dens0
      data.rhou0 = dens0 * velx0
      data.E0 = dens0 * (eint0 + 0.5 * velx0**2)
      dens0 = dens0.reshape((dens.shape[0],1,1)).repeat(3, axis=1)
      velx0 = velx0.reshape((dens.shape[0],1,1)).repeat(3, axis=1)
      eint0 = eint0.reshape((dens.shape[0],1,1)).repeat(3, axis=1)
      zero0 = np.zeros((Nx, Ny, Nz))

      # Make sure the construction works
      ss = SimulationState(data, si)

      # Test the variables
      self.all_nearly_equal(dens, ss.extract("density"))
      self.all_nearly_equal(dens0, ss.extract("mass density", False))
      self.assertEqual(ss.known_variables["density"].zero, "lower bound")

      self.all_nearly_equal(dens*velx, ss.extract("x momentum"))
      self.all_nearly_equal(dens0*velx0,
            ss.extract("x momentum density", False))
      self.assertEqual(ss.known_variables["x momentum"].zero, "center")

      self.all_nearly_equal(dens*vely, ss.extract("y momentum"))
      self.all_nearly_equal(dens0*zero0,
            ss.extract("y momentum density", False))
      self.assertEqual(ss.known_variables["y momentum"].zero, "center")

      self.all_nearly_equal(dens*velz, ss.extract("z momentum"))
      self.all_nearly_equal(dens0*zero0,
            ss.extract("z momentum density", False))
      self.assertEqual(ss.known_variables["z momentum"].zero, "center")

      vmag = np.sqrt(velx**2 + vely**2 + velz**2)
      vmag0 = np.abs(velx0)
      self.all_nearly_equal(dens * vmag, ss.extract("momentum magnitude"))
      self.all_nearly_equal(dens0 * vmag0,
            ss.extract("momentum density magnitude", False))
      self.assertEqual(ss.known_variables["momentum magnitude"].zero,
            "lower bound")

      self.all_nearly_equal(velx, ss.extract("x velocity"))
      self.all_nearly_equal(velx0, ss.extract("specific x momentum", False))
      self.assertEqual(ss.known_variables["x velocity"].zero, "center")

      self.all_nearly_equal(vely, ss.extract("y velocity"))
      self.all_nearly_equal(zero0, ss.extract("specific y momentum", False))
      self.assertEqual(ss.known_variables["y velocity"].zero, "center")

      self.all_nearly_equal(velz, ss.extract("z velocity"))
      self.all_nearly_equal(zero0, ss.extract("specific z momentum", False))
      self.assertEqual(ss.known_variables["z velocity"].zero, "center")

      self.all_nearly_equal(vmag, ss.extract("velocity magnitude"))
      self.all_nearly_equal(vmag0,
            ss.extract("specific momentum magnitude", False))
      self.assertEqual(ss.known_variables["velocity magnitude"].zero,
            "lower bound")

      # TODO : vorticity (x, y, z, magnitude)
      self.assertEqual(ss.known_variables["x vorticity"].zero, "center")
      self.assertEqual(ss.known_variables["y vorticity"].zero, "center")
      self.assertEqual(ss.known_variables["z vorticity"].zero, "center")
      self.assertEqual(ss.known_variables["vorticity magnitude"].zero,
            "lower bound")

      # TODO : specific vorticity (x, y, z, magnitude)
      self.assertEqual(ss.known_variables["specific x vorticity"].zero,
            "center")
      self.assertEqual(ss.known_variables["specific y vorticity"].zero,
            "center")
      self.assertEqual(ss.known_variables["specific z vorticity"].zero,
            "center")
      self.assertEqual(ss.known_variables["specific vorticity magnitude"].zero,
            "lower bound")

      self.all_nearly_equal(0.5 * dens * vmag**2,
            ss.extract("kinetic energy density"))
      self.all_nearly_equal(0.5 * dens0 * vmag0**2,
            ss.extract("kinetic energy density", False))
      self.assertEqual(ss.known_variables["kinetic energy density"].zero,
            "lower bound")

      self.all_nearly_equal(0.5 * vmag**2,
            ss.extract("specific kinetic energy"))
      self.all_nearly_equal(0.5 * vmag0**2,
            ss.extract("specific kinetic energy", False))
      self.assertEqual(ss.known_variables["specific kinetic energy"].zero,
            "lower bound")

      self.all_nearly_equal(dens * (eint + 0.5 * vmag**2),
            ss.extract("total energy density"))
      self.all_nearly_equal(dens0 * (eint0 + 0.5 * vmag0**2),
            ss.extract("total energy density", False))
      self.assertEqual(ss.known_variables["total energy density"].zero,
            "lower bound")

      self.all_nearly_equal(eint + 0.5 * vmag**2,
            ss.extract("specific total energy"))
      self.all_nearly_equal(eint0 + 0.5 * vmag0**2,
            ss.extract("specific total energy", False))
      self.assertEqual(ss.known_variables["specific total energy"].zero,
            "lower bound")

      self.all_nearly_equal(dens * eint, ss.extract("internal energy density"))
      self.all_nearly_equal(dens0 * eint0,
            ss.extract("internal energy density", False))
      self.assertEqual(ss.known_variables["internal energy density"].zero,
            "lower bound")

      self.all_nearly_equal(eint, ss.extract("specific internal energy"))
      self.all_nearly_equal(eint0,
            ss.extract("specific internal energy", False))
      self.assertEqual(ss.known_variables["specific internal energy"].zero,
            "lower bound")

      pres = dens * eint * (ss.params.gamma - 1)
      pres0 = dens0 * eint0 * (ss.params.gamma - 1)
      self.all_nearly_equal(pres, ss.extract("pressure"))
      self.all_nearly_equal(pres0, ss.extract("pressure", False))
      self.assertEqual(ss.known_variables["pressure"].zero, "lower bound")

      entr = (np.log((pres/ss.params.pres_up) *
            (dens/ss.params.dens_up)**(-ss.params.gamma)) /
            (ss.params.gamma - 1))
      entr0 = (np.log((pres0/ss.params.pres_up) *
            (dens0/ss.params.dens_up)**(-ss.params.gamma)) /
            (ss.params.gamma - 1))
      self.all_nearly_equal(entr, ss.extract("entropy"))
      self.all_nearly_equal(entr0, ss.extract("specific entropy", False))
      self.assertIsNone(ss.known_variables["entropy"].zero)

      self.all_nearly_equal(dens * entr, ss.extract("entropy density"))
      self.all_nearly_equal(dens0 * entr0,
            ss.extract("entropy density", False))
      self.assertIsNone(ss.known_variables["entropy density"].zero)

      csnd = np.sqrt(ss.params.gamma * pres / dens)
      csnd0 = np.sqrt(ss.params.gamma * pres0 / dens0)
      self.all_nearly_equal(csnd, ss.extract("sound speed"))
      self.all_nearly_equal(csnd0, ss.extract("sound speed", False))
      self.assertEqual(ss.known_variables["sound speed"].zero, "lower bound")

      self.all_nearly_equal(vmag / csnd, ss.extract("Mach number"))
      self.all_nearly_equal(vmag0 / csnd0, ss.extract("Mach number", False))
      self.assertEqual(ss.known_variables["Mach number"].zero, "lower bound")

      self.all_nearly_equal(eint + pres / dens,
            ss.extract("specific enthalpy"))
      self.all_nearly_equal(eint0 + pres0 / dens0,
            ss.extract("specific enthalpy"))
      self.assertEqual(ss.known_variables["specific enthalpy"].zero,
            "lower bound")

      self.all_nearly_equal(dens * eint + pres,
            ss.extract("enthalpy density"))
      self.all_nearly_equal(dens0 * eint0 + pres0,
            ss.extract("enthalpy density", False))
      self.assertEqual(ss.known_variables["enthalpy density"].zero,
            "lower bound")

      # TODO : convective growth rate
      self.assertEqual(ss.known_variables["convective growth rate"].zero,
            "lower bound")

   def all_nearly_equal(self, npa1, npa2):
      self.assertEqual(len(npa1.shape), len(npa2.shape),
            msg="Mismatching shapes: ({0}, {1}).".format(len(npa1.shape),
            len(npa2.shape)))
      s1 = 1
      s2 = 1
      for i in xrange(len(npa1.shape)):
         s1 *= npa1.shape[i]
         s2 *= npa2.shape[i]
         self.assertEqual(npa1.shape[i], npa2.shape[i],
               msg="Mismatching dimension {0}: ({1}, {2}).".format(i,
               npa1.shape[i], npa2.shape[i]))
      tmp1 = npa1.reshape((s1))
      tmp2 = npa2.reshape((s2))
      for i in xrange(s1):
         self.assertAlmostEqual(tmp1[i], tmp2[i],
               msg="Mismatching value: ({0}, {1}).".format(tmp1[i], tmp2[i]))

if __name__ == "__main__":
   unittest.main()
