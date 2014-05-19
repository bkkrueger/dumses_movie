"""
This module performs a series of tests on the SimulationState object.
It will not be exhaustive, but it should cover basic usage and catch
the most glaring errors.
"""

from __future__ import division
import unittest
import numpy as np
import sys
sys.path.append("..")
import site_setup
from SimulationData import SimulationInput, SimulationError, SimulationState

class SimulationStateTest(unittest.TestCase):
   """
   Test the SimulationState object.
   """

   #===========================================================================

   def test_construction(self):
      """
      Test the construction of a SimulationState.
      """

      si = SimulationInput("input_allvalues")

      # Make a fake "DumsesData"-like object
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

   #===========================================================================

   def test_uninitialized(self):
      """
      Test thet an uninitialized SimulationData raises an error when
      asked to construct a variable. 
      """
      si = SimulationInput()
      ss = SimulationState()

      self.assertRaises(SimulationError, ss.extract, "density")

   #===========================================================================

   def test_variables(self):
      """
      Test that all the variables extract.
      """
      si = SimulationInput("input_allvalues")

      # Make a fake "DumsesData"-like object
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

      ss = SimulationState(data, si)

      # Make sure the extraction works, and also check that the results are the
      # right size
      for v in ss.known_variables:
         self.assertSequenceEqual(ss.extract(v).shape, (Nx, Ny, Nz))

   #===========================================================================

   def test_calculation(self):
      """
      Ensure the correct relations among the variables.
      """

      # Build a SimulationState
      si = SimulationInput("input_allvalues")
      class A(object):
         pass
      data = A()
      Nx = 7
      Ny = 5
      Nz = 3
      data.time = 15.3
      x = (np.arange(Nx, dtype=float) - Nx//2) * 0.75
      y = (np.arange(Ny, dtype=float) - Ny//2) * 0.75
      z = (np.arange(Nz, dtype=float) - Nz//2) * 0.75
      data.x = x
      data.y = y
      data.z = z
      dens = np.random.rand(Nx, Ny, Nz) + 0.1
      eint = np.random.rand(Nx, Ny, Nz) + 0.1
      velx = np.random.rand(Nx, Ny, Nz) - 0.5
      vely = np.random.rand(Nx, Ny, Nz) - 0.5
      velz = np.random.rand(Nx, Ny, Nz) - 0.5
      data.rho = dens
      data.rhou = np.empty((Nx, Ny, Nz, 3))
      data.rhou[...,0] = dens * velx
      data.rhou[...,1] = dens * vely
      data.rhou[...,2] = dens * velz
      data.E = dens * (eint + 0.5 * (velx**2 + vely**2 + velz**2))
      dens0 = np.random.rand(Nx) + 0.1
      eint0 = np.random.rand(Nx) + 0.1
      velx0 = np.random.rand(Nx) - 0.5
      zero0 = np.zeros_like(velx0)
      data.rho0 = dens0
      data.rhou0 = dens0 * velx0
      data.E0 = dens0 * (eint0 + 0.5 * velx0**2)
      dens0 = dens0.reshape((dens.shape[0],1,1)).repeat(Ny, axis=1).repeat(Nz,
            axis=2)
      velx0 = velx0.reshape((dens.shape[0],1,1)).repeat(Ny, axis=1).repeat(Nz,
            axis=2)
      eint0 = eint0.reshape((dens.shape[0],1,1)).repeat(Ny, axis=1).repeat(Nz,
            axis=2)
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

      # boundary conditions: x-bounds are zero-gradient, y & z are periodic
      dvy_dx = np.empty_like(vely)
      dvy_dx[ 0  ,:,:] = 0.0
      dvy_dx[1:-1,:,:] = ((vely[2:,:,:] - vely[:-2,:,:]) /
            (x[2:,None,None] - x[:-2,None,None]))
      dvy_dx[ -1 ,:,:] = 0.0

      dvz_dx = np.empty_like(velz)
      dvz_dx[ 0  ,:,:] = 0.0
      dvz_dx[1:-1,:,:] = ((velz[2:,:,:] - velz[:-2,:,:]) /
            (x[2:,None,None] - x[:-2,None,None]))
      dvz_dx[ -1 ,:,:] = 0.0

      dvx_dy = np.empty_like(velx)
      dvx_dy[:, 0  ,:] = (velx[:, 1,:] - velx[:,-1 ,:]) / ((y[ 1] - y[ 0 ])*2)
      dvx_dy[:,1:-1,:] = ((velx[:,2:,:] - velx[:,:-2,:]) /
            (y[None,2:,None] - y[None,:-2,None]))
      dvx_dy[:, -1 ,:] = (velx[:, 0,:] - velx[:,-2 ,:]) / ((y[-1] - y[-2 ])*2)

      dvz_dy = np.empty_like(velz)
      dvz_dy[:, 0  ,:] = (velz[:, 1,:] - velz[:,-1 ,:]) / ((y[ 1] - y[ 0 ])*2)
      dvz_dy[:,1:-1,:] = ((velz[:,2:,:] - velz[:,:-2,:]) /
            (y[None,2:,None] - y[None,:-2,None]))
      dvz_dy[:, -1 ,:] = (velz[:, 0,:] - velz[:,-2 ,:]) / ((y[-1] - y[-2 ])*2)

      dvx_dz = np.empty_like(velx)
      dvx_dz[:,:, 0  ] = (velx[:,:, 1] - velx[:,:,-1 ]) / ((z[ 1] - z[ 0 ])*2)
      dvx_dz[:,:,1:-1] = ((velx[:,:,2:] - velx[:,:,:-2]) /
            (z[None,None,2:] - z[None,None,:-2]))
      dvx_dz[:,:, -1 ] = (velx[:,:, 0] - velx[:,:,-2 ]) / ((z[-1] - z[-2 ])*2)

      dvy_dz = np.empty_like(vely)
      dvy_dz[:,:, 0  ] = (vely[:,:, 1] - vely[:,:,-1 ]) / ((z[ 1] - z[ 0 ])*2)
      dvy_dz[:,:,1:-1] = ((vely[:,:,2:] - vely[:,:,:-2]) /
            (z[None,None,2:] - z[None,None:-2]))
      dvy_dz[:,:, -1 ] = (vely[:,:, 0] - vely[:,:,-2 ]) / ((z[-1] - z[-2 ])*2)

      vrtx = dvz_dy - dvy_dz
      self.all_nearly_equal(vrtx, ss.extract("x vorticity"))
      self.assertEqual(ss.known_variables["x vorticity"].zero, "center")

      vrty = dvx_dz - dvz_dx
      self.all_nearly_equal(vrty, ss.extract("y vorticity"))
      self.assertEqual(ss.known_variables["y vorticity"].zero, "center")

      vrtz = dvy_dx - dvx_dy
      self.all_nearly_equal(vrtz, ss.extract("z vorticity"))
      self.assertEqual(ss.known_variables["z vorticity"].zero, "center")

      wmag = np.sqrt(vrtx**2 + vrty**2 + vrtz**2)
      self.all_nearly_equal(wmag, ss.extract("vorticity magnitude"))
      self.assertEqual(ss.known_variables["vorticity magnitude"].zero,
            "lower bound")

      svrx = vrtx / dens
      self.all_nearly_equal(svrx, ss.extract("specific x vorticity"))
      self.assertEqual(ss.known_variables["specific x vorticity"].zero,
            "center")

      svry = vrty / dens
      self.all_nearly_equal(svry, ss.extract("specific y vorticity"))
      self.assertEqual(ss.known_variables["specific y vorticity"].zero,
            "center")

      svrz = vrtz / dens
      self.all_nearly_equal(svrz, ss.extract("specific z vorticity"))
      self.assertEqual(ss.known_variables["specific z vorticity"].zero,
            "center")

      swmg = wmag / dens
      self.all_nearly_equal(swmg, ss.extract("specific vorticity magnitude"))
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

      self.all_nearly_equal(dens * eint + pres,
            ss.extract("enthalpy density"))
      self.all_nearly_equal(dens0 * eint0 + pres0,
            ss.extract("enthalpy density", False))
      self.assertEqual(ss.known_variables["enthalpy density"].zero,
            "lower bound")

      self.all_nearly_equal(eint + pres / dens,
            ss.extract("specific enthalpy"))
      self.all_nearly_equal(eint0 + pres0 / dens0,
            ss.extract("specific enthalpy", False))
      self.assertEqual(ss.known_variables["specific enthalpy"].zero,
            "lower bound")

      grad_s = np.empty((Nx,Ny,Nz,3))
      grad_s[ 0  ,:,:,0] = 0.0
      grad_s[1:-1,:,:,0] = ((entr[2:,:,:] - entr[:-2,:,:]) /
            (x[2:,None,None] - x[:-2,None,None]))
      grad_s[ -1 ,:,:,0] = 0.0
      grad_s[:, 0  ,:,1] = (entr[:, 1,:] - entr[:,-1 ,:]) / ((y[ 1] - y[ 0 ])*2)
      grad_s[:,1:-1,:,1] = ((entr[:,2:,:] - entr[:,:-2,:]) /
            (y[None,2:,None] - y[None,:-2,None]))
      grad_s[:, -1 ,:,1] = (entr[:, 0,:] - entr[:,-2 ,:]) / ((y[-1] - y[-2 ])*2)
      grad_s[:,:, 0  ,2] = (entr[:,:, 1] - entr[:,:,-1 ]) / ((z[ 1] - z[ 0 ])*2)
      grad_s[:,:,1:-1,2] = ((entr[:,:,2:] - entr[:,:,:-2]) /
            (z[None,None,2:] - z[None,None:-2]))
      grad_s[:,:, -1 ,2] = (entr[:,:, 0] - entr[:,:,-2 ]) / ((z[-1] - z[-2 ])*2)

      grav = ss.params.gravity(ss.x, ss.y, ss.z)
      g_grad_s = np.sum(grav * grad_s, axis=-1)
      wbv2 = g_grad_s * (ss.params.gamma - 1) / ss.params.gamma
      wbv2[wbv2 > 0] = 0
      w_bv = np.sqrt(-wbv2)
      self.all_nearly_equal(w_bv, ss.extract("Brunt-Vaisala frequency"))
      self.assertEqual(ss.known_variables["convective growth rate"].zero,
            "lower bound")

   #===========================================================================

   def all_nearly_equal(self, npa1, npa2):
      """
      Support function to make sure the data arrays are equal.
      """
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

#==============================================================================

if __name__ == "__main__":
   unittest.main()
