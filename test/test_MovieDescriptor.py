"""
This module performs a series of tests on the MovieDescriptor object.
It will not be exhaustive, but it should cover basic usage and catch
the most glaring errors.
"""

import unittest
import copy
import numpy as np
import sys
sys.path.append("..")
import site_setup
from Descriptors import MovieDescriptor, MaskDescriptor, DescriptorError
from SimulationData import SimulationInput, SimulationState

class MovieDescriptorTest(unittest.TestCase):
   """
   Test the MovieDescriptor object.
   """

   def test_construct(self):
      """
      Can we construct a MovieDescriptor?
      """

      # This should construct a valid movie
      params = {"title" : "title",
                "stub" : "stub",
                "path" : "path",
                "image_type" : "png",
                "variable" : "density",
                "mode" : "profile : contrast",
                "window_x_lo" : -1.2,
                "window_x_hi" :  1.2,
                "window_y_lo" : -1,
                "window_y_hi" :  1,
                "window_z_lo" : -1,
                "window_z_hi" :  1,
                "time_lo" : 13,
                "time_hi" : 73,
                "value_lo" : 1.5,
                "value_hi" : 12.5,
                "make_movie" : True,
                "fps" : 60,
                "movie_type" : "avi",
                "final_pause" : "16.5s",
                "masks" : [],
                "mask_method" : "force low",
                "xlines" : [-1, -0.5, 0, 0.5, 1]
                }
      movie = MovieDescriptor(params, [])

      # Use defaults where allowed
      p = copy.deepcopy(params)
      del p["title"]
      del p["path"]
      del p["mode"]
      del p["window_x_lo"]
      del p["window_x_hi"]
      del p["window_y_lo"]
      del p["window_y_hi"]
      del p["window_z_lo"]
      del p["window_z_hi"]
      del p["time_lo"]
      del p["time_hi"]
      del p["value_lo"]
      del p["value_hi"]
      del p["fps"]
      del p["movie_type"]
      del p["final_pause"]
      del p["masks"]
      del p["mask_method"]
      movie = MovieDescriptor(params, [])

      # Require a stub
      p = copy.deepcopy(params)
      del p["stub"]
      self.assertRaises(KeyError, MovieDescriptor, p, [])

      # Require image type
      p = copy.deepcopy(params)
      del p["image_type"]
      self.assertRaises(KeyError, MovieDescriptor, p, [])

      # Require variable
      p = copy.deepcopy(params)
      del p["variable"]
      self.assertRaises(KeyError, MovieDescriptor, p, [])

      # Require window_x_lo < window_x_hi
      p = copy.deepcopy(params)
      p["window_x_lo"] = p["window_x_hi"]
      self.assertRaises(DescriptorError, MovieDescriptor, p, [])

      # Require window_y_lo < window_y_hi
      p = copy.deepcopy(params)
      p["window_y_lo"] = p["window_y_hi"]
      self.assertRaises(DescriptorError, MovieDescriptor, p, [])

      # Require window_z_lo < window_z_hi
      p = copy.deepcopy(params)
      p["window_z_lo"] = p["window_z_hi"]
      self.assertRaises(DescriptorError, MovieDescriptor, p, [])

      # Require time_lo < time_hi
      p = copy.deepcopy(params)
      p["time_lo"] = p["time_hi"]
      self.assertRaises(DescriptorError, MovieDescriptor, p, [])

      # Require value_lo < value_hi
      p = copy.deepcopy(params)
      p["value_lo"] = p["value_hi"]
      self.assertRaises(DescriptorError, MovieDescriptor, p, [])

      # Allow make_movie to take default
      p = copy.deepcopy(params)
      del p["make_movie"]
      movie = MovieDescriptor(p, [])

      # Require fps only if make_movie is true
      p = copy.deepcopy(params)
      del p["fps"]
      del p["final_pause"]
      p["make_movie"] = False
      movie = MovieDescriptor(p, [])
      p["make_movie"] = True
      self.assertRaises(KeyError, MovieDescriptor, p, [])

      # Require movie_type only if make_movie is true
      p = copy.deepcopy(params)
      del p["movie_type"]
      p["make_movie"] = False
      movie = MovieDescriptor(p, [])
      p["make_movie"] = True
      self.assertRaises(KeyError, MovieDescriptor, p, [])

      # Require fps if final pause specified in seconds
      p = copy.deepcopy(params)
      del p["fps"]
      p["make_movie"] = False
      p["final_pause"] = "16f"
      movie = MovieDescriptor(p, [])   # prove it works up to here
      p["final_pause"] = "13.6s"
      self.assertRaises(DescriptorError, MovieDescriptor, p, [])

      # Try bad value for mode
      p = copy.deepcopy(params)
      p["mode"] = "bad mode"
      self.assertRaises(DescriptorError, MovieDescriptor, p, [])

      # Try bad values for window limits
      p = copy.deepcopy(params)
      p["window_x_lo"] = "yes"
      self.assertRaises(ValueError, MovieDescriptor, p, [])
      p = copy.deepcopy(params)
      p["window_x_hi"] = "yes"
      self.assertRaises(ValueError, MovieDescriptor, p, [])
      p = copy.deepcopy(params)
      p["window_y_lo"] = "yes"
      self.assertRaises(ValueError, MovieDescriptor, p, [])
      p = copy.deepcopy(params)
      p["window_y_hi"] = "yes"
      self.assertRaises(ValueError, MovieDescriptor, p, [])
      p = copy.deepcopy(params)
      p["window_z_lo"] = "yes"
      self.assertRaises(ValueError, MovieDescriptor, p, [])
      p = copy.deepcopy(params)
      p["window_z_hi"] = "yes"
      self.assertRaises(ValueError, MovieDescriptor, p, [])

      # Try bad values for time limits
      p = copy.deepcopy(params)
      p["time_lo"] = "yes"
      self.assertRaises(ValueError, MovieDescriptor, p, [])
      p = copy.deepcopy(params)
      p["time_hi"] = "yes"
      self.assertRaises(ValueError, MovieDescriptor, p, [])

      # Try bad values for value limits
      p = copy.deepcopy(params)
      p["value_lo"] = "yes"
      self.assertRaises(ValueError, MovieDescriptor, p, [])
      p = copy.deepcopy(params)
      p["value_hi"] = "yes"
      self.assertRaises(ValueError, MovieDescriptor, p, [])

      # Try bad value for fps
      p = copy.deepcopy(params)
      p["fps"] = "what?"
      self.assertRaises(ValueError, MovieDescriptor, p, [])

      # Try bad values for final pause
      p = copy.deepcopy(params)
      p["final_pause"] = "what?"
      self.assertRaises(DescriptorError, MovieDescriptor, p, [])
      p["final_pause"] = "13"
      self.assertRaises(DescriptorError, MovieDescriptor, p, [])
      p["final_pause"] = "13.5f"
      self.assertRaises(DescriptorError, MovieDescriptor, p, [])

      # Try bad value for mask method
      p = copy.deepcopy(params)
      p["mask_method"] = "what?"
      self.assertRaises(DescriptorError, MovieDescriptor, p, [])

      # Try bad value for xlines
      p = copy.deepcopy(params)
      p["xlines"] = [0, 0.0, "zero"]
      self.assertRaises(ValueError, MovieDescriptor, p, [])

   def test_build_paths(self):
      # This should construct a valid movie
      params = {"title" : "title",
                "stub" : "stub",
                "path" : "path",
                "image_type" : "png",
                "variable" : "density",
                "mode" : "profile : contrast",
                "window_x_lo" : -1.2,
                "window_x_hi" :  1.2,
                "window_y_lo" : -1,
                "window_y_hi" :  1,
                "window_z_lo" : -1,
                "window_z_hi" :  1,
                "time_lo" : 13,
                "time_hi" : 73,
                "value_lo" : 1.5,
                "value_hi" : 12.5,
                "make_movie" : True,
                "fps" : 60,
                "movie_type" : "avi",
                "final_pause" : "16.5s",
                "masks" : [],
                "mask_method" : "force low",
                "xlines" : [-1, -0.5, 0, 0.5, 1]
                }
      data_dir = "/data_directory"

      # Default path
      del params["path"]
      movie = MovieDescriptor(params, [])
      mp, fp = movie.build_paths(data_dir)
      self.assertEqual(mp, data_dir + "/")
      self.assertEqual(fp, data_dir + "/frames/")

      # Absolute path
      abs_path = "/absolute_path"
      params["path"] = abs_path
      movie = MovieDescriptor(params, [])
      mp, fp = movie.build_paths(data_dir)
      self.assertEqual(mp, abs_path + "/")
      self.assertEqual(fp, abs_path + "/frames/")

      # Relative path
      rel_path = "relative_path"
      params["path"] = rel_path
      movie = MovieDescriptor(params, [])
      mp, fp = movie.build_paths(data_dir)
      self.assertEqual(mp, data_dir + "/" + rel_path + "/")
      self.assertEqual(fp, data_dir + "/" + rel_path + "/frames/")


   def test_frame_data_3D(self):
      """
      Test the frame_data_3D routine, which extracts some data and may
      do some processing.
      """

      params = {"title" : "title",
                "stub" : "stub",
                "path" : "path",
                "image_type" : "png",
                "variable" : "density",
                "mode" : "profile : contrast",
                "window_x_lo" : -1.2,
                "window_x_hi" :  1.2,
                "window_y_lo" : -0.5,
                "window_y_hi" :  1,
                "window_z_lo" : -1,
                "window_z_hi" :  1,
                "time_lo" : 13,
                "time_hi" : 73,
                "value_lo" : 1.5,
                "value_hi" : 12.5,
                "make_movie" : True,
                "fps" : 60,
                "movie_type" : "avi",
                "final_pause" : "16.5s",
                "masks" : [],
                "mask_method" : "force low",
                "xlines" : [-1, -0.5, 0, 0.5, 1]
                }
      movie = MovieDescriptor(params, [])

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
      ss = SimulationState(data, si)

      # Window
      mdata, x, y, z = movie.frame_data_3D(ss)
      self.assertEqual(mdata.shape, (3,2,1))

      # Mask
      ttt = 1.2
      mask = MaskDescriptor({"variable" : "density", "threshold" : ttt,
         "operator" : "<="})
      params["variable"] = "density"
      params["masks"] = ["mask"]
      params["mask_method"] = "force low"
      movie = MovieDescriptor(params, {"mask" : mask})
      mdata, x, y, z = movie.frame_data_3D(ss)
      for i in xrange(mdata.shape[0]):
         for j in xrange(mdata.shape[1]):
            for k in xrange(mdata.shape[2]):
               d = dens[i,j,k]
               if d <= ttt:
                  d = 0
               self.assertEqual(mdata[i,j,k], d)

   # Building the frames is best evaluated by looking at the resulting images
   # generated by test runs, rather than by unittest.  See the sample data set
   # and the reference images.

if __name__ == "__main__":
   unittest.main()
