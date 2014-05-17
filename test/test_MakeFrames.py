"""
This module performs a series of tests on the MakeFrames module.
It will not be exhaustive, but it should cover basic usage and catch
the most glaring errors.
"""

import unittest
import copy
import glob
import numpy as np
import sys
sys.path.append("..")
import site_setup
from Descriptors import MovieDescriptor
from MakeFrames import make_all_frames

class MakeFramesTest(unittest.TestCase):
   """
   Test the MakeFrames module.
   """

   def setUp(self):
      """
      Pre-testing setup.
      """

      # Replace the existing draw_frame with a stub routine
      def draw_frame_stub(self, state, path, number, state0):
         """
         A simplified stub for MovieDescriptor.draw_frame, which merely logs
         the name of the file that it would make into a garbage variable that
         doesn't normally exist in a MovieDescriptor.
         """
         mp, fp = self.build_paths(path)
         image_file_name = "{p}{s}_{n:06d}.{e}".format(p=fp, s=self.stub,
               n=number, e=self.image_type)
         try:
            self.temporary_storage.append(image_file_name)
         except AttributeError:
            self.temporary_storage = list([image_file_name])

         try:
            self.found_state0 = self.found_state0 or (state0 is not None)
         except AttributeError:
            self.found_state0 = (state0 is not None)

      MovieDescriptor.draw_frame = draw_frame_stub

   def test_make_all_frames(self):
      """
      Test the make_all_frames routine, which encapsulates some
      checking/filtering and the main loop.
      """

      data_dir = "/".join(("", "Users", "bkrueger", "research", "CCSNe",
         "heating_layer", "results", "sample", "data"))

      # Make a list of data file names
      data_list = glob.glob("/".join((data_dir, "output_*")))

      # Make a list of MaskDescriptors
      mask_list = dict()
      mask_list["mask1"] = {"variable" : "density",
            "mode" : "pseudocolor: perturbation",
            "operator" : "<=",
            "threshold" : 3.2}
      mask_list["mask2"] = {"variable" : "entropy",
            "mode" : "pseudocolor: contrast",
            "operator" : ">",
            "threshold" : 0.0}
      mask_list["mask3"] = {"variable" : "velocity magnitude",
            "mode" : "pseudocolor: full state",
            "operator" : "<",
            "threshold" : 0.09}

      # Make a list of MovieDescriptors
      # Testing different things:
      # -- different stub, same path
      # -- same stub, different path
      # -- same stub, same path
      movie_list = dict()
      movie_list["movie1"] = MovieDescriptor(
         {"stub" : "movie1",
         "path" : "path1",
         "image_type" : "png",
         "variable" : "enthalpy",
         "make_movie" : False
         }, mask_list)
      movie_list["movie2"] = MovieDescriptor(
         {"stub" : "movie1",
         "path" : "path2",
         "image_type" : "png",
         "variable" : "enthalpy",
         "make_movie" : False,
         "masks" : ["mask1", "mask2"]
         }, mask_list)
      movie_list["movie3"] = MovieDescriptor(
         {"stub" : "movie2",
         "path" : "path1",
         "image_type" : "png",
         "variable" : "enthalpy",
         "make_movie" : False,
         "time_lo" : 99,
         "time_hi" : 201
         }, mask_list)
      movie_list["movie4"] = MovieDescriptor(
         {"stub" : "movie1",
         "path" : "path1",
         "image_type" : "png",
         "variable" : "enthalpy",
         "make_movie" : False,
         "masks" : ["mask1", "mask2"]
         }, mask_list)

      # Add the attributes that the draw_frames stub will use
      for movie in movie_list.values():
         movie.temporary_storage = list()
         movie.found_state0 = False

      # Call the function
      encode_locations = set()
      make_all_frames(data_list, movie_list, encode_locations)

      # Make sure the correct frames were "drawn" (note: the sample data set is
      # known to have output files from 0 to 120, inclusive, with steps of
      # about 2.5s between each)

      # movie 1 collides with movie 4, so it is excluded
      movie = movie_list["movie1"]
      movie.assumed_list = list()
      self.assertSequenceEqual(movie.assumed_list, movie.temporary_storage)
      self.assertFalse(movie.found_state0)

      # movie 2
      movie = movie_list["movie2"]
      movie.assumed_list = list()
      mp, fp = movie.build_paths(data_dir)
      for n in xrange(121):
         image_file_name = "{p}{s}_{n:06d}.{e}".format(
               p=fp, s=movie.stub, n=n, e=movie.image_type)
         movie.assumed_list.append(image_file_name)
      self.assertSequenceEqual(movie.assumed_list, movie.temporary_storage)
      self.assertFalse(movie.found_state0)

      # movie 3 only includes t = 99..201s (files 40..80)
      movie = movie_list["movie3"]
      movie.assumed_list = list()
      mp, fp = movie.build_paths(data_dir)
      for n in xrange(40,81):
         image_file_name = "{p}{s}_{n:06d}.{e}".format(
               p=fp, s=movie.stub, n=n, e=movie.image_type)
         movie.assumed_list.append(image_file_name)
      self.assertSequenceEqual(movie.assumed_list, movie.temporary_storage)
      self.assertFalse(movie.found_state0)

      # movie 4 collides with movie 1, so it is excluded
      movie = movie_list["movie4"]
      movie.assumed_list = list()
      self.assertSequenceEqual(movie.assumed_list, movie.temporary_storage)
      self.assertFalse(movie.found_state0)

      # This pass will test that the initial state is being used for profiles
      movie_list = dict()
      movie_list["movie1"] = MovieDescriptor(
         {"stub" : "movie1",
         "path" : "path1",
         "image_type" : "png",
         "variable" : "enthalpy",
         "mode" : "profile : full state",
         "make_movie" : False
         }, mask_list)

      encode_locations = set()
      make_all_frames(data_list, movie_list, encode_locations)

      # movie 1
      movie = movie_list["movie1"]
      movie.assumed_list = list()
      mp, fp = movie.build_paths(data_dir)
      for n in xrange(121):
         image_file_name = "{p}{s}_{n:06d}.{e}".format(
               p=fp, s=movie.stub, n=n, e=movie.image_type)
         movie.assumed_list.append(image_file_name)
      self.assertSequenceEqual(movie.assumed_list, movie.temporary_storage)
      self.assertTrue(movie.found_state0)

if __name__ == "__main__":
   unittest.main()
