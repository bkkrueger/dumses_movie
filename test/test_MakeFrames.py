"""
This module performs a series of tests on the MakeFrames object.
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

      def draw_frame_stub(self, state, path, number, state0):
         """
         A simplified stub for MovieDescriptor.draw_frame.
         """
         print "CALLED DRAW_FRAME"
         print "   state  : {0}".format(state)
         print "   path   : {0}".format(path)
         print "   number : {0}".format(number)
         print "   state0 : {0}".format(state0)
         # TODO This should be updated to do some assertions

      MovieDescriptor.draw_frame = draw_frame_stub

   def test_make_all_frames(self):
      """
      Some tests
      """
      # TODO Fix this comment

      # Make a list of data file names
      data_list = glob.glob("../../results/sample/data/output_*")

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
      movie_list = dict()
      movie_list["movie1"] = MovieDescriptor(
         {"stub" : "movie1",
         "image_type" : "png",
         "variable" : "enthalpy",
         "make_movie" : False
         }, mask_list)
      movie_list["movie2"] = MovieDescriptor(
         {"stub" : "movie2",
         "image_type" : "png",
         "variable" : "enthalpy",
         "make_movie" : False,
         "masks" : ["mask1", "mask2"]
         }, mask_list)

      # Call the function
      encode_locations = set()
      make_all_frames(data_list, movie_list, encode_locations)

   # TODO : Test some things that should work and some things that shouldn't

if __name__ == "__main__":
   unittest.main()
