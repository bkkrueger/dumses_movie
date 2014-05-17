"""
This module performs a series of tests on the Driver module.
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
import MakeFrames
from Driver import main

# TODO : Don't forget to clean up the comments/docstrings for Driver.py.

class DriverTest(unittest.TestCase):
   """
   Test the Driver module.
   """

   def setUp(self):
      """
      Pre-testing setup.
      """

      # TODO

      # Replace the existing draw_frame with a stub routine
      def make_all_frames_stub(output_list, movies, encode_locations):
         """
         A simplified stub for MakeFrames.make_all_frames, which merely logs
         the relevant information into a garbage variable that doesn't normally
         exist in MakeFrames.
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

   def test_process_command_line(self):
      """
      Test that the command-line parser works as expected
      """

      # TODO
      # Note : The process_command_line function takes either None (defaults to
      #        letting parse_args choose for itself) or a list specifying the
      #        arguments.  It returns the directory holding the data, a list of
      #        descriptor files, and a list of movie names.
      # Note : Also ensure that the directory has a trailing slash to make
      #        DumsesData happy.

   def test_build_descriptors(self):
      """
      Test the function that builds the mask and movie descriptors
      """
      # TODO
      # Note : The build_descriptors function takes a list of descriptor files
      #        and a list of movie names.  It reads the descriptor files and
      #        builds the MaskDescriptors and MovieDescriptors.  If the list of
      #        movie names is None, it returns those descriptors.  If the list
      #        of movie names is actually a list, then it filters out
      #        MovieDescriptors that are not on the list (the MovieDescriptors
      #        have a name in the descriptor files and those names are given in
      #        the list of movie names); it also warns the user if requested
      #        movies are missing.

   def test_get_data_list(self):
      """
      Test the function that constructs the list of data files.
      """

      # TODO
      # Note : The get_data_list function takes a string specifying the
      #        directory to search for data files, gets all the "output_*"
      #        files in that directory, and divides them as evenly as possible
      #        among the processors.

   def test_assign_movies(self):
      """
      Test the function that assigns movies to encode among the processors.
      """

      # TODO
      # Note : The assign_movies function takes a set of encode locations,
      #        merges the sets among the different processors onto processors
      #        zero, then redistributes the encode locations as evenly as
      #        possible among the processors, returning a list of encode
      #        locations.

   def test_encode_movies(self):
      """
      Test the function that encodes the frames into movies.
      """

      # TODO
      # Note : The encode_movies function takes a list of encode locations and
      #        encodes the specified movies.  It returns nothing.

   def test_driver(self):
      """
      Test the overall functioning of the driver.
      """

      # TODO

if __name__ == "__main__":
   unittest.main()
