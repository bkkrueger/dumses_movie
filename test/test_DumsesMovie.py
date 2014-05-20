"""
This module performs a series of tests on the DumsesMovie module.
"""

import unittest
import glob
import sys
sys.path.append("..")
import site_setup
import DumsesMovie

class DumsesMovieTest(unittest.TestCase):
   """
   Test the DumsesMovie module.
   """

   #===========================================================================

   def test_process_command_line(self):
      """
      Test that the command-line parser works as expected
      """

      d = "/path/subpath/data_dir"
      df = ["dfile1.toml", "dfile2.toml"]
      m = ["movie1", "movie2", "movie3"]
      arg_list = " ".join(("--dir", d,
         "--dfile", " ".join(df),
         "--movies", " ".join(m))).split()

      directory, dfiles, movies = DumsesMovie.process_command_line(arg_list)

      if d[-1] != "/":
         d += "/"
      self.assertEqual(directory, d)
      self.assertSequenceEqual(dfiles, df)
      self.assertSequenceEqual(movies, m)

   #===========================================================================

   def test_build_descriptors(self):
      """
      Test the function that builds the mask and movie descriptors
      """

      masks, movies = DumsesMovie.build_descriptors(
            ["descriptors.toml"], ["entr"])

      self.assertSequenceEqual(masks.keys(), ["buoyant"])
      self.assertSequenceEqual(movies.keys(), ["entr"])

   #===========================================================================

   def test_get_data_list(self):
      """
      Test the function that constructs the list of data files.
      """
      data_dir = "/".join(("", "Users", "bkrueger", "research", "CCSNe",
         "heating_layer", "results", "sample", "data", ""))
      full_list = sorted(glob.glob(data_dir + "output_*"))
      partial_list = DumsesMovie.get_data_list(data_dir)
      self.assertSequenceEqual(partial_list,
            full_list[DumsesMovie.ProcID::DumsesMovie.NProcs])

   #===========================================================================

   def test_assign_movies(self):
      """
      Test the function that assigns movies to encode among the processors.
      """

      set_in = set()
      if DumsesMovie.ProcID > 0:
         set_in.add(DumsesMovie.ProcID-1)
      set_in.add(DumsesMovie.ProcID)
      if DumsesMovie.ProcID < DumsesMovie.NProcs-1:
         set_in.add(DumsesMovie.ProcID+1)
      list_out = DumsesMovie.assign_movies(set_in)
      self.assertSequenceEqual(list_out, [DumsesMovie.ProcID])

#==============================================================================

if __name__ == "__main__":
   unittest.main()
