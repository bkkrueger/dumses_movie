"""
This module performs a series of tests on the EncodeInfo class.
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
from MovieLoops import EncodeInfo

class EncodeInfoTest(unittest.TestCase):
   """
   Test the EncodeInfo class.
   """

   def test_creation(self):
      """
      Test the creation of EncodeInfo.
      """

      name = "/path/subpath/movie.ext"
      regex = "/path/subpath/frames/frame_*.ext"
      fps = 60
      ei = EncodeInfo(name, regex, fps)

      self.assertEqual(ei.movie_name, name)
      self.assertEqual(ei.frame_regex, regex)
      self.assertEqual(ei.fps, fps)

   def test_equality(self):
      """
      Test equality and inequalty.
      """

      ei1 = EncodeInfo("/path/subpath/movie.ext",
            "/path/subpath/frames/frame_*.ext", 60)
      self.assertEqual(ei1.movie_name, ei1.movie_name)
      self.assertEqual(ei1.frame_regex, ei1.frame_regex)
      self.assertEqual(ei1.fps, ei1.fps)
      self.assertTrue(ei1==ei1)
      self.assertFalse(ei1!=ei1)

      ei2 = EncodeInfo("/path/subpath/movie.ext",
            "/path/subpath/frames/frame_*.ext", 15)
      self.assertEqual(ei1.movie_name, ei2.movie_name)
      self.assertEqual(ei1.frame_regex, ei2.frame_regex)
      self.assertNotEqual(ei1.fps, ei2.fps)
      self.assertFalse(ei1==ei2)
      self.assertTrue(ei1!=ei2)

      ei3 = EncodeInfo("/path/subpath/movie.ext",
            "/path/subpath/frames2/frame_*.ext", 60)
      self.assertEqual(ei1.movie_name, ei3.movie_name)
      self.assertNotEqual(ei1.frame_regex, ei3.frame_regex)
      self.assertEqual(ei1.fps, ei3.fps)
      self.assertFalse(ei1==ei3)
      self.assertTrue(ei1!=ei3)

      ei4 = EncodeInfo("/path/subpath/movie2.ext",
            "/path/subpath/frames/frame_*.ext", 60)
      self.assertNotEqual(ei1.movie_name, ei4.movie_name)
      self.assertEqual(ei1.frame_regex, ei4.frame_regex)
      self.assertEqual(ei1.fps, ei4.fps)
      self.assertFalse(ei1==ei4)
      self.assertTrue(ei1!=ei4)

      ei5 = EncodeInfo("/path/subpath/movie.ext",
            "/path/subpath/frames/frame_*.ext", 60)
      self.assertEqual(ei1.movie_name, ei5.movie_name)
      self.assertEqual(ei1.frame_regex, ei5.frame_regex)
      self.assertEqual(ei1.fps, ei5.fps)
      self.assertTrue(ei1==ei5)
      self.assertFalse(ei1!=ei5)

      self.assertFalse(ei1=="no")
      self.assertTrue(ei1!="no")

   def test_set(self):
      """
      Test EncodeInfo behavior in sets
      """

      ei1 = EncodeInfo("/path/subpath/movie.ext",
            "/path/subpath/frames/frame_*.ext", 60)
      ei2 = EncodeInfo("/path/subpath/movie.ext",
            "/path/subpath/frames/frame_*.ext", 15)
      ei3 = EncodeInfo("/path/subpath/movie.ext",
            "/path/subpath/frames2/frame_*.ext", 60)
      ei4 = EncodeInfo("/path/subpath/movie2.ext",
            "/path/subpath/frames/frame_*.ext", 60)
      ei5 = EncodeInfo("/path/subpath/movie.ext",
            "/path/subpath/frames/frame_*.ext", 60)

      s = set()
      self.assertEqual(len(s), 0)
      s.add(ei1)
      self.assertEqual(len(s), 1)
      s.add(ei1)
      self.assertEqual(len(s), 1)
      s.add(ei2)
      self.assertEqual(len(s), 2)
      s.add(ei3)
      self.assertEqual(len(s), 3)
      s.add(ei4)
      self.assertEqual(len(s), 4)
      s.add(ei5)  # Equal to 1
      self.assertEqual(len(s), 4)

if __name__ == "__main__":
   unittest.main()
