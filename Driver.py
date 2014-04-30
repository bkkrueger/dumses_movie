"""
Driver for making movies.

This parses the command line, does some housekeeping, compiles the appropriate
descriptors from input files, handles parallelization, and does other
"outermost loop" tasks of this sort to set up what's necessary to run the movie
maker.
"""

# TODO : Develop an extensive testing suite

from __future__ import division

#==============================================================================
#------------------------------------------------------------------------------
# Version 5 of DumPy is not normally in the path on this system, so it needs to
# be explicitly added.  Because this is a for development as part of a local
# test, DumPy_v05 is not set up site-wide, so using the site package and a PTH
# file is not the best plan.  Thus I just use sys.path to add the location of
# the test version of DumPy_v05.
# TODO : take this out to test on other systems
import sys
sys.path.append("/Users/bkrueger/research/CCSNe/heating_layer/results/sample/dumses/visu")
#------------------------------------------------------------------------------
#==============================================================================

import argparse as ap
import glob
import tomlpython as toml
import pypar
import warnings

import Descriptors as Desc
from MakeFrames import make_all_frames

#==============================================================================

def pretty(d, indent=0):
   """
   Pretty-printing for nested dictionaries; convenient for testing.
   """
   for key, value in d.iteritems():
      print '   ' * indent + str(key)
      if isinstance(value, dict):
         pretty(value, indent+1)
      else:
         print '   ' * (indent+1) + str(value)

#==============================================================================

if __name__ == "__main__":

   # Parallel setup - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
   NProcs = pypar.size()
   ProcID = pypar.rank()

   # Parse command-line arguments - - - - - - - - - - - - - - - - - - - - - - -
   prog_desc = "Generate a set of movies from DUMSES data files."
   parser = ap.ArgumentParser(description=prog_desc)

   # Specify the directory containing the outputs
   parser.add_argument("--dir", metavar="DIR", default="./",
         help="directory containing DUMSES outputs", required=False,
         dest="directory")

   # Specify the file describing the movies
   parser.add_argument("--mfile", metavar="MFILE", default="movies.toml",
         help="file describing movies to be made", required=False,
         dest="mfile")

   # Specify which movies to generate
   parser.add_argument("--movies", nargs="+",
         help="list of IDs of movies to create", dest="movie_list")

   args = parser.parse_args()

   # Because some components assume that the directory path ends with "/"
   if args.directory[-1] is not "/":
      args.directory += "/"

   # Set up - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   # Generate the list of DUMSES outputs to plot.
   #output_list = sorted(glob.glob(args.directory + "output_*"))
   if ProcID == 0:
      full_list = sorted(glob.glob(args.directory + "output_*"))
      if NProcs > 1:
         ihi = len(full_list) // NProcs
         output_list = full_list[:ihi]
         for i in xrange(1, NProcs):
            ilo = len(full_list) * i // NProcs
            ihi = len(full_list) * (i+1) // NProcs
            pypar.send(full_list[ilo:ihi], destination=i)
      else:
         output_list = full_list
   else:
      output_list = pypar.receive(0)

   # Get the movie and mask descriptions
   with open(args.mfile) as tomlfile:
      descriptions = toml.parse(tomlfile)

   # Construct the mask descriptors
   masks = {}
   for name, mask_string in descriptions["masks"].items():
      masks[name] = Desc.MaskDescriptor(mask_string)

   # Construct the movie descriptors
   # -- If the user specified a list of movies, only make those movies.  Warn
   #    the user if a requested movies is missing, but continue anyway.
   # -- If the user did not specify a list of movies, make all movies described
   #    in the movies file.
   movies = {}
   for name, movie_string in descriptions["movies"].items():
      if args.movie_list is None or name in args.movie_list:
         movies[name] = Desc.MovieDescriptor(movie_string, masks)
   if args.movie_list is not None:
      for m in args.movie_list:
         if m not in movies:
            msg = '"'.join(('Movie ', m, ' not found in file ',
               args.mfile, '.'))
            warnings.warn(msg, UserWarning)

   # Summarize masks and movies
   if ProcID == 0:
      print "="*79
      print "found {n} masks:".format(n=len(masks))
      for name, mask in masks.items():
         print "   {0}:".format(name), mask
      print "making {n} movies:".format(n=len(movies))
      for name, movie in movies.items():
         print "   {0}:".format(name), movie
      print "="*79

   # Call the primary loop
   make_all_frames(output_list, movies)

   # Encode the movies
   # TODO : For every movie that calls for this step, encode the individual
   # frame images into a movie.

   pypar.barrier()
   pypar.finalize()

