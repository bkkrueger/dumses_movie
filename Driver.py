"""
Driver for making movies

This parses the command line, does some housekeeping, compiling the appropriate
descriptors from input files, any parallelization, and "outermost loop" tasks
of this sort to set up what's necessary to run the movie maker.
"""

# TODO : Rewrite the header comment

# TODO : Develop an extensive testing suite

# TODO : Check all class docstrings and add methods list

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
import warnings

import Descriptors as Desc
from MakeFrames import make_all_frames

#==============================================================================

def pretty(d, indent=0):
   """
   Pretty-printing for nested dictionaries.

   A convenient helper routine for testing.
   """
   for key, value in d.iteritems():
      print '   ' * indent + str(key)
      if isinstance(value, dict):
         pretty(value, indent+1)
      else:
         print '   ' * (indent+1) + str(value)

#==============================================================================

if __name__ == "__main__":

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
   output_list = sorted(glob.glob(args.directory + "output_*"))
   # TODO : Future development will parallelize by having a master processor
   #        get a complete list of the DUMSES outputs, then distribute that
   #        list among all processors (including itself), and then each
   #        processor will proceed as usual from there.

   # Get the movie descriptions
   with open(args.mfile) as tomlfile:
      descriptions = toml.parse(tomlfile)

   # Construct the descriptors
   masks = {}
   for m in descriptions["masks"]:
      masks[m] = Desc.MaskDescriptor(descriptions["masks"][m])

   movies = {}
   for m in descriptions["movies"]:
      if args.movie_list is None:
         movies[m] = Desc.MovieDescriptor(descriptions["movies"][m], masks)
      else:
         if m in args.movie_list:
            movies[m] = Desc.MovieDescriptor(descriptions["movies"][m], masks)
   if args.movie_list is not None:
      for m in args.movie_list:
         if m not in movies:
            msg = '"'.join(('Movie ', m, ' not found in file ',
               args.mfile, '.'))
            warnings.warn(msg, UserWarning)

   # Call the primary loop
   make_all_frames(output_list, movies)

   # Encode the movies
   # TODO : For every movie that calls for this step, encode the individual
   # frame images into a movie.

