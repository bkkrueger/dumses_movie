"""
Driver for making movies.

This parses the command line, does some housekeeping, compiles the appropriate
descriptors from input files, handles parallelization, and does other
"outermost loop" tasks of this sort to set up what's necessary to run the movie
maker.
"""

# TODO : Develop an extensive testing suite

# Python insists that __future__ import(s) should come first
from __future__ import division

# Load site_setup (if it exists) to handle any site-specific details
try:
   import site_setup
except ImportError:
   pass

import argparse as ap
import glob
import tomlpython as toml
import warnings

import Descriptors as Desc
from MakeFrames import make_all_frames

#==============================================================================
# Note on PyPar:
#    I haven't found much documentation or other information online about this,
# but it appears that calling certain PyPar functions causes a segmentation
# fault if this script is launched directly (i.e., "python Driver.py" and not
# "mpirun -np # python Driver.py").  Fortunately, the size and rank functions
# work either way.  Because I want this script to be able to run on a single
# processor (both for users who don't need parallelization and for users who
# have not installed PyPar), I handle the importing of PyPar as a special case,
# and all calls to PyPar functions are wrapped in a check that more than one
# processor exists.
try:
   import pypar
except ImportError:
   NProcs = 1
   ProcID = 0
   warnings.warn("Unable to import PyPar.  Running in single-processor mode.",
         UserWarning)
else:
   NProcs = pypar.size()
   ProcID = pypar.rank()
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

   # Parse command-line arguments - - - - - - - - - - - - - - - - - - - - - - -
   prog_desc = "Generate a set of movies from DUMSES data files."
   parser = ap.ArgumentParser(description=prog_desc)

   # Specify the directory containing the outputs
   parser.add_argument("--dir", metavar="DIR", default="./",
         help="directory containing DUMSES outputs", required=False,
         dest="directory")

   # Specify the file describing the movies
   parser.add_argument("--dfile", metavar="DFILE", nargs="+",
         default=["descriptors.toml"], required=False, dest="dfile_list",
         help="file(s) describing movies and masks")

   # Specify which movies to generate
   parser.add_argument("--movies", metavar="MOVIE", nargs="+",
         dest="movie_list", help="list of IDs of movies to create")

   args = parser.parse_args()

   # Because some components assume that the directory path ends with "/"
   if args.directory[-1] is not "/":
      args.directory += "/"

   # Set up - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   # Generate the list of DUMSES outputs to plot.
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
   movie_descriptions = {}
   movie_blacklist = set()
   mask_descriptions = {}
   mask_blacklist = set()
   for dfile in args.dfile_list:
      with open(dfile) as tomlfile:
         d = toml.parse(tomlfile)
         # Make sure movies are not redundant
         for k, v in d.get("movies", {}).items():
            if k in movie_descriptions and k not in movie_blacklist:
               if ProcID == 0:
                  msg = '"'.join(('Movie ', k, ' is defined multiple times.  ',
                     'It will be removed from the list of movies.'))
                  warnings.warn(msg, UserWarning)
               movie_blacklist.add(k)
            else:
               movie_descriptions[k] = v
         # Make sure masks are not redundant
         for k, v in d.get("masks", {}).items():
            if k in mask_descriptions and k not in mask_blacklist:
               if ProcID == 0:
                  msg = '"'.join(('Mask ', k, ' is defined multiple times.  ',
                     'It will be removed from the list of masks.'))
                  warnings.warn(msg, UserWarning)
               mask_blacklist.add(k)
            else:
               mask_descriptions[k] = v
   for k in movie_blacklist:     # Remove redundant movies
      movie_descriptions.pop(k)
   for k in mask_blacklist:      # Remove redundant masks
      mask_descriptions.pop(k)

   # Construct the mask descriptors
   masks = {}
   for name, mask_string in mask_descriptions.items():
      masks[name] = Desc.MaskDescriptor(mask_string)

   # Construct the movie descriptors
   # -- If the user specified a list of movies, only make those movies.  Warn
   #    the user if a requested movies is missing, but continue anyway.
   # -- If the user did not specify a list of movies, make all movies described
   #    in the movies file.
   movies = {}
   for name, movie_string in movie_descriptions.items():
      if args.movie_list is None or name in args.movie_list:
         movies[name] = Desc.MovieDescriptor(movie_string, masks)
   if args.movie_list is not None:
      for m in args.movie_list:
         if m not in movies:
            msg = '"'.join(('Requested movie ', m, ' not found.'))
            warnings.warn(msg, UserWarning)

   # Summarize masks and movies
   if ProcID == 0:
      print "="*79
      print "found {n} masks:".format(n=len(masks))
      for name, mask in sorted(masks.items()):
         print "   {0}:".format(name), mask
      print "making {n} movies:".format(n=len(movies))
      for name, movie in sorted(movies.items()):
         print "   {0}:".format(name), movie
      print "="*79

   # Summarize outputs list
   if NProcs > 1:
      pypar.barrier()
   for i in xrange(NProcs):
      if i == ProcID:
         print "Processor {p} has {n} output files:".format(p=ProcID,
               n=len(output_list))
         for o in output_list:
            print "    {0}".format(o)
      if NProcs > 1:
         pypar.barrier()

   # Call the primary loop
   make_all_frames(output_list, movies)

   # Encode the movies
   # TODO : For every movie that calls for this step, encode the individual
   # frame images into a movie.

   if NProcs > 1:
      pypar.barrier()
      pypar.finalize()

