"""
Driver for making movies.

This parses the command line, does some housekeeping, compiles the
appropriate descriptors from input files, handles parallelization, and
does other "outermost loop" tasks of this sort to set up what's
necessary to run the movie maker.
"""

# Python insists that __future__ import(s) should come first
from __future__ import division

# Load site_setup (if it exists) to handle any site-specific details
try:
   import site_setup
except ImportError:
   pass

import argparse as ap
import glob
import os
import tomlpython as toml
import subprocess
import sys
import warnings

import Descriptors as Desc
from MovieLoops import draw_frames_loop, encode_movies_loop

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

def process_command_line(argv=None):
   """
   Parse command-line arguments.

   Arguments:
      argv (list) : list of arguments to parse; if None, then parse sys.argv

   Returns:
      directory (string) : name of data directory
      dfile_list (list of strings) : list of descriptor file names
      movie_list (list of strings) : list of names of movies to make

   Exceptions:
      whatever arises from functions called by this routine
   """

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

   if argv is None:
      args = parser.parse_args()
   else:
      args = parser.parse_args(argv)

   # Because some components assume that the directory path ends with "/"
   if args.directory[-1] is not "/":
      args.directory += "/"

   return args.directory, args.dfile_list, args.movie_list

#==============================================================================

def build_descriptors(dfile_list, movie_list):
   """
   Build the mask and movie descriptors.

   Arguments:
      dfile_list (list of strings) : list of descriptor file names
      movie_list (list of strings) : list of names of movies to make

   Returns:
      masks (dict) : {name (from descriptor file) : MaskDescriptor}
      movies (dict) : {name (from descriptor file) : MovieDescriptor}

   Exceptions:
      whatever arises from functions called by this routine
   """

   # Get the movie and mask descriptions
   movie_descriptions = {}
   movie_blacklist = set()
   mask_descriptions = {}
   mask_blacklist = set()
   for dfile in dfile_list:
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
      if movie_list is None or name in movie_list:
         movies[name] = Desc.MovieDescriptor(movie_string, masks)
   if movie_list is not None:
      for m in movie_list:
         if m not in movies:
            msg = '"'.join(('Requested movie ', m, ' not found.'))
            warnings.warn(msg, UserWarning)

   return masks, movies

#==============================================================================

def get_data_list(directory):
   """
   Get the list of data files to visualize and distribute across processors.

   Arguments:
      directory (string) : directory to search for DUMSES output files

   Returns:
      output_list (list of strings) : names of all output files to visualize

   Exceptions:
      whatever arises from functions called by this routine
   """

   if ProcID == 0:
      # Find all the data files in the directory
      full_list = sorted(glob.glob(directory + "output_*"))

      # Share them among all the processors
      if NProcs > 1:
         Nout = len(full_list)
         if NProcs > Nout:
            msg = " ".join(("More processors than output files; processors",
               "{0} through {1} will do no work.".format(Nout, NProcs-1)))
            warnings.warn(msg, UserWarning)
         output_list = full_list[::NProcs]
         for i in xrange(1, NProcs):
            pypar.send(full_list[i::NProcs], destination=i)
      else:
         output_list = full_list
   else:
      output_list = pypar.receive(0)

   return output_list

#==============================================================================

def assign_movies(encode_locations):
   """
   Redistribute the movies to be encoded among the processors available.

   Arguments:
      encode_locations (set) : set of movies to encode

   Returns:
      encode_list (list) : list of movies to be encoded by local processor

   Exceptions:
      whatever arises from functions called by this routine
   """

   # Merge the set of movies to processor zero
   # -- The structure of the code suggests that, in many cases, all processors
   #    will have identical copies of encode_locations prior to this step.
   #    However, that is not guaranteed.  Consider, for example, running this
   #    script on 10 processors with 10 data files, and one of the movies has a
   #    time window that excludes the first 5 data files: processors 0..4 will
   #    not know about that movie unless informed by processors 5..9.
   # -- Make sure they send in-order.  PyPar documentation is a bit thin on the
   #    ground, so I have not yet determined whether or not PyPar's send and
   #    receive are nonblocking.  Other methods (e.g. gather) could also be
   #    used, and would arguably improve performance, but until I find better
   #    PyPar documentation, this will suffice.  And this is absolutely not
   #    even close to a performance-critical operation; the vast majority of
   #    the computing time is spent generating the frames, while communication
   #    is almost negligible.  Hence why I've spent my time on other
   #    components, rather than looking for better PyPar documentation.
   if NProcs > 1:
      pypar.barrier()
   if ProcID == 0:
      for i in xrange(1, NProcs):
         temp_set = pypar.receive(i)
         encode_locations.update(temp_set)
         pypar.barrier()
   else:
      for i in xrange(1, NProcs):
         if ProcID == i:
            pypar.send(encode_locations, destination=0)
         pypar.barrier()

   # Split the movies across the processors
   if NProcs > 1:
      pypar.barrier()
   if ProcID == 0:
      full_list = sorted(list(encode_locations))
      if NProcs > 1:
         Nenc = len(full_list)
         if NProcs > Nenc:
            msg = " ".join(("More processors than movies to encode;"
               "processors {0} through {1} will do no work.".format(
               Nenc, NProcs-1)))
            warnings.warn(msg, UserWarning)
         encode_list = full_list[::NProcs]
         for i in xrange(1, NProcs):
            pypar.send(full_list[i::NProcs], destination=i)
      else:
         encode_list = full_list
   else:
      encode_list = pypar.receive(0)

   return encode_list

#==============================================================================

def main():
   """
   Main driver function.
   """

   # Set up - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   # Parse the command-line arguments
   directory, dfile_list, movie_list = process_command_line(sys.argv[1:])

   # Build the MaskDescriptors and MovieDescriptors
   masks, movies = build_descriptors(dfile_list, movie_list)

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
      sys.stdout.flush()

   # Get the list of data files
   output_list = get_data_list(directory)

   # Summarize outputs list
   if NProcs > 1:
      pypar.barrier()
   for i in xrange(NProcs):
      if i == ProcID:
         print "Processor {p} has {n} output files:".format(p=ProcID,
               n=len(output_list))
         for o in output_list:
            print "    {0}".format(o)
         sys.stdout.flush()
      if NProcs > 1:
         pypar.barrier()
   if ProcID == 0:
      print "="*79

   # Call the primary loop - - - - - - - - - - - - - - - - - - - - - - - - - -
   if len(movies) == 0:
      if ProcID == 0:
         sys.stdout.write("No movies to generate.\n")
   else:
      encode_locations = set()
      draw_frames_loop(output_list, movies, encode_locations)

   # Encode the frames to movies - - - - - - - - - - - - - - - - - - - - - - -

   encode_list = assign_movies(encode_locations)

   # Summarize movies to encode
   if ProcID == 0:
      print "="*79
   if NProcs > 1:
      pypar.barrier()
   for i in xrange(NProcs):
      if i == ProcID:
         print "Processor {p} has {n} movies to encode:".format(p=ProcID,
               n=len(encode_list))
         for o in encode_list:
            print "    {0}".format(o[0])
         sys.stdout.flush()
      if NProcs > 1:
         pypar.barrier()

   encode_movies_loop(encode_list)

   if NProcs > 1:
      pypar.barrier()
      pypar.finalize()

#==============================================================================

if __name__ == "__main__":
   main()

