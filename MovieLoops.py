"""
Make the frames.

This module manages the loops and control logic for making the frames.

Attributes:
   draw_frames_loop : draws all frames for lists of outputs and movies
   encode_movies_loop : encodes all movies for list of frame sources
"""

import errno
import os
import subprocess
import sys
import warnings

from dumpy_v05.data.rd_dumses import DumsesData
from Descriptors import DescriptorError
from SimulationData import SimulationError, SimulationInput, SimulationState

#------------------------------------------------------------------------------
#==============================================================================
class EncodeInfo(object):
   """
   The information necessary to encode a movie.

   This class contains the information necessary to encode a movie.
   This is separated out from the MovieDescriptor, because it is an
   intermediate step in the movie creation process and not part of the
   movie description.  This is a basic convenience class to encapsulate
   the information that needs to be communicated between the
   frame-drawing loop and the movie-encoding loop, with some routines
   for comparison and hashing to enable the use in sets.

   Attributes:
      movie_name (string) : file name for the movie (with full path)
      frame_regex (string) : the regex for the frames (wth full path)
      fps (int) : the frame rate in frames per second
   """

   #===========================================================================
   def __repr__(self):
      """
      Supply a detailed representation of the class.
      """
      string = "Encoding information for movie {0}".format(self.movie_name)
      return string

   #===========================================================================
   def __str__(self):
      """
      Supply a simple representation of the class.
      """
      string = "Encoding information for movie {0}".format(self.movie_name)
      return string

   #===========================================================================
   def __init__(self, name, regex, fps):
      """
      Construct the class.

      Arguments:
         name (str) : the name of the movie (with full path)
         regex (str) : the regular expression for the frames (with full path)
         fps (int) : the frame rate in frames per second
      """

      self.movie_name = name
      self.frame_regex = regex
      self.fps = fps

   #===========================================================================
   def __hash__(self):
      """
      Compute hash value.

      Since tuples and all the components of EncodeInfo are all individually
      hashable, be lazy and put all the components into a tuple and use
      Python's built-in hashing routines to take care of this.
      """

      return hash((self.movie_name, self.frame_regex, self.fps))

   #===========================================================================
   def __eq__(self, other):
      """
      Test equality.
      """

      try:
         if self.movie_name != other.movie_name:
            return False
         elif self.frame_regex != other.frame_regex:
            return False
         elif self.fps != other.fps:
            return False
         else:
            return True
      except AttributeError:
         return False

   #===========================================================================
   def __ne__(self, other):
      """
      Test inequality.
      """
      return not self == other

# End class EncodeInfo
#==============================================================================
#------------------------------------------------------------------------------



#==============================================================================

def draw_frames_loop(data_list, movie_list, encode_detail):
   """
   Makes all frames based on lists of data sources and of movies.

   Arguments:
      data_list (list of strings) : list of filenames of DUMSES outputs
      movie_list (list of MovieDescriptors) : list of movies to make
      encode_detail (set) : set of EncodeInfo (filled by this function)

   Returns:
      None

   Exceptions:
      whatever arises from functions called by this routine
   """

   # Set up checks - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   # Ensure all movies have unique stubs
   list_of_pathstubs = []  # pathstubs we've seen
   blacklist = set()       # pathstubs to remove

   for movie in movie_list.values():
      if movie.pathstub in list_of_pathstubs:
         if movie.pathstub not in blacklist: # Only warn once per repeated stub
            msg = '"'.join(("Multiple movies with path/stub ", movie.pathstub,
               ".  These will be excluded."))
            warnings.warn(msg, UserWarning)
         blacklist.add(movie.pathstub)
      else:
         list_of_pathstubs.append(movie.pathstub)

   movie_list = {m : movie_list[m] for m in movie_list
         if movie_list[m].pathstub not in blacklist}

   # Split the loops because these steps apply only to movies we actually
   # intend to generate
   save_initial_state = False
   for movie in movie_list.values():

      # Do any of the surviving movies need the initial state?
      if movie.mode.dimension == 1:
         save_initial_state = True

      # If any of the movies have absolute paths (allowed in general for the
      # MovieDescriptor, because it is not tied to this control loop) warn the
      # user that multiple data series may overwrite or interleave.
      if movie.path is not None and os.path.isabs(movie.path):
         msg = " ".join(("Absolute paths may cause overwriting and/or",
            "interleaving if you have multiple data series."))
         warnings.warn(msg, UserWarning)

   # Sorting makes it more likely that you'll hit the initial state before the
   # other states in the same series, and less likely that you'll switch
   # between series and then back, thus potentially saving the need to dump and
   # reload different initial states all the time in a poorly-ordered list.
   data_list.sort()

   state0 = None
   series_name = None

   # Master loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   # Loop over each output file first, because loading the files is a
   # significant bottleneck.
   for output_name in data_list:

      # Parse the output_name
      # -- DumsesData expects the path (ending in a slash) and the number
      # -- We will assume that a series consists of files stored in the same
      #    directory, and we assume that a series shares a common inputs file
      #    (also in the same directory).
      path = os.path.dirname(os.path.realpath(output_name)) + "/"
      filename = os.path.basename(output_name)
      try:
         number = int(filename.rpartition("_")[2])
      except ValueError:
         msg = "".join(('Unable to parse file name "', output_name,
            '".  This file will be skipped.'))
         warnings.warn(msg, UserWarning)
         continue

      # Update the inputs file if necessary
      if series_name != path:
         try:
            si = SimulationInput(path + "input")
         except IOError:
            msg = "".join((
               'Unable to open inputs file corresponding to data file "',
               output_name, '".  This file will be skipped.'))
            warnings.warn(msg, UserWarning)
            continue

      # Open the DUMSES output file and repackage into a SimulationState
      try:
         data = DumsesData(number, filedir=path)
      except IOError:
         msg = "".join(('Unable to open data file "', output_name,
            '".  This file will be skipped.'))
         warnings.warn(msg, UserWarning)
         continue
      try:
         state = SimulationState(data, si)
      except SimulationError:
         msg = "".join(('Unable to construct state from data file "',
            output_name, '".  This file will be skipped.'))
         warnings.warn(msg, UserWarning)
         continue

      # Update the (t == 0) entry if necessary
      if save_initial_state and series_name != path:
         if state.t == 0:
            # If we just opened the t == 0 file, then everything is easy
            state0 = state
            series_name = path
         else:
            # We need to open the t == 0 file for this series
            try:
               data0 = DumsesData(0, filedir=path)
               state0 = SimulationState(data0, si)
               series_name = path
            except IOError:
               # Problem opening file; skip all the t == 0 stuff when plotting
               msg = "Unable to open initial state for series " + path + \
                     "; no initial profiles will be plotted."
               warnings.warn(msg, UserWarning)
               state0 = None
               series_name = None

      # Loop over all the movies
      for movie in movie_list.values():
         # Only make frames in the time range
         if movie.time_limits.test(state.t):
            sys.stdout.write("Making frame {s}_{n:06d} (t = {t}).\n".format(
                  t=state.t, s=movie.stub, n=number))
            try:
               movie.draw_frame(state, path, number, state0)
               # Save the movie and frame information for encoding the frame
               # images into a single movie.
               if movie.make_movie:
                  movie_name, frame_regex = movie.build_paths(path)
                  movie_name += ".".join((movie.stub, movie.movie_type))
                  frame_regex += "*.".join((movie.stub, movie.image_type))
                  encode_detail.add(
                        EncodeInfo(movie_name, frame_regex, movie.fps))
            except (DescriptorError, SimulationError) as err:
               # Well that didn't work... guess we'll move along
               msg = "".join(("While generating movie ", str(movie),
                  " with data file ", output_name,
                  ", the following error occurred:\n   ", str(err),
                  "\nThis frame will not be drawn."))
               warnings.warn(msg, UserWarning)
         else:
            sys.stdout.write("Skipping frame {s}_{n:06d} (t = {t}).\n".format(
                  t=state.t, s=movie.stub, n=number))

#==============================================================================

def encode_movies_loop(encode_list):
   """
   Encodes all movies specified in the supplied list.

   Arguments:
      encode_list (list of EncodeInfo) : list of movies to encode

   Returns:
      None

   Exceptions:
      whatever arises from functions called by this routine
   """

   # Encode the movies
   for ei in encode_list:
      # Log mencoder's output into a related logfile
      split = os.path.splitext(ei.movie_name)
      if len(split[1]) > 1:
         split1 = "_" + split[1][1:]
      else:
         split1 = ""
      log_name = "".join((split[0], split1, ".log"))
      with open(log_name, 'w') as log_file:
         # Prepare the mencoder command
         command = ['mencoder', "mf://"+ei.frame_regex,
                    '-mf', ":".join(("type=png", "fps={0}".format(ei.fps))),
                    '-ovc', 'lavc',
                    '-lavcopts', 'vcodec=mpeg4',
                    '-oac', 'copy',
                    '-o', ei.movie_name]
         sys.stdout.write("Making movie {0}.\n".format(ei.movie_name))
         try:
            subprocess.call(command, stdout=log_file, stderr=log_file)
         except OSError as e:
            if e.errno == os.errno.ENOENT:
               # Can't find mencoder; warn the user and kill the loop
               msg = "Cannot locate mencoder to encode frames into movies."
               warnings.warn(msg, UserWarning)
               break
            else:
               raise

#==============================================================================

