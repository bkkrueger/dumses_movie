"""
Make the frames.

This module manages the loops and control logic for making the frames.

Attributes:
   draw_frames_loop : draws all frames for lists of outputs and movies
   encode_movies_loop : encodes all movies for list of frame sources
"""

import copy
import errno
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os
import sys
import warnings

from dumpy_v05.data.rd_dumses import DumsesData
from Descriptors import DescriptorError
from SimulationData import SimulationError, SimulationInput, SimulationState

#==============================================================================

def draw_frames_loop(data_list, movie_list, encode_locations):
   """
   Makes all frames based on lists of data sources and of movies.

   Arguments:
      data_list (list of strings) : list of filenames of DUMSES outputs
      movie_list (list of MovieDescriptors) : list of movies to make
      encode_locations (set) : set of movies to make (filled by this function)

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
      # MovieDescriptor, because it is not tied to this control loop, warn the
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

   # Loop over each output file first, because loading the files seems to be
   # one of the bottlenecks.
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
                  encode_locations.add((movie_name, frame_regex, movie.fps))
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
      encode_list : list of movies to encode

   Returns:
      None

   Exceptions:
      whatever arises from functions called by this routine
   """

   # Encode the movies
   if ProcID == 0:
      print "="*79
   for movie_name, frame_regex, fps in encode_list:
      log_name = os.path.splitext(movie_name)[0] + ".log"
      with open(log_name, 'w') as log_file:
         command = ['mencoder', "mf://"+frame_regex,
                    '-mf', ":".join(("type=png", "fps={0}".format(fps))),
                    '-ovc', 'lavc',
                    '-lavcopts', 'vcodec=mpeg4',
                    '-oac', 'copy',
                    '-o', movie_name]
         sys.stdout.write("Making movie {0}.\n".format(movie_name))
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

