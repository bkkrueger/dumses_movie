"""
Make the frames.

This module manages the loops and control logic for making the frames.  It is
generally entered through the make_all_frames routine, which takes a list of
DUMSES output file names and a list of MovieDescriptors, then loops as
appropriate to make the frames.  It can also be entered one level farther down,
to the make_single_frame routine, which is called by the make_all_frames loops,
and takes the name of a single DUMSES output and a single MovieDescriptor, then
draws the appropriate frame.

Attributes:
   make_all_frames : makes all frames for lists of outputs and movies
   make_single_frame : makes a single frame for a given output and movie
   panel_2D_pseudocolor : draws a 2D pseudocolor plot on a given axes
"""

# TODO : Idea: Allow specification of colorbar limits as function of time.
#        Currently I allow a static colorbar minimum/maximum, but can this be
#        extended in some way to be time-dependent?  In theory, yes: a list of
#        triplets giving a time, a minimum, and a maximum (where min and max
#        can be None).  The simplest implementation would then be to change the
#        colorbar limits in a piecewise constant manner at those times, but it
#        would look better to do a piecewise linear interpolation.  Higher
#        order may or may not add any benefit.  I will have to think about: (a)
#        is such a feature useful; (b) is there a better way to implement it?

import copy
import errno
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os
import warnings

from dumpy_v05.data.rd_dumses import DumsesData
import SimulationData as SD

#==============================================================================

def make_all_frames(data_list, movie_list):
   """
   Makes all the frames based on a list of data sources and a list of movies.

   Arguments:
      data_list (list of strings) : list of DUMSES outputs
      movie_list (list of MovieDescriptors) : list of movies to make
   """

   # Basic sanity checks - - - - - - - - - - - - - - - - - - - - - - - - - - -

   list_of_stubs = []
   s_blacklist = set()
   v_blacklist = set()
   m_blacklist = set()
   for movie in movie_list.values():

      # Check #1 : All movies have unique stubs
      if movie.stub in list_of_stubs:
         if movie.stub not in s_blacklist: # Only warn once per repeated stub
            msg = "".join(("Multiple movies with stub ", movie.stub,
               ".  These will be excluded."))
            warnings.warn(msg, UserWarning)
         s_blacklist.add(movie.stub)
      else:
         list_of_stubs.append(movie.stub)

      # Check #2 : All movies are requesting valid variables
      if movie.variable not in SD.SimulationState.known_variables:
         msg = "".join(('Invalid variable "', movie.variable,
            '" in movie with stub "', movie.stub,
            '".  This will be excluded.'))
         warnings.warn(msg, UserWarning)
         v_blacklist.add(movie.variable)

      # Check #3 : All movies are requesting valid modes
      # Note : We allow any modes known to SimulationState, but also the
      #        "profile" mode, which will be constructed here.
      profile_mode = False
      if movie.mode not in SD.SimulationState.known_modes:
         if movie.mode in ["profile: full state", "profile: perturbation",
               "profile: contrast"]:
            profile_mode = True
         else:
            msg = "".join(('Invalid mode"', movie.mode,
               '" in movie with stub "', movie.stub,
               '".  This will be excluded.'))
            warnings.warn(msg, UserWarning)
            m_blacklist.add(movie.mode)

   # Remove any disallowed MovieDescriptors
   movie_list = {m : movie_list[m] for m in movie_list
         if movie_list[m].stub not in s_blacklist}
   movie_list = {m : movie_list[m] for m in movie_list
         if movie_list[m].variable not in v_blacklist}
   movie_list = {m : movie_list[m] for m in movie_list
         if movie_list[m].mode not in m_blacklist}

   # Sorting makes it more likely that you'll hit the initial state before the
   # other states in the same series.
   data_list.sort()

   # If needed, save the initial state for the current series
   state0 = None
   if profile_mode:
      series_name = None

   # Master loop - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   # -- Loop over each output file first, because we only want to load each
   #    file once.
   for output_name in data_list:

      # Parse the output_name
      # -- The name is in the format "/subpath/subpath/output_######", with
      #    a possible trailing "/".
      # -- DumsesData expects the path (all the "subpath" bits, ending in a
      #    "/") and the number (the ###### as an integer).
      # -- For genericness, we will assume that not all output files in the
      #    list share a common inputs file.
      # -- For simplicity, We will assume that the inputs file is in the same
      #    directory as the output.
      if output_name[-1] == "/":
         output_name = output_name[:-1]
      partition = output_name.rpartition('/')
      path = os.path.realpath(''.join(partition[0:2]))
      if path[-1] != "/":
         path += "/"
      number = int(partition[2].rpartition("_")[2])

      # Open the input and output files to construct the state
      data = DumsesData(number, filedir=path)
      si = SD.SimulationInput(path + "input")
      state = SD.SimulationState(data, si)

      # Update the (t == 0) entry if necessary
      if profile_mode:
         if series_name != path:
            if state.t == 0:
               state0 = state
               series_name = path
            else:
               try:
                  data0 = DumsesData(0, filedir=path)
                  state0 = SD.SimulationState(data0, si)
                  series_name = path
               except:
                  msg = "Unable to open initial state for series " + path + \
                        "; no initial profiles will be plotted."
                  warnings.warn(msg, UserWarning)
                  state0 = None
                  series_name = None

      for movie in movie_list.values():
         if movie.within_time_limits(state.t):
            make_single_frame(state, state0, movie, path, number)

#==============================================================================

def make_single_frame(state, state0, movie, directory, number):
   """
   Make a single frame based on the current state and a movie description.

   In order to be generic for future development, this will be a fairly
   lightweight wrapper.  It creates a new figure, creates a subplot that fills
   the figure, and passes the axes to a routine that does 2D plotting.  Future
   development will allow multiple subplots (for example, a side-by-side
   comparison of vorticity and entropy, with a trajectory of the mean value of
   Foglizzo's chi value), which would require extracting multiple subplots and
   calling the appropriate drawing routines.

   Note on file locations: If the MovieDescriptor has an absolute path, the
   image will be stored there.  If the MovieDescriptor has a relative path, it
   will be assumed to be relative to the directory where the simulation data
   file is stored.  If the MovieDescriptor has no path (i.e., is None), the
   image will be stored in the same directory where the simulation data is
   found.

   Arguments:
      state (SimulationState) : the current state of the simulation
      movie (MovieDescriptor) : the description of the movie to be made
      directory (string) : the directory where the simulation data is stored
   """

   print "Making frame {s}{n:06d} (t = {t}).".format(t=state.t, s=movie.stub,
         n=number)

   # Generate and partition the figures
   plt.ioff()
   fig = plt.figure()
   ax1 = fig.add_axes([0, 0.05, 1, 0.9])

   # Fill the panels (Axes)
   if movie.mode in ["profile: full state", "profile: perturbation",
         "profile: contrast"]:
      panel_profile(state, state0, movie, ax1)
   else:
      panel_2D_pseudocolor(state, movie, ax1)

   # Verify that the directory exists
   if movie.path is None:        # Default to directory
      path = directory
   elif movie.path[0] == "/":    # Absolute path specified
      path = movie.path
   else:                         # Relative path from directory
      if directory[-1] != "/":
         directory += "/"
      path = directory + movie.path
   if path[-1] != "/":           # Make sure there's a trailing "/"
      path += "/"
   try:
      os.makedirs(path)
   except OSError as e:
      if e.errno != errno.EEXIST:
         raise

   # Save the figures
   image_file_name = "{p}{s}{n:06d}.{e}".format(p=path, s=movie.stub, n=number,
         e=movie.image_type)
   fig.savefig(image_file_name, bbox_inches="tight")
   fig.clear()

#==============================================================================

def panel_2D_pseudocolor(state, movie, axes):
   """
   Draw a 2D pseudocolor plot as one panel of a figure.

   This routine assumes that we are plotting a slice in the xy-plane at a
   z-index of zero.  This can be modified later if I need to be worrying about
   3D data.

   Arguments:
      state (SimulationState) : the current state of the simulation
      movie (MovieDescriptor) : the description of the movie to be made
      axes (matplotlib.axes.Axes) : the axes to draw the pseudocolor plot on
   """

   # Get the data to be drawn
   data, x, y, z, vlo, vhi = movie.frame_data(state)

   # Since we are drawing in 2D, reduce to the assumed geometry
   # Also, transform to the orientation expected by imshow
   plot_data = np.flipud(data[:,:,0])
   xx = y[0,:,0]
   yy = x[:,0,0]
   y_grid = np.abs(x[:,:,0].repeat(y.shape[1],1))

   # Construct the desired colormap
   my_cmap = copy.copy(cm.seismic)
   my_cmap.set_bad(color="k", alpha=1)

   # Split the axes to have space for the colormap
   divider = make_axes_locatable(axes)
   cbar_axes = divider.append_axes("right", size="20%", pad=0.1)

   # Draw the plot
   image = axes.imshow(plot_data,
         extent=[xx.min(), xx.max(), yy.min(), yy.max()],
         aspect=1.0, cmap=my_cmap, vmin=vlo, vmax=vhi)
   axes.contour(xx, yy, y_grid, levels=[state.params.layer_width],
         colors="#00FF00")
   cbar = plt.colorbar(image, cax=cbar_axes)

   # Construct the title
   axes.set_title(" ".join((movie.variable, movie.mode)))

   # Compute ticks (x, y, color)
   ntickx = 6
   nticky = int(round(ntickx * (yy.max() - yy.min()) / (xx.max() - xx.min())))
   axes.xaxis.set_major_locator(MaxNLocator(ntickx))
   axes.yaxis.set_major_locator(MaxNLocator(nticky))
   cbar.locator = MaxNLocator(nbins=10)

#==============================================================================

def panel_profile(state, state0, movie, ax_prof):
   """
   Draw the profiles as one panel of a figure.

   This routine assumes that we are plotting a profile along the x-axis, and
   computes the mean and standard deviation along that axis.

   Arguments:
      state (SimulationState) : the current state of the simulation
      state0 (SimulationState) : the initial state of the simulation
      movie (MovieDescriptor) : the description of the movie to be made
      axes (matplotlib.axes.Axes) : the axes to draw the profile plot on
   """

   # Get the data to be drawn
   movie_copy = copy.copy(movie)

   movie_copy.mode = "base state"
   base, x, y, z, vlo, vhi = movie_copy.frame_data(state)
   base = base[:,0,0]

   movie_copy.mode = "full state"
   if state0 is not None:
      init, x, y, z, vlo, vhi = movie_copy.frame_data(state0)
      view = init.reshape((init.shape[0], init.shape[1]*init.shape[2]))
      std0 = np.std(view, axis=1)
      min0 = np.amin(view, axis=1)
      max0 = np.amax(view, axis=1)

   data, x, y, z, vlo, vhi = movie_copy.frame_data(state)
   view = data.reshape((data.shape[0], data.shape[1]*data.shape[2]))
   mean = np.mean(view, axis=1)
   stdv = np.std(view, axis=1)

   x = x[:,0,0]

   # Split the axes in half
   divider = make_axes_locatable(ax_prof)
   ax_stdv = divider.append_axes("bottom", size="40%", pad=0.1)

   # Draw the top plot
   if state0 is not None:
      ax_prof.plot(x, min0, lw=2, linestyle="dashed", color='grey')
      ax_prof.plot(x, max0, lw=2, linestyle="dashed", color='grey')
      ax_prof.fill_between(x, base-std0, base+std0,
            facecolor='grey', alpha=0.5)
   ax_prof.plot(x, base, lw=2, label="steady state", color='grey')
   ax_prof.fill_between(x, mean-stdv, mean+stdv, facecolor='blue', alpha=0.5)
   ax_prof.plot(x, mean, lw=2, label="mean profile", color='blue')

   # Draw the bottom plot
   if state0 is not None:
      ax_stdv.fill_between(x, std0, facecolor='grey', alpha=0.5)
      ax_stdv.plot(x, std0, lw=2, color='grey')
   ax_stdv.fill_between(x, stdv, facecolor='blue', alpha=0.5)
   ax_stdv.plot(x, stdv, lw=2, label="standard deviation", color='blue')

   # Draw lines marking the source layer
   ax_prof.axvline(x=-state.params.layer_width,color='black')
   ax_prof.axvline(x= state.params.layer_width,color='black')

   # Construct the title
   ax_prof.set_title(" ".join((movie.variable, movie.mode)))

   # Label axes and set ticks
   plt.setp(ax_prof.get_xticklabels(), visible=False)
   ax_prof.set_ylim((vlo, vhi))
   ax_prof.set_ylabel("profile")
   ax_stdv.set_xlabel("position")
   ax_stdv.set_ylabel("standard deviation")

   # Set the plotting limits
   bounds_padding = 0.075

   bounds_type = state.bounds_type(movie.variable, "full state")
   if bounds_type == "non-negative":
      if movie.value_limits.hi is None:
         bounds_max = max(base.max(), mean.max(), (mean+stdv).max())
         if state0 is not None:
            bounds_max = max(bounds_max, max0.max(), (base+std0).max())
         if bounds_max == 0.0:
            bounds_max = 1.0
         else:
            bounds_max = (1.0 + bounds_padding) * bounds_max
      else:
         bounds_max = movie.value_limits.hi
      if movie.value_limits.lo is None:
         bounds_min = 0.0
      else:
         bounds_min = movie.value_limits.lo

   elif bounds_type == "symmetric":
      if movie.value_limits.lo is None and movie.value_limits.hi is None:
         bounds = max(np.abs(base).max(), np.abs(mean).max(),
               np.abs(mean-stdv).max())
         bounds = max(bounds, np.abs(base).max(), np.abs(mean).max(),
               np.abs(mean+stdv).max())
         if state0 is not None:
            bounds = max(bounds, np.abs(min0).max(), np.abs(base-std0).max())
            bounds = max(bounds, np.abs(max0).max(), np.abs(base+std0).max())
         if bounds == 0.0:
            bounds_max = 1.0
         else:
            bounds_max = (1.0 + bounds_padding) * bounds
         bounds_min = -bounds_max
      elif movie.value_limits.lo is not None:
         bounds_max = abs(movie.value_limits.lo)
         bounds_min = -bounds_max
      elif movie.value_limits.hi is not None:
         bounds_max = abs(movie.value_limits.hi)
         bounds_min = -bounds_max
      else:
         bounds_max = movie.value_limits.hi
         bounds_min = movie.value_limits.lo

   elif bounds_type == "general":
      bounds_min = min(base.min(), mean.min(), (mean-stdv).min())
      if state0 is not None:
         bounds_min = min(bounds_min, min0.min(), (base-std0).min())
      bounds_max = max(base.max(), mean.max(), (mean+stdv).max())
      if state0 is not None:
         bounds_max = max(bounds_max, max0.max(), (base+std0).max())
      span = bounds_max - bounds_min
      if span == 0.0:
         bounds_max = 1.0
         bounds_min = -1.0
      else:
         bounds_min -= bounds_padding * span
         bounds_max += bounds_padding * span
      if movie.value_limits.lo is not None:
         bounds_min = movie.value_limits.lo
      if movie.value_limits.hi is not None:
         bounds_max = movie.value_limits.hi
   else:
      # TODO
      pass
   ax_prof.set_ylim([bounds_min, bounds_max])

   bounds_max = stdv.max()
   if state0 is not None:
      bounds_max = max(bounds_max, std0.max())
   if bounds_max == 0.0:
      bounds_max = 1.0
   else:
      bounds_max = (1.0 + bounds_padding) * bounds_max
   ax_stdv.set_ylim([0.0, bounds_max])

