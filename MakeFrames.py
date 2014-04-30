"""
Make the frames.

This module manages the loops and control logic for making the frames.  It is
generally entered through the make_all_frames routine, which takes a list of
DUMSES output file names and a list of MovieDescriptors, then loops as
appropriate to make the frames.  It can also be entered one level farther down,
to the make_single_frame routine (which is called by the make_all_frames
loops), which takes the name of a single DUMSES output and a single
MovieDescriptor, then draws the appropriate frame.  The make_single_frame
routine uses several support routines to draw the panels for each image.  The
panel_2D_pseudocolor routine takes a state, movie description, and axes, and
draws a two-dimensional psuedocolor plot  on the axes using the data in the
state and the description of the movie.  The panel_profile function is
analogous, but draws a one-dimensional profile (along the x axis due to its
special nature).

Attributes:
   make_all_frames : makes all frames for lists of outputs and movies
   make_single_frame : makes a single frame for a given output and movie
   panel_2D_pseudocolor : draws a 2D pseudocolor plot on a given axes
   panel_panel : draws a 1D profile plot on a given axes
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
import Descriptors as Desc
import SimulationData as SD

#==============================================================================

def make_all_frames(data_list, movie_list):
   """
   Makes all the frames based on a list of data sources and a list of movies.

   Arguments:
      data_list (list of strings) : list of filenames of DUMSES outputs
      movie_list (list of MovieDescriptors) : list of movies to make

   Returns:
      None

   Exceptions:
      whatever arises from functions called by this routine
   """

   # Set up checks - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

   # Ensure all movies have unique stubs
   list_of_stubs = []   # stubs we've seen
   s_blacklist = set()  # stubs to remove

   for movie in movie_list.values():
      if movie.stub in list_of_stubs:
         if movie.stub not in s_blacklist: # Only warn once per repeated stub
            msg = "".join(("Multiple movies with stub ", movie.stub,
               ".  These will be excluded."))
            warnings.warn(msg, UserWarning)
         s_blacklist.add(movie.stub)
      else:
         list_of_stubs.append(movie.stub)

   movie_list = {m : movie_list[m] for m in movie_list
         if movie_list[m].stub not in s_blacklist}

   # Split the loops because these steps apply only to movies we actually
   # intend to generate
   profile_mode = False
   for movie in movie_list.values():

      # Do any of the surviving movies need the initial state?
      if movie.mode.dimension == 1:
         profile_mode = True

      # If any of the movies have absolute paths (allowed in general for the
      # MovieDescriptor, because it is not tied to this set of drawing
      # routines), warn the user that multiple data series may overwrite or
      # interleave.
      if os.path.isabs(movie.path):
         msg = " ".join(("Absolute paths may cause overwriting and/or",
            "interleaving if you have multiple data series."))
         warnings.warn(msg, UserWarning)

   # Sorting makes it more likely that you'll hit the initial state before the
   # other states in the same series, thus potentially saving the need to dump
   # and reload different initial states all the time in a poorly-ordered list.
   data_list.sort()

   # If needed, save the initial state for the current series
   state0 = None
   if profile_mode:
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
            si = SD.SimulationInput(path + "input")
         except IOError:
            msg = "".join((
               'Unable to open inputs file corresponding to data file "',
               output_name, '".  This file will be skipped.'))
            warnings.warn(msg, UserWarning)
            continue

      # Open the DUMSES output file
      try:
         data = DumsesData(number, filedir=path)
      except IOError:
         msg = "".join(('Unable to open data file "', output_name,
            '".  This file will be skipped.'))
         warnings.warn(msg, UserWarning)
         continue

      # Repackage into a SimulationState object
      state = SD.SimulationState(data, si)

      # Update the (t == 0) entry if necessary
      if profile_mode and series_name != path:
         if state.t == 0:
            # If we just opened the t == 0 file, then everything is easy
            state0 = state
            series_name = path
         else:
            # We need to open the t == 0 file for this series
            try:
               data0 = DumsesData(0, filedir=path)
               state0 = SD.SimulationState(data0, si)
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
         if movie.within_time_limits(state.t):
            try:
               make_single_frame(state, movie, path, number, state0)
            except (Desc.DescriptorError, SD.SimulationError) as err:
               # Well that didn't work... guess we'll move on
               msg = "".join(("While generating movie ", str(movie),
                  " with data file ", output_name,
                  ", the following error occurred:\n   ", str(err),
                  "\nThis frame will not be drawn."))
               warnings.warn(msg, UserWarning)

#==============================================================================

def make_single_frame(state, movie, directory, number, state0=None):
   """
   Draw a single frame of the movie

   In order to be generic for future development, this will be a fairly
   lightweight wrapper.  It creates a new figure, creates a subplot that fills
   the figure, and passes the axes to a routine that does 2D plotting.  Future
   development will allow multiple subplots (for example, a side-by-side
   comparison of vorticity and entropy, with a trajectory of the mean value of
   Foglizzo's chi value), which would require extracting multiple subplots and
   calling the appropriate drawing routines.

   Note on file locations: If the MovieDescriptor has an absolute path, the
   movie will be stored there.  If the MovieDescriptor has a relative path, it
   will be assumed to be relative to the directory where the simulation data
   file is stored.  If the MovieDescriptor has no path (i.e., is None), the
   movie will be stored in the same directory where the simulation data is
   found.  The frame images will be stored in a subdirectory called "frames"
   relative to the directory where the movie will be stored.

   Arguments:
      state (SimulationState) : the current state of the simulation
      movie (MovieDescriptor) : the description of the movie to be made
      directory (string) : the directory where the simulation data is stored
      number (int) : the number of the current data file
      state0 (SimulationState) : the initial state of the simulation (or None)

   Returns:
      None

   Exceptions:
      whatever arises from functions called by this routine
   """

   print "Making frame {s}{n:06d} (t = {t}).".format(t=state.t, s=movie.stub,
         n=number)

   # Generate and partition the figures
   plt.ioff()
   fig = plt.figure()
   ax1 = fig.add_axes([0, 0.05, 1, 0.9]) # left, bottom, width, height

   # Fill the panels (Axes)
   if movie.mode.dimension == 1:
      panel_profile(state, movie, ax1, state0)
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
   path += "frames/"             # Subdirectory "frames"
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
   z-index of zero.  This comes from the structure of data in DUMSES for a 2D
   simulation, and can be modified later if I generate 3D data to visualize.

   Arguments:
      state (SimulationState) : the current state of the simulation
      movie (MovieDescriptor) : the description of the movie to be made
      axes (matplotlib.axes.Axes) : the axes to draw the pseudocolor plot on

   Returns:
      None

   Exceptions:
      whatever arises from functions called by this routine
   """

   # Get the data to be drawn
   data, x, y, z, vlo, vhi = movie.frame_data(state)

   # Since we are drawing in 2D, reduce to the assumed geometry
   # --> Also, transform to the orientation expected by imshow
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
   axes.set_title(" ".join((movie.variable, str(movie.mode))))

   # Compute ticks (x, y, color)
   ntickx = 6
   nticky = int(round(ntickx * (yy.max() - yy.min()) / (xx.max() - xx.min())))
   axes.xaxis.set_major_locator(MaxNLocator(ntickx))
   axes.yaxis.set_major_locator(MaxNLocator(nticky))
   cbar.locator = MaxNLocator(nbins=10)

#==============================================================================

def panel_profile(state, movie, ax_prof, state0=None):
   """
   Draw the profiles as one panel of a figure.

   This routine assumes that we are plotting a profile along the x-axis, and
   computes the mean and standard deviation along that axis.

   The figure is split into two panels:
      the upper panel:
         the mean state +/- the standard deviation
         the unperturbed base state +/- the initial standard deviation
         the initial maximum and minimum states
      the lower panel:
         the standard deviation
         the initial standard deviation
   The profiles derived from the initial state will not be plotted if None is
   supplied for state0.

   Arguments:
      state (SimulationState) : the current state of the simulation
      movie (MovieDescriptor) : the description of the movie to be made
      axes (matplotlib.axes.Axes) : the axes to draw the profile plot on
      state0 (SimulationState) : the initial state of the simulation

   Returns:
      None

   Exceptions:
      whatever arises from functions called by this routine
   """

   pad = 0.075

   # Get the data to be drawn
   data, x, y, z, vlo, vhi = movie.frame_data(state, pad_bounds=pad)
   view = data.reshape((data.shape[0], data.shape[1]*data.shape[2]))
   mean = np.mean(view, axis=1)
   stdv = np.std(view, axis=1)

   if state0 is not None:
      init, x2, y2, z2, vlo2, vhi2 = movie.frame_data(state0)
      view = init.reshape((init.shape[0], init.shape[1]*init.shape[2]))
      std0 = np.std(view, axis=1)
      min0 = np.amin(view, axis=1)
      max0 = np.amax(view, axis=1)

   if movie.mode.transform == "none":
      movie_copy = copy.deepcopy(movie)
      movie_copy.mode.absolute = False
      movie_copy.mode.transform = "none"
      movie_copy.mode.reference = "base"
      data, x2, y2, z2, vlo2, vhi2 = movie_copy.frame_data(state)
      base = data[:,0,0]
   else:
      base = np.zeros_like(mean)

   x = x[:,0,0]

   # Split the axes in half
   divider = make_axes_locatable(ax_prof)
   ax_stdv = divider.append_axes("bottom", size="40%", pad=0.1)

   # TODO : If the data is a vertical or horizontal line (all x- or
   #        y-coordinates are the same), then a warning appears.  Look into
   #        whether or not that can be silenced elegantly without hiding
   #        anything important.

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
   ax_prof.set_title(" ".join((movie.variable, str(movie.mode))))

   # Label axes and set ticks
   plt.setp(ax_prof.get_xticklabels(), visible=False)
   ax_prof.set_ylabel("profile")
   ax_stdv.set_xlabel("position")
   ax_stdv.set_ylabel("standard deviation")

   # Set the plotting limits
   bounds_max = stdv.max()
   if state0 is not None:
      bounds_max = max(bounds_max, std0.max())
   if bounds_max == 0.0:
      bounds_max = 1.0
   else:
      bounds_max = (1.0 + pad) * bounds_max
   ax_stdv.set_ylim([0.0, bounds_max])
   ax_prof.set_ylim([vlo, vhi])

