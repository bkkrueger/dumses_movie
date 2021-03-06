"""
Specification of movie descriptors.

The MovieDescriptor object describes the details necessary to construct
a movie from data.  It describes the design of the images for each
individual frame, and includes information relating to how those images
will be encoded together into a movie.

The MaskDescriptor object describes the details of a mask.  This is a
component used by MovieDescriptors to describe how to mask cells that
are not "interesting enough" to be plotted.

The ModeDescriptor object describes the details of how to process a
variable for a movie.  It allows specification of profiles or 2D
pseudocolor plots, with a selection of transformations that may be
applied.

The MovieLimits object describes the lower and upper limits on some
quantity.  The MovieWindow class packages three of these together to
describe limits in 3D coordinates.  These two are essentially for
packaging more than functionality.

The DescriptorError class describes an error raised by one of the
Descriptors.  The errors raised by them tend to be related to illegal
values (or illegal combinations of values), so they can be logically
grouped together.  Giving them their own class means that calling
routines can catch errors raised by the Descriptors separately from
other errors (e.g. a ValueError in something used by one of the
Descriptors, etc).
"""

import copy
import errno
import math
import matplotlib.cm as cm
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.ticker import MaxNLocator, FuncFormatter
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np
import os
import warnings

#------------------------------------------------------------------------------
#==============================================================================

class DescriptorError(StandardError):
   """
   Signals that a descriptor encountered an error.
   """

   def __init__(self, *args, **kwargs):
      """
      Initialize the exception.

      Because StandardError explicitly states that it does not use any
      keyword arguments, that gives the freedom to use any keyword
      arguments to customize the initialization of DescriptorError.
      Based on the arguments supplied, the initialization selects
      between the following modes:
      -- invalid=True:
         A very common message: 'Invalid <parameter>: "<value>".',
         where <parameter> and <value> are given as the first and
         second positional arguments respectively.  The name of the
         parameter should be a string, while the value may be any
         quantity that can be cast into a string.
      -- join=True:
         Sometimes the message is the concatenation of a collection of
         substrings.  In this case all positional arguments are cast to
         strings and concatenated in order.
      -- message:
         If there is only one positional argument and no special
         keyword arguments, the positional argument is saved as the
         message.
      -- other:
         In order to mimic the (apparent) behavior of StandardError, if
         multiple positional arguments are supplied without any special
         keyword arguments, then no message is stored.
      """

      # Save the arguments
      self.args = args
      self.kwargs = kwargs

      # Check that no more than one mutually-exclusive keywords are given
      exclusive = ["invalid", "join"]
      sum_exclusive = sum(1 for kw in exclusive if kwargs.get(kw, False))
      if sum_exclusive > 1:
         msg = "".join(("DescriptorError does not accept more than one of [",
               ", ".join(exclusive), "]."))
         raise TypeError(msg)

      # Construct the message
      if kwargs.get("invalid", False):
         self.message = 'Invalid {p}: "{v}".'.format(p=args[0], v=str(args[1]))
      elif kwargs.get("join", False):
         self.message = "".join(str(a) for a in args)
      elif len(args) == 1:
         self.message = args[0]
      else:
         self.message = ""

# End of DescriptorError class
#==============================================================================
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#==============================================================================

class MovieLimits(object):
   """
   Pair of low/high limits on some quantity.
   """

   def __init__(self):
      """
      Initialize an empty pair of limits.
      """
      self.lo = None
      self.hi = None

   def __repr__(self):
      """
      String representation.
      """
      if self.lo is None:
         lo = "*"
      else:
         lo = self.lo
      if self.hi is None:
         hi = "*"
      else:
         hi = self.hi
      return "Limits: [{l},{h}]".format(l=lo, h=hi)

   def __str__(self):
      """
      String representation.
      """
      if self.lo is None:
         lo = "*"
      else:
         lo = self.lo
      if self.hi is None:
         hi = "*"
      else:
         hi = self.hi
      return "[{l},{h}]".format(l=lo, h=hi)

   def valid(self):
      """
      Verify that the lower limit is below the upper limit.
      """
      if None not in [self.lo, self.hi]:
         return self.lo < self.hi
      else:
         return True

   def test(self, a):
      """
      Test whether the supplied argument is in the limits.

      Note: This may be called with a scalar or a NumPy array.  Using
      NumPy's ones_like and logical_and handles both cases.
      """
      if self.lo is None:
         lo = np.ones_like(a, dtype=bool)
      else:
         lo = (self.lo <= a)
      if self.hi is None:
         hi = np.ones_like(a, dtype=bool)
      else:
         hi = (a <= self.hi)
      return np.logical_and(lo, hi)

# End of MovieLimits class
#==============================================================================
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#==============================================================================

class MovieWindow(object):
   """
   Viewing limits on a plot, to restrict the image to a certain region.
   """

   def __init__(self):
      """
      Initialize an empty window."
      """
      self.x = MovieLimits()
      self.y = MovieLimits()
      self.z = MovieLimits()

   def __repr__(self):
      """
      String representation.
      """
      return "   \n".join(("Movie window:",
         repr(self.x), repr(self.y), repr(self.z)))

   def __str__(self):
      """
      String representation.
      """
      return " ".join((str(self.x), str(self.y), str(self.z)))

# End of MovieWindow class
#==============================================================================
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#==============================================================================

class MaskDescriptor(object):
   """
   A description of a mask to be applied to the data prior to plotting.

   This class describes a data mask.  It specifies a variable, mode,
   threshold, and operator.  As an example to explain, if the variable
   is entropy, the mode is perturbation, the threshold is zero, and the
   operator is less than, the applied mask will hide all cells where
   the entropy perturbation is less than zero.  The operators are
   supplied by a string and the permitted values are "<=", "<", ">",
   and ">=".

   Attributes:
      _variable (string) : the variable used for computing the mask
      _mode (ModeDescriptor) : the mode for the variable
      _threshold (float) : the cutoff value for the mask
      _operator (string) : which side of the threshold to hide
   """

   #===========================================================================
   def __init__(self, dictionary):
      """
      Construct the class from a dictionary.
      """

      self.construct(dictionary)

   #===========================================================================
   def construct(self, dictionary):
      """
      Initialize a mask that is fully described.

      Note: This constructor permits a default value for mode of
      "pseudocolor: full state".

      Arguments:
         dictionary (dict) : a dictionary with the needed information

      Returns:
         None

      Exceptions:
         DescriptorError for invalid options or combinations of options
      """

      # Get the values and do conversions before changing the internal storage
      v = dictionary["variable"]

      m = ModeDescriptor(dictionary.get("mode", "pseudocolor: full state"))
      if m.dimension == 1:
         raise DescriptorError("Masks cannot be defined based on a profile.")

      try:
         t = float(dictionary["threshold"])
      except ValueError:
         raise DescriptorError(
               "Mask threshold must be a floating-point number.")

      o = dictionary["operator"]
      if o not in ["<", "<=", ">", ">="]:
         raise DescriptorError("operator for constructing a mask", o,
               invalid=True)

      # Now that we've confirmed the options are valid, we store them
      self._variable = v
      self._mode = m
      self._threshold = t
      self._operator = o

   #===========================================================================
   def __repr__(self):
      """
      Defines the string representation of the instance.
      """

      return "".join(('MaskDescriptor: "', self._variable, " (", self._mode,
         ") ", self._operator, " ", str(self._threshold), '"'))

   #===========================================================================
   def __str__(self):
      """
      Defines the string representation of the instance.
      """

      return ''.join((self._variable, " (", str(self._mode), ") ",
         self._operator, " ", str(self._threshold)))

   #===========================================================================
   def apply(self, state):
      """
      Compute the described mask based on the supplied data.

      The mask is the object that hides part of what it is applied to,
      so the mask will be true for cells that are to be hidden and
      false for cells that are to be displayed.

      Arguments:
         state (SimulationState) : the current state of the simulation

      Returns:
         an array specifying the mask

      Exceptions:
         DescriptorError for invalid options or combinations of options
      """

      var = self._mode.extract(self._variable, state)

      if self._operator == "<":
         mask = (var < self._threshold)
      elif self._operator == "<=":
         mask = (var <= self._threshold)
      elif self._operator == ">":
         mask = (var > self._threshold)
      elif self._operator == ">=":
         mask = (var >= self._threshold)
      else:
         raise DescriptorError("operator for applying a mask", self._operator,
               invalid=True)

      return mask

# End of MaskDescriptor class
#==============================================================================
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#==============================================================================

class MovieDescriptor(object):
   """
   A description of a movie to be constructed.

   This class describes a movie to be made from some data set.  It
   includes information regarding the layout and design of the images
   to be used in the movie, as well as information regarding the
   encoding of the individual frames into a movie.

   Attributes:
      title (string) : title for top of frames (or None for default)
      stub (string) : the stub of the output files (images and movie)
      path (string) : the path to where the movie will be saved
      pathstub (string) : property to merge path and stub
      image_type (string) : the extension for the image files
      variable (string) : the variable to be plotted
      mode (ModeDescriptor) : the mode for computing the variable
      window (MovieWindow) : the section of the data to be plotted
      time_limits (MovieLimits) : the time range to be plotted
      value_limits (MovieLimits) : the limits on the colorbar
      fps (float) : the frame rate for the movie in frames per second
      movie_type (string) : the extension for the movie file
      make_movie (bool) : generate the movie after making the frames
      final_pause (int) : the number of extra copies of the final frame
      masks (list of MaskDescriptors) : the masks to be applied
      mask_method (string) : how to apply the masks
      xlines (list of floats) : list of x-positions to mark with lines
      colormap (string) : select the colormap
   """

   #===========================================================================
   @property
   def pathstub(self):
      if self.path is None:
         pathstub = self.stub
      elif self.path[-1] == "/":
         pathstub = "".join((self.path, self.stub))
      else:
         pathstub = "/".join((self.path, self.stub))
      return pathstub

   #===========================================================================
   def __repr__(self):
      """
      String representation.
      """

      return "".join(('MovieDescriptor: "', self.variable, " (",
         str(self.mode), ')"'))

   #===========================================================================
   def __str__(self):
      """
      String representation.
      """

      return "".join((self.variable, " (", str(self.mode), ')'))

   #===========================================================================
   def __init__(self, dictionary, mask_list):
      """
      Initialize the instance from a dictionary.
      """

      self.construct(dictionary, mask_list)

   #===========================================================================
   def construct(self, dictionary, mask_list):
      """
      Initialize the instance from a dictionary.

      Arguments:
         dictionary (dict) : a dictionary describing the parameters
         mask_list (list of MaskDescriptors) : the available masks

      Returns:
         an array specifying the mask

      Exceptions:
         DescriptorError for invalid options or combinations of options
      """

      # Use temporaries to store the values until all the processing is done.
      # That way, if a bad value is supplied, no values will be changed and you
      # won't end up with an awkward half-way blend of the old values and the
      # new values (which may be inconsistent).

      # Because the float built-in rejects None as an input, I use the
      # oh-so-cleverly-named flote (don't judge me!) function that passes None
      # through unchanged, but tries to cast anything else to a floating-point
      # number.  That lets me use dict's get() method to supply None as a
      # default value and make the code a bit shorter/cleaner.
      def flote(x):
         """
         Same as float built-in, but allows None to pass unchanged.
         """
         if x is None:
            return None
         else:
            return float(x)

      title = dictionary.get("title", None)
      stub = dictionary["stub"]
      path = dictionary.get("path", None)
      image_type = dictionary["image_type"]
      variable = dictionary["variable"]
      mode = ModeDescriptor(dictionary.get("mode", "pseudocolor: full state"))

      window = MovieWindow()
      window.x.lo = flote(dictionary.get("window_x_lo", None))
      window.x.hi = flote(dictionary.get("window_x_hi", None))
      if not window.x.valid():
         raise DescriptorError("Window lower limit must be below upper limit.")
      window.y.lo = flote(dictionary.get("window_y_lo", None))
      window.y.hi = flote(dictionary.get("window_y_hi", None))
      if not window.y.valid():
         raise DescriptorError("Window lower limit must be below upper limit.")
      window.z.lo = flote(dictionary.get("window_z_lo", None))
      window.z.hi = flote(dictionary.get("window_z_hi", None))
      if not window.z.valid():
         raise DescriptorError("Window lower limit must be below upper limit.")

      time_limits = MovieLimits()
      time_limits.lo = flote(dictionary.get("time_lo", None))
      time_limits.hi = flote(dictionary.get("time_hi", None))
      if not time_limits.valid():
         raise DescriptorError("Time lower limit must be below upper limit.")

      value_limits = MovieLimits()
      value_limits.lo = flote(dictionary.get("value_lo", None))
      value_limits.hi = flote(dictionary.get("value_hi", None))
      if not value_limits.valid():
         raise DescriptorError("Value lower limit must be below upper limit.")

      # If we are merging the frames into a movie, the frame rate (fps) and
      # movie type (movie_type) parameters are required.  However, if we are
      # not actually making the movie, they are only read in for informational
      # purposes, which means they are not required and can be None.
      make_movie = bool(dictionary.get("make_movie", True))
      if make_movie:
         fps = float(dictionary["fps"])
         movie_type = dictionary["movie_type"]
      else:
         fps = flote(dictionary.get("fps", None))
         movie_type = dictionary.get("movie_type", None)

      # The final pause can be specified as <N>f or <F>s
      # -- <N>f says to pause for <N> frames, where <N> is an integer (e.g. 12f
      #    pauses for 12 frames)
      # -- <F>s says to pause for <F> seconds, where <F> is a floating-point
      #    number (e.g. 1.7s pauses for at least 1.7 seconds; if the frame rate
      #    is 12 frames per second, then that gives 20.4 frames, and the number
      #    is always raised to the next-higher integer, so 21 frames)
      pause = dictionary.get("final_pause", "0f")
      pause_length = pause[:-1]
      unit = pause[-1]
      if unit == "f":
         try:
            final_pause = int(pause_length)
         except ValueError:
            raise DescriptorError("specification of final pause", pause,
                  invalid=True)
      elif unit == "s":
         if fps is None:
            raise DescriptorError("Cannot specify final pause in ",
               "seconds without a frame rate.", join=True)
         try:
            final_pause = int(math.ceil(float(pause_length) * fps))
         except ValueError:
            raise DescriptorError("specification of final pause", pause,
                  invalid=True)
      else:
         raise DescriptorError("specification of final pause", pause,
               invalid=True)

      # Add all the masks; skip any that are missing, but warn the user.
      mask_names = dictionary.get("masks", [])
      masks = []
      for name in mask_names:
         try:
            masks.append(mask_list[name])
         except KeyError:
            msg = "".join(('Mask "', name,
               '" is undefined and will be skipped.'))
            warnings.warn(msg, UserWarning)

      # Valid values for the mask method are "explicit" (explicitly masks
      # values using NumPy's masked arrays), "force low" (force the masked
      # cells to be equal to the lowest non-masked value), "force high" (force
      # the masked cells to be equal to the highest non-masked value in the
      # array), or a floating-point number (force the masked cells to be equal
      # to the specified value).
      mask_method = dictionary.get("mask_method", "explicit")
      if mask_method not in ["explicit", "force low", "force high"]:
         try:
            mask_method = float(mask_method)
         except ValueError:
            raise DescriptorError("mask method", mask_method, invalid=True)

      xlines = list(float(x) for x in dictionary.get("xlines", []))

      # Valid values for the colormap are "old", "new", "divergent",
      # "sequential".  The "old" colormap uses the seismic colormap with green
      # contours and black masks.  The "new" colormap actually selects between
      # two colormaps: the divergent and sequential.  The divergent colormap is
      # used for data where there is a central zero value (e.g., velocities or
      # perturbations), while the sequential colormap is used for all other
      # cases (zero is a lower-bound or zero has no special meaning).  The
      # "divergent" and "sequential" options exist in case the user truly wants
      # one of them even when the standard behavior would choose the other.
      # The default is to use the "new" colormap.
      colormap = dictionary.get("colormap", "new")
      if colormap not in ["old", "new", "sequential", "divergent"]:
         raise DescriptorError("colormap", colormap, invalid=True)

      # Now that all the values have been verified, save them
      self.title = title
      self.stub = stub
      self.path = path
      self.image_type = image_type
      self.variable = variable
      self.mode = mode
      self.window = window
      self.time_limits = time_limits
      self.value_limits = value_limits
      self.make_movie = make_movie
      self.fps = fps
      self.movie_type = movie_type
      self.final_pause = final_pause
      self.masks = masks
      self.mask_method = mask_method
      self.xlines = xlines
      self.colormap = colormap

   #===========================================================================
   def frame_data_3D(self, state):
      """
      Compute the data for a single frame of the described movie.

      Collapsing the data to a lower dimensionality (2D slice, 1D
      profile, etc) will be handled elsewhere.  This routine ignores
      the dimensionality of the requested mode and returns the full 3D
      state.  It is up to the plotting routine to collapse as needed.
      The reason for this is that the profile plots don't extract a
      single line from the data, but instead 6 lines: two can be
      processed from the 3D equivalent, one can be processed from the
      3D equivalent of the base state, and three can be processed from
      the 3D equivalent of the initial state.  Thus to plot the
      profile, the plotter will call frame_data_3D three times with
      variations on the original mode or with different states, and
      process the three resulting 3D arrays into the six desired 1D
      lines that are plotted.

      Arguments:
         state (SimulationState) : the current state of the simulation

      Returns:
         variable (numpy.ndarray) : the masked data array
         x (numpy.ndarray) : the x-coordinates
         y (numpy.ndarray) : the y-coordinates
         z (numpy.ndarray) : the z-coordinates
         vlo (float) : the lowest value for the plot range
         vhi (float) : the highest value for the plot range

      Exceptions:
         DescriptorError for invalid options or combinations of options
      """

      # Get the data for the desired variable
      var = self.mode.extract(self.variable, state)

      # Apply the limits from the MovieWindow
      x_idx = self.window.x.test(state.x[:,0,0])
      x = state.x[x_idx,0:1,0:1]

      y_idx = self.window.y.test(state.y[0,:,0])
      y = state.y[0:1,y_idx,0:1]

      z_idx = self.window.z.test(state.z[0,0,:])
      z = state.z[0:1,0:1,z_idx]

      var = var[x_idx,:,:]
      var = var[:,y_idx,:]
      var = var[:,:,z_idx]

      # Get the bounds type (may be used in masking)
      bounds_type = state.known_variables[self.variable].zero

      # Construct the overall mask (union of all masks)
      cumulative_mask = np.zeros_like(var, dtype=bool)
      for mask in self.masks:
         tmp = mask.apply(state)
         tmp = tmp[x_idx,:,:]
         tmp = tmp[:,y_idx,:]
         tmp = tmp[:,:,z_idx]
         cumulative_mask = np.logical_or(cumulative_mask, tmp)

      # Apply the mask
      if self.mask_method == "explicit":
         var = np.ma.masked_where(cumulative_mask, var)
      elif self.mask_method == "force low":
         if (self.mode.transform in ["perturbation", "contrast"] or
               bounds_type in ["lower bound", "center"]):
            # The "force low" mask forces the data to the lowest value.  In a
            # general context, that is simply the minimum value of the data
            # array.  However, if we know that zero holds a special meaning, we
            # can come up with a better value.  For the "lower bound" case, we
            # know that zero is the lowest possible value, so forcing to the
            # lowest value becomes forcing to the lowest possible value.  For
            # the "center" case, the magnitude and the sign both have meaning,
            # so we are forcing not to the lowest value but to the
            # lowest-magnitude value.
            maskval = 0.0
         else:
            maskval = np.amin(var[np.logical_not(cumulative_mask)])
         var[cumulative_mask] = maskval
      elif self.mask_method == "force high":
         maskval = np.amax(var[np.logical_not(cumulative_mask)])
         var[cumulative_mask] = maskval
      else:
         var[cumulative_mask] = self.mask_method

      return var, x, y, z

   #===========================================================================
   def build_paths(self, data_dir):
      """
      Build the paths to the movie and frame directories.
      """

      if self.path is None:
         # Default to directory
         movie_path = data_dir
      elif self.path[0] == "/":
         # Absolute path specified --- use it
         movie_path = self.path
      else:
         # Relative path from data_dir assumed
         if data_dir[-1] != "/":
            movie_path = data_dir + "/" + self.path
         else:
            movie_path = data_dir + self.path
      if movie_path[-1] != "/":
         movie_path += "/"
      frame_path = movie_path + "frames/"

      return movie_path, frame_path

   #===========================================================================

   def draw_frame(self, state, data_dir, data_num, state0=None):
      """
      Draw a single frame of the movie

      In order to be generic for future development, this will be a
      fairly lightweight wrapper.  It creates a new figure, creates a
      subplot that fills the figure, and passes the axes to a routine
      for plotting.  Future development will allow multiple subplots,
      which would require extracting multiple subplots and calling the
      appropriate drawing routines.

      Arguments:
         state (SimulationState) : the current state of the simulation
         data_dir (string) : directory where the simulation data is
         data_num (int) : the number of the current data file
         state0 (SimulationState) : the initial state of the simulation

      Returns:
         None

      Exceptions:
         whatever arises from functions called by this routine
      """

      # Generate and partition the figures
      plt.ioff()
      fig = plt.figure()
      ax1 = fig.add_axes([0, 0.05, 1, 0.9]) # left, bottom, width, height

      # Fill the panels (Axes)
      if self.mode.dimension == 1:
         self.panel_profile(state, ax1, state0)
      else:
         self.panel_2D_pseudocolor(state, ax1)

      # Verify that the directory exists
      movie_path, frame_path = self.build_paths(data_dir)
      try:
         os.makedirs(frame_path)
      except OSError as e:
         if e.errno != errno.EEXIST:
            raise

      # Save the figures
      image_file_name = "{p}{s}_{n:06d}.{e}".format(p=frame_path, s=self.stub,
            n=data_num, e=self.image_type)
      fig.savefig(image_file_name, bbox_inches="tight")
      fig.clear()

   #===========================================================================

   def panel_2D_pseudocolor(self, state, axes):
      """
      Draw a 2D pseudocolor plot as one panel of a figure.

      This routine assumes that we are plotting a slice in the xy-plane
      at a z-index of zero.  This comes from the structure of data in
      DUMSES for a 2D simulation, and can be modified later if I
      generate 3D data.

      Arguments:
         state (SimulationState) : the current state of the simulation
         axes (matplotlib.axes.Axes) : the axes to draw on

      Returns:
         None

      Exceptions:
         whatever arises from functions called by this routine
      """

      # Get the data to be drawn
      data, x, y, z = self.frame_data_3D(state)

      # Since we are drawing in 2D, reduce to the assumed geometry
      # --> Also, transform to the orientation expected by imshow
      plot_data = np.flipud(data[:,:,0])
      xx = y[0,:,0]
      yy = x[:,0,0]

      # Split the axes to have space for the colormap
      divider = make_axes_locatable(axes)
      cbar_axes = divider.append_axes("right", size="20%", pad=0.1)

      # Compute value limits
      if self.mode.absolute:
         # Absolute value by definition means no negative numbers
         bounds_type = "lower bound"
      elif self.mode.transform in ["perturbation", "contrast"]:
         # Non-absolute value perturbations/contrasts vary around zero
         bounds_type = "center"
      else:
         # Otherwise, just use the variable's own default behavior
         bounds_type = state.known_variables[self.variable].zero

      if bounds_type == "lower bound":
         vlo = 0.0
         vhi = data.max()
      elif bounds_type == "center":
         vhi = np.abs(data).max()
         vlo = -vhi
      elif bounds_type is None:
         vlo = data.min()
         vhi = data.max()
      else:
         raise DescriptorError("bounds type", bounds_type, invalid=True)
      if self.value_limits.lo is not None:
         vlo = self.value_limits.lo
      if self.value_limits.hi is not None:
         vhi = self.value_limits.hi

      # Select colormap
      if self.colormap == "new":
         if bounds_type == "center":
            colormap_choice = "divergent"
         else:
            colormap_choice = "sequential"
      else:
         colormap_choice = self.colormap

      if colormap_choice == "old":
         cdict = {'red'   : ((0.00, 0.0, 0.0),
                             (0.25, 0.0, 0.0),
                             (0.50, 1.0, 1.0),
                             (0.75, 1.0, 1.0),
                             (1.00, 0.5, 0.5)),
                  'green' : ((0.00, 0.0, 0.0),
                             (0.25, 0.0, 0.0),
                             (0.50, 1.0, 1.0),
                             (0.75, 0.0, 0.0),
                             (1.00, 0.0, 0.0)),
                  'blue'  : ((0.00, 0.3, 0.3),
                             (0.25, 1.0, 1.0),
                             (0.50, 1.0, 1.0),
                             (0.75, 0.0, 0.0),
                             (1.00, 0.0, 0.0))
                 }
         c_contour = "#00FF00"
         c_mask    = "#000000"
      elif colormap_choice == "divergent":
         cdict = {'red'   : ((0.00, 0.0, 0.0),
                             (0.30, 0.3, 0.3),
                             (0.50, 1.0, 1.0),
                             (0.70, 1.0, 1.0),
                             (1.00, 0.5, 0.5)),
                  'green' : ((0.00, 0.0, 0.0),
                             (0.30, 0.3, 0.3),
                             (0.50, 1.0, 1.0),
                             (0.70, 0.3, 0.3),
                             (1.00, 0.0, 0.0)),
                  'blue'  : ((0.00, 0.7, 0.7),
                             (0.30, 1.0, 1.0),
                             (0.50, 1.0, 1.0),
                             (0.70, 0.3, 0.3),
                             (1.00, 0.0, 0.0))
                 }
         c_contour = "#000000"
         c_mask    = "#2F2F2F"
      elif colormap_choice == "sequential":
         cdict = {'red'   : ((0.00, 0.0, 0.0),
                             (0.30, 0.0, 0.0),
                             (0.70, 1.0, 1.0),
                             (1.00, 0.5, 0.5)),
                  'green' : ((0.00, 0.0, 0.0),
                             (0.30, 0.0, 0.0),
                             (0.70, 0.0, 0.0),
                             (1.00, 0.0, 0.0)),
                  'blue'  : ((0.00, 0.7, 0.7),
                             (0.30, 1.0, 1.0),
                             (0.70, 0.0, 0.0),
                             (1.00, 0.0, 0.0))
                 }
         c_contour = "#000000"
         c_mask    = "#2F2F2F"
      else:
         raise DescriptorError("colormap", colormap_choice, invalid=True)
      my_cmap = LinearSegmentedColormap('custom', cdict)
      my_cmap.set_bad(color=c_mask, alpha=1)

      # Draw the plot
      image = axes.imshow(plot_data,
            extent=[xx.min(), xx.max(), yy.min(), yy.max()],
            aspect=1.0, cmap=my_cmap, vmin=vlo, vmax=vhi)
      axes.contour(xx, yy, x[:,:,0].repeat(y.shape[1],1),
            levels=state.params.layer_width*np.array((-1.0, 1.0)),
            colors=c_contour, linestyles="solid")
      if len(self.xlines) > 0:
         axes.contour(xx, yy, x[:,:,0].repeat(y.shape[1],1),
               levels=self.xlines, colors=c_contour,
               linestyles="dashed")
      def y_format(x, pos):
         return "{0:>9.2e}".format(x)
      formatter = FuncFormatter(y_format)
      cbar = plt.colorbar(image, cax=cbar_axes, format=formatter)

      # Construct the title
      if self.title is None:
         title = "{v} {m}\n{t:10.3e}s".format(
               v=self.variable, m=self.mode, t=state.t)
         #axes.set_title(" ".join((self.variable, str(self.mode))))
      elif self.title == "":
         title = "{t:10.3e}s".format(t=state.t)
      else:
         title = "{s}\n{t:10.3e}s".format(s=self.title, t=state.t)
         #axes.set_title(self.title)
      axes.set_title(title)

      # Compute ticks (x, y, color)
      ntickx = 6
      nticky = int(round(ntickx * (yy.max()-yy.min()) / (xx.max()-xx.min())))
      axes.xaxis.set_major_locator(MaxNLocator(ntickx))
      axes.yaxis.set_major_locator(MaxNLocator(nticky))
      cbar.locator = MaxNLocator(nbins=10)

   #===========================================================================

   def panel_profile(self, state, ax_prof, state0=None):
      """
      Draw the profiles as one panel of a figure.

      This routine assumes that we are plotting a profile along the
      x-axis, and computes the mean and standard deviation along that
      axis.

      The figure is split into two panels:
         the upper panel:
            the mean state +/- the standard deviation
            the unperturbed base state +/- initial standard deviation
            the initial maximum and minimum states
         the lower panel:
            the standard deviation
            the initial standard deviation
      The profiles derived from the initial state will not be plotted
      if None is supplied for state0.

      Arguments:
         state (SimulationState) : the current state of the simulation
         axes (matplotlib.axes.Axes) : the axes to draw on
         state0 (SimulationState) : the initial state of the simulation

      Returns:
         None

      Exceptions:
         whatever arises from functions called by this routine
      """

      # Get the data to be drawn
      data, x, y, z = self.frame_data_3D(state)
      view = data.reshape((data.shape[0], data.shape[1]*data.shape[2]))
      mean = np.mean(view, axis=1)
      stdv = np.std(view, axis=1)

      if state0 is not None:
         init, x, y, z = self.frame_data_3D(state0)
         view = init.reshape((init.shape[0], init.shape[1]*init.shape[2]))
         std0 = np.std(view, axis=1)
         min0 = np.amin(view, axis=1)
         max0 = np.amax(view, axis=1)

      if self.mode.transform == "none":
         original_mode = self.mode  # Save the mode so as to not muck it up
         self.mode = copy.deepcopy(original_mode)  # Tweak a copy
         self.mode.absolute = False
         self.mode.transform = "none"
         self.mode.reference = "base"
         data, x, y, z = self.frame_data_3D(state)
         self.mode = original_mode  # Restore the original mode
         base = data[:,0,0]
      else:
         base = np.zeros_like(mean)

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
      ax_prof.fill_between(x, mean-stdv, mean+stdv,
            facecolor='blue', alpha=0.5)
      ax_prof.plot(x, mean, lw=2, label="mean profile", color='blue')

      # Draw the bottom plot
      if state0 is not None:
         ax_stdv.fill_between(x, std0, facecolor='grey', alpha=0.5)
         ax_stdv.plot(x, std0, lw=2, color='grey')
      ax_stdv.fill_between(x, stdv, facecolor='blue', alpha=0.5)
      ax_stdv.plot(x, stdv, lw=2, label="standard deviation", color='blue')

      # Draw lines marking the source layer
      ax_prof.axvline(x=-state.params.layer_width, color='black')
      ax_prof.axvline(x= state.params.layer_width, color='black')
      ax_stdv.axvline(x=-state.params.layer_width, color='black')
      ax_stdv.axvline(x= state.params.layer_width, color='black')

      # Draw other lines
      for ax in [ax_prof, ax_stdv]:
         for xpos in self.xlines:
            ax.axvline(x=xpos, color="black", linestyle="--")

      # Construct the title
      if self.title is None:
         title = "{v} {m}\n{t:10.3e}s".format(
               v=self.variable, m=self.mode, t=state.t)
         #axes.set_title(" ".join((self.variable, str(self.mode))))
      elif self.title == "":
         title = "{t:10.3e}s".format(t=state.t)
      else:
         title = "{s}\n{t:10.3e}s".format(s=self.title, t=state.t)
         #axes.set_title(self.title)
      ax_prof.set_title(title)

      # Label axes and set ticks
      plt.setp(ax_prof.get_xticklabels(), visible=False)
      ax_prof.set_ylabel("profile")
      ax_stdv.set_xlabel("position")
      ax_stdv.set_ylabel("standard deviation")

      # Set the plotting limits (profiles panel)
      #pad = 0.075
      if self.mode.absolute:
         # Absolute value by definition means no negative numbers
         bounds_type = "lower bound"
      elif self.mode.transform in ["perturbation", "contrast"]:
         # Non-absolute value perturbations/contrasts vary around zero
         bounds_type = "center"
      else:
         # Otherwise, just use the variable's own default behavior
         bounds_type = state.known_variables[self.variable].zero

      if bounds_type == "lower bound":
         vlo = 0.0
         vhi = (mean+stdv).max()
      elif bounds_type == "center":
         vhi = max(np.abs(mean+stdv).max(), np.abs(mean-stdv).max())
         vlo = -vhi
      elif bounds_type is None:
         vlo = (mean-stdv).min()
         vhi = (mean+stdv).max()
      else:
         raise DescriptorError("bounds type", bounds_type, invalid=True)
      if self.value_limits.lo is not None:
         vlo = self.value_limits.lo
      if self.value_limits.hi is not None:
         vhi = self.value_limits.hi
      ax_prof.set_ylim([vlo, vhi])
      ax_prof.set_xlim([x.min(), x.max()])

      # Set the plotting limits (deviations panel)
      vhi = stdv.max()
      if vhi == 0.0:
         vhi = 1.0
      ax_stdv.set_ylim([0.0, vhi])
      ax_stdv.set_xlim([x.min(), x.max()])

      # Adjust tick formatting
      def y_format(x, pos):
         return "{0:>9.2e}".format(x)
      formatter = FuncFormatter(y_format)
      ax_prof.yaxis.set_major_formatter(formatter)
      ax_stdv.yaxis.set_major_formatter(formatter)

# End of MovieDescriptor class
#==============================================================================
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#==============================================================================

class ModeDescriptor(object):
   """
   The mode to use when plotting a quantity.

   This object has several different initialization routines: default,
   from a descriptive string, or from the specific properties of the
   mode.  The mode properties are visible to users, but filling them in
   with the initialization routine performs consistency checks.

   The attributes in a bit more detail:
   -- The dimension can take the value 1 (for a profile) or 2 (for a
      pseudocolor).  This can be extended to 3 when I have a chance to
      implement 3D visualizations.  When the dimension is 1, only 3
      modes are allowed: full state, or the perturbation or contrast
      relative to the base state.
   -- The transform may be "none" (the reference state), "perturbation"
      (subtract the reference state from the full state), or "contrast"
      (subtract and divide by the reference state).
   -- The reference state may be "full" (only allowed in conjunction
      with the "none" transformation to specify that the full state
      should be plotted), "base" (the steady base state for the
      unperturbed initial conditions), or "mean" (the horizontal
      average of the state).
   -- The absolute value operator may only be applied to pseudocolor
      plots of the perturbation or the contrast.

   Attributes:
      dimension (int) : the dimensionality of the result
      absolute (bool) : apply the absolute value operator or not
      transform (str) : name of the transformation (if any) to apply
      reference (str) : the reference state to use
   """

   #===========================================================================
   def __init__(self, *args, **kwargs):
      """
      Initialize the descriptor.
      """

      if len(args) == 0:
         self.__init_default()
      elif len(args) == 1:
         self.__init_from_string(args[0])
      elif len(args) == 4:
         self.__init_from_properties(*args)
      else:
         msg = "__init__ takes 0, 1, or 4 positions arguments ({n} given}"
         raise TypeError(msg.format(n=len(args)))

   #===========================================================================
   def __init_default(self):
      """
      Set to the default mode: full state pseudocolor plot.
      """
      self.dimension = 2
      self.absolute = False
      self.transform = "none"
      self.reference = "full"

   #===========================================================================
   def __init_from_string(self, string):
      """
      Initialize from a descriptive string.

      The string may be None, in which case the default occurs.

      Arguments:
         string (string) : string describing the mode, or can be None

      Returns:
         None

      Exceptions:
         DescriptorError for invalid options or combinations of options
      """

      if string is None:
         self.__init_default()
      else:
         # Normalize spacing
         mode = ' '.join(string.split())

         # Split at the colon to divide the dimension notation from the rest of
         # the details.
         dim, junk, detail = mode.partition(":")
         dim = dim.strip()
         detail = detail.strip()

         if dim == "profile":
            # 1D profiles
            self.dimension = 1
            self.absolute = False
            if detail == "full state":
               self.transform = "none"
               self.reference = "full"
            elif detail in ["perturbation", "contrast"]:
               self.transform = detail
               self.reference = "base"
            else:
               raise DescriptorError("profile", detail, invalid=True)
         elif dim == "pseudocolor":
            # 2D pseudocolor plots
            self.dimension = 2
            d1, junk, d2 = detail.rpartition(" ")
            if d2 == "state" and d1 in ["full", "base", "mean"]:
               self.absolute = False
               self.transform = "none"
               self.reference = d1
            else:
               if "perturbation" in detail:
                  self.transform = "perturbation"
               elif "contrast" in detail:
                  self.transform = "contrast"
               else:
                  raise DescriptorError("pseudocolor", detail, invalid=True)
               a, f = detail.split(self.transform)
               a = a.strip()
               f = f.strip()
               if a == "":
                  self.absolute = False
               elif a == "absolute":
                  self.absolute = True
               else:
                  raise DescriptorError("modifier", a, invalid=True)
               if f in ["", "from base"]:
                  self.reference = "base"
               elif f == "from mean":
                  self.reference = "mean"
               else:
                  raise DescriptorError("reference state", f, invalid=True)
         else:
            raise DescriptorError("dimension", dim, invalid=True)

   #===========================================================================
   def __init_from_properties(self, d, a, t, r):
      """
      Initialize from mode properties.

      Arguments:
         d (int) : dimensionality
         a (bool) : absolute value or not
         t (string) : transformation to be applied
         r (string) : reference state

      Returns:
         None

      Exceptions:
         DescriptorError for invalid options or combinations of options
      """

      if a not in [True, False]:
         raise DescriptorError("The absolute flag must be true or false.")

      if d not in [1, 2]:
         raise DescriptorError("The dimension must be one or two.")

      if a and d == 1:
         raise DescriptorError("Cannot apply absolute value to a profile.")

      if a and t not in ["perturbation", "contrast"]:
         msg = "Can only apply absolute value to a perturbation or contrast."
         raise DescriptorError(msg)

      if d == 1:
         allowed = {"full" : ["none"],
                    "base" : ["perturbation", "contrast"]}
      elif d == 2:
         allowed = {"full" : ["none"],
                    "base" : ["none", "perturbation", "contrast"],
                    "mean" : ["none", "perturbation", "contrast"]}
      else:
         raise DescriptorError("dimension", d, invalid=True)

      if t in allowed.get(r, []):
         self.dimension = d
         self.absolute = a
         self.transform = t
         self.reference = r
      else:
         raise DescriptorError("reference state/transformation pair",
               "{r}, {t}".format(r=r, t=t), invalid=True)

   #===========================================================================
   def __repr__(self):
      """
      Supply a detailed representation of the mode.

      Arguments:
         None

      Returns:
         None

      Exceptions:
         DescriptorError for invalid options or combinations of options
      """

      if self.dimension == 1:
         if self.absolute:
            raise DescriptorError("Cannot apply absolute value to a profile.")
         if self.reference == "full" and self.transform == "none":
            return "profile: full state"
         elif (self.reference == "base" and
               self.transform in ["perturbation", "contrast"]):
            return "profile: " + self.transform
         else:
            raise DescriptorError("reference state/transformation pair",
                  "{r}, {t}".format(r=r, t=t), invalid=True)
      elif self.dimension == 2:
         if self.transform == "none":
            if self.absolute:
               raise DescriptorError("Can only apply absolute value to a ",
                     "perturbation or contrast.", join=True)
            if self.reference in ["full", "base", "mean"]:
               return "pseudocolor: " + self.reference + " state"
            else:
               raise DescriptorError("reference state", self.reference,
                     invalid=True)
         elif self.transform in ["perturbation", "contrast"]:
            if self.reference in ["base", "mean"]:
               if self.absolute:
                  detail_list = ["absolute"]
               else:
                  detail_list = []
               detail_list.extend([self.transform, "from", self.reference])
               return "pseudocolor: " + " ".join(detail_list)
            else:
               raise DescriptorError("reference state", self.reference,
                     invalid=True)
         else:
            raise DescriptorError("transformation", self.transform,
                  invalid=True)
      else:
         raise DescriptorError("dimension", d, invalid=True)

   #===========================================================================
   def extract(self, variable, state):
      """
      Extract the data from the state (with appropriate operations).

      This serves as a wrapper around the SimulationState's extract
      method in order to handle complicated modes (SimulationState is
      only "aware" of the full-state and base-state modes).

      Dimensionality is ignored.  As discussed in the MovieDescriptor,
      all extraction is done on the full 3D state, and the plotting
      routine is responsible for understanding enough about the mode to
      process the data into its final form.

      Arguments:
         variable (string) : the name of the variable to be extracted
         state (SimulationState) : the state of the simulation

      Returns:
         the variable array

      Exceptions:
         DescriptorError for invalid options or combinations of options
      """

      if self.transform == "none":
         if self.absolute:
            raise DescriptorError("Can only apply absolute value to a ",
                  "perturbation or contrast.", join=True)
         if self.reference == "full":
            var = state.extract(variable)
         elif self.reference == "base":
            var = state.extract(variable, False)
         elif state.reference == "mean":
            var = state.extract(variable)
            shape = var.shape
            var = np.mean(var.reshape((shape[0], shape[1]*shape[2])), axis=1)
            var = var[:,None,None]
            var = var.repeat(shape[1], axis=1).repeat(shape[2], axis=2)
         else:
            raise DescriptorError("reference state", self.reference,
                  invalid=True)
      elif self.transform in ["perturbation", "contrast"]:
         full = state.extract(variable)
         if self.reference == "base":
            ref = state.extract(variable, False)
         elif self.reference == "mean":
            ref = state.extract(variable)
            shape = ref.shape
            ref = np.mean(ref.reshape((shape[0], shape[1]*shape[2])), axis=1)
            ref = ref[:,None,None]
            ref = ref.repeat(shape[1], axis=1).repeat(shape[2], axis=2)
         else:
            raise DescriptorError("reference state", self.reference,
                  invalid=True)
         var = full - ref
         if self.transform == "contrast":
            var = var / ref
         if self.absolute:
            var = np.abs(var)
      else:
         raise DescriptorError("transformation", self.transform, invalid=True)

      return var

# End class ModeDescriptor
#==============================================================================
#------------------------------------------------------------------------------

