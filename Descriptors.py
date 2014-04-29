"""
Specification of movie descriptors.

The MovieDescriptor object describes the details necessary to construct a movie
from data.  It describes the design of the images for each individual frame,
and includes information relating to how those images will be encoded together
into a movie.

The MaskDescriptor object describes the details of a mask.  This is a component
used by MovieDescriptors to describe how to mask cells that are not
"interesting enough" to be plotted.

The ModeDescriptor object describes the details of how to process a variable
for a movie.  It allows specification of profiles or 2D pseudocolor plots, with
a selection of transformations that may be applied.

The MovieWindow class describes a viewing window (limits in the coordinates)
for plotting.  It uses a MovieLimits object for each coordinate direction.
These two essentially exist for packaging more than functionality.

The DescriptorError class describes an error raised by one of the Descriptors.
The errors raised by them tend to be related to illegal values (or illegal
combinations of values), so they can be logically grouped together.  Giving
them their own class means that calling routines can catch errors raised by the
Descriptors separately from other errors (e.g. a ValueError in something used
by one of the Descriptors, etc).

Attributes:
   MovieDescriptor (class) : a description of the movie to be made
   MaskDescriptor (class) : a description of a mask to be applied
   ModeDescriptor (class) : a description of a plotting mode
   MovieLimits (class) : low/high limits of some quantity
   MovieWindow (class) : low/high limits in 3 dimensions for a viewing space
   DescriptorError (class) : exception for Descriptors
"""

import math
import numpy as np
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

      Because StandardError explicitly states that it does not use any keyword
      arguments, that gives the freedom to use any keyword arguments to
      customize the initialization of DescriptorError.  Based on the arguments
      supplied, the initialization selects between the following modes:
      -- invalid=True:
         A very common message: 'Invalid <parameter>: "<value>".', where
         <parameter> and <value> are given as the first and second positional
         arguments respectively.
      -- join=True:
         It is not unusual to have to construct the message string by
         concatenating a collection of substrings, so in this case all
         positional arguments are merged.
      -- message:
         If there is only one positional argument and no special keyword
         arguments, the positional argument is saved as the message.
      -- other:
         In order to mimic the (apparent) behavior of StandardError, if
         multiple positional arguments are supplied without any special keyword
         arguments, then no message is stored.
      """

      # Save the arguments
      self.args = args
      self.kwargs = kwargs

      # Check that no more than one mutually-exclusive keywords are given
      exclusive_keywords = ["invalid", "join"]
      sum_exclusive = sum(1 for kw in exclusive_keywords if kwargs[kw])
      if sum_exclusive > 1:
         msg = "".join(("DescriptorError does not accept more than one of [",
               ", ".join(exclusive_keywords), "]."))
         raise TypeError(msg)

      # Construct the message
      if kwargs.get("invalid", False):
         self.message = 'Invalid {p}: "{v}".'.format(p=args[0], v=str(args[1]))
      elif kwargs.get("join", False):
         self.message = "".join(args)
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
      return "Limits: [" + str(self.lo) + "," + str(self.hi) + "]"

   __str__ == __repr__

# End of MovieLimits class
#==============================================================================
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#==============================================================================

class MovieWindow(object):
   """
   Viewing limits on a plot, to restrict the image to only a certain region.
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

# End of MovieWindow class
#==============================================================================
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#==============================================================================

class MaskDescriptor(object):
   """
   A description of a mask to be applied to the data prior to plotting.

   This class describes a data mask.  It specifies a variable, mode, threshold,
   and operator.  As an example to explain, if the variable is entropy, the
   mode is perturbation, the threshold is zero, and the operator is less than,
   the applied mask will hide all cells where the entropy perturbation is less
   than zero.  The operators are supplied by a string and the permitted values
   are "<=", "<", ">", and ">=".

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

      self.construct(*args)

   #===========================================================================
   def construct(self, dictionary):
      """
      Initialize a mask that is fully described.

      Note: This constructor permits a default value for mode of
      "pseudocolor: full state".

      Arguments:
         dictionary (dict) : a dictionary with the needed information
      """

      # Get the values and do conversions before changing the internal storage
      v = dictionary["variable"]
      m = ModeDescriptor(dictionary.get("mode", "pseudocolor: full state"))
      if m.dimension == 1:
         raise DescriptorError("Masks cannot be defined based on a profile.")
      t = float(dictionary["threshold"])
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

      return 'MaskDescriptor: "' + ' '.join((self._variable, self._mode,
         self._operator, str(self._threshold))) + '"'

   #===========================================================================
   def __str__(self):
      """
      Defines the string representation of the instance.
      """

      return ' '.join((self._variable, self._mode, self._operator,
         str(self._threshold)))

   #===========================================================================
   def apply(self, state):
      """
      Compute the described mask based on the supplied data.

      The mask is the object that hides part of what it is applied to, so the
      mask will be true for cells that are to be hidden and false for cells
      that are to be displayed.

      Arguments:
         state (SimulationState) : the current state of the simulation
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

   This class describes a movie to be made from some data set.  It includes
   information regarding the layout and design of the images to be used in the
   movie, as well as information regarding the encoding of the individual
   frames into a movie.

   Attributes:
      stub (string) : the stub of the output files (images and movie)
      path (string) : the path to where the files will be saved
      image_type (string) : the extension for the image files
      variable (string) : the variable to be plotted
      mode (ModeDescriptor) : the mode for computing the variable
      window (MovieWindow) : the section of the data to be plotted
      time_limits (MovieLimits) : the time range to be plotted
      value_limits (MovieLimits) : the limits on the colorbar
      fps (float) : the frames per second of the movie
      movie_type (string) : the extension for the movie file
      make_movie (bool) : flag to generate the movie after making the frames
      final_pause (int) : the number of extra copies of the final frame
      masks (list of MaskDescriptors) : the masks to be applied (if any)
      mask_method (string) : how to apply the masks
   """

   #===========================================================================
   def __repr__(self):
      """
      String representation.
      """

      return " ".join(("Movie of the", self.variable, self.mode))

   __str__ == __repr__

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

      stub = dictionary["stub"]
      path = dictionary["path"]
      image_type = dictionary["image_type"]
      variable = dictionary["variable"]
      mode = ModeDescriptor(dictionary.get("mode", "pseudocolor: full state"))

      window = MovieWindow()
      window.x.lo = flote(dictionary.get("window_x_lo", None))
      window.x.hi = flote(dictionary.get("window_x_hi", None))
      window.y.lo = flote(dictionary.get("window_y_lo", None))
      window.y.hi = flote(dictionary.get("window_y_hi", None))
      window.z.lo = flote(dictionary.get("window_z_lo", None))
      window.z.hi = flote(dictionary.get("window_z_hi", None))

      time_limits = MovieLimits()
      time_limits.lo = flote(dictionary.get("time_lo", None))
      time_limits.hi = flote(dictionary.get("time_hi", None))

      value_limits = MovieLimits()
      value_limits.lo = flote(dictionary.get("value_lo", None))
      value_limits.hi = flote(dictionary.get("value_hi", None))

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
            pause_length = int(pause_length)
         except:
            raise DescriptorError("specification of final pause", pause.
                  invalid=True)
      elif unit == "s":
         if fps is None:
            raise DescriptorError("Cannot specify final pause in ",
               "seconds without a frame rate.", join=True)
         try:
            pause_length = int(math.ceil(float(pause_length) * fps))
         except:
            raise DescriptorError("specification of final pause", pause,
                  invalid=True)
      else:
         raise DescriptorError("specification of final pause", pause,
               invalid=True)
      final_pause = pause_length

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

      # Value values for the mask method are "explicit" (explicitly masks
      # values using NumPy's masked arrays), "force-lo" (force the masked cells
      # to be equal to the lowest non-masked value in the array), "force-hi"
      # (force the masked cells to be equal to the highest non-masked value in
      # the array), or a floating-point number (force the masked cells to be
      # equal to the specified value).
      mask_method = dictionary.get("mask_method", "explicit")
      if mask_method not in ["explicit", "force-lo", "force-hi"]:
         try:
            mask_method = float(mask_method)
         except ValueError:
            raise DescriptorError("mask method", mask_method, invalid=True)

      # Now that all the values have been verified, save them
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

   #===========================================================================
   def within_time_limits(self, time):
      """
      Check if the given time is within the specified time limits.
      """

      if self.time_limits.lo is not None:
         if time < self.time_limits.lo:
            return False

      if self.time_limits.hi is not None:
         if time > self.time_limits.hi:
            return False

      return True

   #===========================================================================
   def frame_data(self, state, pad_bounds=0.0):
      """
      Compute the data for a single frame of the described movie.

      Collapsing the data to a lower dimensionality (2D slice, 1D profile, etc)
      will be handled elsewhere.  This routine ignores the dimensionality of
      the requested mode and returns the full 3D state.  It is up to the
      plotting routine to collapse as needed.  The reason for this is that the
      profile plots don't extract a single line from the data, but instead 6
      lines: two can be processed from the 3D equivalent, one can be processed
      from the 3D equivalent of the base state, and three can be processed from
      the 3D equivalent of the initial state.  Thus to plot the profile, the
      plotter will call frame_data three times with variations on the original
      mode or with different states, and process the three resulting 3D arrays
      into the six desired 1D lines that are plotted.

      Arguments:
         state (SimulationState) : the current state of the simulation
         pad_bounds (float) : extra padding around value bounds

      Returns (in order):
         variable (numpy.ndarray) : the masked data array
         x (numpy.ndarray) : the x-coordinates
         y (numpy.ndarray) : the y-coordinates
         z (numpy.ndarray) : the z-coordinates
      """

      # Get the variable
      var = self.mode.extract(self.variable, state)

      # Apply the limits from the MovieWindow
      x = state.x[:,0,0]
      idx = np.ones_like(x, dtype=bool)
      if self.window.x.lo is not None:
         idx = np.logical_and(idx, x >= self.window.x.lo)
      if self.window.x.hi is not None:
         idx = np.logical_and(idx, x <= self.window.x.hi)
      var = var[idx,:,:]
      x = x[idx]
      x = x.reshape((len(x),1,1))

      y = state.y[0,:,0]
      idx = np.ones_like(y, dtype=bool)
      if self.window.y.lo is not None:
         idx = np.logical_and(idx, y >= self.window.y.lo)
      if self.window.y.hi is not None:
         idx = np.logical_and(idx, y <= self.window.y.hi)
      var = var[:,idx,:]
      y = y[idx]
      y = y.reshape((1,len(y),1))

      z = state.z[0,0,:]
      idx = np.ones_like(z, dtype=bool)
      if self.window.z.lo is not None:
         idx = np.logical_and(idx, z >= self.window.z.lo)
      if self.window.z.hi is not None:
         idx = np.logical_and(idx, z <= self.window.z.hi)
      var = var[:,:,idx]
      z = z[idx]
      z = z.reshape((1,1,len(z)))

      # Construct the overall mask (union of all masks)
      cumulative_mask = np.zeros_like(var, dtype=bool)
      for mask in self.masks:
         cumulative_mask = np.logical_or(cumulative_mask, mask.apply(state))

      # Apply the mask
      if self.mask_method == "explicit":
         var = np.ma.masked_where(cumulative_mask, var)
      elif self.mask_method == "force-lo":
         maskval = np.amin(var[np.logical_not(cumulative_mask)])
         var[cumulative_mask] = maskval
      elif self.mask_method == "force-hi":
         maskval = np.amax(var[np.logical_not(cumulative_mask)])
         var[cumulative_mask] = maskval
      else:
         var[cumulative_mask] = self.mask_method

      # Compute value limits
      bounds_type = state.known_variables[self.variable].zero
      if bounds_type == "lower bound":
         vlo = 0.0
         vhi = (1.0 + pad_bounds) * var.max()
      elif bounds_type == "center":
         vhi = (1.0 + 2.0*pad_bounds) * np.abs(var).max()
         vlo = -vhi
      elif bounds_type is None:
         vlo = var.min()
         vhi = var.max()
         pad = pad_bounds * (vhi - vlo)
         vhi += pad
         vlo -= pad
      else:
         raise DescriptorError("bounds type", bounds_type, invalid=True)
      if self.value_limits.lo is not None:
         vlo = self.value_limits.lo
      if self.value_limits.hi is not None:
         vhi = self.value_limits.hi

      return var, x, y, z, vlo, vhi

# End of MovieDescriptor class
#==============================================================================
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#==============================================================================

class ModeDescriptor(object):
   """
   The mode to use when plotting a quantity.

   This object has several different initialization routines: default, from a
   descriptive string, or from the specific properties of the mode.  The mode
   properties are visible to users, but filling them in with the initialization
   routine performs consistency checks.

   The attributes in a bit more detail:
   -- The dimension can take the value 1 (for a profile) or 2 (for a
      pseudocolor).  This can be extended to 3 when I have a chance to
      implement 3D visualizations.  When the dimension is 1, only 3 modes are
      allowed: full state, or the perturbation or contrast relative to the base
      state.
   -- The transform may be "none" (the full state), "perturbation" (subtract
      the reference state), or "contrast" (subtract and divide by the reference
      state).
   -- The reference state may be "full" (only allowed in conjunction with the
      "none" transformation to specify that the full state should be plotted),
      "base" (the steady base state for the unperturbed initial conditions), or
      "mean" (the horizontal average of the state).
   -- The absolute value operator may only be applied to pseudocolor plots of
      the perturbation or the contrast.

   Attributes:
      dimension (int) : the dimensionality of the result
      absolute (bool) : apply the absolute value operator or not
      transform (str) : name of the transformation (if any) to apply
      reference (str) : the reference state to use
   """

   #===========================================================================
   def __init__(self, *args, **kwargs):
      """
      Select the appropriate initialization method.
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
      """

      if string is None:
         self.__init_default()
      else:
         # Normalize spacing
         mode = ' '.join(string.split())

         dim, detail = mode.split(":")
         dim = dim.strip()
         detail = detail.strip()

         if dim == "profile":
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
               if f == "" or f == "from base":
                  self.reference = "base"
               elif f == "from mean":
                  self.reference = "mean"
               else:
                  raise DescriptorError("reference state", f, invalid=True)

   #===========================================================================
   def __init_from_properties(self, d, a, t, r):
      """
      Initialize from mode properties.
      """

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
      Extract the data from the state and apply the appropriate operations.

      This serves as a wrapper around the SimulationState's extract method in
      order to handle complicated modes (SimulationState is only "aware" of the
      full-state and base-state modes).

      Dimensionality is ignored.  As discussed in the MovieDescriptor, all
      extraction is done on the full 3D state, and the plotting routine is
      responsible for understanding enough about the mode to process the data
      into its final form.
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
            var = var.reshape([shape[0], 1, 1])
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
            ref = ref.reshape([shape[0], 1, 1])
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
