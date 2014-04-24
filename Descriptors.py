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

The InvalidParameterError states that a specified parameter has been assigned
an invalid value.

The MovieWindow class describes a viewing window (limits in the coordinates)
for plotting.  It uses a MovieLimits object for each coordinate direction.
These two essentially exist for packaging more than functionality.

Attributes:
   MovieDescriptor (class) : a description of the movie to be made
   MaskDescriptor (class) : a description of a mask to be applied
   ModeDescriptor (class) : a description of a plotting mode
   InvalidParameterError (class) : an exception for invalid parameter values
   MovieLimits (class) : low/high limits of some quantity
   MovieWindow (class) : low/high limits in 3 dimensions for a viewing space
"""

import math
import numpy as np
import warnings

#==============================================================================

class InvalidParameterError(KeyError):
   """
   Signals an invalid choice of parameter.

   Attributes:
      param (string) : the parameter that was assigned an invalid value
      value (string) : the invalid value
   """

   param = None
   value = None

   def __init__(self, p, v):
      """
      Initialize with a string specifying parameter and a stringifiable value.
      """
      self.param = p
      self.value = str(v)

   def __repr__(self):
      """
      String representation.
      """
      return '"'.join(('Invalid parameter: ', self.param,
         ' may not take the value ', self.value, '.'

   __str__ = __repr__  # To mask parent class implementations of __str__

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
      __initialized (bool) : flag stating if instance has been initialized
      _variable (string) : the variable used for computing the mask
      _mode (ModeDescriptor) : the mode for the variable
      _threshold (float) : the cutoff value for the mask
      _operator (string) : which side of the threshold to hide
   """

   #===========================================================================
   def __init__(self):
      """
      A bare initialization method for a dummy mask.
      """

      self.__clear()

   #===========================================================================
   def __init__(self, dictionary):
      """
      Initializes a mask that is fully described.

      Note: This constructor permits a default value for mode of
      "pseudocolor: full state".

      Arguments:
         dictionary (dict) : a dictionary with the needed information
      """

      self.__construct(dictionary)

   #===========================================================================
   def clear(self):
      """
      Flush the data to give an uninitialized instance.
      """

      self.__initialized = False

      self._variable = None
      self._mode = None
      self._threshold = None
      self._operator = None

   __clear = clear # private copy to preserve __init__ behavior

   #===========================================================================
   def construct(self, dictionary):
      """
      Initialize a mask that is fully described.

      Note: This constructor permits a default value for mode of
      "pseudocolor: full state".

      Arguments:
         dictionary (dict) : a dictionary with the needed information
      """

      self._variable = dictionary["variable"]
      mode_name = dictionary.get("mode", "pseudocolor: full state")
      self._mode = ModeDescriptor(mode_name)
      self._threshold = float(dictionary["threshold"])
      self._operator = dictionary["operator"]

      self.__initialized = True

   __construct = construct # private copy to preserve __init__ behavior

   #===========================================================================

   def __repr__(self):
      """
      Defines the string representation of the instance.
      """

      return 'MaskDescriptor: "' + ' '.join((self._variable, self._mode,
         self._operator, str(self._threshold))) + '"'

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

      if not self.__initialized:
         raise UninitializedObjectError()

      var = state.extract(self._variable, self._mode)

      if self._operator == "<":
         mask = (var < self._threshold)
      elif self._operator == "<=":
         mask = (var <= self._threshold)
      elif self._operator == ">":
         mask = (var > self._threshold)
      elif self._operator == ">=":
         mask = (var >= self._threshold)
      else:
         raise InvalidParameterError("operator", self._operator)

      return mask

# End of MaskDescriptor class
#==============================================================================



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
   """

   #===========================================================================
   def __repr__(self):
      """
      String representation.
      """

      return " ".join(("Movie of the", self.variable, self.mode))

   #===========================================================================
   def __init__(self):
      """
      A bare initialization method for a dummy movie.
      """

      self.__clear()

   #===========================================================================
   def clear(self):
      """
      Reset to default values and clear values that have no default.
      """

      self.stub = None
      self.path = None
      self.image_type = None
      self.variable = None
      self.mode = ModeDescriptor("pseudocolor: full state")
      self.window = MovieWindow()
      self.time_limits = MovieLimits()
      self.value_limits = MovieLimits()
      self.fps = None
      self.movie_type = None
      self.make_movie = True
      self.final_pause = 0
      self.masks = []

   __clear = clear # private copy to preserve __init__ behavior

   #===========================================================================
   def __init__(self, dictionary, mask_list):
      """
      Initialize the instance from a dictionary."
      """

      self.__construct(dictionary, mask_list)

   #===========================================================================
   def construct(self, dictionary, mask_list):
      """
      Initialize the instance from a dictionary."
      """

      self.stub = dictionary["stub"]
      self.path = dictionary["path"]
      self.image_type = dictionary["image_type"]
      self.variable = dictionary["variable"]
      mode_name = dictionary.get("mode", "pseudocolor: full state")
      self.mode = ModeDescriptor(mode_name)

      self.window = MovieWindow()
      try:
         self.window.x.lo = float(dictionary["window_x_lo"])
      except KeyError:
         self.window.x.lo = None
      try:
         self.window.x.hi = float(dictionary["window_x_hi"])
      except KeyError:
         self.window.x.hi = None
      try:
         self.window.y.lo = float(dictionary["window_y_lo"])
      except KeyError:
         self.window.y.lo = None
      try:
         self.window.y.hi = float(dictionary["window_y_hi"])
      except KeyError:
         self.window.y.hi = None
      try:
         self.window.z.lo = float(dictionary["window_z_lo"])
      except KeyError:
         self.window.z.lo = None
      try:
         self.window.z.hi = float(dictionary["window_z_hi"])
      except KeyError:
         self.window.z.hi = None

      self.time_limits = MovieLimits()
      try:
         self.time_limits.lo = float(dictionary["time_lo"])
      except KeyError:
         self.time_limits.lo = None
      try:
         self.time_limits.hi = float(dictionary["time_hi"])
      except KeyError:
         self.time_limits.hi = None

      self.value_limits = MovieLimits()
      try:
         self.value_limits.lo = float(dictionary["value_lo"])
      except KeyError:
         self.value_limits.lo = None
      try:
         self.value_limits.hi = float(dictionary["value_hi"])
      except KeyError:
         self.value_limits.hi = None

      # TODO : This current logic looks like it will require FPS no matter
      #        what.  I should change this so that if make_movie is false, then
      #        options like fps that only apply to the movie will be ignored.
      self.fps = float(dictionary["fps"])
      self.movie_type = dictionary["movie_type"]
      self.make_movie = bool(dictionary.get("make_movie", True))

      pause = dictionary.get("final_pause", "0f")
      pause_length = pause[:-1]
      unit = pause[-1]
      if unit == "f":
         try:
            pause_length = int(pause_length)
         except:
            raise ValueError("Invalid spefication of final pause.")
      elif unit == "s":
         try:
            pause_length = int(math.ceil(float(pause_length) * self.fps))
         except:
            raise ValueError("Invalid spefication of final pause.")
      else:
         raise ValueError("Invalid spefication of final pause.")
      self.final_pause = pause_length

      mask_names = dictionary.get("masks", [])
      self.masks = []
      for name in mask_names:
         try:
            self.masks.append(mask_list[name])
         except KeyError:
            msg = "".join(('Mask "', name,
               '" is undefined and will be skipped.'))
            warnings.warn(msg, UserWarning)

   __construct = construct # private copy to preserve __init__ behavior

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
   def frame_data(self, state):
      """
      Compute the data for a single frame of the described movie.

      Arguments:
         state (SimulationState) : the current state of the simulation

      Returns (in order):
         variable (numpy.ndarray) : the masked data array
         x (numpy.ndarray) : the x-coordinates
         y (numpy.ndarray) : the y-coordinates
         z (numpy.ndarray) : the z-coordinates
      """

      # TODO : Need to update this to reflect new ModeDescriptor

      # Make sure we're actually computing a variable
      if self.variable is None:
         raise UninitializedObjectError()

      # Get the variable
      var = state.extract(self.variable, self.mode)

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

      # Apply the masks
      cumulative_mask = np.zeros_like(var, dtype=bool)
      for mask in self.masks:
         cumulative_mask = np.logical_or(cumulative_mask, mask.apply(state))

      var = np.ma.masked_where(cumulative_mask, var)

      # Compute value limits
      bounds_type = state.bounds_type(self.variable, self.mode)
      if bounds_type == "non-negative":
         vlo = 0.0
         vhi = var.max()
      elif bounds_type == "symmetric":
         vhi = np.abs(var).max()
         vlo = -vhi
      elif bounds_type == "general":
         vlo = var.min()
         vhi = var.max()
      else:
         pass  # TODO
      if self.value_limits.lo is not None:
         vlo = self.value_limits.lo
      if self.value_limits.hi is not None:
         vhi = self.value_limits.hi

      return var, x, y, z, vlo, vhi

# End of MovieDescriptor class
#==============================================================================



#==============================================================================

class ModeDescriptor(object):
   """
   The mode to use when plotting a quantity.

   This object has several different initialization routines: default, from a
   descriptive string, or from the specific properties of the mode.  The mode
   properties are visible to users, but filling them in with the initialization
   routine performs consistency checks.

   The attributes in a bit more detail:  The dimension can take the value 1
   (for a profile) or 2 (for a pseudocolor).  The transform may be "none" (the
   full state), "perturbation" (subtract the reference state), or "contrast"
   (subtract and divide by the reference state).  The reference state may be
   "full" (only allowed in conjunction with the "none" transformation to
   specify that the full state should be plotted), "base" (the steady base
   state for the initial conditions), or "mean" (the horizontal average of the
   state).  Some restrictions apply:  The absolute value operator may only be
   applied to pseudocolor plots of the perturbation or contrast.  The profile
   only allows the full state ("none" transformation), or the perturbation or
   contrast (only relative to the base state).

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
         raise InvalidParameterError("nargs", len(args))

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
         self.__default()
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
               raise InvalidParameterError("mode", mode)
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
                  raise InvalidParameterError("mode", mode)
               a, f = detail.split(self.transform)
               a = a.strip()
               f = f.strip()
               if a == "":
                  self.absolute = False
               elif a == "absolute":
                  self.absolute = True
               else:
                  raise InvalidParameterError("mode", mode)
               if f == "" or f == "from base":
                  self.reference = "base"
               elif f == "from mean":
                  self.reference = "mean"
               else:
                  raise InvalidParameterError("mode", mode)

   #===========================================================================
   def __init_from_properties(self, d, a, t, r):
      """
      Initialize from mode properties.
      """

      if a and (d != 2 or t not in ["perturbation", "contrast"]):
         raise InvalidParameterError("absolute", a)

      if d == 1:
         allowed = {"full" : ["none"],
                    "base" : ["perturbation", "contrast"]}
      elif d == 2:
         allowed = {"full" : ["none"],
                    "base" : ["none", "perturbation", "contrast"],
                    "mean" : ["none", "perturbation", "contrast"]}
      else:
         raise InvalidParameterError("dimension", d)

      if t in allowed.get(r, []):
         self.dimension = d
         self.absolute = a
         self.transform = t
         self.reference = r
      else:
         raise InvalidParameterError("transform", t)

   #===========================================================================
   def __repr__(self):
      """
      Supply a detailed representation of the mode.
      """

      # TODO : This was originally written for compactness.  I need to rewrite
      #        it to give clear error messages (from InvalidParameterError).
      if self.dimension == 1:
         if self.absolute:
            raise InvalidParameterError("absolute", self.absolute)
         if self.reference == "full" and self.transform == "none":
            return "profile: full state"
         else:
            if (self.reference == "base" and
                  self.transform in ["perturbation", "contrast"]):
               return "profile: " + self.transform
            else:
               raise InvalidParameterError("transform", self.transform)
      elif self.dimension == 2:
         if self.transform == "none":
            if self.absolute:
               raise InvalidParameterError("absolute", self.absolute)
            if self.reference in ["full", "base", "mean"]:
               return "pseudocolor: " + self.reference + " state"
            else:
               raise InvalidParameterError("reference", self.reference)
         elif self.transform in ["perturbation", "contrast"]:
            if self.reference in ["base", "mean"]:
               if self.absolute:
                  detail_list = ["absolute"]
               else:
                  detail_list = []
               detail_list.extend([self.transform, "from", self.reference])
               return "pseudocolor: " + " ".join(detail_list)
            else:
               raise InvalidParameterError("reference", self.reference)
         else:
            raise InvalidParameterError("transform", self.transform)
      else:
         raise InvalidParameterError("dimension", self.dimension)

#==============================================================================
