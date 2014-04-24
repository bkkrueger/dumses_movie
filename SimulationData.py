"""
Specification of simulation descriptions

This module provides the specification of five classes:
   SimulationState
   SimulationInput
   InvalidVariableError
   InvalidModeError

The SimulationInput object describes the setup of the simulation.  This
includes parameters such as the upstream conditions used to construct the
initial state, the adiabatic index for the equation of state, and the
description of the sources.

The SimulationState object packages the data from the current state of a
simulation.  It also provides a mechanism for extracting a number of different
quantities in a convenient way.

The InvalidVariableError and InvalidModeError are exceptions inheriting from
KeyError signaling invalid choices of variable and mode respectively.

Attributes:
   InvalidVariableError (class) : an error for invalid variable requests
   InvalidModeError (class) : an error for invalid mode requests
   SimulationInput (class) : a description of the simulation
   SimulationState (class) : the current state of the simulation
"""

import numpy as np

from dumpy_v05.data.rd_dumses import DumsesData
from UninitializedObject import UninitializedObjectError

#==============================================================================

class InvalidVariableError(KeyError):
   """
   Signals an invalid choice of variable.

   Attributes:
      invalid_variable (string) : name of the invalid variable requested
   """

   invalid_variable = None

   def __init__(self, string):
      """
      Initialize with a string naming the invalid variable.
      """
      self.invalid_variable = string

   def __repr__(self):
      """
      Supply a detailed representation of the exception.
      """
      return '"'.join(("Invalid choice of variable: ",
         str(self.invalid_variable), "."))

#==============================================================================

class InvalidModeError(KeyError):
   """
   Signals an invalid choice of mode.

   Attributes:
      invalid_mode (string) : name of the invalid mode requested
   """
   invalid_mode = None

   def __init__(self, string):
      """
      Initialize with a string naming the invalid mode.
      """
      self.invalid_mode == string

   def __repr__(self):
      """
      Supply a detailed representation of the exception.
      """
      return '"'.join(("Invalid choice of mode: ",
         str(self.invalid_mode), "."))

#==============================================================================

class SimulationInput(object):
   """
   A description of the simulation.

   This class describes the simulation based on the values extracted from the
   inputs file.  It includes methods to initialize from a file and to construct
   the source functions.

   Some attributes are carries for computation purposes (such as _heat_coef and
   _grav_coef), while others for informational purposes (such as _Kheat and
   _Kgrav).  Only the necessary attributes are "public", while the others are
   "hidden".

   Attributes:
      __initialized (bool) : flag to specify whether the object is initialized
      gamma (float) : adiabatic index
      pres_up (float) : scaling value for pressure (used especially in entropy)
      dens_up (float) : scaling value for density (used especially in entropy)
      _csnd_up (float) : scaling value for sound speed
      _mach_up (float) : scaling value for Mach number
      layer_width (float) : width parameter for source layer
      _shape (string) : name of the shape function for the sources
      _heat_coef (float) : coefficient for heating
      _grav_coef (float) : coefficient for gravity
      _Kheat (float) : dimensionless heating strength factor
      _Kgrav (float) : dimensionless gravity strength factor
   """

   #===========================================================================
   def __repr__(self):
      """
      Supply a detailed representation of the class.
      """
      if self.__initialized:
         string = "   \n".join(("Simulation inputs:",
            "gamma = {0}".format(self.gamma),
            "shape = {0}".format(self._shape),
            "Kheat = {0}".format(self._Kheat),
            "Kgrav = {0}".format(self._Kgrav),
            "upstream conditions:",
            "   dens = {0}".format(self.dens_up),
            "   csnd = {0}".format(self._csnd_up),
            "   Mach = {0}".format(self._mach_up)
            ))
      else:
         string = "Simulation inputs (uninitialized)"
      return string

   #===========================================================================
   def __init__(self):
      """
      A bare initialization method for an empty class.
      """

      self.__clear()

   #===========================================================================
   def __init__(self, inputs_file_name):
      """
      Initialize the class from an inputs file.

      Arguments:
         inputs_file_name (string) : the name of the inputs file to parse
      """

      self.__construct(inputs_file_name)

   #===========================================================================
   def clear(self):
      """
      Flush the data to give an uninitialized instance.
      """

      self.__initialized = False

      self.gamma = None

      self.pres_up = None
      self.dens_up = None
      self._csnd_up = None
      self._mach_up = None

      self.layer_width = None
      self._shape = None

      self._heat_coef = None
      self._grav_coef = None
      self._Kheat = None
      self._Kgrav = None

   __clear = clear # private copy to preserve __init__ behavior

   #===========================================================================
   def construct(self, inputs_file_name):
      """
      Initialize the class from an inputs file.

      Arguments:
         inputs_file_name (string) : the name of the inputs file to parse
      """

      # Data to extract from the parameters file, using defaults as appropriate
      # to match DUMSES code
      inputs_dict = {"gamma" : str(4.0 / 3.0),
                     "csnd_up" : str(1.0),
                     "dens_up" : str(1.0),
                     "mach_up" : str(0.1),
                     "layer_limit" : str(1.0),
                     "layer_shape" : 'trapezoid',
                     "Kheat" : str(1.0e-2),
                     "Kgrav" : str(1.0)}

      # Loop over inputs file and extract values
      with open(inputs_file_name) as inputsfile:
         for line in inputsfile: # Loop over lines in inputs file
            line = line.lstrip() # Remove leading spaces
            if len(line) == 0:   # If the line is blank, skip it
               continue
            if line[0] == '!':   # If the line is a comment, skip it
               continue
            for parameter in inputs_dict.keys():   # Check if line contains any
               if parameter in line:               #   of parameters we want
                  line = line.rstrip()             # Remove trailing "\n"
                  # Get everything after "=" and before (optional) "/"
                  inputs_dict[parameter] = line.split("=")[1].split("/")[0]

      # Conversions
      for key in inputs_dict.keys():
         # Handle special cases
         if key == "layer_shape":
            temp = inputs_dict[key].strip()
            if ((temp[0] == "'" and temp[-1] == "'") or
                (temp[0] == '"' and temp[-1] == '"')):
               inputs_dict[key] = temp[1:-1].strip()
            else:
               inputs_dict[key] = temp
               
         else:
            # Default case: conversion to float (with accounting for Fortran's
            # scientific notation for double-precision constants)
            inputs_dict[key] = float(inputs_dict[key].replace("d","e"))

      # Construct desired quantities
      cs = inputs_dict["csnd_up"] # I'm sick of typing 'inputs_dict["csnd_up"]'

      self.gamma = inputs_dict["gamma"]

      self._csnd_up = cs
      self._mach_up = inputs_dict["mach_up"]
      self.dens_up = inputs_dict["dens_up"]
      self.pres_up = self.dens_up * cs**2 / self.gamma

      self.layer_width = inputs_dict["layer_limit"]
      self._shape = inputs_dict["layer_shape"]

      self._Kheat = inputs_dict["Kheat"]
      self._heat_coef = self._Kheat * \
            (self._mach_up * cs**3) / (self.layer_width * self.gamma)
      self._Kgrav = inputs_dict["Kgrav"]
      self._grav_coef = - self._Kgrav *  cs**2 / self.layer_width

      self.__initialized = True

   __construct = construct # private copy to preserve __init__ behavior

   #===========================================================================
   def _shape_function(self, x, y, z):
      """
      Return the shape function.

      Arguments:
         x (numpy.ndarray) : the x coordinates
         y (numpy.ndarray) : the y coordinates
         z (numpy.ndarray) : the z coordinates
      """

      if not self.__initialized:
         raise UninitializedObjectError()

      XX = np.abs(x / self.layer_width)

      if self._shape == 'trapezoid':
         f = np.maximum(0.0, np.minimum(2.0 * (1.0 - XX), 1.0))
      elif self._shape == 'square':
         f = np.zeros_like(x)
         f[XX <= 1.0] = 1.0
      else:
         # Default to no source: While this behavior may allow some careless
         # usage to slip through, it mimics the shape function code in DUMSES.
         f = np.zeros_like(x)

      return f

   #===========================================================================
   def gravity(self, x, y, z):
      """
      Return the gravitational acceleration

      Arguments:
         x (numpy.ndarray) : the x coordinates
         y (numpy.ndarray) : the y coordinates
         z (numpy.ndarray) : the z coordinates
      """

      if not self.__initialized:
         raise UninitializedObjectError()

      grav = np.zeros((x.shape[0], y.shape[1], z.shape[2], 3))
      grav[...,0] = self._grav_coef * self._shape_function(x,y,z)

      return grav

   #===========================================================================
   def heating(self, x, y, z, density):
      """
      Return the heating rate (source for total energy density per time).

      Arguments:
         x (numpy.ndarray) : the x coordinates
         y (numpy.ndarray) : the y coordinates
         z (numpy.ndarray) : the z coordinates
         density (numpy.ndarray) : the density
      """

      if not self.__initialized:
         raise UninitializedObjectError()

      return (self._heat_coef * density * self._shape_function(x, y, z))

# End class SimulationInput
#==============================================================================












#==============================================================================

class SimulationState(object):
   """
   The current state of a simulation, with methods to extract quantities.

   This class stores the current state of a simulation, which includes the full
   and base states of mass density, momentum density, total energy density, as
   well as the time and axes for the state.  This much is similar to the
   DumsesData class, to which this is closely related (and a method exists to
   construct the SimulationState from a DumsesData because of this close
   relation).  The primary difference is that the SimulationState will also
   have a routine that accepts the name of a variable and of a mode, which will
   use private methods to compute the desired variable in the specified mode.
   Having a separate class also allows me to tweak the internals to fit my
   needs without changing the standard DumsesData class (e.g. if I find time to
   update this to be parallel so that I can visualize large simulations that
   would fill the memory on a single core).

   Attributes:
      __initialized (bool) : flag stating if instance has been initialized
      _dens (numpy.ndarray) : mass density full state array
      _momv (numpy.ndarray) : momentum density full state array
      _Ener (numpy.ndarray) : total energy density full state array
      _dens0 (numpy.ndarray) : mass density base state array
      _momv0 (numpy.ndarray) : momentum density base state array
      _Ener0 (numpy.ndarray) : total energy density base state array
      x (numpy.ndarray) : x coordinates array
      y (numpy.ndarray) : y coordinates array
      z (numpy.ndarray) : z coordinates array
      t (float) : time
      params (SimulationInputs) : important parameters from inputs file
      _variable_functions (dict) : known variables and associated functions
      known_variables (list of strings) : list of the names of known variables
      _mode_codes (dict) : known modes and associated codes
      known_modes (list of strings) : list of the names of known modes
   """

   #===========================================================================
   def __repr__(self):
      """
      Supply a detailed representation of the class.
      """
      if self.__initialized:
         string = "   \n".join(("Simulation data:",
            "size = ({Nx}, {Ny}, {Nz})".format(
               Nx=self.x.size[0], Ny=self.y.size[1], Nz=self.z.size[2]),
            "time = {0}".format(self.t)
            ))
      else:
         string = "Simulation data (uninitialized)"
      return string

   #===========================================================================
   def __init__(self):
      """
      A bare initialization method for an empty class.
      """

      self.__clear()

   #===========================================================================
   def __init__(self, dumpy, input_parameters):
      """
      Initialize the SimulationState from the supplied information.

      Arguments:
         dumpy (DumsesData) : to be used as a basis to construct this object
         input_parameters (SimulationInputs) : parameters from inputs file
      """

      self.__construct(dumpy, input_parameters)

   #===========================================================================
   def clear(self):
      """
      Flush the data to give an uninitialized instance.
      """

      self.__initialized = False

      self._dens = None
      self._momv = None
      self._Ener = None

      self._dens0 = None
      self._momv0 = None
      self._Ener0 = None

      self.x = None
      self.y = None
      self.z = None

      self.t = None

      self.params = None

   __clear = clear # private copy to preserve __init__ behavior

   #===========================================================================
   def construct(self, dumpy, input_parameters):
      """
      Construct the SimulationState from the supplied information.

      Arguments:
         dumpy (DumsesData) : to be used as a basis to construct this object
         input_parameters (SimulationInputs) : parameters from inputs file
      """

      self.params = input_parameters

      self.t = dumpy.time

      self.x = dumpy.x.reshape((len(dumpy.x),1,1))
      self.y = dumpy.y.reshape((1,len(dumpy.y),1))
      self.z = dumpy.z.reshape((1,1,len(dumpy.z)))

      self._dens = dumpy.rho
      self._momv = dumpy.rhou
      self._Ener = dumpy.E

      temp = dumpy.rho0.reshape((len(dumpy.x),1,1))
      temp = temp.repeat(len(dumpy.y), axis=1)
      temp = temp.repeat(len(dumpy.z), axis=2)
      self._dens0 = temp

      temp = dumpy.rhou0.reshape((len(dumpy.x),1,1))
      temp = temp.repeat(len(dumpy.y), axis=1)
      temp = temp.repeat(len(dumpy.z), axis=2)
      self._momv0 = temp

      temp = dumpy.E0.reshape((len(dumpy.x),1,1))
      temp = temp.repeat(len(dumpy.y), axis=1)
      temp = temp.repeat(len(dumpy.z), axis=2)
      self._Ener0 = temp

      self.__initialized = True

   __construct = construct # private copy to preserve __init__ behavior

   #===========================================================================
   # TODO : Is SimulationData the best place for this?  There should be a
   # better way to implement this that does not require an instance of
   # SimulationData, but still ties it to the lists of known variables and
   # modes.  Perhaps a public dictionary would be more effective.
   def bounds_type(self, variable_name, mode):
      """
      Specify the structure of the bounds.

      Most physical quantities consider zero to be a special value.  In some
      cases, zero is a lower bound.  An example of this would be density, where
      a negative density has no physical meaning, so that physically density
      may only be zero or positive.  We refer to these quantities as
      "non-negative".  In other cases, zero is a central value.  An example of
      this would be x-velocity, where we want to preserve zero in the center of
      the range to clearly demonstrate both magnitude and direction.  We refer
      to these quantities as "symmetric".  Other quantities make no reference
      to zero as a significant value.  An example of this would be entropy,
      where an arbitrary shift may be applied as only entropy differences are
      physical.  We refer to these quantities as "general".  This routine
      determines, based on the variable and mode, which type of bounds is
      appropriate for plotting.

      Arguments:
         variable_name (string) : name of the variable
         mode (ModeDescriptor) : mode for extracting the data
      """

      if variable_name not in self.known_variables:
         raise InvalidVariableError(variable_name)

      if mode.transform == "none":
         # Get the appropriate variable function
         # TODO : This can be better implemented, but for the moment it is a
         #        functional (if inelegant) placeholder.
         if variable_name == "density":
            return "non-negative"
         elif variable_name == "x momentum":
            return "symmetric"
         elif variable_name == "y momentum":
            return "symmetric"
         elif variable_name == "z momentum":
            return "symmetric"
         elif variable_name == "momentum magnitude":
            return "non-negative"
         elif variable_name == "x velocity":
            return "symmetric"
         elif variable_name == "y velocity":
            return "symmetric"
         elif variable_name == "z velocity":
            return "symmetric"
         elif variable_name == "velocity magnitude":
            return "non-negative"
         elif variable_name == "x vorticity":
            return "symmetric"
         elif variable_name == "y vorticity":
            return "symmetric"
         elif variable_name == "z vorticity":
            return "symmetric"
         elif variable_name == "vorticity magnitude":
            return "non-negative"
         elif variable_name == "specific x vorticity":
            return "symmetric"
         elif variable_name == "specific y vorticity":
            return "symmetric"
         elif variable_name == "specific z vorticity":
            return "symmetric"
         elif variable_name == "specific vorticity magnitude":
            return "non-negative"
         elif variable_name == "kinetic energy density":
            return "non-negative"
         elif variable_name == "specific kinetic energy":
            return "non-negative"
         elif variable_name == "total energy density":
            return "non-negative"
         elif variable_name == "specific total energy":
            return "non-negative"
         elif variable_name == "internal energy density":
            return "non-negative"
         elif variable_name == "specific internal energy":
            return "non-negative"
         elif variable_name == "pressure":
            return "non-negative"
         elif variable_name == "specific entropy":
            return "general"
         elif variable_name == "entropy density":
            return "general"
         elif variable_name == "sound speed":
            return "non-negative"
         elif variable_name == "Mach number":
            return "non-negative"
         elif variable_name == "specific enthalpy":
            return "non-negative"
         elif variable_name == "enthalpy density":
            return "non-negative"
         elif variable_name == "convective growth rate":
            return "non-negative"
         else:
            raise InvalidVariableError(variable_name)
         #try:
         #   func = self._variable_functions[variable_name]
         #except KeyError, ke:
         #   raise InvalidVariableError(str(ke))
      elif mode.transform in ["perturbation", "contrast"]:
         if mode.absolute:
            return "non-negative"
         else:
            return "symmetric"
      else:
         raise InvalidModeError(mode_name)

   #===========================================================================
   def extract(self, variable_name, mode):
      """
      Extract the desired variable in the desired mode.

      Arguments:
         variable_name (string) : name of the variable to be extracted
         mode (ModeDescriptor) : mode for extracting the data
      """

      if not self.__initialized:
         raise UninitializedObjectError()

      # Get the appropriate variable function
      # TODO : This can be better implemented, but for the moment it is a
      #        functional (if inelegant) placeholder.  It would be better to
      #        use a dictionary to allow hashing, but I had some issues on the
      #        first pass and decided to save it for later.  See the notes
      #        below (by the known_variables attribute) for some comments.
      if variable_name == "density":
         func = self._func_density
      elif variable_name == "x momentum":
         func = self._func_x_momentum
      elif variable_name == "y momentum":
         func = self._func_y_momentum
      elif variable_name == "z momentum":
         func = self._func_z_momentum
      elif variable_name == "momentum magnitude":
         func = self._func_momentum_magnitude
      elif variable_name == "x velocity":
         func = self._func_x_velocity
      elif variable_name == "y velocity":
         func = self._func_y_velocity
      elif variable_name == "z velocity":
         func = self._func_z_velocity
      elif variable_name == "velocity magnitude":
         func = self._func_velocity_magnitude
      elif variable_name == "x vorticity":
         func = self._func_x_vorticity
      elif variable_name == "y vorticity":
         func = self._func_y_vorticity
      elif variable_name == "z vorticity":
         func = self._func_z_vorticity
      elif variable_name == "vorticity magnitude":
         func = self._func_vorticity_magnitude
      elif variable_name == "specific x vorticity":
         func = self._func_x_specific_vorticity
      elif variable_name == "specific y vorticity":
         func = self._func_y_specific_vorticity
      elif variable_name == "specific z vorticity":
         func = self._func_z_specific_vorticity
      elif variable_name == "specific vorticity magnitude":
         func = self._func_specific_vorticity_magnitude
      elif variable_name == "kinetic energy density":
         func = self._func_kinetic_energy_density
      elif variable_name == "specific kinetic energy":
         func = self._func_kinetic_energy_specific
      elif variable_name == "total energy density":
         func = self._func_total_energy_density
      elif variable_name == "specific total energy":
         func = self._func_total_energy_specific
      elif variable_name == "internal energy density":
         func = self._func_internal_energy_density
      elif variable_name == "specific internal energy":
         func = self._func_internal_energy_specific
      elif variable_name == "pressure":
         func = self._func_pressure
      elif variable_name == "specific entropy":
         func = self._func_entropy_specific
      elif variable_name == "entropy density":
         func = self._func_entropy_density
      elif variable_name == "sound speed":
         func = self._func_sound_speed
      elif variable_name == "Mach number":
         func = self._func_Mach_number
      elif variable_name == "specific enthalpy":
         func = self._func_enthalpy_specific
      elif variable_name == "enthalpy density":
         func = self._func_enthalpy_density
      elif variable_name == "convective growth rate":
         func = self._func_convective_growth_rate
      else:
         raise InvalidVariableError(variable_name)
      #try:
      #   func = self._variable_functions[variable_name]
      #except KeyError, ke:
      #   raise InvalidVariableError(str(ke))

      # The only time we don't need





      # Get the appropriate mode code
      try:
         code = self._mode_codes[mode_name]
      except KeyError, ke:
         raise InvalidModeError(str(ke))

      # Compute the correct mode
      if self._mode_codes[mode_name] == 0:
         # Mode 0: full state
         var = func('full state')
      else:
         # Any other mode: compare against some reference state
         if code > 0:
            # Positive mode: compare against the base state
            ref = func('base state')
         else: # code < 0
            # Negative mode: compare against the mean state
            ref = func('full state')
            shape = ref.shape
            ref = np.mean(ref.reshape((shape[0], shape[1]*shape[2])), axis=1)
            newshape = np.ones_like(shape)
            newshape[0] = shape[0]
            ref = ref.reshape(newshape)
            ref = ref.repeat(shape[1], axis=1).repeat(shape[2], axis=2)
         # Now that we've selected the correct reference state, only the
         # absolute value of the code matters
         abscode = abs(code)
         if abscode == 1:
            # |Mode| == 1: reference state
            var = ref
         else:
            # Otherwise: we are looking at the perturbation relative to the
            # reference state
            var = func('full state') - ref
            # |Mode| == 2: perturbation
            if abscode == 3:
               # |Mode| == 3: absolute value of the perturbation
               var = np.abs(var)
            elif abscode > 3:
               # Otherwise: we are looking at the contrast (the perturbation
               # scaled against the reference state)
               var = var / ref
               # |Mode| == 4: contrast (nothing more needed)
               # |Mode| == 5: absolute value of the contrast
               if abscode == 5:
                  var = np.abs(var)

      return var

   #===========================================================================
   def _x_derivative(self, q):
      """
      The derivative of quantity q with respect to x.

      NOTE: I assume that the boundary conditions in the x direction are
      zero-gradient.  I could write this to be flexible by having the
      SimulationInputs class parse the bval_in and bval_out entries in the
      inputs file, but if a user-defined boundary value is used then there is
      no way to know what to do.  Thus I will for now simply be lazy and assume
      the boundary conditions that I use in all of my own simulations and leave
      it to other users to read the documentation and adjust accordingly for
      their own uses.

      Arguments:
         q (numpy.ndarray) : the quantity to be differentiated
      """

      if q.shape[0] == 1:
         dqdx = np.zeros_like(q)
      else:
         dqdx = np.empty_like(q)
         dqdx[0,...] = 0.0
         if q.shape[0] >= 3:
            dqdx[1:-1,...] = (q[2:,...] - q[:-2,...]) / \
                             (self.x[2:,...] - self.x[:-2,...])
         dqdx[-1,...] = 0.0

      return dqdx

   #===========================================================================
   def _y_derivative(self, q):
      """
      The derivative of quantity q with respect to y.

      NOTE: I assume that the boundary conditions in the y direction are
      periodic.  I could write this to be flexible by having the
      SimulationInputs class parse the bval_in and bval_out entries in the
      inputs file, but if a user-defined boundary value is used then there is
      no way to know what to do.  Thus I will for now simply be lazy and assume
      the boundary conditions that I use in all of my own simulations and leave
      it to other users to read the documentation and adjust accordingly for
      their own uses.

      Arguments:
         q (numpy.ndarray) : the quantity to be differentiated
      """

      if q.shape[1] == 1:
         dqdy = np.zeros_like(q)
      else:
         dqdy = np.empty_like(q)
         dqdy[:,0,...] = (q[:,1,...] - q[:,-1,...]) / \
                         (2.0 * (self.y[:,1,...] - self.y[:,0,...]))
         if q.shape[1] >= 3:
            dqdy[:,1:-1,...] = (q[:,2:,...] - q[:,:-2,...]) / \
                               (self.y[:,2:,...] - self.y[:,:-2,...])
         dqdy[:,-1,...] = (q[:,0,...] - q[:,-2,...]) / \
                          (2.0 * (self.y[:,-1,...] - self.y[:,-2,...]))

      return dqdy

   #===========================================================================
   def _z_derivative(self, q):
      """
      The derivative of quantity q with respect to z.

      NOTE: I assume that the boundary conditions in the z direction are
      periodic.  I could write this to be flexible by having the
      SimulationInputs class parse the bval_in and bval_out entries in the
      inputs file, but if a user-defined boundary value is used then there is
      no way to know what to do.  Thus I will for now simply be lazy and assume
      the boundary conditions that I use in all of my own simulations and leave
      it to other users to read the documentation and adjust accordingly for
      their own uses.

      Arguments:
         q (numpy.ndarray) : the quantity to be differentiated
      """

      if q.shape[2] == 1:
         dqdz = np.zeros_like(q)
      else:
         dqdz = np.empty_like(q)
         dqdz[:,:,0,...] = (q[:,:,1,...] - q[:,:,-1,...]) / \
                           (2.0 * (self.z[:,:,1,...] - self.z[:,:,0,...]))
         if q.shape[2] >= 3:
            dqdz[:,:,1:-1,...] = (q[:,:,2:,...] - q[:,:,:-2,...]) / \
                                 (self.z[:,:,2:,...] - self.z[:,:,:-2,...])
         dqdz[:,:,-1,...] = (q[:,:,0,...] - q[:,:,-2,...]) / \
                            (2.0 * (self.z[:,:,-1,...] - self.z[:,:,-2,...]))

      return dqdz

   #===========================================================================
   def _func_density(self, full_or_base):
      """
      Extract the density in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      if full_or_base == "full state":
         return self._dens
      elif full_or_base == "base state":
         return self._dens0
      else:
         raise InvalidModeError(full_or_base)

   #===========================================================================
   def _func_x_momentum(self, full_or_base):
      """
      Extract the x-momentum in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      if full_or_base == "full state":
         return self._momv[:,:,:,0]
      elif full_or_base == "base state":
         return self._momv0
      else:
         raise InvalidModeError(full_or_base)

   #===========================================================================
   def _func_y_momentum(self, full_or_base):
      """
      Extract the y-momentum in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      if full_or_base == "full state":
         return self._momv[:,:,:,1]
      elif full_or_base == "base state":
         return np.zeros_like(self._momv0)
      else:
         raise InvalidModeError(full_or_base)

   #===========================================================================
   def _func_z_momentum(self, full_or_base):
      """
      Extract the z-momentum in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      if full_or_base == "full state":
         return self._momv[:,:,:,2]
      elif full_or_base == "base state":
         return np.zeros_like(self._momv0)
      else:
         raise InvalidModeError(full_or_base)

   #===========================================================================
   def _func_momentum_magnitude(self, full_or_base):
      """
      Extract the magnitude of the momentum in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      if full_or_base == "full state":
         return np.sqrt(np.sum(self._momv**2, axis=-1))
      elif full_or_base == "base state":
         return np.abs(self._momv0)
      else:
         raise InvalidModeError(full_or_base)

   #===========================================================================
   def _func_x_velocity(self, full_or_base):
      """
      Extract the x-velocity in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      return self._func_x_momentum(full_or_base) / \
             self._func_density(full_or_base)

   #===========================================================================
   def _func_y_velocity(self, full_or_base):
      """
      Extract the y-velocity in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      return self._func_y_momentum(full_or_base) / \
             self._func_density(full_or_base)

   #===========================================================================
   def _func_z_velocity(self, full_or_base):
      """
      Extract the z-velocity in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      return self._func_z_momentum(full_or_base) / \
             self._func_density(full_or_base)

   #===========================================================================
   def _func_velocity_magnitude(self, full_or_base):
      """
      Extract the magnitude of the velocity in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      return self._func_momentum_magnitude(full_or_base) / \
             self._func_density(full_or_base)

   #===========================================================================
   def _func_x_vorticity(self, full_or_base):
      """
      Extract the x-vorticity in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      vely = self._func_y_velocity(full_or_base)
      velz = self._func_z_velocity(full_or_base)
      vrtx = self._y_derivative(velz) - self._z_derivative(vely)
      return vrtx

   #===========================================================================
   def _func_y_vorticity(self, full_or_base):
      """
      Extract the y-vorticity in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      velx = self._func_x_velocity(full_or_base)
      velz = self._func_z_velocity(full_or_base)
      vrty = self._z_derivative(velx) - self._x_derivative(velz)
      return vrty

   #===========================================================================
   def _func_z_vorticity(self, full_or_base):
      """
      Extract the z-vorticity in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      velx = self._func_x_velocity(full_or_base)
      vely = self._func_y_velocity(full_or_base)
      vrtz = self._x_derivative(vely) - self._y_derivative(velx)
      return vrtz

   #===========================================================================
   def _func_vorticity_magnitude(self, full_or_base):
      """
      Extract the magnitude of the vorticity in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      vrtx = self._func_x_vorticity(full_or_base)
      vrty = self._func_y_vorticity(full_or_base)
      vrtz = self._func_z_vorticity(full_or_base)
      magnitude = np.sqrt(vrtx**2 + vrty**2 + vrtz**2)
      return magnitude

   #===========================================================================
   def _func_x_specific_vorticity(self, full_or_base):
      """
      Extract the specific x-vorticity in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      vrtx = self._func_x_vorticity(full_or_base)
      dens = self._func_density(full_or_base)
      return vrtx / dens

   #===========================================================================
   def _func_y_specific_vorticity(self, full_or_base):
      """
      Extract the specific y-vorticity in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      vrty = self._func_y_vorticity(full_or_base)
      dens = self._func_density(full_or_base)
      return vrty / dens

   #===========================================================================
   def _func_z_specific_vorticity(self, full_or_base):
      """
      Extract the specific z-vorticity in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      vrtz = self._func_z_vorticity(full_or_base)
      dens = self._func_density(full_or_base)
      return vrtz / dens

   #===========================================================================
   def _func_specific_vorticity_magnitude(self, full_or_base):
      """
      Extract the magnitude of the specific vorticity in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      vrtx = self._func_x_vorticity(full_or_base)
      vrty = self._func_y_vorticity(full_or_base)
      vrtz = self._func_z_vorticity(full_or_base)
      magnitude = np.sqrt(vrtx**2 + vrty**2 + vrtz**2)
      dens = self._func_density(full_or_base)
      return magnitude / dens

   #===========================================================================
   def _func_kinetic_energy_density(self, full_or_base):
      """
      Extract the kinetic energy density in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      momentum = self._func_momentum_magnitude(full_or_base)
      kinetic = 0.5 * momentum**2 / self._func_density(full_or_base)
      return kinetic

   #===========================================================================
   def _func_kinetic_energy_specific(self, full_or_base):
      """
      Extract the specific kinetic energy in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      momentum = self._func_momentum_magnitude(full_or_base)
      kinetic = 0.5 * (momentum / self._func_density(full_or_base))**2
      return kinetic

   #===========================================================================
   def _func_total_energy_density(self, full_or_base):
      """
      Extract the total energy density in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      if full_or_base == "full state":
         return self._Ener
      elif full_or_base == "base state":
         return self._Ener0
      else:
         raise InvalidModeError(full_or_base)

   #===========================================================================
   def _func_total_energy_specific(self, full_or_base):
      """
      Extract the specific total energy in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      Ener = self._func_total_energy_density(full_or_base) / \
             self._func_density(full_or_base)

      return Ener

   #===========================================================================
   def _func_internal_energy_density(self, full_or_base):
      """
      Extract the internal energy density in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      return self._func_total_energy_density(full_or_base) - \
             self._func_kinetic_energy_density(full_or_base)

   #===========================================================================
   def _func_internal_energy_specific(self, full_or_base):
      """
      Extract the specific internal energy in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      return self._func_total_energy_specific(full_or_base) - \
             self._func_kinetic_energy_specific(full_or_base)

   #===========================================================================
   def _func_pressure(self, full_or_base):
      """
      Extract the pressure in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      return self._func_internal_energy_density(full_or_base) * \
             (self.params.gamma - 1.0)

   #===========================================================================
   def _func_entropy_specific(self, full_or_base):
      """
      Extract the specific entropy in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      PP = self._func_pressure(full_or_base) / self.params.pres_up
      dd = self._func_density(full_or_base) / self.params.dens_up
      Y = self.params.gamma
      return np.log(PP * dd**(-Y)) / (Y - 1.0)

   #===========================================================================
   def _func_entropy_density(self, full_or_base):
      """
      Extract the entropy density in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      return self._func_entropy_specific(full_or_base) * \
             self._func_density(full_or_base)

   #===========================================================================
   def _func_sound_speed(self, full_or_base):
      """
      Extract the speed of sound in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      c2 = self.params.gamma * self._func_pressure(full_or_base) / \
            self._func_density(full_or_base)
      return np.sqrt(c2)

   #===========================================================================
   def _func_Mach_number(self, full_or_base):
      """
      Extract the Mach number in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      return self._func_velocity_magnitude(full_or_base) / \
             self._func_sound_speed(full_or_base)

   #===========================================================================
   def _func_enthalpy_specific(self, full_or_base):
      """
      Extract the specific enthalpy in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      u = self._func_internal_energy_specific(full_or_base)
      pdv = self._func_pressure(full_or_base) / \
            self._func_density(full_or_base)
      return u + pdv

   #===========================================================================
   def _func_enthalpy_density(self, full_or_base):
      """
      Extract the enthalpy density in the desired mode.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      u = self._func_internal_energy_density(full_or_base)
      pdv = self._func_pressure(full_or_base)
      return u + pdv

   #===========================================================================
   def _func_convective_growth_rate(self, full_or_base):
      """
      Extract the convective growth rate in the desired mode.

      NOTE: The convective growth rate is the imaginary part of the Brunt
      Vaisala frequency.

      NOTE: It is not immediately clear what the meaning of the "base state"
      value is for this quantity.  The square of the Brunt-Vaisala frequency is
      proportional to the dot product of the gravitational acceleration and the
      specific entropy gradient.  I choose the most direct definition, which is
      that the "base state" value of the square of the Brunt-Vaisala frequency
      is proportional to the dot product of the gravitational acceleration and
      the gradient of the base state specific entropy.

      Arguments:
         full_or_base (string) : name of the mode to use
      """

      g = self.params.gravity(self.x, self.y, self.z)
      s = self._func_entropy_specific(full_or_base)
      grad_s_x = self._x_derivative(s)
      grad_s_y = self._y_derivative(s)
      grad_s_z = self._z_derivative(s)
      g_dot_grad_s = g[...,0] * grad_s_x + \
                     g[...,1] * grad_s_y + \
                     g[...,2] * grad_s_z
      kk = (self.params.gamma - 1.0) / self.params.gamma
      BV_squared = kk * g_dot_grad_s
      BV_squared[BV_squared > 0] = 0.0
      return np.sqrt(-BV_squared)

   #===========================================================================

   # NOTE : The method commented out had issues getting the bindings correct.
   # I have to think more about how to do this properly.  For now, I just have
   # to duplicate the list of "known" variables.
   ## This is the list of known quantities and the corresponding dictionary
   ## pointing to the associated functions.  Since we are processing strings
   ## anyway, instead of having short, "variable-like" names, we use longer,
   ## more-descriptive names.
   #_variable_functions = \
   #   {"density" : self._func_density,
   #    "x momentum" : self._func_x_momentum,
   #    "y momentum" : self._func_y_momentum,
   #    "z momentum" : self._func_z_momentum,
   #    "momentum magnitude" : self._func_momentum_magnitude,
   #    "x velocity" : self._func_x_velocity,
   #    "y velocity" : self._func_y_velocity,
   #    "z velocity" : self._func_z_velocity,
   #    "velocity magnitude" : self._func_velocity_magnitude,
   #    "x vorticity" : self._func_x_vorticity,
   #    "y vorticity" : self._func_y_vorticity,
   #    "z vorticity" : self._func_z_vorticity,
   #    "vorticity magnitude" : self._func_vorticity_magnitude,
   #    "specific x vorticity" : self._func_x_specific_vorticity,
   #    "specific y vorticity" : self._func_y_specific_vorticity,
   #    "specific z vorticity" : self._func_z_specific_vorticity,
   #    "specific vorticity magnitude" : self._func_specific_vorticity_magnitude,
   #    "kinetic energy density" : self._func_kinetic_energy_density,
   #    "specific kinetic energy" : self._func_kinetic_energy_specific,
   #    "total energy density" : self._func_total_energy_density,
   #    "specific total energy" : self._func_total_energy_specific,
   #    "internal energy density" : self._func_internal_energy_density,
   #    "specific internal energy" : self._func_internal_energy_specific,
   #    "pressure" : self._func_pressure,
   #    "specific entropy" : self._func_entropy_specific,
   #    "entropy density" : self._func_entropy_density,
   #    "sound speed" : self._func_sound_speed,
   #    "Mach number" : self._func_Mach_number,
   #    "specific enthalpy" : self._func_enthalpy_specific,
   #    "enthalpy density" : self._func_enthalpy_density,
   #    "convective growth rate" : self._func_convective_growth_rate}
   #known_variables = tuple(_variable_functions.keys())
   known_variables = ("density",
       "x momentum",
       "y momentum",
       "z momentum",
       "momentum magnitude",
       "x velocity",
       "y velocity",
       "z velocity",
       "velocity magnitude",
       "x vorticity",
       "y vorticity",
       "z vorticity",
       "vorticity magnitude",
       "specific x vorticity",
       "specific y vorticity",
       "specific z vorticity",
       "specific vorticity magnitude",
       "kinetic energy density",
       "specific kinetic energy",
       "total energy density",
       "specific total energy",
       "internal energy density",
       "specific internal energy",
       "pressure",
       "specific entropy",
       "entropy density",
       "sound speed",
       "Mach number",
       "specific enthalpy",
       "enthalpy density",
       "convective growth rate")

   # This is the list of known modes.  Since we are processing strings anyway,
   # instead of having short, "variable-like" names, we use longer,
   # more-descriptive names.
   _mode_codes = \
      {"full state" : 0,
       "base state" : 1,
       "perturbation" : 2,
       "absolute perturbation" : 3,
       "contrast" : 4,
       "absolute contrast" : 5,
       "mean state" : -1,
       "residual" : -2,
       "absolute residual" : -3,
       "fractional residual" : -4,
       "absolute fractional residual" : -5}
   known_modes = tuple(_mode_codes.keys())

# End class SimulationState
#==============================================================================

