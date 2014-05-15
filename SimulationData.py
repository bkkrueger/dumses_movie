"""
Specification of simulation descriptions

The SimulationError object describes errors that arise from classes
defined here.  This allows calling routines to distinguish errors from
this body of code from other, more generic errors.

The SimulationInput object describes the setup of the simulation.  This
includes parameters such as the upstream conditions used to construct
the initial state, the adiabatic index for the equation of state, and
the description of the sources.

The SimulationVariable object describes a variable that can be
extracted from the data.  This includes a function to compute the
variable, the default name of the variable, and how the variable
relates to zero ("lower bound", "center", or None).

The SimulationState object packages the data from the current state of
a simulation.  It also provides a mechanism for extracting a number of
different quantities in a convenient way.
"""

import numpy as np

from dumpy_v05.data.rd_dumses import DumsesData

#------------------------------------------------------------------------------
#==============================================================================
class SimulationError(StandardError):
   """
   Signals that a simulation data object encountered an error.
   """
   pass

# End class SimulationError
#==============================================================================
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#==============================================================================
class SimulationInput(object):
   """
   A description of the simulation.

   This class describes the simulation based on the values extracted
   from the inputs file.  It includes methods to initialize from a file
   and to construct the source functions.

   Attributes:
      gamma (float) : adiabatic index
      shape (string) : name of the shape function for the sources
      layer_width (float) : width parameter for source layer
      dens_up (float) : scaling value for density
      pres_up (float) : scaling value for pressure (unsettable)
      csnd_up (float) : scaling value for sound speed
      mach_up (float) : scaling value for Mach number
      K_heat (float) : dimensionless heating strength factor
      K_grav (float) : dimensionless gravity strength factor
   """

   #===========================================================================
   def __repr__(self):
      """
      Supply a detailed representation of the class.
      """
      string = "   \n".join(("Simulation inputs:",
         "gamma = {0}".format(self.gamma),
         "source layer:",
         "   shape = {0}".format(self.shape),
         "   width = {0}".format(self.layer_width),
         "   Kheat = {0}".format(self.K_heat),
         "   Kgrav = {0}".format(self.K_grav),
         "upstream conditions:",
         "   dens = {0}".format(self.dens_up),
         "   csnd = {0}".format(self.csnd_up),
         "   Mach = {0}".format(self.mach_up)
         ))
      return string

   #===========================================================================
   def __str__(self):
      """
      Supply a simple representation of the class.
      """
      return "SimulationInput"

   #===========================================================================
   def __init__(self, *args, **kwargs):
      """
      Construct the class by selecting the appropriate initializer.
      """

      if len(args) == 0:
         self.clear()
      elif len(args) == 1:
         self.construct(*args)
      else:
         msg = "__init__ takes at most 1 positional argument ({n} given)"
         raise TypeError(msg.format(n=len(args)))

   #===========================================================================
   def clear(self):
      """
      Flush the data to give an uninitialized instance.
      """

      self.gamma = None

      self.dens_up = None
      self.csnd_up = None
      self.mach_up = None

      self.layer_width = None
      self.shape = None

      self.K_heat = None
      self.K_grav = None

   #===========================================================================
   def construct(self, inputs_file_name):
      """
      Initialize the class from an inputs file.

      Arguments:
         inputs_file_name (string) : name of the inputs file to parse

      Returns:
         None

      Exceptions:
         whatever arises from functions called by this routine
      """

      # Data to extract from the parameters file, using defaults as appropriate
      # to match DUMSES code
      #    By default all values should be cast to strings.  This is because
      # the input is read from a Fortran name list file, and Fortran allows
      # scientific notation with a "d" instead of an "e" to specify
      # double-precision, so I run a replace on all values, which requires that
      # the values be strings.  The exception is values that require special
      # handling anyway, so that they are not run through the default process.
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

      # Store quantities now that they've been constructed/verified
      self.gamma = inputs_dict["gamma"]

      self.csnd_up = inputs_dict["csnd_up"]
      self.mach_up = inputs_dict["mach_up"]
      self.dens_up = inputs_dict["dens_up"]

      self.layer_width = inputs_dict["layer_limit"]
      self.shape = inputs_dict["layer_shape"]

      self.K_heat = inputs_dict["Kheat"]
      self.K_grav = inputs_dict["Kgrav"]

   #===========================================================================
   @property
   def pres_up(self):
      if None in [self.dens_up, self.csnd_up, self.gamma]:
         msg = " ".join(("Cannot compute upstream pressure scale without",
            "defining the adiabatic index and the upstream scales for density",
            "and sound speed."))
         raise SimulationError(msg)
      return self.dens_up * self.csnd_up**2 / self.gamma

   #===========================================================================
   def shape_function(self, x, y, z):
      """
      Compute the shape function.

      Arguments:
         x (numpy.ndarray) : the x coordinates
         y (numpy.ndarray) : the y coordinates
         z (numpy.ndarray) : the z coordinates

      Returns:
         an array specifying the shape function on the grid

      Exceptions:
         SimulationError flagging an uninitialized object
      """

      if None in [self.layer_width, self.shape]:
         msg = " ".join(("Cannot compute shape function without defining the",
            "shape name and width."))

      XX = np.abs(x / self.layer_width)

      if self.shape == "trapezoid":
         f = np.clip(2 * (1.0 - XX), 0.0, 1.0)
      elif self.shape == "square":
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
      Compute the gravitational acceleration

      Arguments:
         x (numpy.ndarray) : the x coordinates
         y (numpy.ndarray) : the y coordinates
         z (numpy.ndarray) : the z coordinates

      Returns:
         an array specifying the gravitational acceleration on the grid

      Exceptions:
         SimulationError flagging an uninitialized object
      """

      if None in [self.K_grav, self.csnd_up, self.layer_width]:
         msg = " ".join(("Cannot compute gravity without defining gravity",
            "scaling constant, upstream sound speed scaling,",
            "and shape width."))
         raise SimulationError(msg)

      grav_coef = - self.K_grav *  self.csnd_up**2 / self.layer_width
      grav = np.zeros((x.shape[0], y.shape[1], z.shape[2], 3), dtype=float)
      grav[...,0] = grav_coef * self.shape_function(x,y,z)

      return grav

   #===========================================================================
   def heating(self, x, y, z, density):
      """
      Compute the heating rate (total energy density per time source).

      Arguments:
         x (numpy.ndarray) : the x coordinates
         y (numpy.ndarray) : the y coordinates
         z (numpy.ndarray) : the z coordinates
         density (numpy.ndarray) : the density

      Returns:
         an array specifying the heating rate on the grid

      Exceptions:
         SimulationError flagging an uninitialized object
      """

      if None in [self.K_heat, self.mach_up, self.csnd_up, self.layer_width,
            self.gamma]:
         msg = " ".join(("Cannot compute heating without defining heating",
            "scaling constant, shape width, adiabatic index, and upstream",
            "scaling of sound speed and Mach number."))
         raise SimulationError(msg)

      heat_coef = (self.K_heat * (self.mach_up * self.csnd_up**3) /
            (self.layer_width * self.gamma))
      return (heat_coef * density * self.shape_function(x, y, z))

# End class SimulationInput
#==============================================================================
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#==============================================================================
class SimulationVariable(object):
   """
   A variable that may be computed from the simulation data.

   The SimulationState "knows" a collection of variables.  To "know" a
   variable, it must have one or more names (SimulationState maintains
   a complete list, including aliases, although although
   SimulationVariable will hold a default name) and two attributes
   (held by SimulationVariable): a function to compute the variable and
   a parameter describing how the variable relates to zero (generally
   one of "lower bound", "center", or None).  SimulationState is the
   object that knows the variables, and SimulationVariable is only a
   convenience object to package that information together.  Thus the
   SimulationVariable object is mostly dummy information that will be
   filled in by the SimulationState when it generates its collection of
   known variables.

   Attributes:
      name (str) : the name of the variable
      zero (str) : how the variable relates to zero
   """

   #===========================================================================
   def __repr__(self):
      """
      Supply a detailed representation of the class.
      """
      return "SimulationVariable (" + self.name + ")"

   #===========================================================================
   def __str__(self):
      """
      Supply a simple representation of the class.
      """
      return self.name

   #===========================================================================
   def __init__(self, name, function, zero):
      """
      Initialize the variable.
      """
      self.name = name
      self.compute = function
      self.zero = zero

# End class SimulationVariable
#==============================================================================
#------------------------------------------------------------------------------



#------------------------------------------------------------------------------
#==============================================================================
class SimulationState(object):
   """
   The current state of a simulation.

   This class stores the current state of a simulation, which includes
   the full and base states of mass density, momentum density, total
   energy density, as well as the time and the coordinate axes.  This
   much is similar to the DumsesData class, to which this is closely
   related (and a method exists to construct the SimulationState from a
   DumsesData because of this close relation).  The primary difference
   is that the SimulationState will also have a routine that accepts
   the name of a variable and the name of a state (full state or base
   state), which will use private methods to compute the desired
   variable in the specified state.  Having a separate class also
   allows me to tweak the internals to fit my needs without changing
   the standard DumsesData class (e.g. the names of variables, or if I
   find time to update this to be parallel so that I can visualize
   large simulations that would fill the memory on a single core).

   An important note on DumsesData: It is not actually necessary to supply a
   DumsesData (think duck typing).  Instead the object passed in for a
   DumsesData requires the following:
   -- t: A floating-point value for the time
   -- x: A 1D array of cell-centered x-coordinates (length Nx)
   -- y: A 1D array of cell-centered y-coordinates (length Ny)
   -- z: A 1D array of cell-centered z-coordinates (length Nz)
   -- rho: A 3D array for mass density (Nx x Ny x Nz)
   -- rhou: A 4D array for momentum density vector (Nx x Ny x Nz x 3)
   -- E: A 3D array for total energy density (Nx x Ny x Nz)
   -- rho0: A 1D array for base state mass density (length Nx)
   -- rhou0: A 1D array for base state momentum density (length Ny)
   -- E0: A 1D array for base state total energy density (length Nz)

   Attributes:
      __initialized (bool) : is the object initialized
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
      known_variables (list) : list of known variables
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
   def __str__(self):
      """
      Supply a simple representation of the class.
      """
      if self.__initialized:
         string = "Simulation data at t = {0:13.6e}".format(self.t)
      else:
         string = "Simulation data (uninitialized)"
      return string

   #===========================================================================
   def __init__(self, *args, **kwargs):
      """
      Select the appropriate initialization method.
      """

      if len(args) == 0:
         self.clear()
      elif len(args) == 2:
         self.construct(*args)
      else:
         msg = "__init__ takes either 0 or 2 positional arguments ({n} given)"
         raise TypeError(msg.format(n=len(args)))

      # Generate the collection of known variables with appropriate function
      # bindings.  Some variables will have multiple aliases, so they will be
      # constructed once with a default name and saved under multiple keys.
      self.known_variables = {}

      var = SimulationVariable("density", self._func_density, "lower bound")
      names = ["density", "mass density"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("x momentum", self._func_x_momentum, "center")
      names = ["x momentum", "x momentum density"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("y momentum", self._func_y_momentum, "center")
      names = ["y momentum", "y momentum density"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("z momentum", self._func_z_momentum, "center")
      names = ["z momentum", "z momentum density"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("momentum magnitude",
            self._func_momentum_magnitude, "lower bound")
      names = ["momentum magnitude", "momentum density magnitude"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("x velocity", self._func_x_velocity, "center")
      names = ["x velocity", "specific x momentum"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("y velocity", self._func_y_velocity, "center")
      names = ["y velocity", "specific y momentum"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("z velocity", self._func_z_velocity, "center")
      names = ["z velocity", "specific z momentum"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable( "velocity magnitude",
            self._func_velocity_magnitude, "lower bound")
      names = ["velocity magnitude", "specific momentum magnitude"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("x vorticity", self._func_x_vorticity, "center")
      names = ["x vorticity"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("y vorticity", self._func_y_vorticity, "center")
      names = ["y vorticity"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("z vorticity", self._func_z_vorticity, "center")
      names = ["z vorticity"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("vorticity magnitude",
            self._func_vorticity_magnitude, "lower bound")
      names = ["vorticity magnitude"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("specific x vorticity",
            self._func_x_specific_vorticity, "center")
      names = ["specific x vorticity"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("specific y vorticity",
            self._func_y_specific_vorticity, "center")
      names = ["specific y vorticity"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("specific z vorticity",
            self._func_z_specific_vorticity, "center")
      names = ["specific z vorticity"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("specific vorticity magnitude",
            self._func_specific_vorticity_magnitude, "lower bound")
      names = ["specific vorticity magnitude"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("kinetic energy density",
            self._func_kinetic_energy_density, "lower bound")
      names = ["kinetic energy density"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("specific kinetic energy",
            self._func_kinetic_energy_specific, "lower bound")
      names = ["specific kinetic energy"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("total energy density",
            self._func_total_energy_density, "lower bound")
      names = ["total energy density"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("specific total energy",
            self._func_total_energy_specific, "lower bound")
      names = ["specific total energy"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("internal energy density",
            self._func_internal_energy_density, "lower bound")
      names = ["internal energy density"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("specific internal energy",
            self._func_internal_energy_specific, "lower bound")
      names = ["specific internal energy"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("pressure", self._func_pressure, "lower bound")
      names = ["pressure"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("entropy", self._func_entropy_specific, None)
      names = ["entropy", "specific entropy"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("entropy density",
            self._func_entropy_density, None)
      names = ["entropy density"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("sound speed",
            self._func_sound_speed, "lower bound")
      names = ["sound speed"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("Mach number",
            self._func_Mach_number, "lower bound")
      names = ["Mach number"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("specific enthalpy",
            self._func_enthalpy_specific, "lower bound")
      names = ["specific enthalpy"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("enthalpy density",
            self._func_enthalpy_density, "lower bound")
      names = ["enthalpy density"]
      self.known_variables.update(dict.fromkeys(names, var))

      var = SimulationVariable("convective growth rate",
            self._func_convective_growth_rate, "lower bound")
      names = ["convective growth rate", "Brunt-Vaisala frequency"]
      self.known_variables.update(dict.fromkeys(names, var))

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

   #===========================================================================
   def construct(self, dumpy, input_parameters):
      """
      Construct the SimulationState from the supplied information.

      Arguments:
         dumpy (DumsesData) : used as a basis to construct this object
         input_parameters (SimulationInputs) : input file parameters

      Returns:
         None

      Exceptions:
         whatever arises from functions called by this routine
      """

      # Extract values
      params = input_parameters

      t = dumpy.time

      x = dumpy.x[:,None,None]
      y = dumpy.y[None,:,None]
      z = dumpy.z[None,None,:]

      dens = dumpy.rho
      momv = dumpy.rhou
      Ener = dumpy.E

      temp = dumpy.rho0[:,None,None]
      temp = temp.repeat(len(dumpy.y), axis=1)
      temp = temp.repeat(len(dumpy.z), axis=2)
      dens0 = temp

      temp = dumpy.rhou0[:,None,None]
      temp = temp.repeat(len(dumpy.y), axis=1)
      temp = temp.repeat(len(dumpy.z), axis=2)
      momv0 = temp

      temp = dumpy.E0[:,None,None]
      temp = temp.repeat(len(dumpy.y), axis=1)
      temp = temp.repeat(len(dumpy.z), axis=2)
      Ener0 = temp

      # Save values
      self.params = params

      self.t = t

      self.x = x
      self.y = y
      self.z = z

      self._dens = dens
      self._momv = momv
      self._Ener = Ener

      self._dens0 = dens0
      self._momv0 = momv0
      self._Ener0 = Ener0

      self.__initialized = True

   #===========================================================================
   def extract(self, variable_name, use_full_state=True):
      """
      Extract the desired variable.

      Arguments:
         variable_name (string) : name of the variable to be extracted
         use_full_state (bool) : full (True) or base (False) state

      Returns:
         array holding the data for the requested variable

      Exceptions:
         SimulationError for uninitialized object or unknown variable
      """

      if not self.__initialized:
         msg = " ".join(("Cannot extract data from an uninitialized",
            "SimulationState object."))
         raise SimulationError(msg)

      # Get the appropriate variable function
      try:
         func = self.known_variables[variable_name].compute
      except KeyError:
         msg = 'Unrecognized variable: "{v}"'.format(v=variable_name)
         raise SimulationError(msg)

      var = func(use_full_state)

      return var

   #===========================================================================
   def _x_derivative(self, q):
      """
      The derivative of quantity q with respect to x.

      NOTE: I assume that the boundary conditions in the x direction
      are zero-gradient.  I could write this to be flexible by having
      the SimulationInputs class parse the bval_in and bval_out entries
      in the inputs file, but if a user-defined boundary value is used
      then there is no way to know what to do.  Thus I will for now
      simply be lazy and assume the boundary conditions that I use in
      all of my own simulations and leave it to other users to read the
      documentation and adjust accordingly for their own uses.

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

      NOTE: I assume that the boundary conditions in the y direction
      are periodic.  I could write this to be flexible by having the
      SimulationInputs class parse the bval_in and bval_out entries in
      the inputs file, but if a user-defined boundary value is used
      then there is no way to know what to do.  Thus I will for now
      simply be lazy and assume the boundary conditions that I use in
      all of my own simulations and leave it to other users to read the
      documentation and adjust accordingly for their own uses.

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

      NOTE: I assume that the boundary conditions in the z direction
      are periodic.  I could write this to be flexible by having the
      SimulationInputs class parse the bval_in and bval_out entries in
      the inputs file, but if a user-defined boundary value is used
      then there is no way to know what to do.  Thus I will for now
      simply be lazy and assume the boundary conditions that I use in
      all of my own simulations and leave it to other users to read the
      documentation and adjust accordingly for their own uses.

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
   def _func_density(self, use_full_state):
      """
      Extract the density.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      if use_full_state:
         return self._dens
      else:
         return self._dens0

   #===========================================================================
   def _func_x_momentum(self, use_full_state):
      """
      Extract the x-momentum.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      if use_full_state:
         return self._momv[:,:,:,0]
      else:
         return self._momv0

   #===========================================================================
   def _func_y_momentum(self, use_full_state):
      """
      Extract the y-momentum.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      if use_full_state:
         return self._momv[:,:,:,1]
      else:
         return np.zeros_like(self._momv0)

   #===========================================================================
   def _func_z_momentum(self, use_full_state):
      """
      Extract the z-momentum.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      if use_full_state:
         return self._momv[:,:,:,2]
      else:
         return np.zeros_like(self._momv0)

   #===========================================================================
   def _func_momentum_magnitude(self, use_full_state):
      """
      Extract the magnitude of the momentum.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      if use_full_state:
         return np.sqrt(np.sum(self._momv**2, axis=-1))
      else:
         return np.abs(self._momv0)

   #===========================================================================
   def _func_x_velocity(self, use_full_state):
      """
      Extract the x-velocity.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      return self._func_x_momentum(use_full_state) / \
             self._func_density(use_full_state)

   #===========================================================================
   def _func_y_velocity(self, use_full_state):
      """
      Extract the y-velocity.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      return self._func_y_momentum(use_full_state) / \
             self._func_density(use_full_state)

   #===========================================================================
   def _func_z_velocity(self, use_full_state):
      """
      Extract the z-velocity.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      return self._func_z_momentum(use_full_state) / \
             self._func_density(use_full_state)

   #===========================================================================
   def _func_velocity_magnitude(self, use_full_state):
      """
      Extract the magnitude of the velocity.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      return self._func_momentum_magnitude(use_full_state) / \
             self._func_density(use_full_state)

   #===========================================================================
   def _func_x_vorticity(self, use_full_state):
      """
      Extract the x-vorticity.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      vely = self._func_y_velocity(use_full_state)
      velz = self._func_z_velocity(use_full_state)
      vrtx = self._y_derivative(velz) - self._z_derivative(vely)
      return vrtx

   #===========================================================================
   def _func_y_vorticity(self, use_full_state):
      """
      Extract the y-vorticity.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      velx = self._func_x_velocity(use_full_state)
      velz = self._func_z_velocity(use_full_state)
      vrty = self._z_derivative(velx) - self._x_derivative(velz)
      return vrty

   #===========================================================================
   def _func_z_vorticity(self, use_full_state):
      """
      Extract the z-vorticity.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      velx = self._func_x_velocity(use_full_state)
      vely = self._func_y_velocity(use_full_state)
      vrtz = self._x_derivative(vely) - self._y_derivative(velx)
      return vrtz

   #===========================================================================
   def _func_vorticity_magnitude(self, use_full_state):
      """
      Extract the magnitude of the vorticity.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      vrtx = self._func_x_vorticity(use_full_state)
      vrty = self._func_y_vorticity(use_full_state)
      vrtz = self._func_z_vorticity(use_full_state)
      magnitude = np.sqrt(vrtx**2 + vrty**2 + vrtz**2)
      return magnitude

   #===========================================================================
   def _func_x_specific_vorticity(self, use_full_state):
      """
      Extract the specific x-vorticity.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      vrtx = self._func_x_vorticity(use_full_state)
      dens = self._func_density(use_full_state)
      return vrtx / dens

   #===========================================================================
   def _func_y_specific_vorticity(self, use_full_state):
      """
      Extract the specific y-vorticity.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      vrty = self._func_y_vorticity(use_full_state)
      dens = self._func_density(use_full_state)
      return vrty / dens

   #===========================================================================
   def _func_z_specific_vorticity(self, use_full_state):
      """
      Extract the specific z-vorticity.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      vrtz = self._func_z_vorticity(use_full_state)
      dens = self._func_density(use_full_state)
      return vrtz / dens

   #===========================================================================
   def _func_specific_vorticity_magnitude(self, use_full_state):
      """
      Extract the magnitude of the specific vorticity.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      vrtx = self._func_x_vorticity(use_full_state)
      vrty = self._func_y_vorticity(use_full_state)
      vrtz = self._func_z_vorticity(use_full_state)
      magnitude = np.sqrt(vrtx**2 + vrty**2 + vrtz**2)
      dens = self._func_density(use_full_state)
      return magnitude / dens

   #===========================================================================
   def _func_kinetic_energy_density(self, use_full_state):
      """
      Extract the kinetic energy density.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      momentum = self._func_momentum_magnitude(use_full_state)
      kinetic = 0.5 * momentum**2 / self._func_density(use_full_state)
      return kinetic

   #===========================================================================
   def _func_kinetic_energy_specific(self, use_full_state):
      """
      Extract the specific kinetic energy.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      momentum = self._func_momentum_magnitude(use_full_state)
      kinetic = 0.5 * (momentum / self._func_density(use_full_state))**2
      return kinetic

   #===========================================================================
   def _func_total_energy_density(self, use_full_state):
      """
      Extract the total energy density.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      if use_full_state:
         return self._Ener
      else:
         return self._Ener0

   #===========================================================================
   def _func_total_energy_specific(self, use_full_state):
      """
      Extract the specific total energy.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      Ener = self._func_total_energy_density(use_full_state) / \
             self._func_density(use_full_state)

      return Ener

   #===========================================================================
   def _func_internal_energy_density(self, use_full_state):
      """
      Extract the internal energy density.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      return self._func_total_energy_density(use_full_state) - \
             self._func_kinetic_energy_density(use_full_state)

   #===========================================================================
   def _func_internal_energy_specific(self, use_full_state):
      """
      Extract the specific internal energy.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      return self._func_total_energy_specific(use_full_state) - \
             self._func_kinetic_energy_specific(use_full_state)

   #===========================================================================
   def _func_pressure(self, use_full_state):
      """
      Extract the pressure.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      return self._func_internal_energy_density(use_full_state) * \
             (self.params.gamma - 1.0)

   #===========================================================================
   def _func_entropy_specific(self, use_full_state):
      """
      Extract the specific entropy.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      PP = self._func_pressure(use_full_state) / self.params.pres_up
      dd = self._func_density(use_full_state) / self.params.dens_up
      Y = self.params.gamma
      return np.log(PP * dd**(-Y)) / (Y - 1.0)

   #===========================================================================
   def _func_entropy_density(self, use_full_state):
      """
      Extract the entropy density.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      return self._func_entropy_specific(use_full_state) * \
             self._func_density(use_full_state)

   #===========================================================================
   def _func_sound_speed(self, use_full_state):
      """
      Extract the speed of sound.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      c2 = self.params.gamma * self._func_pressure(use_full_state) / \
            self._func_density(use_full_state)
      return np.sqrt(c2)

   #===========================================================================
   def _func_Mach_number(self, use_full_state):
      """
      Extract the Mach number.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      return self._func_velocity_magnitude(use_full_state) / \
             self._func_sound_speed(use_full_state)

   #===========================================================================
   def _func_enthalpy_specific(self, use_full_state):
      """
      Extract the specific enthalpy.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      u = self._func_internal_energy_specific(use_full_state)
      pdv = self._func_pressure(use_full_state) / \
            self._func_density(use_full_state)
      return u + pdv

   #===========================================================================
   def _func_enthalpy_density(self, use_full_state):
      """
      Extract the enthalpy density.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      u = self._func_internal_energy_density(use_full_state)
      pdv = self._func_pressure(use_full_state)
      return u + pdv

   #===========================================================================
   def _func_convective_growth_rate(self, use_full_state):
      """
      Extract the convective growth rate.

      NOTE: The convective growth rate is the imaginary part of the
      Brunt-Vaisala frequency.

      NOTE: It is not immediately clear what the meaning of the "base
      state" value is for this quantity.  The square of the
      Brunt-Vaisala frequency is proportional to the dot product of the
      gravitational acceleration and the specific entropy gradient.  I
      choose the most direct definition, which is that the "base state"
      value of the square of the Brunt-Vaisala frequency is
      proportional to the dot product of the gravitational acceleration
      and the gradient of the base state specific entropy.

      Arguments:
         use_full_state (bool) : full (True) or base (False) state
      """

      g = self.params.gravity(self.x, self.y, self.z)
      s = self._func_entropy_specific(use_full_state)
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

# End class SimulationState
#==============================================================================
#------------------------------------------------------------------------------

