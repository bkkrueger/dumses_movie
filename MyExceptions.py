# TODO : Add custom exceptions here.  Some exceptions will be specific to a
#        certain, self-contained part of the code; those may go with their
#        associated code.  However most of my custom exceptions will be more
#        general-use throughout the code, so they will go here to be packaged
#        into a central location.

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

"""
An exception raised if an object is not initialized.

This module provides the specification of the UninitializedObjectError.  This
is a general-purpose error that I will use in many objects to state that an
object's method has been invoked when the object has not been initialized.

Attributes:
   UninitializedObjectError (class) : error when uninitialized object is called
"""

#==============================================================================

class UninitializedObjectError(Exception):
   """
   Signals that an object has not been initialized.
   """

   def __repr__(self):
      """
      Supply a detailed representation of the exception.
      """
      return "Object is uninitialized."

#==============================================================================

