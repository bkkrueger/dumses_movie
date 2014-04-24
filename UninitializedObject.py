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

