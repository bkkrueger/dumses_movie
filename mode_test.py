#==============================================================================

class PlotMode(object):
   """
   The mode to use when plotting a quantity.
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
         raise StandardError() # TODO

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
               raise StandardError() # TODO
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
                  raise StandardError() # TODO
               a, f = detail.split(self.transform)
               a = a.strip()
               f = f.strip()
               if a == "":
                  self.absolute = False
               elif a == "absolute":
                  self.absolute = True
               else:
                  raise StandardError() # TODO
               if f == "" or f == "from base":
                  self.reference = "base"
               elif f == "from mean":
                  self.reference = "mean"
               else:
                  raise StandardError() # TODO

   #===========================================================================
   def __init_from_properties(self, d, a, t, r):
      """
      Initialize from mode properties.
      """

      if a and (d != 2 or t not in ["perturbation", "contrast"]):
         raise StandardError() # TODO

      if d == 1:
         allowed = {"full" : ["none"],
                    "base" : ["perturbation", "contrast"]}
      elif d == 2:
         allowed = {"full" : ["none"],
                    "base" : ["none", "perturbation", "contrast"],
                    "mean" : ["none", "perturbation", "contrast"]}
      else:
         raise StandardError() # TODO

      if t in allowed.get(r, []):
         self.dimension = d
         self.absolute = a
         self.transform = t
         self.reference = r
      else:
         raise StandardError() # TODO

   #===========================================================================
   def __repr__(self):
      """
      Supply a detailed representation of the mode.
      """

      if self.dimension == 1:
         if self.absolute:
            raise StandardError() # TODO
         if self.reference == "full" and self.transform == "none":
            return "profile: full state"
         else:
            if (self.reference == "base" and
                  self.transform in ["perturbation", "contrast"]):
               return "profile: " + self.transform
            else:
               raise StandardError() # TODO
      elif self.dimension == 2:
         if self.transform == "none":
            if self.absolute:
               raise StandardError() # TODO
            if self.reference in ["full", "base", "mean"]:
               return "pseudocolor: " + self.reference + " state"
            else:
               raise StandardError() # TODO
         elif self.transform in ["perturbation", "contrast"]:
            if self.reference in ["base", "mean"]:
               if self.absolute:
                  detail_list = ["absolute"]
               else:
                  detail_list = []
               detail_list.extend([self.transform, "from", self.reference])
               return "pseudocolor: " + " ".join(detail_list)
            else:
               raise StandardError() # TODO
         else:
            raise StandardError() # TODO
      else:
         raise StandardError() # TODO

#==============================================================================

def parse_mode(mode):
   """
   profile: (full state|perturbation|contrast)
   pseudocolor: ([absolute ](perturbation|contrast)[ from (base|mean)]|
                 (full|base|mean) state)
   """

   # Normalize spacing
   mode = ' '.join(mode.split())

   dimension, detail = mode.split(":")
   dimension = dimension.strip()
   detail = detail.strip()

   if dimension == "profile":
      dimension = 1
      absolute = False
      if detail == "full state":
         transform = "none"
         reference = "full"
      elif detail in ["perturbation", "contrast"]:
         transform = detail
         reference = "base"
      else:
         raise StandardError() # TODO
   elif dimension == "pseudocolor":
      dimension = 2
      if detail[-5:] == "state" and detail[:-6] in ["full", "base", "mean"]:
         absolute = False
         transform = "none"
         reference = detail[:-6]
      else:
         if "perturbation" in detail:
            transform = "perturbation"
         elif "contrast" in detail:
            transform = "contrast"
         else:
            raise StandardError() # TODO
         a, f = detail.split(transform)
         a = a.strip()
         f = f.strip()
         if a == "":
            absolute = False
         elif a == "absolute":
            absolute = True
         else:
            raise StandardError() # TODO
         if f == "" or f == "from base":
            reference = "base"
         elif f == "from mean":
            reference = "mean"
         else:
            raise StandardError() # TODO

   return (dimension, absolute, transform, reference)

#==============================================================================

def mode_name(dimension, absolute, transform, reference):

   if dimension == 1:
      if absolute:
         raise StandardError() # TODO
      if reference == "full" and transform == "none":
         return "profile: full state"
      else:
         if reference == "base" and transform in ["perturbation", "contrast"]:
            return "profile: " + transform
         else:
            raise StandardError() # TODO
   elif dimension == 2:
      if transform == "none":
         if absolute:
            raise StandardError() # TODO
         if reference in ["full", "base", "mean"]:
            return "pseudocolor: " + reference + " state"
         else:
            raise StandardError() # TODO
      elif transform in ["perturbation", "contrast"]:
         if reference in ["base", "mean"]:
            if absolute:
               detail_list = ["absolute"]
            else:
               detail_list = []
            detail_list.extend([transform, "from", reference])
            return "pseudocolor: " + " ".join(detail_list)
         else:
            raise StandardError() # TODO
      else:
         raise StandardError() # TODO
   else:
      raise StandardError() # TODO

#==============================================================================

n = 1
for dim in ["pseudocolor", "profile"]:
   if dim == "pseudocolor":
      detail_list = []
      for det in ["transform", "state"]:
         if det == "transform":
            for a in ["", "absolute"]:
               for t in ["perturbation", "contrast"]:
                  for f in ["", "from base", "from mean"]:
                     detail_list.append(" ".join((a,t,f)))
         else:
            detail_list.extend(["full state", "base state", "mean state"])
   else:
      detail_list = ["full state", "perturbation", "contrast"]
   for detail in detail_list:
      mode = dim + ": " + detail
      mode = ' '.join(mode.split())
      pm = PlotMode(mode)
      d = pm.dimension
      a = pm.absolute
      t = pm.transform
      r = pm.reference
      print "{n:3} | {d:1} {a:1} {t:1} {r:1} | {m}".format(n=n, a=a, d=d,
            t=t[0], r=r[0], m=mode)
      n += 1

print ""
print ""

n = 1
for d in [1, 2]:
   for a in [False, True]:
      for t in ["none", "perturbation", "contrast"]:
         for r in ["full", "base", "mean"]:
            f_str = "{n:3} | {d:1} {a:1} {t:1} {r:1} | {m}"
            try:
               pm = PlotMode(d, a, t, r)
               print f_str.format(n=n, a=a, d=d, t=t[0], r=r[0], m=pm)
               n += 1
            except StandardError:
               print f_str.format(n="   ",a=a, d=d, t=t[0], r=r[0],
                     m="unable to construct name")

