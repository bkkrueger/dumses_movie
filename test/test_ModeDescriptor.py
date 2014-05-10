import unittest
import numpy as np
import sys
sys.path.append("..")
import site_setup
from Descriptors import ModeDescriptor, DescriptorError

class ModeDescriptorTest(unittest.TestCase):

   def test_construct(self):
      mm = ModeDescriptor()

      self.assertEqual(str(mm), "pseudocolor: full state",
            msg="Failure of default initialization.")

      for d in [1, 2, "junk"]:
         for a in [False, True, "junk"]:
            for t in ["none", "perturbation", "contrast", "junk"]:
               for r in ["full", "base", "mean", "junk"]:
                  pc = ["perturbation", "contrast"]
                  if (("junk" in [d, a, t, r]) or
                        (r == "full" and t != "none") or
                        (a and t not in pc) or
                        (a and d == 1) or
                        (d == 1 and r == "mean") or
                        (d == 1 and r == "base" and t not in pc)):
                     print "FAIL :", d, a, t, r
                     self.assertRaises(DescriptorError,
                           mm = ModeDescriptor(d, a, t, r))
                  else:
                     print "TEST :", d, a, t, r
                     mm = ModeDescriptor(d, a, t, r)

   def test_loop(self):

      for d in [1, 2]:
         for a in [False, True]:
            for t in ["none", "perturbation", "contrast"]:
               for r in ["full", "base", "mean"]:
                  if r == "full" and t != "none":
                     self.assertRaises(DescriptorError,
                           mm = ModeDescriptor(d, a, t, r))
                  elif a and t not in ["perturbation", "contrast"]:
                     self.assertRaises(DescriptorError,
                           mm = ModeDescriptor(d, a, t, r))
                  elif a and d != 2:
                     self.assertRaises(DescriptorError,
                           mm = ModeDescriptor(d, a, t, r))
                  else:
                     mm = ModeDescriptor(d, a, t, r)
                     self.assertEqual(d, mm.dimension,
                           msg="Failed to maintain dimension.")
                     self.assertEqual(a, mm.absolute,
                           msg="Failed to maintain absolute.")
                     self.assertEqual(t, mm.transform,
                           msg="Failed to maintain transform.")
                     self.assertEqual(r, mm.reference,
                           msg="Failed to maintain reference.")
                     mm2 = ModeDescriptor(str(mm))
                     self.assertEqual(d, mm2.dimension,
                           msg="Failed to maintain dimension.")
                     self.assertEqual(a, mm2.absolute,
                           msg="Failed to maintain absolute.")
                     self.assertEqual(t, mm2.transform,
                           msg="Failed to maintain transform.")
                     self.assertEqual(r, mm2.reference,
                           msg="Failed to maintain reference.")

if __name__ == "__main__":
   unittest.main()
