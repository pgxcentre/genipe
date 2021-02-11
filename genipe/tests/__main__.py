"""Main for running tests."""


import sys
import unittest

from . import test_suite


result = unittest.TextTestRunner(verbosity=1).run(test_suite)
sys.exit(not result.wasSuccessful())
