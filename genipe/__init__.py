
# This file is part of genipe.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import importlib

try:
    from .version import genipe_version as __version__
except ImportError:
    __version__ = None


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__credits__ = ["Louis-Philippe Lemieux Perreault", "Marc-Andre Legault"]
__license__ = "CC BY-NC 4.0"
__maintainer__ = "Louis-Philippe Lemieux Perreault"
__email__ = "louis-philippe.lemieux.perreault@statgen.org"
__status__ = "Production"


# The chromosomes
autosomes = range(1, 23)
chromosomes_23 = (23, 25)
chromosomes = tuple(autosomes) + chromosomes_23


# Checking once for 'pyfaidx' and 'matplotlib'
HAS_PYFAIDX = importlib.find_loader("pyfaidx") is not None
HAS_MATPLOTLIB = importlib.find_loader("matplotlib") is not None
HAS_DRMAA = importlib.find_loader("drmaa") is not None


def test(verbosity=1):
    """Executes all the tests for genipe.

    Args:
        verbosity (int): the verbosity level

    Just set ``verbosity`` to an integer higher than 1 to have more information
    about the tests.

    """
    import logging
    import unittest
    from .tests import test_suite

    # Disabling the INFO logging
    logging.disable(logging.INFO)

    unittest.TextTestRunner(verbosity=verbosity).run(test_suite)

    # Enabling back the logging
    logging.disable(logging.NOTSET)
