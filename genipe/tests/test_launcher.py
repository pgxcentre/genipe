
# This file is part of genipe.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import os
import shutil
import unittest
from tempfile import TemporaryDirectory

from ..task.launcher import _check_output_files


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["TestLauncher"]


class TestLauncher(unittest.TestCase):

    def setUp(self):
        """Setup the tests."""
        # Creating the temporary directory
        self.output_dir = TemporaryDirectory(prefix="genipe_test_")

    def tearDown(self):
        """Finishes the test."""
        # Deleting the output directory
        self.output_dir.cleanup()

    def test_check_output_files(self):
        """Tests the '_check_output_files' function."""
        # First a list of files were they are all available
        filenames = [
            os.path.join(self.output_dir.name, "file_{}.txt".format(i))
            for i in range(1, 11)
        ]
        for filename in filenames:
            with open(filename, "w"):
                pass
        self.assertTrue(_check_output_files(filenames, "dummy_task_id"))

        # Adding an 'impute2' file
        filenames.append(os.path.join(self.output_dir.name, "test.impute2"))
        with open(filenames[-1], "w"):
            pass
        self.assertTrue(_check_output_files(filenames, "dummy_task_id"))

        # 'Compressing' the file
        shutil.move(filenames[-1], filenames[-1] + ".gz")
        self.assertTrue(_check_output_files(filenames, "dummy_task_id"))

        # Deleting a file
        os.remove(filenames[0])
        self.assertFalse(_check_output_files(filenames, "dummy_task_id"))

    @unittest.skip("Test not implemented")
    def test_check_missing_impute2(self):
        """Tests the '_check_output_files' for missing impute2 file."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_check_impute2_file(self):
        """Tests the '_check_impute2_file' function."""
        self.fail("Test not implemented")
