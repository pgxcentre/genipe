
# This file is part of gwip.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import unittest
from tempfile import TemporaryDirectory


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


class TestDB(unittest.TestCase):

    def setUp(self):
        """Setup the tests."""
        # Creating the temporary directory
        self.output_dir = TemporaryDirectory(prefix="gwip_test_")

    def tearDown(self):
        """Finishes the test."""
        # Deleting the output directory
        self.output_dir.cleanup()

    @unittest.skip("Test not implemented")
    def test_create_task_db(self):
        """Tests the 'create_task_db' function."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_create_db_connection(self):
        """Tests the '_create_db_connection' function."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_check_task_completion(self):
        """Tests the 'check_task_completion' function."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_create_task_entry(self):
        """Tests the 'create_task_entry' function."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_mark_task_completed(self):
        """Tests the 'mark_task_completed' function."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_mark_task_incomplete(self):
        """Tests the 'mark_task_incomplete' function."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_mark_drmaa_task_completed(self):
        """Tests the 'mark_drmaa_task_completed' function."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_get_task_runtime(self):
        """Tests the 'task_runtime' function."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_get_all_runtimes(db_name):
        """Tests the 'get_all_runtimes' function."""
        self.fail("Test not implemented")
