
# This file is part of gwip.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import os
import unittest
from tempfile import TemporaryDirectory

from ..pipeline import *


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


class TestMainPipeline(unittest.TestCase):

    def setUp(self):
        """Setup the tests."""
        # Creating the temporary directory
        self.output_dir = TemporaryDirectory(prefix="gwip_test_")

    def tearDown(self):
        """Finishes the test."""
        # Deleting the output directory
        self.output_dir.cleanup()

    def test_file_sorter(self):
        """Tests the 'file_sorter' function."""
        filenames = ["chr1.1000_100000.impute2", "chr2.10000_1002300.impute2",
                     "chr1.1_100.impute2", "/some/path/chr1.1_10.impute2",
                     "chr1.100000_2000000.impute2.some_extension"]
        expected_filenames = ["/some/path/chr1.1_10.impute2",
                              "chr1.1_100.impute2", "chr1.1000_100000.impute2",
                              "chr1.100000_2000000.impute2.some_extension",
                              "chr2.10000_1002300.impute2"]
        expected_results = [(1, 1000, 100000), (2, 10000, 1002300),
                            (1, 1, 100), (1, 1, 10), (1, 100000, 2000000)]

        # Trying the function
        for filename, expected in zip(filenames, expected_results):
            self.assertEqual(expected, file_sorter(filename))

        # Trying the sort function
        filenames.sort(key=file_sorter)
        self.assertEqual(filenames, expected_filenames)

    def test_get_chromosome_length(self):
        """Tests the 'get_chromosome_length' function."""
        # Tests that all chromosome are there
        expected_chrom = {str(i) for i in range(1, 23)} | {"X"}
        chrom_length = get_chromosome_length(self.output_dir.name)
        self.assertTrue(len(expected_chrom - set(chrom_length.keys())) == 0)
        self.assertTrue(os.path.isfile(os.path.join(self.output_dir.name,
                                                    "chromosome_lengths.txt")))

        # Tests that we correctly read the file
        expected_chrom = {"1": 249250621, "10": 135534747, "11": 135006516,
                          "12": 133851895, "13": 115169878, "14": 107349540,
                          "15": 102531392, "16": 90354753, "17": 81195210,
                          "18": 78077248, "19": 59128983, "2": 243199373,
                          "20": 63025520, "21": 48129895, "22": 51304566,
                          "3": 198022430, "4": 191154276, "5": 180915260,
                          "6": 171115067, "7": 159138663, "8": 146364022,
                          "9": 141213431}

        # Writing the file
        chrom_filename = os.path.join(self.output_dir.name,
                                      "chromosome_lengths.txt")
        with open(chrom_filename, "w") as o_file:
            for k, v in expected_chrom.items():
                print(k, v, sep="\t", file=o_file)

        # Comparing what we got
        chrom_length = get_chromosome_length(self.output_dir.name)
        self.assertEqual(expected_chrom, chrom_length)

        # Removing some chromosomes from the file
        del expected_chrom["9"]
        del expected_chrom["12"]
        with open(chrom_filename, "w") as o_file:
            for k, v in expected_chrom.items():
                print(k, v, sep="\t", file=o_file)

        # Tests that an exception is raised if there is a missing chromosome
        with self.assertRaises(ProgramError) as e:
            get_chromosome_length(self.output_dir.name)
        self.assertEqual("missing chromosomes: 12, 9", e.exception.message)

    @unittest.skip("Test not implemented")
    def test_read_preamble(self):
        """Tests the 'read_preamble' function."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_get_cross_validation_results(self):
        """Tests the 'get_cross_validation_results' function."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_gather_imputation_stats(self):
        """Tests the 'gather_imputation_stats' function."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_gather_execution_time(self):
        """Tests the 'gather_execution_time' function."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_check_args(self):
        """Tests the 'check_args' function."""
        self.fail("Test not implemented")
