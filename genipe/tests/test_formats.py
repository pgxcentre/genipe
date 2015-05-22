
# This file is part of genipe.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import unittest
from tempfile import TemporaryDirectory

import numpy as np

from ..formats import *
from ..error import ProgramError


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["TestFormats"]


class TestFormats(unittest.TestCase):

    def test_matrix_from_line(self):
        """Tests the 'matrix_from_line' function."""
        # The IMPUTE2 line
        input_line = ("1 marker_1 1 A AT 0 0 1 1 0 0 0 1 0 0.9 0.033 0.067 "
                      "0 0 0")

        # The expected results
        expected_info = ["1", "marker_1", "1", "A", "AT"]
        expected_geno = np.array([[0, 0, 1], [1, 0, 0], [0, 1, 0],
                                  [0.9, 0.033, 0.067], [0, 0, 0]], dtype=float)

        # The observed results
        observed_info, observed_geno = matrix_from_line(input_line.split(" "))

        # The marker information should be the same
        self.assertEqual(expected_info, observed_info)

        # The probability matrix should have the same shape
        self.assertEqual(expected_geno.shape, observed_geno.shape)

        # The values in the probability matrix should be the same as expected
        self.assertTrue(np.allclose(expected_geno, observed_geno))

        # An invalid line should raise an exception
        with self.assertRaises(ValueError) as cm:
            matrix_from_line(input_line.split(" ")[:-1])
        error_m = "total size of new array must be unchanged"
        self.assertEqual(error_m, str(cm.exception))

    def test_get_good_probs(self):
        """Tests the 'get_good_probs' function."""
        # The probability matrix
        input_data = np.array([[0, 0, 1], [0, 1, 0], [1, 0, 0], [0.9, 0.1, 0],
                               [0, 0.9, 0.1], [0, 0.1, 0.9], [0.8, 0.2, 0],
                               [0.1, 0.8, 0.1], [0, 0.2, 0.8], [0.5, 0.1, 0.4],
                               [0, 0, 0], [0.05, 0.05, 0.9]], dtype=float)

        # The expected results (for a probability of 0.9)
        expected = np.array([True] * 6 + [False] * 5 + [True], dtype=bool)
        observed = get_good_probs(input_data, 0.9)
        self.assertTrue((expected == observed).all())

        # The expected results (for a probability of 0.8)
        expected = np.array([True] * 9 + [False] * 2 + [True], dtype=bool)
        observed = get_good_probs(input_data, 0.8)
        self.assertTrue((expected == observed).all())

        # The expected results (for a probability of 0)
        expected = np.ones(12, dtype=bool)
        observed = get_good_probs(input_data, 0)
        self.assertTrue((expected == observed).all())

    def test_maf_from_probs(self):
        """Tests the 'maf_from_probs' function."""
        # The probability matrix (10 samples)
        probs = np.array(
            [[0.9, 0.1, 0.0],   # AA, male   (A)
             [0.1, 0.9, 0.0],   # AB, female (AB)
             [0.0, 0.1, 0.9],   # BB, male   (B)*
             [0.9, 0.1, 0.0],   # AA, female (AA)
             [0.0, 0.1, 0.9],   # BB, male   (B)
             [0.9, 0.1, 0.0],   # AA, male   (A)
             [0.9, 0.1, 0.0],   # AA, male   (A)
             [0.1, 0.9, 0.0],   # AB, female (AB)*
             [0.1, 0.9, 0.0],   # AB, female (AB)
             [0.9, 0.1, 0.0]],  # AA, male   (A)
            dtype=float,
        )

        # The reversed probability matrix
        r_probs = np.array([i[::-1] for i in probs], dtype=float)

        # The gender array (10 samples, 6 males, 4 females)
        gender = np.array([1, 2, 1, 2, 1, 1, 1, 2, 2, 1], dtype=int)
        unknown_gender = np.array([1, 2, 0, 2, 1, 1, 1, 0, 2, 1], dtype=int)

        # Checking without gender
        expected_maf = 7 / 20
        expected_minor = "B"
        expected_major = "A"
        observed_results = maf_from_probs(probs, "A", "B")
        observed_maf, observed_minor, observed_major = observed_results
        self.assertAlmostEqual(expected_maf, observed_maf)
        self.assertEqual(expected_minor, observed_minor)
        self.assertEqual(expected_major, observed_major)

        # Reversing the matrix should give the same maf, but different major
        # and minor allele
        expected_maf = 7 / 20
        expected_minor = "A"
        expected_major = "B"
        observed_results = maf_from_probs(r_probs, "A", "B")
        observed_maf, observed_minor, observed_major = observed_results
        self.assertAlmostEqual(expected_maf, observed_maf)
        self.assertEqual(expected_minor, observed_minor)
        self.assertEqual(expected_major, observed_major)

        # Checking with gender
        expected_maf = 5 / 14
        expected_minor = "B"
        expected_major = "A"
        observed_results = maf_from_probs(probs, "A", "B", gender)
        observed_maf, observed_minor, observed_major = observed_results
        self.assertAlmostEqual(expected_maf, observed_maf)
        self.assertEqual(expected_minor, observed_minor)
        self.assertEqual(expected_major, observed_major)

        # Reversing the matrix should give the same maf, but different major
        # and minor allele
        expected_maf = 5 / 14
        expected_minor = "A"
        expected_major = "B"
        observed_results = maf_from_probs(r_probs, "A", "B", gender)
        observed_maf, observed_minor, observed_major = observed_results
        self.assertAlmostEqual(expected_maf, observed_maf)
        self.assertEqual(expected_minor, observed_minor)
        self.assertEqual(expected_major, observed_major)

        # Checking with gender, unknown samples
        expected_maf = 3 / 11
        expected_minor = "B"
        expected_major = "A"
        observed_results = maf_from_probs(probs, "A", "B", unknown_gender)
        observed_maf, observed_minor, observed_major = observed_results
        self.assertAlmostEqual(expected_maf, observed_maf)
        self.assertEqual(expected_minor, observed_minor)
        self.assertEqual(expected_major, observed_major)

        # Reversing the matrix should give the same maf, but different major
        # and minor allele
        expected_maf = 3 / 11
        expected_minor = "A"
        expected_major = "B"
        observed_results = maf_from_probs(r_probs, "A", "B", unknown_gender)
        observed_maf, observed_minor, observed_major = observed_results
        self.assertAlmostEqual(expected_maf, observed_maf)
        self.assertEqual(expected_minor, observed_minor)
        self.assertEqual(expected_major, observed_major)

        # Checking with gender (only males)
        tmp_gender = np.ones(10, dtype=int)
        tmp_probs = np.array(
            [[0.9, 0.1, 0.0],   # AA, male   (A)
             [0.0, 0.1, 0.9],   # BB, male   (B)
             [0.9, 0.1, 0.0],   # AA, female (A)
             [0.9, 0.1, 0.0],   # AA, female (A)
             [0.9, 0.1, 0.0],   # AA, male   (A)
             [0.9, 0.1, 0.0],   # AA, male   (A)
             [0.9, 0.1, 0.0],   # AA, male   (A)
             [0.0, 0.1, 0.9],   # BB, male   (B)
             [0.0, 0.1, 0.9],   # BB, male   (B)
             [0.9, 0.1, 0.0]],  # AA, male   (A)
            dtype=float,
        )
        tmp_r_probs = np.array([i[::-1] for i in tmp_probs], dtype=float)
        expected_maf = 3 / 10
        expected_minor = "B"
        expected_major = "A"
        observed_results = maf_from_probs(tmp_probs, "A", "B", tmp_gender)
        observed_maf, observed_minor, observed_major = observed_results
        self.assertAlmostEqual(expected_maf, observed_maf)
        self.assertEqual(expected_minor, observed_minor)
        self.assertEqual(expected_major, observed_major)

        # Reversing the matrix should give the same maf, but different major
        # and minor allele
        expected_maf = 3 / 10
        expected_minor = "A"
        expected_major = "B"
        observed_results = maf_from_probs(tmp_r_probs, "A", "B", tmp_gender)
        observed_maf, observed_minor, observed_major = observed_results
        self.assertAlmostEqual(expected_maf, observed_maf)
        self.assertEqual(expected_minor, observed_minor)
        self.assertEqual(expected_major, observed_major)

        # Checking with gender (only females)
        tmp_gender += 1
        expected_maf = 7 / 20
        expected_minor = "B"
        expected_major = "A"
        observed_results = maf_from_probs(probs, "A", "B", tmp_gender)
        observed_maf, observed_minor, observed_major = observed_results
        self.assertAlmostEqual(expected_maf, observed_maf)
        self.assertEqual(expected_minor, observed_minor)
        self.assertEqual(expected_major, observed_major)

        # Reversing the matrix should give the same maf, but different major
        # and minor allele
        expected_maf = 7 / 20
        expected_minor = "A"
        expected_major = "B"
        observed_results = maf_from_probs(r_probs, "A", "B", tmp_gender)
        observed_maf, observed_minor, observed_major = observed_results
        self.assertAlmostEqual(expected_maf, observed_maf)
        self.assertEqual(expected_minor, observed_minor)
        self.assertEqual(expected_major, observed_major)

        # Checking empty probabilities (without gender)
        tmp_probs = np.array([], dtype=float)
        expected_maf = "NA"
        expected_minor = "B"
        expected_major = "A"
        observed_results = maf_from_probs(tmp_probs, "A", "B")
        observed_maf, observed_minor, observed_major = observed_results
        self.assertEqual(expected_maf, observed_maf)
        self.assertEqual(expected_minor, observed_minor)
        self.assertEqual(expected_major, observed_major)

        # Checking empty probabilities (with gender)
        expected_maf = "NA"
        expected_minor = "B"
        expected_major = "A"
        observed_results = maf_from_probs(tmp_probs, "A", "B", gender)
        observed_maf, observed_minor, observed_major = observed_results
        self.assertEqual(expected_maf, observed_maf)
        self.assertEqual(expected_minor, observed_minor)
        self.assertEqual(expected_major, observed_major)

        # Checking all unknown gender
        all_unknown = np.zeros(10, dtype=int)
        expected_maf = "NA"
        expected_minor = "B"
        expected_major = "A"
        observed_results = maf_from_probs(probs, "A", "B", all_unknown)
        observed_maf, observed_minor, observed_major = observed_results
        self.assertEqual(expected_maf, observed_maf)
        self.assertEqual(expected_minor, observed_minor)
        self.assertEqual(expected_major, observed_major)

        # An heterozygous male (using gender) should raise an exception
        with self.assertRaises(ProgramError) as cm:
            maf_from_probs(probs, "A", "B", np.ones(10, dtype=int), "marker_1")
        self.assertEqual("marker_1: heterozygous male present",
                         str(cm.exception))

    def test_dosage_from_probs(self):
        """Tests the 'dosage_from_probs' function."""
        # The probability matrices
        homo_probs = np.array([1, 0.9, 0.98, 0.99, 0, 0.1, 0.5, 0.3],
                              dtype=float)
        het_probs = np.array([0, 0.1, 0.01, 0.01, 0.1, 0.2, 0.4, 0.6],
                             dtype=float)

        # The expected results (for different scale)
        for scale in [1, 2]:
            expected = np.array([(1 + (0 / 2)) * scale,
                                 (0.9 + (0.1 / 2)) * scale,
                                 (0.98 + (0.01 / 2)) * scale,
                                 (0.99 + (0.01 / 2)) * scale,
                                 (0 + (0.1 / 2)) * scale,
                                 (0.1 + (0.2 / 2)) * scale,
                                 (0.5 + (0.4 / 2)) * scale,
                                 (0.3 + (0.6 / 2)) * scale], dtype=float)
            observed = dosage_from_probs(homo_probs, het_probs, scale)
            self.assertTrue(np.allclose(expected, observed))

        # This one should cause an exception
        with self.assertRaises(ValueError) as cm:
            dosage_from_probs(homo_probs, het_probs[:-1])
        error_m = ("operands could not be broadcast together with shapes "
                   "(8,) (7,)")
        self.assertEqual(error_m, str(cm.exception).strip())

    def test_hard_calls_from_probs(self):
        """Tests the 'hard_calls_from_probs' function."""
        prob_matrix = np.array([
            [0.01, 0.01, 0.98],
            [0.01, 0.98, 0.01],
            [0.98, 0.01, 0.01],
        ])

        # Checking with single alleles
        a1 = "A"
        a2 = "B"
        expected = np.array(["B B", "A B", "A A"])
        observed = hard_calls_from_probs(a1, a2, prob_matrix)
        self.assertEqual(expected.shape, observed.shape)
        self.assertTrue((expected == observed).all())

        # Checking with indels
        a1 = "ABB"
        a2 = "AB"
        expected = np.array(["AB AB", "ABB AB", "ABB ABB"])
        observed = hard_calls_from_probs(a1, a2, prob_matrix)
        self.assertEqual(expected.shape, observed.shape)
        self.assertTrue((expected == observed).all())
