
# This file is part of genipe.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import os
import logging
import unittest
from tempfile import TemporaryDirectory

from ..tools import impute2_merger


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["TestImpute2Merger"]


class TestImpute2Merger(unittest.TestCase):

    @staticmethod
    def clean_logging_handlers():
        handlers = list(logging.root.handlers)
        for handler in handlers:
            logging.root.removeHandler(handler)

    def setUp(self):
        """Setup the tests."""
        # Creating the temporary directory
        self.output_dir = TemporaryDirectory(prefix="genipe_test_")

        # Creating small input files
        self.filenames = [
            os.path.join(self.output_dir.name, "input_1.impute2"),
            os.path.join(self.output_dir.name, "input_2.impute2"),
            os.path.join(self.output_dir.name, "input_3.impute2"),
            os.path.join(self.output_dir.name, "input_4.impute2"),
            os.path.join(self.output_dir.name, "input_5.impute2"),
            os.path.join(self.output_dir.name, "input_6.impute2"),
            os.path.join(self.output_dir.name, "input_7.impute2"),
        ]

        # The content of the files
        file_content = [
            "--- rs12345 1231415 A G 1 0 0 0.988 0.002 0 0 0.997 0.003",
            "1 rs23456 3214569 T C 0.869 0.130 0 0.903 0.095 0.002 0 0 1",
            "--- rs23457 3214570 T TC 0.869 0.130 0 0 1 0 0 0 1",
            "--- rs23457 3214570 T TC 0.869 0.130 0 0 1 0 0 0 1",
            "--- rs23457 3214570 T TC 0.869 0.130 0 0 1 0 0 0 1",
            "--- . 3214570 T TC 0.869 0.130 0 0 1 0 0 0 1",
            "--- . 3214570 T TC 0.869 0.130 0 0.869 0.130 0 0.869 0.130 0",
        ]

        # The summary content
        summary_content = [
            "-Output file\n"
            " --1 type 0 SNPs\n"
            " --0 type 1 SNPs\n"
            " --1 type 2 SNPs\n"
            " --1 type 3 SNPs\n"
            " --1 total SNPs"
        ] * len(self.filenames)

        # The SNP-wise information content
        info_content = [
            "--- rs12345 1231415 A G 0.006 0.359 0.987 0 -1 -1 -1",
            "1 rs23456 3214569 T C 0.082 0.362 0.866 0 -1 -1 -1",
            "--- rs23457 3214570 T TC 0.126 0.299 0.832 0 -1 -1 -1",
            "--- rs23457 3214570 T TC 0.060 0.300 0.909 0 -1 -1 -1",
            "--- rs23457 3214570 T TC 0.084 0.203 0.854 0 -1 -1 -1",
            "--- . 3214570 T TC 0.371 0.339 0.619 0 -1 -1 -1",
            "--- . 3214570 T TC 0.174 0.589 0.831 0 -1 -1 -1",
        ]

        # Creating the files
        zipped = zip(self.filenames, file_content, summary_content,
                     info_content)
        for filename, content, s_content, i_content in zipped:
            with open(filename, "w") as o_file:
                print(content, file=o_file)
            with open(filename + "_summary", "w") as o_file:
                print(s_content, file=o_file)
            with open(filename + "_info", "w") as o_file:
                print("snp_id", "rs_id", "position", "a0", "a1", "exp_freq_a1",
                      "info", "certainty", "type", "info_type0",
                      "concord_type0", "r2_type0", file=o_file)
                print(i_content, file=o_file)

        # Executing the script for the first time
        prefix1 = os.path.join(self.output_dir.name, "genipe_results_1")
        args = [
            "--chr", "1",
            "--probability", "0.9",
            "--completion", "0.98",
            "--prefix", prefix1,
            "--impute2",
        ]
        args += self.filenames
        impute2_merger.main(args=args)

        # Cleaning the handlers
        TestImpute2Merger.clean_logging_handlers()

        # Executing the script for the second time with different values
        prefix2 = os.path.join(self.output_dir.name, "genipe_results_2")
        args = [
            "--chr", "1",
            "--probability", "0.8",
            "--completion", "0.98",
            "--prefix", prefix2,
            "--info", "0.21",
            "--impute2",
        ]
        args += self.filenames
        impute2_merger.main(args=args)

        # Cleaning the handlers
        TestImpute2Merger.clean_logging_handlers()

        # Executing the script for the second time with different values
        prefix3 = os.path.join(self.output_dir.name, "genipe_results_3")
        args = [
            "--chr", "1",
            "--probability", "0.9",
            "--completion", "0.6",
            "--prefix", prefix3,
            "--info", "0.3",
            "--impute2",
        ]
        args += self.filenames
        impute2_merger.main(args=args)

        # Cleaning the handlers
        TestImpute2Merger.clean_logging_handlers()

        # Saving the prefixes
        self.prefixes = [prefix1, prefix2, prefix3]

    def tearDown(self):
        """Finishes the test."""
        # Deleting the output directory
        self.output_dir.cleanup()

    def test_check_output_files(self):
        """Checks the presence of all the output files."""
        suffixes = [".alleles", ".completion_rates", ".good_sites", ".impute2",
                    ".imputed_sites", ".log", ".maf", ".map", ".impute2_info"]
        for suffix in suffixes:
            for prefix in self.prefixes:
                self.assertTrue(os.path.isfile(prefix + suffix))

    def test_check_alleles(self):
        """Checks the '.alleles' file."""
        expected = (
            "name\ta1\ta2\n"
            "rs12345\tA\tG\n"
            "rs23456\tT\tC\n"
            "rs23457\tT\tTC\n"
            "rs23457_1\tT\tTC\n"
            "rs23457_2\tT\tTC\n"
            "1:3214570\tT\tTC\n"
            "1:3214570_1\tT\tTC\n"
        )

        # Checking the files
        for prefix in self.prefixes:
            observed = None
            with open(prefix + ".alleles", "r") as i_file:
                observed = i_file.read()
            self.assertEqual(expected, observed)

    def test_completion_rates(self):
        """Checks the '.completion_rates' file."""
        all_expected = []
        all_expected.append((
            "name\tnb_missing\tcompletion_rate\n"
            "rs12345\t0\t{}\n"
            "rs23456\t1\t{}\n"
            "rs23457\t1\t{}\n"
            "rs23457_1\t1\t{}\n"
            "rs23457_2\t1\t{}\n"
            "1:3214570\t1\t{}\n"
            "1:3214570_1\t3\t{}\n"
        ).format(1, 2/3, 2/3, 2/3, 2/3, 2/3, 0))
        all_expected.append((
            "name\tnb_missing\tcompletion_rate\n"
            "rs12345\t0\t{}\n"
            "rs23456\t0\t{}\n"
            "rs23457\t0\t{}\n"
            "rs23457_1\t0\t{}\n"
            "rs23457_2\t0\t{}\n"
            "1:3214570\t0\t{}\n"
            "1:3214570_1\t0\t{}\n"
        ).format(1, 1, 1, 1, 1, 1, 1))
        all_expected.append((
            "name\tnb_missing\tcompletion_rate\n"
            "rs12345\t0\t{}\n"
            "rs23456\t1\t{}\n"
            "rs23457\t1\t{}\n"
            "rs23457_1\t1\t{}\n"
            "rs23457_2\t1\t{}\n"
            "1:3214570\t1\t{}\n"
            "1:3214570_1\t3\t{}\n"
        ).format(1, 2/3, 2/3, 2/3, 2/3, 2/3, 0))

        # There should be enough expected values
        if len(self.prefixes) != len(all_expected):
            self.fail("Wrong number of expected values...")

        # Checking the files
        for prefix, expected in zip(self.prefixes, all_expected):
            observed = None
            with open(prefix + ".completion_rates", "r") as i_file:
                observed = i_file.read()

            # Splitting lines
            expected = expected.splitlines()
            observed = observed.splitlines()

            # Should have the same number of lines
            self.assertEqual(len(expected), len(observed))

            # Comparing the results
            for i, (e, o) in enumerate(zip(expected, observed)):
                if i == 0:
                    self.assertEqual(e, o)
                    continue

                # Splitting
                e = e.split("\t")
                o = o.split("\t")

                # Should be the same length
                self.assertEqual(len(e), len(o))

                # Values should be the same
                for j in range(len(e)):
                    if j == 0:
                        self.assertEqual(e[j], o[j])
                    elif j == 1:
                        self.assertEqual(int(e[j]), int(o[j]))
                    elif j == 2:
                        self.assertAlmostEqual(float(e[j]), float(o[j]))
                    else:
                        self.fail("Wrong number of values")

    def test_good_sites(self):
        """Checks the '.good_sites' file."""
        all_expected = []
        all_expected.append(
            "rs12345\n"
        )
        all_expected.append(
            "rs12345\n"
            "rs23456\n"
            "rs23457\n"
            "rs23457_1\n"
            "1:3214570\n"
            "1:3214570_1\n"
        )
        all_expected.append(
            "rs12345\n"
            "rs23456\n"
            "rs23457_1\n"
            "1:3214570\n"
        )

        # There should be enough expected values
        if len(self.prefixes) != len(all_expected):
            self.fail("Wrong number of expected values...")

        for prefix, expected in zip(self.prefixes, all_expected):
            observed = None
            with open(prefix + ".good_sites", "r") as i_file:
                observed = i_file.read()
            self.assertEqual(expected, observed)

    def test_impute2(self):
        """Checks the '.impute2' file."""
        expected = (
            "1 rs12345 1231415 A G 1 0 0 0.988 0.002 0 0 0.997 0.003\n"
            "1 rs23456 3214569 T C 0.869 0.130 0 0.903 0.095 0.002 0 0 1\n"
            "1 rs23457 3214570 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 rs23457_1 3214570 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 rs23457_2 3214570 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 1:3214570 3214570 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 1:3214570_1 3214570 T TC 0.869 0.130 0 0.869 0.130 0 0.869 "
            "0.130 0\n"
        )

        # Checking the files
        for prefix in self.prefixes:
            observed = None
            with open(prefix + ".impute2", "r") as i_file:
                observed = i_file.read()
            self.assertEqual(expected, observed)

    def test_impute2_info(self):
        """Checks the '.impute2_info' file."""
        expected = (
            "chr\tname\tposition\ta0\ta1\texp_freq_a1\tinfo\tcertainty\ttype\t"
            "info_type0\tconcord_type0\tr2_type0\n"
            "1\trs12345\t1231415\tA\tG\t0.006\t0.359\t0.987\t0\t-1\t-1\t-1\n"
            "1\trs23456\t3214569\tT\tC\t0.082\t0.362\t0.866\t0\t-1\t-1\t-1\n"
            "1\trs23457\t3214570\tT\tTC\t0.126\t0.299\t0.832\t0\t-1\t-1\t-1\n"
            "1\trs23457_1\t3214570\tT\tTC\t0.060\t0.300\t0.909\t0\t-1\t-1\t"
            "-1\n"
            "1\trs23457_2\t3214570\tT\tTC\t0.084\t0.203\t0.854\t0\t-1\t-1\t"
            "-1\n"
            "1\t1:3214570\t3214570\tT\tTC\t0.371\t0.339\t0.619\t0\t-1\t-1\t"
            "-1\n"
            "1\t1:3214570_1\t3214570\tT\tTC\t0.174\t0.589\t0.831\t0\t-1\t-1\t"
            "-1\n"
        )

        # Checking the files
        for prefix in self.prefixes:
            observed = None
            with open(prefix + ".impute2_info", "r") as i_file:
                observed = i_file.read()
            self.assertEqual(expected, observed)

    def test_imputed_sites(self):
        """Checks the '.imputed_sites' file."""
        expected = (
            "rs12345\n"
            "rs23457\n"
            "rs23457_1\n"
            "rs23457_2\n"
            "1:3214570\n"
            "1:3214570_1\n"
        )

        # Checking the files
        for prefix in self.prefixes:
            observed = None
            with open(prefix + ".imputed_sites", "r") as i_file:
                observed = i_file.read()
            self.assertEqual(expected, observed)

    def test_maf(self):
        """Checks the '.maf' file."""
        all_expected = []
        all_expected.append((
            "name\tmajor\tminor\tmaf\n"
            "rs12345\tA\tG\t{}\n"
            "rs23456\tT\tC\t{}\n"
            "rs23457\tTC\tT\t{}\n"
            "rs23457_1\tTC\tT\t{}\n"
            "rs23457_2\tTC\tT\t{}\n"
            "1:3214570\tTC\tT\t{}\n"
            "1:3214570_1\tT\tTC\t{}\n"
        ).format(1/6, 0.5, 1/4, 1/4, 1/4, 1/4, "NA"))
        all_expected.append((
            "name\tmajor\tminor\tmaf\n"
            "rs12345\tA\tG\t{}\n"
            "rs23456\tT\tC\t{}\n"
            "rs23457\tT\tTC\t{}\n"
            "rs23457_1\tT\tTC\t{}\n"
            "rs23457_2\tT\tTC\t{}\n"
            "1:3214570\tT\tTC\t{}\n"
            "1:3214570_1\tT\tTC\t{}\n"
        ).format(1/6, 2/6, 0.5, 0.5, 0.5, 0.5, 0))
        all_expected.append((
            "name\tmajor\tminor\tmaf\n"
            "rs12345\tA\tG\t{}\n"
            "rs23456\tT\tC\t{}\n"
            "rs23457\tTC\tT\t{}\n"
            "rs23457_1\tTC\tT\t{}\n"
            "rs23457_2\tTC\tT\t{}\n"
            "1:3214570\tTC\tT\t{}\n"
            "1:3214570_1\tT\tTC\t{}\n"
        ).format(1/6, 0.5, 1/4, 1/4, 1/4, 1/4, "NA"))

        # There should be enough expected values
        if len(self.prefixes) != len(all_expected):
            self.fail("Wrong number of expected values...")

        for prefix, expected in zip(self.prefixes, all_expected):
            observed = None
            with open(prefix + ".maf", "r") as i_file:
                observed = i_file.read()

            # Splitting lines
            expected = expected.splitlines()
            observed = observed.splitlines()

            # Should have the same number of lines
            self.assertEqual(len(expected), len(observed))

            # Comparing the results
            for i, (e, o) in enumerate(zip(expected, observed)):
                if i == 0:
                    self.assertEqual(e, o)
                    continue

                # Splitting
                e = e.split("\t")
                o = o.split("\t")

                # Should be the same length
                self.assertEqual(len(e), len(o))

                # Values should be the same
                self.assertEqual(e[:3], o[:3])

                for j in range(len(e)):
                    if j < 3:
                        self.assertEqual(e[j], o[j])
                    elif j == 3:
                        if e[j] != "NA":
                            self.assertAlmostEqual(float(e[j]), float(o[j]))
                        else:
                            self.assertEqual(e[j], o[j])
                    else:
                        self.fail("Wrong number of values")

    def test_map(self):
        """Checks the '.map' file."""
        expected = (
            "1\trs12345\t0\t1231415\n"
            "1\trs23456\t0\t3214569\n"
            "1\trs23457\t0\t3214570\n"
            "1\trs23457_1\t0\t3214570\n"
            "1\trs23457_2\t0\t3214570\n"
            "1\t1:3214570\t0\t3214570\n"
            "1\t1:3214570_1\t0\t3214570\n"
        )

        # Checking the files
        for prefix in self.prefixes:
            observed = None
            with open(prefix + ".map", "r") as i_file:
                observed = i_file.read()
            self.assertEqual(expected, observed)
