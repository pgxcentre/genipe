
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

from ..tools import impute2_extractor


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["TestImpute2Extractor"]


class TestImpute2Extractor(unittest.TestCase):

    @staticmethod
    def clean_logging_handlers():
        handlers = list(logging.root.handlers)
        for handler in handlers:
            logging.root.removeHandler(handler)

    def setUp(self):
        """Setup the tests."""
        # Creating the temporary directory
        self.output_dir = TemporaryDirectory(prefix="genipe_test_")

        # Creating the input files (impute2, map, maf, completion_rate and
        # impute2_info)
        filename = os.path.join(self.output_dir.name, "genipe.impute2")
        with open(filename, "w") as o_file:
            o_file.write((
                "1 rs12345 1231415 A G 1 0 0 0.988 0.002 0 0 0.997 0.003\n"
                "1 rs23456 3214569 T C 0.869 0.130 0 0.903 0.095 0.002 0 0 1\n"
                "1 rs23457 3214570 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
                "1 rs23457_1 3214571 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
                "1 rs23457_2 3214572 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
                "1 1:3214573 3214573 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
                "1 1:4214570_1 4214570 T TC 0.869 0.130 0 0.869 0.130 0 0.869 "
                "0.130 0\n"
            ))

        filename = os.path.join(self.output_dir.name, "genipe.impute2_info")
        with open(filename, "w") as o_file:
            o_file.write((
                "chr\tname\tposition\ta0\ta1\texp_freq_a1\tinfo\tcertainty\t"
                "type\tinfo_type0\tconcord_type0\tr2_type0\n"
                "1\trs12345\t1231415\tA\tG\t0.006\t0.359\t0.987\t0\t-1\t-1\t"
                "-1\n"
                "1\trs23456\t3214569\tT\tC\t0.082\t0.362\t0.866\t0\t-1\t-1\t"
                "-1\n"
                "1\trs23457\t3214570\tT\tTC\t0.126\t0.299\t0.832\t0\t-1\t-1\t"
                "-1\n"
                "1\trs23457_1\t3214571\tT\tTC\t0.060\t0.300\t0.909\t0\t-1\t"
                "-1\t-1\n"
                "1\trs23457_2\t3214572\tT\tTC\t0.084\t0.203\t0.854\t0\t-1\t"
                "-1\t-1\n"
                "1\t1:3214573\t3214573\tT\tTC\t0.371\t0.339\t0.619\t0\t-1\t"
                "-1\t-1\n"
                "1\t1:4214570_1\t4214570\tT\tTC\t0.174\t0.589\t0.831\t0\t-1\t"
                "-1\t-1\n"
            ))

        filename = os.path.join(self.output_dir.name, "genipe.map")
        with open(filename, "w") as o_file:
            o_file.write((
                "1\trs12345\t0\t1231415\n"
                "1\trs23456\t0\t3214569\n"
                "1\trs23457\t0\t3214570\n"
                "1\trs23457_1\t0\t3214571\n"
                "1\trs23457_2\t0\t3214572\n"
                "1\t1:3214573\t0\t3214573\n"
                "1\t1:4214570_1\t0\t4214570\n"
            ))

        filename = os.path.join(self.output_dir.name, "genipe.maf")
        with open(filename, "w") as o_file:
            o_file.write((
                "name\tmajor\tminor\tmaf\n"
                "rs12345\tA\tG\t{}\n"
                "rs23456\tT\tC\t{}\n"
                "rs23457\tTC\tT\t{}\n"
                "rs23457_1\tTC\tT\t{}\n"
                "rs23457_2\tTC\tT\t{}\n"
                "1:3214573\tTC\tT\t{}\n"
                "1:4214570_1\tT\tTC\t{}\n"
            ).format(1/6, 0.5, 1/4, 1/4, 1/4, 1/4, "NA"))

        filename = os.path.join(self.output_dir.name,
                                "genipe.completion_rates")
        with open(filename, "w") as o_file:
            o_file.write((
                "name\tnb_missing\tcompletion_rate\n"
                "rs12345\t0\t{}\n"
                "rs23456\t1\t{}\n"
                "rs23457\t1\t{}\n"
                "rs23457_1\t1\t{}\n"
                "rs23457_2\t1\t{}\n"
                "1:3214573\t1\t{}\n"
                "1:4214570_1\t3\t{}\n"
            ).format(1, 2/3, 2/3, 2/3, 2/3, 2/3, 0))

        filename = os.path.join(self.output_dir.name, "genipe.sample")
        with open(filename, "w") as o_file:
            o_file.write(
             "ID_1 ID_2 missing father mother sex plink_pheno\n"
             "0 0 0 D D D B\n"
             "f1 s1 0 0 0 0 -9\n"
             "f2 s2 0 0 0 0 -9\n"
             "f3 s3 0 0 0 0 -9\n"
            )

        self.common_args = [
            "--impute2", os.path.join(self.output_dir.name, "genipe.impute2"),
            "--format", "impute2", "dosage", "calls",
            "--out", os.path.join(self.output_dir.name, "results"),
        ]

    def tearDown(self):
        """Finishes the test."""
        # Deleting the output directory
        self.output_dir.cleanup()

    def test_extract(self):
        """Tests the extraction by marker name."""
        # Creating a file with markers to extract
        extract_filename = os.path.join(self.output_dir.name, "to_extract")
        with open(extract_filename, "w") as o_file:
            o_file.write("rs23456\nrs23457_2\n1:4214570_1\n")

        # Executing the script
        args = self.common_args + [
            "--extract", extract_filename,
        ]
        impute2_extractor.main(args=args)
        TestImpute2Extractor.clean_logging_handlers()

        # Testing we have the three output files
        template_name = os.path.join(self.output_dir.name, "results.{ext}")
        for suffix in ("impute2", "dosage", "calls"):
            self.assertTrue(os.path.isfile(template_name.format(ext=suffix)))

        # Checking the impute2 file
        expected = (
            "1 rs23456 3214569 T C 0.869 0.130 0 0.903 0.095 0.002 0 0 1\n"
            "1 rs23457_2 3214572 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 1:4214570_1 4214570 T TC 0.869 0.130 0 0.869 0.130 0 0.869 "
            "0.130 0\n"
        )
        observed = None
        with open(template_name.format(ext="impute2"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the dosage file
        expected = (
            "chrom\tpos\tname\tminor\tmajor\tf1/s1\tf2/s2\tf3/s3\n"
            "1\t3214569\trs23456\tC\tT\tnan\t0.099\t2.0\n"
            "1\t3214572\trs23457_2\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t4214570\t1:4214570_1\tTC\tT\tnan\tnan\tnan\n"
        )
        observed = None
        with open(template_name.format(ext="dosage"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the hard calls file
        expected = (
            "chrom\tname\tcm\tpos\tf1/s1\tf2/s2\tf3/s3\n"
            "1\trs23456\t0\t3214569\t0 0\tT T\tC C\n"
            "1\trs23457_2\t0\t3214572\t0 0\tT TC\tTC TC\n"
            "1\t1:4214570_1\t0\t4214570\t0 0\t0 0\t0 0\n"
        )
        observed = None
        with open(template_name.format(ext="calls"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (impute2_info)
        info_fn = template_name.format(ext="impute2_info")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "chr\tname\tposition\ta0\ta1\texp_freq_a1\tinfo\tcertainty\t"
            "type\tinfo_type0\tconcord_type0\tr2_type0\n"
            "1\trs23456\t3214569\tT\tC\t0.082\t0.362\t0.866\t0\t-1\t-1\t-1\n"
            "1\trs23457_2\t3214572\tT\tTC\t0.084\t0.203\t0.854\t0\t-1\t-1\t"
            "-1\n"
            "1\t1:4214570_1\t4214570\tT\tTC\t0.174\t0.589\t0.831\t0\t-1\t-1\t"
            "-1\n"
        )
        observed = None
        with open(info_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (map)
        map_fn = template_name.format(ext="map")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "1\trs23456\t0\t3214569\n"
            "1\trs23457_2\t0\t3214572\n"
            "1\t1:4214570_1\t0\t4214570\n"
        )
        observed = None
        with open(map_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion file (maf)
        maf_fn = template_name.format(ext="maf")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "name\tmajor\tminor\tmaf\n"
            "rs23456\tT\tC\t{}\n"
            "rs23457_2\tTC\tT\t{}\n"
            "1:4214570_1\tT\tTC\t{}\n"
        ).format(0.5, 1/4, "NA")
        observed = None
        with open(maf_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

    def test_genomic(self):
        """Tests the extraction by genomic location."""
        # Executing the script
        args = self.common_args + [
            "--genomic", "chr1:3214570-3214573",
        ]
        impute2_extractor.main(args=args)
        TestImpute2Extractor.clean_logging_handlers()

        # Testing we have the three output files
        template_name = os.path.join(self.output_dir.name, "results.{ext}")
        for suffix in ("impute2", "dosage", "calls"):
            self.assertTrue(os.path.isfile(template_name.format(ext=suffix)))

        # Checking the impute2 file
        expected = (
            "1 rs23457 3214570 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 rs23457_1 3214571 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 rs23457_2 3214572 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 1:3214573 3214573 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
        )
        observed = None
        with open(template_name.format(ext="impute2"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the dosage file
        expected = (
            "chrom\tpos\tname\tminor\tmajor\tf1/s1\tf2/s2\tf3/s3\n"
            "1\t3214570\trs23457\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t3214571\trs23457_1\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t3214572\trs23457_2\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t3214573\t1:3214573\tT\tTC\tnan\t1.0\t0.0\n"
        )
        observed = None
        with open(template_name.format(ext="dosage"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the hard calls file
        expected = (
            "chrom\tname\tcm\tpos\tf1/s1\tf2/s2\tf3/s3\n"
            "1\trs23457\t0\t3214570\t0 0\tT TC\tTC TC\n"
            "1\trs23457_1\t0\t3214571\t0 0\tT TC\tTC TC\n"
            "1\trs23457_2\t0\t3214572\t0 0\tT TC\tTC TC\n"
            "1\t1:3214573\t0\t3214573\t0 0\tT TC\tTC TC\n"
        )
        observed = None
        with open(template_name.format(ext="calls"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (impute2_info)
        info_fn = template_name.format(ext="impute2_info")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "chr\tname\tposition\ta0\ta1\texp_freq_a1\tinfo\tcertainty\t"
            "type\tinfo_type0\tconcord_type0\tr2_type0\n"
            "1\trs23457\t3214570\tT\tTC\t0.126\t0.299\t0.832\t0\t-1\t-1\t-1\n"
            "1\trs23457_1\t3214571\tT\tTC\t0.060\t0.300\t0.909\t0\t-1\t-1\t"
            "-1\n"
            "1\trs23457_2\t3214572\tT\tTC\t0.084\t0.203\t0.854\t0\t-1\t-1\t"
            "-1\n"
            "1\t1:3214573\t3214573\tT\tTC\t0.371\t0.339\t0.619\t0\t-1\t-1\t"
            "-1\n"
        )
        observed = None
        with open(info_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (map)
        map_fn = template_name.format(ext="map")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "1\trs23457\t0\t3214570\n"
            "1\trs23457_1\t0\t3214571\n"
            "1\trs23457_2\t0\t3214572\n"
            "1\t1:3214573\t0\t3214573\n"
        )
        observed = None
        with open(map_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion file (map)
        maf_fn = template_name.format(ext="maf")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "name\tmajor\tminor\tmaf\n"
            "rs23457\tTC\tT\t{}\n"
            "rs23457_1\tTC\tT\t{}\n"
            "rs23457_2\tTC\tT\t{}\n"
            "1:3214573\tTC\tT\t{}\n"
        ).format(1/4, 1/4, 1/4, 1/4)
        observed = None
        with open(maf_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

    def test_maf(self):
        """Tests the extraction by maf."""
        # Executing the script
        args = self.common_args + [
            "--maf", "0.25",
        ]
        impute2_extractor.main(args=args)
        TestImpute2Extractor.clean_logging_handlers()

        # Testing we have the three output files
        template_name = os.path.join(self.output_dir.name, "results.{ext}")
        for suffix in ("impute2", "dosage", "calls"):
            self.assertTrue(os.path.isfile(template_name.format(ext=suffix)))

        # Checking the impute2 file
        expected = (
            "1 rs23456 3214569 T C 0.869 0.130 0 0.903 0.095 0.002 0 0 1\n"
            "1 rs23457 3214570 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 rs23457_1 3214571 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 rs23457_2 3214572 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 1:3214573 3214573 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
        )
        observed = None
        with open(template_name.format(ext="impute2"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the dosage file
        expected = (
            "chrom\tpos\tname\tminor\tmajor\tf1/s1\tf2/s2\tf3/s3\n"
            "1\t3214569\trs23456\tC\tT\tnan\t0.099\t2.0\n"
            "1\t3214570\trs23457\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t3214571\trs23457_1\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t3214572\trs23457_2\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t3214573\t1:3214573\tT\tTC\tnan\t1.0\t0.0\n"
        )
        observed = None
        with open(template_name.format(ext="dosage"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the hard calls file
        expected = (
            "chrom\tname\tcm\tpos\tf1/s1\tf2/s2\tf3/s3\n"
            "1\trs23456\t0\t3214569\t0 0\tT T\tC C\n"
            "1\trs23457\t0\t3214570\t0 0\tT TC\tTC TC\n"
            "1\trs23457_1\t0\t3214571\t0 0\tT TC\tTC TC\n"
            "1\trs23457_2\t0\t3214572\t0 0\tT TC\tTC TC\n"
            "1\t1:3214573\t0\t3214573\t0 0\tT TC\tTC TC\n"
        )
        observed = None
        with open(template_name.format(ext="calls"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (impute2_info)
        info_fn = template_name.format(ext="impute2_info")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "chr\tname\tposition\ta0\ta1\texp_freq_a1\tinfo\tcertainty\t"
            "type\tinfo_type0\tconcord_type0\tr2_type0\n"
            "1\trs23456\t3214569\tT\tC\t0.082\t0.362\t0.866\t0\t-1\t-1\t"
            "-1\n"
            "1\trs23457\t3214570\tT\tTC\t0.126\t0.299\t0.832\t0\t-1\t-1\t"
            "-1\n"
            "1\trs23457_1\t3214571\tT\tTC\t0.060\t0.300\t0.909\t0\t-1\t"
            "-1\t-1\n"
            "1\trs23457_2\t3214572\tT\tTC\t0.084\t0.203\t0.854\t0\t-1\t"
            "-1\t-1\n"
            "1\t1:3214573\t3214573\tT\tTC\t0.371\t0.339\t0.619\t0\t-1\t"
            "-1\t-1\n"
        )
        observed = None
        with open(info_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (map)
        map_fn = template_name.format(ext="map")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "1\trs23456\t0\t3214569\n"
            "1\trs23457\t0\t3214570\n"
            "1\trs23457_1\t0\t3214571\n"
            "1\trs23457_2\t0\t3214572\n"
            "1\t1:3214573\t0\t3214573\n"
        )
        observed = None
        with open(map_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion file (map)
        maf_fn = template_name.format(ext="maf")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "name\tmajor\tminor\tmaf\n"
            "rs23456\tT\tC\t{}\n"
            "rs23457\tTC\tT\t{}\n"
            "rs23457_1\tTC\tT\t{}\n"
            "rs23457_2\tTC\tT\t{}\n"
            "1:3214573\tTC\tT\t{}\n"
        ).format(0.5, 1/4, 1/4, 1/4, 1/4)
        observed = None
        with open(maf_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

    def test_rate(self):
        """Tests the extraction by completion rate."""
        # Executing the script
        args = self.common_args + [
            "--rate", "0.5",
        ]
        impute2_extractor.main(args=args)
        TestImpute2Extractor.clean_logging_handlers()

        # Testing we have the three output files
        template_name = os.path.join(self.output_dir.name, "results.{ext}")
        for suffix in ("impute2", "dosage", "calls"):
            self.assertTrue(os.path.isfile(template_name.format(ext=suffix)))

        # Checking the impute2 file
        expected = (
            "1 rs12345 1231415 A G 1 0 0 0.988 0.002 0 0 0.997 0.003\n"
            "1 rs23456 3214569 T C 0.869 0.130 0 0.903 0.095 0.002 0 0 1\n"
            "1 rs23457 3214570 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 rs23457_1 3214571 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 rs23457_2 3214572 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 1:3214573 3214573 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
        )
        observed = None
        with open(template_name.format(ext="impute2"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the dosage file
        expected = (
            "chrom\tpos\tname\tminor\tmajor\tf1/s1\tf2/s2\tf3/s3\n"
            "1\t1231415\trs12345\tG\tA\t0.0\t0.002\t1.003\n"
            "1\t3214569\trs23456\tC\tT\tnan\t0.099\t2.0\n"
            "1\t3214570\trs23457\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t3214571\trs23457_1\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t3214572\trs23457_2\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t3214573\t1:3214573\tT\tTC\tnan\t1.0\t0.0\n"
        )
        observed = None
        with open(template_name.format(ext="dosage"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the hard calls file
        expected = (
            "chrom\tname\tcm\tpos\tf1/s1\tf2/s2\tf3/s3\n"
            "1\trs12345\t0\t1231415\tA A\tA A\tA G\n"
            "1\trs23456\t0\t3214569\t0 0\tT T\tC C\n"
            "1\trs23457\t0\t3214570\t0 0\tT TC\tTC TC\n"
            "1\trs23457_1\t0\t3214571\t0 0\tT TC\tTC TC\n"
            "1\trs23457_2\t0\t3214572\t0 0\tT TC\tTC TC\n"
            "1\t1:3214573\t0\t3214573\t0 0\tT TC\tTC TC\n"
        )
        observed = None
        with open(template_name.format(ext="calls"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (impute2_info)
        info_fn = template_name.format(ext="impute2_info")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "chr\tname\tposition\ta0\ta1\texp_freq_a1\tinfo\tcertainty\t"
            "type\tinfo_type0\tconcord_type0\tr2_type0\n"
            "1\trs12345\t1231415\tA\tG\t0.006\t0.359\t0.987\t0\t-1\t-1\t"
            "-1\n"
            "1\trs23456\t3214569\tT\tC\t0.082\t0.362\t0.866\t0\t-1\t-1\t"
            "-1\n"
            "1\trs23457\t3214570\tT\tTC\t0.126\t0.299\t0.832\t0\t-1\t-1\t"
            "-1\n"
            "1\trs23457_1\t3214571\tT\tTC\t0.060\t0.300\t0.909\t0\t-1\t"
            "-1\t-1\n"
            "1\trs23457_2\t3214572\tT\tTC\t0.084\t0.203\t0.854\t0\t-1\t"
            "-1\t-1\n"
            "1\t1:3214573\t3214573\tT\tTC\t0.371\t0.339\t0.619\t0\t-1\t"
            "-1\t-1\n"
        )
        observed = None
        with open(info_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (map)
        map_fn = template_name.format(ext="map")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "1\trs12345\t0\t1231415\n"
            "1\trs23456\t0\t3214569\n"
            "1\trs23457\t0\t3214570\n"
            "1\trs23457_1\t0\t3214571\n"
            "1\trs23457_2\t0\t3214572\n"
            "1\t1:3214573\t0\t3214573\n"
        )
        observed = None
        with open(map_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion file (map)
        maf_fn = template_name.format(ext="maf")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "name\tmajor\tminor\tmaf\n"
            "rs12345\tA\tG\t{}\n"
            "rs23456\tT\tC\t{}\n"
            "rs23457\tTC\tT\t{}\n"
            "rs23457_1\tTC\tT\t{}\n"
            "rs23457_2\tTC\tT\t{}\n"
            "1:3214573\tTC\tT\t{}\n"
        ).format(1/6, 0.5, 1/4, 1/4, 1/4, 1/4)
        observed = None
        with open(maf_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

    def test_info(self):
        """Tests the extraction by information value."""
        # Executing the script
        args = self.common_args + [
            "--info", "0.3",
        ]
        impute2_extractor.main(args=args)
        TestImpute2Extractor.clean_logging_handlers()

        # Testing we have the three output files
        template_name = os.path.join(self.output_dir.name, "results.{ext}")
        for suffix in ("impute2", "dosage", "calls"):
            self.assertTrue(os.path.isfile(template_name.format(ext=suffix)))

        # Checking the impute2 file
        expected = (
            "1 rs12345 1231415 A G 1 0 0 0.988 0.002 0 0 0.997 0.003\n"
            "1 rs23456 3214569 T C 0.869 0.130 0 0.903 0.095 0.002 0 0 1\n"
            "1 rs23457_1 3214571 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 1:3214573 3214573 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 1:4214570_1 4214570 T TC 0.869 0.130 0 0.869 0.130 0 0.869 "
            "0.130 0\n"
        )
        observed = None
        with open(template_name.format(ext="impute2"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the dosage file
        expected = (
            "chrom\tpos\tname\tminor\tmajor\tf1/s1\tf2/s2\tf3/s3\n"
            "1\t1231415\trs12345\tG\tA\t0.0\t0.002\t1.003\n"
            "1\t3214569\trs23456\tC\tT\tnan\t0.099\t2.0\n"
            "1\t3214571\trs23457_1\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t3214573\t1:3214573\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t4214570\t1:4214570_1\tTC\tT\tnan\tnan\tnan\n"
        )
        observed = None
        with open(template_name.format(ext="dosage"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the hard calls file
        expected = (
            "chrom\tname\tcm\tpos\tf1/s1\tf2/s2\tf3/s3\n"
            "1\trs12345\t0\t1231415\tA A\tA A\tA G\n"
            "1\trs23456\t0\t3214569\t0 0\tT T\tC C\n"
            "1\trs23457_1\t0\t3214571\t0 0\tT TC\tTC TC\n"
            "1\t1:3214573\t0\t3214573\t0 0\tT TC\tTC TC\n"
            "1\t1:4214570_1\t0\t4214570\t0 0\t0 0\t0 0\n"
        )
        observed = None
        with open(template_name.format(ext="calls"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (impute2_info)
        info_fn = template_name.format(ext="impute2_info")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "chr\tname\tposition\ta0\ta1\texp_freq_a1\tinfo\tcertainty\t"
            "type\tinfo_type0\tconcord_type0\tr2_type0\n"
            "1\trs12345\t1231415\tA\tG\t0.006\t0.359\t0.987\t0\t-1\t-1\t"
            "-1\n"
            "1\trs23456\t3214569\tT\tC\t0.082\t0.362\t0.866\t0\t-1\t-1\t"
            "-1\n"
            "1\trs23457_1\t3214571\tT\tTC\t0.060\t0.300\t0.909\t0\t-1\t"
            "-1\t-1\n"
            "1\t1:3214573\t3214573\tT\tTC\t0.371\t0.339\t0.619\t0\t-1\t"
            "-1\t-1\n"
            "1\t1:4214570_1\t4214570\tT\tTC\t0.174\t0.589\t0.831\t0\t-1\t"
            "-1\t-1\n"
        )
        observed = None
        with open(info_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (map)
        map_fn = template_name.format(ext="map")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "1\trs12345\t0\t1231415\n"
            "1\trs23456\t0\t3214569\n"
            "1\trs23457_1\t0\t3214571\n"
            "1\t1:3214573\t0\t3214573\n"
            "1\t1:4214570_1\t0\t4214570\n"
        )
        observed = None
        with open(map_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion file (map)
        maf_fn = template_name.format(ext="maf")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "name\tmajor\tminor\tmaf\n"
            "rs12345\tA\tG\t{}\n"
            "rs23456\tT\tC\t{}\n"
            "rs23457_1\tTC\tT\t{}\n"
            "1:3214573\tTC\tT\t{}\n"
            "1:4214570_1\tT\tTC\t{}\n"
        ).format(1/6, 0.5, 1/4, 1/4, "NA")
        observed = None
        with open(maf_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

    def test_genomic_maf(self):
        """Tests the extraction by genomic location and maf."""
        # Executing the script
        args = self.common_args + [
            "--genomic", "chr1:3214569-3214573",
            "--maf", "0.3",
        ]
        impute2_extractor.main(args=args)
        TestImpute2Extractor.clean_logging_handlers()

        # Testing we have the three output files
        template_name = os.path.join(self.output_dir.name, "results.{ext}")
        for suffix in ("impute2", "dosage", "calls"):
            self.assertTrue(os.path.isfile(template_name.format(ext=suffix)))

        # Checking the impute2 file
        expected = (
            "1 rs23456 3214569 T C 0.869 0.130 0 0.903 0.095 0.002 0 0 1\n"
        )
        observed = None
        with open(template_name.format(ext="impute2"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the dosage file
        expected = (
            "chrom\tpos\tname\tminor\tmajor\tf1/s1\tf2/s2\tf3/s3\n"
            "1\t3214569\trs23456\tC\tT\tnan\t0.099\t2.0\n"
        )
        observed = None
        with open(template_name.format(ext="dosage"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the hard calls file
        expected = (
            "chrom\tname\tcm\tpos\tf1/s1\tf2/s2\tf3/s3\n"
            "1\trs23456\t0\t3214569\t0 0\tT T\tC C\n"
        )
        observed = None
        with open(template_name.format(ext="calls"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (impute2_info)
        info_fn = template_name.format(ext="impute2_info")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "chr\tname\tposition\ta0\ta1\texp_freq_a1\tinfo\tcertainty\t"
            "type\tinfo_type0\tconcord_type0\tr2_type0\n"
            "1\trs23456\t3214569\tT\tC\t0.082\t0.362\t0.866\t0\t-1\t-1\t"
            "-1\n"
        )
        observed = None
        with open(info_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (map)
        map_fn = template_name.format(ext="map")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "1\trs23456\t0\t3214569\n"
        )
        observed = None
        with open(map_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion file (map)
        maf_fn = template_name.format(ext="maf")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "name\tmajor\tminor\tmaf\n"
            "rs23456\tT\tC\t{}\n"
        ).format(0.5)
        observed = None
        with open(maf_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

    def test_genomic_rate(self):
        """Tests the extraction by genomic location and completion rate."""
        # Executing the script
        args = self.common_args + [
            "--genomic", "chr1:1231415-3214573",
            "--rate", "0.7",
        ]
        impute2_extractor.main(args=args)
        TestImpute2Extractor.clean_logging_handlers()

        # Testing we have the three output files
        template_name = os.path.join(self.output_dir.name, "results.{ext}")
        for suffix in ("impute2", "dosage", "calls"):
            self.assertTrue(os.path.isfile(template_name.format(ext=suffix)))

        # Checking the impute2 file
        expected = (
            "1 rs12345 1231415 A G 1 0 0 0.988 0.002 0 0 0.997 0.003\n"
        )
        observed = None
        with open(template_name.format(ext="impute2"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the dosage file
        expected = (
            "chrom\tpos\tname\tminor\tmajor\tf1/s1\tf2/s2\tf3/s3\n"
            "1\t1231415\trs12345\tG\tA\t0.0\t0.002\t1.003\n"
        )
        observed = None
        with open(template_name.format(ext="dosage"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the hard calls file
        expected = (
            "chrom\tname\tcm\tpos\tf1/s1\tf2/s2\tf3/s3\n"
            "1\trs12345\t0\t1231415\tA A\tA A\tA G\n"
        )
        observed = None
        with open(template_name.format(ext="calls"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (impute2_info)
        info_fn = template_name.format(ext="impute2_info")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "chr\tname\tposition\ta0\ta1\texp_freq_a1\tinfo\tcertainty\t"
            "type\tinfo_type0\tconcord_type0\tr2_type0\n"
            "1\trs12345\t1231415\tA\tG\t0.006\t0.359\t0.987\t0\t-1\t-1\t"
            "-1\n"
        )
        observed = None
        with open(info_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (map)
        map_fn = template_name.format(ext="map")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "1\trs12345\t0\t1231415\n"
        )
        observed = None
        with open(map_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion file (map)
        maf_fn = template_name.format(ext="maf")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "name\tmajor\tminor\tmaf\n"
            "rs12345\tA\tG\t{}\n"
        ).format(1/6)
        observed = None
        with open(maf_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

    def test_genomic_info(self):
        """Tests the extraction by genomic location and information value."""
        # Executing the script
        args = self.common_args + [
            "--genomic", "chr1:1231415-3214573",
            "--info", "0.28",
        ]
        impute2_extractor.main(args=args)
        TestImpute2Extractor.clean_logging_handlers()

        # Testing we have the three output files
        template_name = os.path.join(self.output_dir.name, "results.{ext}")
        for suffix in ("impute2", "dosage", "calls"):
            self.assertTrue(os.path.isfile(template_name.format(ext=suffix)))

        # Checking the impute2 file
        expected = (
            "1 rs12345 1231415 A G 1 0 0 0.988 0.002 0 0 0.997 0.003\n"
            "1 rs23456 3214569 T C 0.869 0.130 0 0.903 0.095 0.002 0 0 1\n"
            "1 rs23457 3214570 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 rs23457_1 3214571 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 1:3214573 3214573 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
        )
        observed = None
        with open(template_name.format(ext="impute2"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the dosage file
        expected = (
            "chrom\tpos\tname\tminor\tmajor\tf1/s1\tf2/s2\tf3/s3\n"
            "1\t1231415\trs12345\tG\tA\t0.0\t0.002\t1.003\n"
            "1\t3214569\trs23456\tC\tT\tnan\t0.099\t2.0\n"
            "1\t3214570\trs23457\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t3214571\trs23457_1\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t3214573\t1:3214573\tT\tTC\tnan\t1.0\t0.0\n"
        )
        observed = None
        with open(template_name.format(ext="dosage"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the hard calls file
        expected = (
            "chrom\tname\tcm\tpos\tf1/s1\tf2/s2\tf3/s3\n"
            "1\trs12345\t0\t1231415\tA A\tA A\tA G\n"
            "1\trs23456\t0\t3214569\t0 0\tT T\tC C\n"
            "1\trs23457\t0\t3214570\t0 0\tT TC\tTC TC\n"
            "1\trs23457_1\t0\t3214571\t0 0\tT TC\tTC TC\n"
            "1\t1:3214573\t0\t3214573\t0 0\tT TC\tTC TC\n"
        )
        observed = None
        with open(template_name.format(ext="calls"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (impute2_info)
        info_fn = template_name.format(ext="impute2_info")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "chr\tname\tposition\ta0\ta1\texp_freq_a1\tinfo\tcertainty\t"
            "type\tinfo_type0\tconcord_type0\tr2_type0\n"
            "1\trs12345\t1231415\tA\tG\t0.006\t0.359\t0.987\t0\t-1\t-1\t"
            "-1\n"
            "1\trs23456\t3214569\tT\tC\t0.082\t0.362\t0.866\t0\t-1\t-1\t"
            "-1\n"
            "1\trs23457\t3214570\tT\tTC\t0.126\t0.299\t0.832\t0\t-1\t-1\t"
            "-1\n"
            "1\trs23457_1\t3214571\tT\tTC\t0.060\t0.300\t0.909\t0\t-1\t"
            "-1\t-1\n"
            "1\t1:3214573\t3214573\tT\tTC\t0.371\t0.339\t0.619\t0\t-1\t"
            "-1\t-1\n"
        )
        observed = None
        with open(info_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (map)
        map_fn = template_name.format(ext="map")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "1\trs12345\t0\t1231415\n"
            "1\trs23456\t0\t3214569\n"
            "1\trs23457\t0\t3214570\n"
            "1\trs23457_1\t0\t3214571\n"
            "1\t1:3214573\t0\t3214573\n"
        )
        observed = None
        with open(map_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion file (map)
        maf_fn = template_name.format(ext="maf")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "name\tmajor\tminor\tmaf\n"
            "rs12345\tA\tG\t{}\n"
            "rs23456\tT\tC\t{}\n"
            "rs23457\tTC\tT\t{}\n"
            "rs23457_1\tTC\tT\t{}\n"
            "1:3214573\tTC\tT\t{}\n"
        ).format(1/6, 0.5, 1/4, 1/4, 1/4)
        observed = None
        with open(maf_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

    def test_maf_rate(self):
        """Tests the extraction by maf and completion rate."""
        # Executing the script
        args = self.common_args + [
            "--maf", "0.2",
            "--rate", "0.6",
        ]
        impute2_extractor.main(args=args)
        TestImpute2Extractor.clean_logging_handlers()

        # Testing we have the three output files
        template_name = os.path.join(self.output_dir.name, "results.{ext}")
        for suffix in ("impute2", "dosage", "calls"):
            self.assertTrue(os.path.isfile(template_name.format(ext=suffix)))

        # Checking the impute2 file
        expected = (
            "1 rs23456 3214569 T C 0.869 0.130 0 0.903 0.095 0.002 0 0 1\n"
            "1 rs23457 3214570 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 rs23457_1 3214571 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 rs23457_2 3214572 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 1:3214573 3214573 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
        )
        observed = None
        with open(template_name.format(ext="impute2"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the dosage file
        expected = (
            "chrom\tpos\tname\tminor\tmajor\tf1/s1\tf2/s2\tf3/s3\n"
            "1\t3214569\trs23456\tC\tT\tnan\t0.099\t2.0\n"
            "1\t3214570\trs23457\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t3214571\trs23457_1\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t3214572\trs23457_2\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t3214573\t1:3214573\tT\tTC\tnan\t1.0\t0.0\n"
        )
        observed = None
        with open(template_name.format(ext="dosage"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the hard calls file
        expected = (
            "chrom\tname\tcm\tpos\tf1/s1\tf2/s2\tf3/s3\n"
            "1\trs23456\t0\t3214569\t0 0\tT T\tC C\n"
            "1\trs23457\t0\t3214570\t0 0\tT TC\tTC TC\n"
            "1\trs23457_1\t0\t3214571\t0 0\tT TC\tTC TC\n"
            "1\trs23457_2\t0\t3214572\t0 0\tT TC\tTC TC\n"
            "1\t1:3214573\t0\t3214573\t0 0\tT TC\tTC TC\n"
        )
        observed = None
        with open(template_name.format(ext="calls"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (impute2_info)
        info_fn = template_name.format(ext="impute2_info")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "chr\tname\tposition\ta0\ta1\texp_freq_a1\tinfo\tcertainty\t"
            "type\tinfo_type0\tconcord_type0\tr2_type0\n"
            "1\trs23456\t3214569\tT\tC\t0.082\t0.362\t0.866\t0\t-1\t-1\t"
            "-1\n"
            "1\trs23457\t3214570\tT\tTC\t0.126\t0.299\t0.832\t0\t-1\t-1\t"
            "-1\n"
            "1\trs23457_1\t3214571\tT\tTC\t0.060\t0.300\t0.909\t0\t-1\t"
            "-1\t-1\n"
            "1\trs23457_2\t3214572\tT\tTC\t0.084\t0.203\t0.854\t0\t-1\t"
            "-1\t-1\n"
            "1\t1:3214573\t3214573\tT\tTC\t0.371\t0.339\t0.619\t0\t-1\t"
            "-1\t-1\n"
        )
        observed = None
        with open(info_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (map)
        map_fn = template_name.format(ext="map")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "1\trs23456\t0\t3214569\n"
            "1\trs23457\t0\t3214570\n"
            "1\trs23457_1\t0\t3214571\n"
            "1\trs23457_2\t0\t3214572\n"
            "1\t1:3214573\t0\t3214573\n"
        )
        observed = None
        with open(map_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion file (map)
        maf_fn = template_name.format(ext="maf")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "name\tmajor\tminor\tmaf\n"
            "rs23456\tT\tC\t{}\n"
            "rs23457\tTC\tT\t{}\n"
            "rs23457_1\tTC\tT\t{}\n"
            "rs23457_2\tTC\tT\t{}\n"
            "1:3214573\tTC\tT\t{}\n"
        ).format(0.5, 1/4, 1/4, 1/4, 1/4)
        observed = None
        with open(maf_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

    def test_maf_info(self):
        """Tests the extraction by maf and information value."""
        # Executing the script
        args = self.common_args + [
            "--maf", "0.2",
            "--info", "0.3",
        ]
        impute2_extractor.main(args=args)
        TestImpute2Extractor.clean_logging_handlers()

        # Testing we have the three output files
        template_name = os.path.join(self.output_dir.name, "results.{ext}")
        for suffix in ("impute2", "dosage", "calls"):
            self.assertTrue(os.path.isfile(template_name.format(ext=suffix)))

        # Checking the impute2 file
        expected = (
            "1 rs23456 3214569 T C 0.869 0.130 0 0.903 0.095 0.002 0 0 1\n"
            "1 rs23457_1 3214571 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 1:3214573 3214573 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
        )
        observed = None
        with open(template_name.format(ext="impute2"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the dosage file
        expected = (
            "chrom\tpos\tname\tminor\tmajor\tf1/s1\tf2/s2\tf3/s3\n"
            "1\t3214569\trs23456\tC\tT\tnan\t0.099\t2.0\n"
            "1\t3214571\trs23457_1\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t3214573\t1:3214573\tT\tTC\tnan\t1.0\t0.0\n"
        )
        observed = None
        with open(template_name.format(ext="dosage"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the hard calls file
        expected = (
            "chrom\tname\tcm\tpos\tf1/s1\tf2/s2\tf3/s3\n"
            "1\trs23456\t0\t3214569\t0 0\tT T\tC C\n"
            "1\trs23457_1\t0\t3214571\t0 0\tT TC\tTC TC\n"
            "1\t1:3214573\t0\t3214573\t0 0\tT TC\tTC TC\n"
        )
        observed = None
        with open(template_name.format(ext="calls"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (impute2_info)
        info_fn = template_name.format(ext="impute2_info")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "chr\tname\tposition\ta0\ta1\texp_freq_a1\tinfo\tcertainty\t"
            "type\tinfo_type0\tconcord_type0\tr2_type0\n"
            "1\trs23456\t3214569\tT\tC\t0.082\t0.362\t0.866\t0\t-1\t-1\t"
            "-1\n"
            "1\trs23457_1\t3214571\tT\tTC\t0.060\t0.300\t0.909\t0\t-1\t"
            "-1\t-1\n"
            "1\t1:3214573\t3214573\tT\tTC\t0.371\t0.339\t0.619\t0\t-1\t"
            "-1\t-1\n"
        )
        observed = None
        with open(info_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (map)
        map_fn = template_name.format(ext="map")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "1\trs23456\t0\t3214569\n"
            "1\trs23457_1\t0\t3214571\n"
            "1\t1:3214573\t0\t3214573\n"
        )
        observed = None
        with open(map_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion file (map)
        maf_fn = template_name.format(ext="maf")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "name\tmajor\tminor\tmaf\n"
            "rs23456\tT\tC\t{}\n"
            "rs23457_1\tTC\tT\t{}\n"
            "1:3214573\tTC\tT\t{}\n"
        ).format(0.5, 1/4, 1/4)
        observed = None
        with open(maf_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

    def test_rate_info(self):
        """Tests the extraction by completion rate and information value."""
        # Executing the script
        args = self.common_args + [
            "--info", "0.35",
            "--rate", "0.6",
        ]
        impute2_extractor.main(args=args)
        TestImpute2Extractor.clean_logging_handlers()

        # Testing we have the three output files
        template_name = os.path.join(self.output_dir.name, "results.{ext}")
        for suffix in ("impute2", "dosage", "calls"):
            self.assertTrue(os.path.isfile(template_name.format(ext=suffix)))

        # Checking the impute2 file
        expected = (
            "1 rs12345 1231415 A G 1 0 0 0.988 0.002 0 0 0.997 0.003\n"
            "1 rs23456 3214569 T C 0.869 0.130 0 0.903 0.095 0.002 0 0 1\n"
        )
        observed = None
        with open(template_name.format(ext="impute2"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the dosage file
        expected = (
            "chrom\tpos\tname\tminor\tmajor\tf1/s1\tf2/s2\tf3/s3\n"
            "1\t1231415\trs12345\tG\tA\t0.0\t0.002\t1.003\n"
            "1\t3214569\trs23456\tC\tT\tnan\t0.099\t2.0\n"
        )
        observed = None
        with open(template_name.format(ext="dosage"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the hard calls file
        expected = (
            "chrom\tname\tcm\tpos\tf1/s1\tf2/s2\tf3/s3\n"
            "1\trs12345\t0\t1231415\tA A\tA A\tA G\n"
            "1\trs23456\t0\t3214569\t0 0\tT T\tC C\n"
        )
        observed = None
        with open(template_name.format(ext="calls"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (impute2_info)
        info_fn = template_name.format(ext="impute2_info")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "chr\tname\tposition\ta0\ta1\texp_freq_a1\tinfo\tcertainty\t"
            "type\tinfo_type0\tconcord_type0\tr2_type0\n"
            "1\trs12345\t1231415\tA\tG\t0.006\t0.359\t0.987\t0\t-1\t-1\t"
            "-1\n"
            "1\trs23456\t3214569\tT\tC\t0.082\t0.362\t0.866\t0\t-1\t-1\t"
            "-1\n"
        )
        observed = None
        with open(info_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (map)
        map_fn = template_name.format(ext="map")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "1\trs12345\t0\t1231415\n"
            "1\trs23456\t0\t3214569\n"
        )
        observed = None
        with open(map_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion file (map)
        maf_fn = template_name.format(ext="maf")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "name\tmajor\tminor\tmaf\n"
            "rs12345\tA\tG\t{}\n"
            "rs23456\tT\tC\t{}\n"
        ).format(1/6, 0.5)
        observed = None
        with open(maf_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

    def test_genomic_maf_rate(self):
        """Tests the extraction by genomic location, maf and rate."""
        # Executing the script
        args = self.common_args + [
            "--genomic", "chr1:3214568-4514570",
            "--maf", "0.01",
            "--rate", "0.6",
        ]
        impute2_extractor.main(args=args)
        TestImpute2Extractor.clean_logging_handlers()

        # Testing we have the three output files
        template_name = os.path.join(self.output_dir.name, "results.{ext}")
        for suffix in ("impute2", "dosage", "calls"):
            self.assertTrue(os.path.isfile(template_name.format(ext=suffix)))

        # Checking the impute2 file
        expected = (
            "1 rs23456 3214569 T C 0.869 0.130 0 0.903 0.095 0.002 0 0 1\n"
            "1 rs23457 3214570 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 rs23457_1 3214571 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 rs23457_2 3214572 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
            "1 1:3214573 3214573 T TC 0.869 0.130 0 0 1 0 0 0 1\n"
        )
        observed = None
        with open(template_name.format(ext="impute2"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the dosage file
        expected = (
            "chrom\tpos\tname\tminor\tmajor\tf1/s1\tf2/s2\tf3/s3\n"
            "1\t3214569\trs23456\tC\tT\tnan\t0.099\t2.0\n"
            "1\t3214570\trs23457\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t3214571\trs23457_1\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t3214572\trs23457_2\tT\tTC\tnan\t1.0\t0.0\n"
            "1\t3214573\t1:3214573\tT\tTC\tnan\t1.0\t0.0\n"
        )
        observed = None
        with open(template_name.format(ext="dosage"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the hard calls file
        expected = (
            "chrom\tname\tcm\tpos\tf1/s1\tf2/s2\tf3/s3\n"
            "1\trs23456\t0\t3214569\t0 0\tT T\tC C\n"
            "1\trs23457\t0\t3214570\t0 0\tT TC\tTC TC\n"
            "1\trs23457_1\t0\t3214571\t0 0\tT TC\tTC TC\n"
            "1\trs23457_2\t0\t3214572\t0 0\tT TC\tTC TC\n"
            "1\t1:3214573\t0\t3214573\t0 0\tT TC\tTC TC\n"
        )
        observed = None
        with open(template_name.format(ext="calls"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (impute2_info)
        info_fn = template_name.format(ext="impute2_info")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "chr\tname\tposition\ta0\ta1\texp_freq_a1\tinfo\tcertainty\t"
            "type\tinfo_type0\tconcord_type0\tr2_type0\n"
            "1\trs23456\t3214569\tT\tC\t0.082\t0.362\t0.866\t0\t-1\t-1\t"
            "-1\n"
            "1\trs23457\t3214570\tT\tTC\t0.126\t0.299\t0.832\t0\t-1\t-1\t"
            "-1\n"
            "1\trs23457_1\t3214571\tT\tTC\t0.060\t0.300\t0.909\t0\t-1\t"
            "-1\t-1\n"
            "1\trs23457_2\t3214572\tT\tTC\t0.084\t0.203\t0.854\t0\t-1\t"
            "-1\t-1\n"
            "1\t1:3214573\t3214573\tT\tTC\t0.371\t0.339\t0.619\t0\t-1\t"
            "-1\t-1\n"
        )
        observed = None
        with open(info_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (map)
        map_fn = template_name.format(ext="map")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "1\trs23456\t0\t3214569\n"
            "1\trs23457\t0\t3214570\n"
            "1\trs23457_1\t0\t3214571\n"
            "1\trs23457_2\t0\t3214572\n"
            "1\t1:3214573\t0\t3214573\n"
        )
        observed = None
        with open(map_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion file (map)
        maf_fn = template_name.format(ext="maf")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "name\tmajor\tminor\tmaf\n"
            "rs23456\tT\tC\t{}\n"
            "rs23457\tTC\tT\t{}\n"
            "rs23457_1\tTC\tT\t{}\n"
            "rs23457_2\tTC\tT\t{}\n"
            "1:3214573\tTC\tT\t{}\n"
        ).format(0.5, 1/4, 1/4, 1/4, 1/4)
        observed = None
        with open(maf_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

    def test_genomic_maf_info(self):
        """Tests the extraction by genomic location, maf and information."""
        # Executing the script
        args = self.common_args + [
            "--genomic", "chr1:3214568-4514570",
            "--maf", "0.01",
            "--info", "0.35",
        ]
        impute2_extractor.main(args=args)
        TestImpute2Extractor.clean_logging_handlers()

        # Testing we have the three output files
        template_name = os.path.join(self.output_dir.name, "results.{ext}")
        for suffix in ("impute2", "dosage", "calls"):
            self.assertTrue(os.path.isfile(template_name.format(ext=suffix)))

        # Checking the impute2 file
        expected = (
            "1 rs23456 3214569 T C 0.869 0.130 0 0.903 0.095 0.002 0 0 1\n"
        )
        observed = None
        with open(template_name.format(ext="impute2"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the dosage file
        expected = (
            "chrom\tpos\tname\tminor\tmajor\tf1/s1\tf2/s2\tf3/s3\n"
            "1\t3214569\trs23456\tC\tT\tnan\t0.099\t2.0\n"
        )
        observed = None
        with open(template_name.format(ext="dosage"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the hard calls file
        expected = (
            "chrom\tname\tcm\tpos\tf1/s1\tf2/s2\tf3/s3\n"
            "1\trs23456\t0\t3214569\t0 0\tT T\tC C\n"
        )
        observed = None
        with open(template_name.format(ext="calls"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (impute2_info)
        info_fn = template_name.format(ext="impute2_info")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "chr\tname\tposition\ta0\ta1\texp_freq_a1\tinfo\tcertainty\t"
            "type\tinfo_type0\tconcord_type0\tr2_type0\n"
            "1\trs23456\t3214569\tT\tC\t0.082\t0.362\t0.866\t0\t-1\t-1\t"
            "-1\n"
        )
        observed = None
        with open(info_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (map)
        map_fn = template_name.format(ext="map")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "1\trs23456\t0\t3214569\n"
        )
        observed = None
        with open(map_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion file (map)
        maf_fn = template_name.format(ext="maf")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "name\tmajor\tminor\tmaf\n"
            "rs23456\tT\tC\t{}\n"
        ).format(0.5)
        observed = None
        with open(maf_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

    def test_maf_rate_info(self):
        """Tests the extraction by maf, completion rate and information."""
        # Executing the script
        args = self.common_args + [
            "--rate", "0.7",
            "--maf", "0.01",
            "--info", "0.35",
        ]
        impute2_extractor.main(args=args)
        TestImpute2Extractor.clean_logging_handlers()

        # Testing we have the three output files
        template_name = os.path.join(self.output_dir.name, "results.{ext}")
        for suffix in ("impute2", "dosage", "calls"):
            self.assertTrue(os.path.isfile(template_name.format(ext=suffix)))

        # Checking the impute2 file
        expected = (
            "1 rs12345 1231415 A G 1 0 0 0.988 0.002 0 0 0.997 0.003\n"
        )
        observed = None
        with open(template_name.format(ext="impute2"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the dosage file
        expected = (
            "chrom\tpos\tname\tminor\tmajor\tf1/s1\tf2/s2\tf3/s3\n"
            "1\t1231415\trs12345\tG\tA\t0.0\t0.002\t1.003\n"
        )
        observed = None
        with open(template_name.format(ext="dosage"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the hard calls file
        expected = (
            "chrom\tname\tcm\tpos\tf1/s1\tf2/s2\tf3/s3\n"
            "1\trs12345\t0\t1231415\tA A\tA A\tA G\n"
        )
        observed = None
        with open(template_name.format(ext="calls"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (impute2_info)
        info_fn = template_name.format(ext="impute2_info")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "chr\tname\tposition\ta0\ta1\texp_freq_a1\tinfo\tcertainty\t"
            "type\tinfo_type0\tconcord_type0\tr2_type0\n"
            "1\trs12345\t1231415\tA\tG\t0.006\t0.359\t0.987\t0\t-1\t-1\t"
            "-1\n"
        )
        observed = None
        with open(info_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (map)
        map_fn = template_name.format(ext="map")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "1\trs12345\t0\t1231415\n"
        )
        observed = None
        with open(map_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion file (map)
        maf_fn = template_name.format(ext="maf")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "name\tmajor\tminor\tmaf\n"
            "rs12345\tA\tG\t{}\n"
        ).format(1/6)
        observed = None
        with open(maf_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

    def test_genomic_maf_rate_info(self):
        """Tests the extraction by genomic location, maf, rate and info."""
        # Executing the script
        args = self.common_args + [
            "--genomic", "chr1:1231415-3214572",
            "--rate", "0.6",
            "--maf", "0.2",
            "--info", "0.35",
        ]
        impute2_extractor.main(args=args)
        TestImpute2Extractor.clean_logging_handlers()

        # Testing we have the three output files
        template_name = os.path.join(self.output_dir.name, "results.{ext}")
        for suffix in ("impute2", "dosage", "calls"):
            self.assertTrue(os.path.isfile(template_name.format(ext=suffix)))

        # Checking the impute2 file
        expected = (
            "1 rs23456 3214569 T C 0.869 0.130 0 0.903 0.095 0.002 0 0 1\n"
        )
        observed = None
        with open(template_name.format(ext="impute2"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the dosage file
        expected = (
            "chrom\tpos\tname\tminor\tmajor\tf1/s1\tf2/s2\tf3/s3\n"
            "1\t3214569\trs23456\tC\tT\tnan\t0.099\t2.0\n"
        )
        observed = None
        with open(template_name.format(ext="dosage"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checking the hard calls file
        expected = (
            "chrom\tname\tcm\tpos\tf1/s1\tf2/s2\tf3/s3\n"
            "1\trs23456\t0\t3214569\t0 0\tT T\tC C\n"
        )
        observed = None
        with open(template_name.format(ext="calls"), "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (impute2_info)
        info_fn = template_name.format(ext="impute2_info")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "chr\tname\tposition\ta0\ta1\texp_freq_a1\tinfo\tcertainty\t"
            "type\tinfo_type0\tconcord_type0\tr2_type0\n"
            "1\trs23456\t3214569\tT\tC\t0.082\t0.362\t0.866\t0\t-1\t-1\t"
            "-1\n"
        )
        observed = None
        with open(info_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion files (map)
        map_fn = template_name.format(ext="map")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "1\trs23456\t0\t3214569\n"
        )
        observed = None
        with open(map_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)

        # Checks companion file (map)
        maf_fn = template_name.format(ext="maf")
        self.assertTrue(os.path.isfile(info_fn))
        expected = (
            "name\tmajor\tminor\tmaf\n"
            "rs23456\tT\tC\t{}\n"
        ).format(0.5)
        observed = None
        with open(maf_fn, "r") as i_file:
            observed = i_file.read()
        self.assertEqual(expected, observed)
