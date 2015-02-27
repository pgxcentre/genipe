
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
from .. import chromosomes
from ..error import ProgramError


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

    def test_gather_maf_stats(self):
        """Tests the 'gather_maf_stats' function."""
        # Creating one file per chromosome with two markers in each
        header = ["name", "maf"]
        content = [
            ["marker_1", "0.0"],
            ["marker_2", "0.2"],
            ["marker_3", "NA"],
            ["marker_4", "0.3"],
            ["marker_5", "9.4e-05"],
            ["marker_6", "0.04"],
            ["marker_7", "0.01"],
            ["marker_8", "0.003"],
            ["marker_9", "0.015"],
            ["marker_10", "0.020"],
            ["marker_11", "0.005"],
            ["marker_12", "0.004"],
            ["marker_13", "0.05"],
            ["marker_14", "0.055"],
            ["marker_15", "0.001"],
            ["marker_16", "0.004"],
            ["marker_17", "0.056"],
            ["marker_18", "0.0123"],
            ["marker_19", "0.005"],
            ["marker_20", "0.012"],
            ["marker_21", "0.001"],
            ["marker_22", "0.0316"],
            ["marker_23", "0.432"],
            ["marker_24", "0.03423"],
            ["marker_25", "0.00514"],
            ["marker_26", "0.01004"],
            ["marker_27", "0.011"],
            ["marker_28", "0.051"],
            ["marker_29", "0.048"],
            ["marker_30", "0.0484"],
            ["marker_31", "0.0871"],
            ["marker_32", "0.5"],
            ["marker_33", "0.006"],
            ["marker_34", "0.06"],
            ["marker_35", "0.08"],
            ["marker_36", "0.0784"],
            ["marker_37", "0.0984"],
            ["marker_38", "0.19444"],
            ["marker_39", "1.5e-04"],
            ["marker_40", "NA"],
            ["marker_41", "4.87e-07"],
            ["marker_42", "5.4e-08"],
            ["marker_43", "0.394"],
            ["marker_44", "0.004"],
        ]
        good_sites = ["marker_1", "marker_2", "marker_5", "marker_6",
                      "marker_7", "marker_8", "marker_9", "marker_10",
                      "marker_11", "marker_12", "marker_13", "marker_14",
                      "marker_15", "marker_16", "marker_17", "marker_18",
                      "marker_19", "marker_20", "marker_21", "marker_22",
                      "marker_23", "marker_24", "marker_25", "marker_26",
                      "marker_27", "marker_28", "marker_29", "marker_30",
                      "marker_31", "marker_32", "marker_33", "marker_34",
                      "marker_35", "marker_36", "marker_37", "marker_38",
                      "marker_39", "marker_41", "marker_43", "marker_44"]

        # The PDF generated
        frequency_pie = ""
        if HAS_MATPLOTLIB:
            frequency_pie = os.path.join(self.output_dir.name,
                                         "frequency_pie.pdf")
        # The expected results
        nb_sites = len(good_sites)
        expected_results = {
            "nb_marker_with_maf":   str(nb_sites),
            "nb_maf_geq_01":        "26",
            "pct_maf_geq_01":       "{:.1f}".format(26 / nb_sites * 100),
            "nb_maf_geq_05":        "14",
            "pct_maf_geq_05":       "{:.1f}".format(14 / nb_sites * 100),
            "nb_maf_lt_05":         "26",
            "pct_maf_lt_05":        "{:.1f}".format(26 / nb_sites * 100),
            "nb_maf_lt_01":         "14",
            "pct_maf_lt_01":        "{:.1f}".format(14 / nb_sites * 100),
            "nb_maf_geq_01_lt_05":  "12",
            "pct_maf_geq_01_lt_05": "{:.1f}".format(12 / nb_sites * 100),
            "nb_maf_nan":           "0",
            "frequency_pie":        frequency_pie,
        }

        # Creating the files for the test
        filename_template = os.path.join(self.output_dir.name, "chr{chrom}",
                                         "final_impute2",
                                         "chr{chrom}.imputed.{suffix}")
        for chrom in chromosomes:
            # Getting the name of the file
            maf_filename = filename_template.format(chrom=chrom, suffix="maf")
            good_sites_filename = filename_template.format(chrom=chrom,
                                                           suffix="good_sites")

            # Getting the directory and create it
            dirname = os.path.dirname(maf_filename)
            if not os.path.isdir(dirname):
                os.makedirs(dirname)
            self.assertTrue(os.path.isdir(dirname))

            # Creating the content of the maf file
            with open(maf_filename, "w") as o_file:
                print(*header, sep="\t", file=o_file)
                for i in range(2):
                    print(*content.pop(), sep="\t", file=o_file)

            # Creating the content of the good sites file
            with open(good_sites_filename, "w") as o_file:
                print(*good_sites, sep="\n", file=o_file)

            # Checking the files were created
            self.assertTrue(os.path.isfile(maf_filename))
            self.assertTrue(os.path.isfile(good_sites_filename))

        # Checking we passed all the content (MAF)
        self.assertEqual(0, len(content))

        # Executing the command (getting the observed data)
        observed = gather_maf_stats(self.output_dir.name)

        # Checking the observed results
        self.assertEqual(len(expected_results), len(observed))
        for expected_key, expected_value in expected_results.items():
            self.assertTrue(expected_key in observed)
            self.assertEqual(expected_value, observed[expected_key])

        # If matplotlib is installed, checking we have a figure (and not
        # otherwise)
        if HAS_MATPLOTLIB:
            self.assertTrue(os.path.isfile(frequency_pie))
        else:
            self.assertFalse(os.path.isfile(frequency_pie))

        # Testing an invalid entry
        changed_filename = filename_template.format(chrom=1, suffix="maf")
        with open(changed_filename, "w") as o_file:
            print(*header, sep="\t", file=o_file)
            print("marker_1", "0.6", sep="\t", file=o_file)

        # This should raise an exception
        with self.assertRaises(ProgramError) as cm:
            gather_maf_stats(self.output_dir.name)
        self.assertEqual("{}: {}: invalid MAF".format("marker_1",
                                                      round(0.6, 3)),
                         str(cm.exception))

        # Testing an invalid entry
        changed_filename = filename_template.format(chrom=1, suffix="maf")
        with open(changed_filename, "w") as o_file:
            print(*header, sep="\t", file=o_file)
            print("marker_1", "-0.01", sep="\t", file=o_file)

        # This should raise an exception
        with self.assertRaises(ProgramError) as cm:
            gather_maf_stats(self.output_dir.name)
        self.assertEqual("{}: {}: invalid MAF".format("marker_1",
                                                      round(-0.01, 3)),
                         str(cm.exception))

        # Testing a good site with NA MAF
        changed_filename = filename_template.format(chrom=1, suffix="maf")
        with open(changed_filename, "w") as o_file:
            print(*header, sep="\t", file=o_file)
            print("marker_1", "NA", sep="\t", file=o_file)

        # This should issue a warning
        with self._my_compatibility_assertLogs(level="WARNING") as cm:
            gather_maf_stats(self.output_dir.name)
        log_m = "WARNING:root:chr1: good sites with invalid MAF (NaN)"
        self.assertEqual(1, len(cm.output))
        self.assertEqual(log_m, cm.output[0])

        # Clearing the good sites file to see if we have a warning
        for chrom in chromosomes:
            filename = filename_template.format(chrom=chrom,
                                                suffix="good_sites")
            with open(filename, "w") as o_file:
                pass

        # This should issue a warning
        with self._my_compatibility_assertLogs(level="WARNING") as cm:
            gather_maf_stats(self.output_dir.name)
        log_m = ("WARNING:root:There were no marker with MAF (something went "
                 "wrong)")
        self.assertEqual(1, len(cm.output))
        self.assertEqual(log_m, cm.output[0])

        # Deleting a good sites file, and checking we have an error
        removed_filename = filename_template.format(chrom=1,
                                                    suffix="good_sites")
        os.remove(removed_filename)
        self.assertFalse(os.path.isfile(removed_filename))

        # This should raise an exception
        with self.assertRaises(ProgramError) as cm:
            gather_maf_stats(self.output_dir.name)
        self.assertEqual("{}: no such file".format(removed_filename),
                         str(cm.exception))

        # Deleting a MAF sites file, and checking we have an error
        removed_filename = filename_template.format(chrom=1, suffix="maf")
        os.remove(removed_filename)
        self.assertFalse(os.path.isfile(removed_filename))

        # This should raise an exception
        with self.assertRaises(ProgramError) as cm:
            gather_maf_stats(self.output_dir.name)
        self.assertEqual("{}: no such file".format(removed_filename),
                         str(cm.exception))

    @unittest.skip("Test not implemented")
    def test_gather_execution_time(self):
        """Tests the 'gather_execution_time' function."""
        self.fail("Test not implemented")

    def test_check_args(self):
        """Tests the 'check_args' function."""
        # Creating a dummy object that 'check_args' can work with
        class Dummy(object):
            pass
        args = Dummy()

        # Setting all the attributes (and creating the files if required)
        # bfile
        bfile = os.path.join(self.output_dir.name, "input_file")
        for extension in [".bed", ".bim", ".fam"]:
            with open(bfile + extension, "w") as o_file:
                pass
        args.bfile = bfile

        # thread
        args.thread = 1

        # shapeit_thread
        args.shapeit_thread = 1

        # hap_template, legend_template and map_template
        hap_template = os.path.join(self.output_dir.name, "chr{chrom}.hap.gz")
        leg_template = os.path.join(self.output_dir.name, "chr{chrom}.leg.gz")
        map_template = os.path.join(self.output_dir.name, "chr{chrom}.map")
        for filename in [hap_template, leg_template, map_template]:
            for chrom in chromosomes:
                with open(filename.format(chrom=chrom), "w") as o_file:
                    pass
        args.hap_template = hap_template
        args.legend_template = leg_template
        args.map_template = map_template

        # sample_file
        sample_file = os.path.join(self.output_dir.name, "sample.txt")
        with open(sample_file, "w") as o_file:
            pass
        args.sample_file = sample_file

        # shapeit_bin
        shapeit_bin = os.path.join(self.output_dir.name, "shapeit")
        with open(shapeit_bin, "w") as o_file:
            pass
        args.shapeit_bin = shapeit_bin

        # impute2_bin
        impute2_bin = os.path.join(self.output_dir.name, "impute2")
        with open(impute2_bin, "w") as o_file:
            pass
        args.impute2_bin = impute2_bin

        # plink_bin
        plink_bin = os.path.join(self.output_dir.name, "plink")
        with open(plink_bin, "w") as o_file:
            pass
        args.plink_bin = plink_bin

        # segment_length
        args.segment_length = 5e6

        # preamble
        preamble = os.path.join(self.output_dir.name, "preamble.txt")
        with open(preamble, "w") as o_file:
            pass
        args.preamble = preamble

        # use_drmaa
        args.use_drmaa = True

        # args.drmaa_config
        drmaa_config = os.path.join(self.output_dir.name, "drmaa_config.txt")
        with open(drmaa_config, "w") as o_file:
            pass
        args.drmaa_config = drmaa_config

        # Testing begins
        # Everything should work
        self.assertTrue(check_args(args))

        # Now, setting an invalid thread number (negative or 0)
        for i in range(-1, 1):
            args.thread = i
            with self.assertRaises(ProgramError) as cm:
                check_args(args)
            self.assertEqual("thread should be one or more", str(cm.exception))
        args.thread = 2

        # Now, setting an invalid shapeit_thread number (negative or 0)
        original_value = args.shapeit_thread
        for i in range(-1, 1):
            args.shapeit_thread = i
            with self.assertRaises(ProgramError) as cm:
                check_args(args)
            self.assertEqual("thread should be one or more", str(cm.exception))
        args.shapeit_thread = 2

        # Deleting each of the required plink file
        for extension in [".bed", ".bim", ".fam"]:
            os.remove(args.bfile + extension)
            self.assertFalse(os.path.isfile(args.bfile + extension))
            with self.assertRaises(ProgramError) as cm:
                check_args(args)
            self.assertEqual("{}: no such file".format(args.bfile + extension),
                             str(cm.exception))
            with open(args.bfile + extension, "w") as o_file:
                pass

        # Deleting each of the template, legend or map file
        for template in [args.hap_template, args.legend_template,
                         args.map_template]:
            for chrom in chromosomes:
                filename = template.format(chrom=chrom)
                os.remove(filename)
                self.assertFalse(os.path.isfile(filename))
                with self.assertRaises(ProgramError) as cm:
                    check_args(args)
                self.assertEqual("{}: no such file".format(filename),
                                 str(cm.exception))
                with open(filename, "w") as o_file:
                    pass

        # Deleting the sample file
        os.remove(args.sample_file)
        self.assertFalse(os.path.isfile(args.sample_file))
        with self.assertRaises(ProgramError) as cm:
            check_args(args)
        self.assertEqual("{}: no such file".format(args.sample_file),
                         str(cm.exception))
        with open(args.sample_file, "w") as o_file:
            pass

        # Deleting the shapeit dummy binary
        os.remove(args.shapeit_bin)
        self.assertFalse(os.path.isfile(args.shapeit_bin))
        with self.assertRaises(ProgramError) as cm:
            check_args(args)
        self.assertEqual("{}: no such file".format(args.shapeit_bin),
                         str(cm.exception))
        with open(args.shapeit_bin, "w") as o_file:
            pass

        # Setting the shapeit binary to None (if not in the path)
        original_value = args.shapeit_bin
        args.shapeit_bin = None
        if which("shapeit") is None:
            with self.assertRaises(ProgramError) as cm:
                check_args(args)
            self.assertEqual("shapeit: not in the path (use --shapeit-bin)",
                             str(cm.exception))
        else:
            self.assertTrue(check_args(args))
        args.shapeit_bin = original_value

        # Deleting the impute2 dummy binary
        os.remove(args.impute2_bin)
        self.assertFalse(os.path.isfile(args.impute2_bin))
        with self.assertRaises(ProgramError) as cm:
            check_args(args)
        self.assertEqual("{}: no such file".format(args.impute2_bin),
                         str(cm.exception))
        with open(args.impute2_bin, "w") as o_file:
            pass

        # Setting the shapeit binary to None (if not in the path)
        original_value = args.impute2_bin
        args.impute2_bin = None
        if which("shapeit") is None:
            with self.assertRaises(ProgramError) as cm:
                check_args(args)
            self.assertEqual("impute2: not in the path (use --impute2-bin)",
                             str(cm.exception))
        else:
            self.assertTrue(check_args(args))
        args.impute2_bin = original_value

        # Deleting the plink dummy binary
        os.remove(args.plink_bin)
        self.assertFalse(os.path.isfile(args.plink_bin))
        with self.assertRaises(ProgramError) as cm:
            check_args(args)
        self.assertEqual("{}: no such file".format(args.plink_bin),
                         str(cm.exception))
        with open(args.plink_bin, "w") as o_file:
            pass

        # Setting the shapeit binary to None (if not in the path)
        original_value = args.plink_bin
        args.plink_bin = None
        if which("shapeit") is None:
            with self.assertRaises(ProgramError) as cm:
                check_args(args)
            self.assertEqual("plink: not in the path (use --plink-bin)",
                             str(cm.exception))
        else:
            self.assertTrue(check_args(args))
        args.plink_bin = original_value

        # Modifying the segment length
        original_value = args.segment_length
        for value in [-1, 0, 1e3 - 1, 5e6 + 1]:
            args.segment_length = value
            if value == 5e6 + 1:
                # Too big
                with self._my_compatibility_assertLogs(level="WARNING") as cm:
                    check_args(args)
                log_m = ("WARNING:root:segment length ({:g} bp) is more "
                         "than 5Mb")
                self.assertEqual(1, len(cm.output))
                self.assertEqual(log_m.format(value), cm.output[0])

            elif value == 1e3 - 1:
                # Too small
                with self._my_compatibility_assertLogs(level="WARNING") as cm:
                    check_args(args)
                log_m = "WARNING:root:segment length ({:g} bp) is too small"
                self.assertEqual(1, len(cm.output))
                self.assertEqual(log_m.format(value), cm.output[0])

            else:
                # Invalid
                with self.assertRaises(ProgramError) as cm:
                    check_args(args)
                error_m = "{}: invalid segment length"
                self.assertEqual(error_m.format(value), str(cm.exception))
        args.segment_length = original_value

        # Deleting the preamble file
        os.remove(args.preamble)
        self.assertFalse(os.path.isfile(args.preamble))
        with self.assertRaises(ProgramError) as cm:
            check_args(args)
        self.assertEqual("{}: no such file".format(args.preamble),
                         str(cm.exception))

        # Setting the preamble to None should do the trick
        args.preamble = None
        self.assertTrue(check_args(args))

        # Deleting the drmaa config file
        os.remove(args.drmaa_config)
        self.assertFalse(os.path.isfile(args.drmaa_config))
        with self.assertRaises(ProgramError) as cm:
            check_args(args)
        self.assertEqual("{}: no such file".format(args.drmaa_config),
                         str(cm.exception))

        # Setting to None should do the tick
        original_value = args.drmaa_config
        args.drmaa_config = None
        self.assertTrue(
            check_args(args) and not os.path.isfile(original_value)
        )
        args.drmaa_config = original_value

        # Or, setting use_drmaa to false
        args.use_drmaa = False
        self.assertTrue(
            check_args(args) and not os.path.isfile(args.drmaa_config)
        )
        with open(args.drmaa_config, "w") as o_file:
            pass

        # Final check
        self.assertTrue(check_args(args))

    def _my_compatibility_assertLogs(self, logger=None, level=None):
        """Compatibility 'assertLogs' function for Python 3.3."""
        if hasattr(self, "assertLogs"):
            return self.assertLogs(logger, level)

        else:
            from .python_3_3_compatibility import Python_3_4_AssertLogsContext
            return Python_3_4_AssertLogsContext(self, logger, level)
