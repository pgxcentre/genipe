
# This file is part of genipe.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import os
import unittest
from random import randint
from tempfile import TemporaryDirectory

from .. import autosomes
from ..pipeline.cli import *
from ..error import GenipeError

if HAS_PYFAIDX:
    import pyfaidx


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["TestMainPipeline"]


class TestArguments(unittest.TestCase):

    def setUp(self):
        """Setup the tests."""
        # Creating the temporary directory
        self.output_dir = TemporaryDirectory(prefix="genipe_test_")

        # Creating a complete argument object
        class Dummy(object):
            pass
        self.args = Dummy()

        # The required chromosomes
        self.args.required_chrom = (1, 2, 3, 6, 9, 23, 25)

        # Setting all the attributes (and creating the files if required)
        # bfile
        bfile = os.path.join(self.output_dir.name, "input_file")
        for extension in [".bed", ".bim", ".fam"]:
            with open(bfile + extension, "w") as o_file:
                pass
        self.args.bfile = bfile

        # thread
        self.args.thread = 1

        # shapeit_thread
        self.args.shapeit_thread = 1

        # hap_template, legend_template and map_template
        hap_template = os.path.join(self.output_dir.name, "chr{chrom}.hap.gz")
        leg_template = os.path.join(self.output_dir.name, "chr{chrom}.leg.gz")
        map_template = os.path.join(self.output_dir.name, "chr{chrom}.map")
        for filename in [hap_template, leg_template, map_template]:
            for chrom in self.args.required_chrom + ("23_PAR1", "23_PAR2"):
                with open(filename.format(chrom=chrom), "w") as o_file:
                    pass
        self.args.hap_template = hap_template
        self.args.legend_template = leg_template
        self.args.map_template = map_template

        # hap, legend and map for chromosome 23 (non pseudo-autosomal region)
        self.args.hap_chr23 = hap_template.format(chrom=23)
        self.args.legend_chr23 = leg_template.format(chrom=23)
        self.args.map_chr23 = map_template.format(chrom=23)

        # hap, legend and map for chromosome 23 (first pseudo-autosomal region)
        self.args.hap_par1 = hap_template.format(chrom="23_PAR1")
        self.args.legend_par1 = leg_template.format(chrom="23_PAR1")
        self.args.map_par1 = map_template.format(chrom="23_PAR1")

        # hap, legend and map for chromosome 23 (first pseudo-autosomal region)
        self.args.hap_par2 = hap_template.format(chrom="23_PAR2")
        self.args.legend_par2 = leg_template.format(chrom="23_PAR2")
        self.args.map_par2 = map_template.format(chrom="23_PAR2")

        # sample_file
        sample_file = os.path.join(self.output_dir.name, "sample.txt")
        with open(sample_file, "w") as o_file:
            pass
        self.args.sample_file = sample_file

        # shapeit_bin
        shapeit_bin = os.path.join(self.output_dir.name, "shapeit")
        with open(shapeit_bin, "w") as o_file:
            pass
        self.args.shapeit_bin = shapeit_bin

        # impute2_bin
        impute2_bin = os.path.join(self.output_dir.name, "impute2")
        with open(impute2_bin, "w") as o_file:
            pass
        self.args.impute2_bin = impute2_bin

        # plink_bin
        plink_bin = os.path.join(self.output_dir.name, "plink")
        with open(plink_bin, "w") as o_file:
            pass
        self.args.plink_bin = plink_bin

        # segment_length
        self.args.segment_length = 5e6

        # preamble
        preamble = os.path.join(self.output_dir.name, "preamble.txt")
        with open(preamble, "w") as o_file:
            pass
        self.args.preamble = preamble

        # use_drmaa
        self.args.use_drmaa = True

        # drmaa_config
        drmaa_config = os.path.join(self.output_dir.name, "drmaa_config.txt")
        with open(drmaa_config, "w") as o_file:
            pass
        self.args.drmaa_config = drmaa_config

        # reference
        reference = os.path.join(self.output_dir.name, "h_ref.fasta")
        for filename in [reference, reference + ".fai"]:
            with open(filename, "w") as o_file:
                pass
        self.args.reference = reference

        # bgzip
        self.args.bgzip = False

    def tearDown(self):
        """Finishes the test."""
        # Deleting the output directory
        self.output_dir.cleanup()

    def test_normal(self):
        """Tests under normal utilization"""
        self.assertTrue(check_args(self.args))

    def test_invalid_thread(self):
        """Tests with an invalid thread number."""
        for i in range(-1, 1):
            self.args.thread = i
            with self.assertRaises(GenipeError) as cm:
                check_args(self.args)
            self.assertEqual("thread should be one or more", str(cm.exception))

    def test_invalid_shapeit_thread(self):
        """Tests with an invalid shapeit thread number."""
        # Now, setting an invalid shapeit_thread number (negative or 0)
        for i in range(-1, 1):
            self.args.shapeit_thread = i
            with self.assertRaises(GenipeError) as cm:
                check_args(self.args)
            self.assertEqual("thread should be one or more", str(cm.exception))

    def test_missing_plink_file(self):
        """Tests with missing Plink file."""
        # Deleting each of the required plink file
        for ext in (".bed", ".bim", ".fam"):
            os.remove(self.args.bfile + ext)
            self.assertFalse(os.path.isfile(self.args.bfile + ext))

            # Checking the exception is raised
            with self.assertRaises(GenipeError) as cm:
                check_args(self.args)
            self.assertEqual("{}: no such file".format(self.args.bfile + ext),
                             str(cm.exception))

            # Creating back the file
            with open(self.args.bfile + ext, "w") as o_file:
                pass
            self.assertTrue(os.path.isfile(self.args.bfile + ext))

    def test_missing_template_autosomes(self):
        """Tests with missing template file (autosome)."""
        # Deleting each of the template files
        for chrom in self.args.required_chrom:
            if chrom not in autosomes:
                continue

            # Deleting the file
            filename = self.args.legend_template.format(chrom=chrom)
            os.remove(filename)
            self.assertFalse(os.path.isfile(filename))

            # Checking the exception is raised
            with self.assertRaises(GenipeError) as cm:
                check_args(self.args)
            self.assertEqual("{}: no such file".format(filename),
                             str(cm.exception))

            # Creating back the file
            with open(filename, "w") as o_file:
                pass
            self.assertTrue(os.path.isfile(filename))

    def test_missing_template_autosomes_no_autosomes(self):
        """Tests with missing template file (autosome, but none required)."""
        # Asking just for sexual chromosomes
        self.args.required_chrom = [
            chrom
            for chrom in self.args.required_chrom
            if chrom not in autosomes
        ]

        # Deleting
        for chrom in self.args.required_chrom:
            if chrom not in autosomes:
                continue

            filename = self.args.legend_template.format(chrom=chrom)
            os.remove(filename)
            self.assertFalse(os.path.isfile(filename))

        # Checking everything works fine
        self.assertTrue(check_args(self.args))

    def test_missing_template_chr23(self):
        """Tests with missing template file (non pseudo-autosomal region)."""
        # Deleting the legend file
        os.remove(self.args.legend_chr23)
        self.assertFalse(os.path.isfile(self.args.legend_chr23))

        # Checking the exception is raised
        with self.assertRaises(GenipeError) as cm:
            check_args(self.args)
        self.assertEqual("{}: no such file".format(self.args.legend_chr23),
                         str(cm.exception))

    def test_missing_template_chr23_no_chr23(self):
        """Tests with missing template file (chr23, but not asked)."""
        self.args.required_chrom = [
            chrom
            for chrom in self.args.required_chrom
            if chrom != 23
        ]

        # Deleting the legend file
        os.remove(self.args.legend_chr23)
        self.assertFalse(os.path.isfile(self.args.legend_chr23))

        # Checking everything works
        self.assertTrue(check_args(self.args))

    def test_missing_template_chr25_par1(self):
        """Tests with missing template file (pseudo-autosomal region 1)."""
        # Deleting the legend file
        os.remove(self.args.legend_par1)
        self.assertFalse(os.path.isfile(self.args.legend_par1))

        # Checking the exception is raised
        with self.assertRaises(GenipeError) as cm:
            check_args(self.args)
        self.assertEqual("{}: no such file".format(self.args.legend_par1),
                         str(cm.exception))

    def test_missing_template_chr25_par1_no_chr25_par1(self):
        """Tests with missing template file (PAR1, but not asked)."""
        self.args.required_chrom = [
            chrom
            for chrom in self.args.required_chrom
            if chrom != 25
        ]

        # Deleting the legend file
        os.remove(self.args.legend_par1)
        self.assertFalse(os.path.isfile(self.args.legend_par1))

        # Checking everything works
        self.assertTrue(check_args(self.args))

    def test_missing_template_chr25_par2(self):
        """Tests with missing template file (pseudo-autosomal region 2)."""
        # Deleting the legend file
        os.remove(self.args.legend_par2)
        self.assertFalse(os.path.isfile(self.args.legend_par2))

        # Checking the exception is raised
        with self.assertRaises(GenipeError) as cm:
            check_args(self.args)
        self.assertEqual("{}: no such file".format(self.args.legend_par2),
                         str(cm.exception))

    def test_missing_template_chr25_par2_no_chr25_par2(self):
        """Tests with missing template file (PAR2, but not asked)."""
        self.args.required_chrom = [
            chrom
            for chrom in self.args.required_chrom
            if chrom != 25
        ]

        # Deleting the legend file
        os.remove(self.args.legend_par2)
        self.assertFalse(os.path.isfile(self.args.legend_par2))

        # Checking everything works
        self.assertTrue(check_args(self.args))

    def test_missing_map_autosomes(self):
        """Tests with missing map file (autosome)."""
        # Deleting each of the template, legend or map file
        for chrom in self.args.required_chrom:
            if chrom not in autosomes:
                continue

            # Deleting the file
            filename = self.args.map_template.format(chrom=chrom)
            os.remove(filename)
            self.assertFalse(os.path.isfile(filename))

            # Checking the exception is raised
            with self.assertRaises(GenipeError) as cm:
                check_args(self.args)
            self.assertEqual("{}: no such file".format(filename),
                             str(cm.exception))

            # Creating back the file
            with open(filename, "w") as o_file:
                pass
            self.assertTrue(os.path.isfile(filename))

    def test_missing_map_autosomes_no_autosomes(self):
        """Tests with missing map file (autosome), but none required."""
        # Deleting each of the template, legend or map file
        self.args.required_chrom = [
            chrom
            for chrom in self.args.required_chrom
            if chrom not in autosomes
        ]

        # Deleting
        for chrom in self.args.required_chrom:
            if chrom not in autosomes:
                continue

            filename = self.args.map_template.format(chrom=chrom)
            os.remove(filename)
            self.assertFalse(os.path.isfile(filename))

        # Checking everything works
        self.assertTrue(check_args(self.args))

    def test_missing_hap_autosomes(self):
        """Tests with missing haplotype file (autosome)."""
        # Deleting each of the template, legend or map file
        for chrom in self.args.required_chrom:
            if chrom not in autosomes:
                continue

            # Deleting the file
            filename = self.args.hap_template.format(chrom=chrom)
            os.remove(filename)
            self.assertFalse(os.path.isfile(filename))

            # Checking the exception is raised
            with self.assertRaises(GenipeError) as cm:
                check_args(self.args)
            self.assertEqual("{}: no such file".format(filename),
                             str(cm.exception))

            # Creating back the file
            with open(filename, "w") as o_file:
                pass
            self.assertTrue(os.path.isfile(filename))

    def test_missing_hap_autosomes_no_autosomes(self):
        """Tests with missing haplotype file (autosome), none required."""
        self.args.required_chrom = [
            chrom
            for chrom in self.args.required_chrom
            if chrom not in autosomes
        ]

        # Deleting
        for chrom in self.args.required_chrom:
            if chrom not in autosomes:
                continue

            # Deleting the file
            filename = self.args.hap_template.format(chrom=chrom)
            os.remove(filename)
            self.assertFalse(os.path.isfile(filename))

        # Checking everything works
        self.assertTrue(check_args(self.args))

    def test_missing_sample_file(self):
        """Tests with missing sample file."""
        # Deleting the sample file
        os.remove(self.args.sample_file)
        self.assertFalse(os.path.isfile(self.args.sample_file))

        # Checking the exception is raised
        with self.assertRaises(GenipeError) as cm:
            check_args(self.args)
        self.assertEqual("{}: no such file".format(self.args.sample_file),
                         str(cm.exception))

    def test_missing_shapeit_binary(self):
        """Tests with missing shapeit binary file."""
        # Deleting the shapeit dummy binary
        os.remove(self.args.shapeit_bin)
        self.assertFalse(os.path.isfile(self.args.shapeit_bin))

        # Checking the exception is raised
        with self.assertRaises(GenipeError) as cm:
            check_args(self.args)
        self.assertEqual("{}: no such file".format(self.args.shapeit_bin),
                         str(cm.exception))

    def test_shapeit_in_path(self):
        """Tests what happens if shapeit needs to be in the path."""
        # Setting the shapeit binary to None (if not in the path)
        self.args.shapeit_bin = None

        if which("shapeit") is None:
            with self.assertRaises(GenipeError) as cm:
                check_args(self.args)
            self.assertEqual("shapeit: not in the path (use --shapeit-bin)",
                             str(cm.exception))

        else:
            self.assertTrue(check_args(self.args))

    def test_missing_impute2_binary(self):
        """Tests with missing impute2 dummy binary."""
        # Deleting the impute2 dummy binary
        os.remove(self.args.impute2_bin)
        self.assertFalse(os.path.isfile(self.args.impute2_bin))

        # Checking the exception is raised
        with self.assertRaises(GenipeError) as cm:
            check_args(self.args)
        self.assertEqual("{}: no such file".format(self.args.impute2_bin),
                         str(cm.exception))

    def test_impute2_in_path(self):
        """Tests what happens if impute2 needs to be in the path."""
        # Setting the impute2 binary to None (if not in the path)
        self.args.impute2_bin = None

        if which("impute2") is None:
            with self.assertRaises(GenipeError) as cm:
                check_args(self.args)
            self.assertEqual("impute2: not in the path (use --impute2-bin)",
                             str(cm.exception))

        else:
            self.assertTrue(check_args(self.args))

    def test_missing_plink_binary(self):
        """Tests with missing plink dummy binary."""
        # Deleting the plink dummy binary
        os.remove(self.args.plink_bin)
        self.assertFalse(os.path.isfile(self.args.plink_bin))

        # Checking the exception is raised
        with self.assertRaises(GenipeError) as cm:
            check_args(self.args)
        self.assertEqual("{}: no such file".format(self.args.plink_bin),
                         str(cm.exception))

    def test_plink_in_path(self):
        """Tests what happens if plink needs to be in the path."""
        # Setting the plink binary to None (if not in the path)
        self.args.plink_bin = None

        if which("plink") is None:
            with self.assertRaises(GenipeError) as cm:
                check_args(self.args)
            self.assertEqual("plink: not in the path (use --plink-bin)",
                             str(cm.exception))

        else:
            self.assertTrue(check_args(self.args))

    def test_bgzip_in_path(self):
        """Tests what happens if bgzip needs to be in the path."""
        # If bgzip is not in the path, it should raise an error
        self.args.bgzip = True

        if which("bgzip") is None:
            with self.assertRaises(GenipeError) as cm:
                check_args(self.args)
            self.assertEqual("bgzip: no installed", str(cm.exception))

        else:
            self.assertTrue(check_args(self.args))

    def test_segment_length_large(self):
        """Tests different invalid segment length (too large)."""
        self.args.segment_length = 5e6 + 1

        # Checking the warning is logged
        with self._my_compatibility_assertLogs(level="WARNING") as cm:
            check_args(self.args)
        log_m = ["WARNING:root:segment length (5e+06 bp) is more than 5Mb"]
        self.assertEqual(log_m, cm.output)

    def test_segment_length_small(self):
        """Tests different invalid segment length (too small)."""
        self.args.segment_length = 1e3 - 1

        # Checking the warning is logged
        with self._my_compatibility_assertLogs(level="WARNING") as cm:
            check_args(self.args)
        log_m = ["WARNING:root:segment length (999 bp) is too small"]
        self.assertEqual(log_m, cm.output)

    def test_segment_length_invalid(self):
        """Tests different invalid segment length."""
        # Checking the exception is raised
        self.args.segment_length = -1
        with self.assertRaises(GenipeError) as cm:
            check_args(self.args)
        self.assertEqual("-1: invalid segment length", str(cm.exception))

        # Checking the exception is raised
        self.args.segment_length = 0
        with self.assertRaises(GenipeError) as cm:
            check_args(self.args)
        self.assertEqual("0: invalid segment length", str(cm.exception))

    def test_missing_preamble(self):
        """Tests with a missing preamble file."""
        # Deleting the preamble file
        os.remove(self.args.preamble)
        self.assertFalse(os.path.isfile(self.args.preamble))

        # Checking the exception is raised
        with self.assertRaises(GenipeError) as cm:
            check_args(self.args)
        self.assertEqual("{}: no such file".format(self.args.preamble),
                         str(cm.exception))

    def test_missing_preamble_none(self):
        """Tests with a missing preamble file, but with a None setting."""
        # Deleting the preamble file
        os.remove(self.args.preamble)
        self.assertFalse(os.path.isfile(self.args.preamble))

        # Setting the preamble to None should do the trick
        self.args.preamble = None
        self.assertTrue(check_args(self.args))

    def test_missing_drmaa_config(self):
        """Tests with a missing DRMAA config file."""
        # Deleting the drmaa config file
        os.remove(self.args.drmaa_config)
        self.assertFalse(os.path.isfile(self.args.drmaa_config))

        # Checking the exception is raised
        with self.assertRaises(GenipeError) as cm:
            check_args(self.args)
        self.assertEqual("{}: no such file".format(self.args.drmaa_config),
                         str(cm.exception))

    def test_missing_drmaa_config_none(self):
        """Tests with a missing DRMAA config file, but with a None setting."""
        # Removing the DRMAA config file
        os.remove(self.args.drmaa_config)
        self.assertFalse(os.path.isfile(self.args.drmaa_config))
        self.args.drmaa_config = None

        # Checking the exception is raised
        with self.assertRaises(GenipeError) as cm:
            check_args(self.args)
        self.assertEqual(
            "DRMAA configuration file was not provided (--drmaa-config), but "
            "DRMAA is used (--use-drmaa)",
            str(cm.exception)
        )

    def test_missing_drmaa_config_none_false(self):
        """Tests with a missing drmaa config file, but None and False."""
        # Changing the configuration
        os.remove(self.args.drmaa_config)
        self.assertFalse(os.path.isfile(self.args.drmaa_config))
        self.args.drmaa_config = None
        self.args.use_drmaa = False

        # Checking everything is fine
        self.assertTrue(check_args(self.args))

    @unittest.skipIf(not HAS_PYFAIDX, "no module 'pyfaidx'")
    def test_missing_reference_index_file(self):
        """Tests with a missing index for the reference file."""
        os.remove(self.args.reference + ".fai")
        self.assertFalse(os.path.isfile(self.args.reference + ".fai"))

        # Checking the exception is raised
        with self.assertRaises(GenipeError) as cm:
            check_args(self.args)
        self.assertEqual(
            "{}: should be indexed using FAIDX".format(
                self.args.reference,
            ),
            str(cm.exception),
        )

    @unittest.skipIf(not HAS_PYFAIDX, "no module 'pyfaidx'")
    def test_missing_reference_index_file_none(self):
        """Tests with a missing index for the reference file (with None)."""
        os.remove(self.args.reference + ".fai")
        self.assertFalse(os.path.isfile(self.args.reference + ".fai"))
        self.args.reference = None

        # Checking everything went fine
        self.assertTrue(check_args(self.args))

    @unittest.skipIf(not HAS_PYFAIDX, "no module 'pyfaidx'")
    def test_missing_reference_file(self):
        """Tests with a missing reference file."""
        # Removing the reference file should raise an exception
        os.remove(self.args.reference)
        self.assertFalse(os.path.isfile(self.args.reference))

        # Checking the exception is raised
        with self.assertRaises(GenipeError) as cm:
            check_args(self.args)
        self.assertEqual("{}: no such file".format(self.args.reference),
                         str(cm.exception))

    @unittest.skipIf(not HAS_PYFAIDX, "no module 'pyfaidx'")
    def test_missing_reference_file_none(self):
        """Tests with a missing reference file (with None)."""
        # Removing the reference file should raise an exception
        os.remove(self.args.reference)
        self.assertFalse(os.path.isfile(self.args.reference))
        self.args.reference = None

        # Everything is fine
        self.assertTrue(check_args(self.args))

    def _my_compatibility_assertLogs(self, logger=None, level=None):
        """Compatibility 'assertLogs' function for Python 3.3."""
        if hasattr(self, "assertLogs"):
            return self.assertLogs(logger, level)

        else:
            from .python_3_3_compatibility import Python_3_4_AssertLogsContext
            return Python_3_4_AssertLogsContext(self, logger, level)
