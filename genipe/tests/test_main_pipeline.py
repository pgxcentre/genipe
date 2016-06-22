
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
from ..pipeline import cli
from ..error import GenipeError

if cli.HAS_PYFAIDX:
    import pyfaidx


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["TestMainPipeline"]


class TestMainPipeline(unittest.TestCase):

    def setUp(self):
        """Setup the tests."""
        # Creating the temporary directory
        self.output_dir = TemporaryDirectory(prefix="genipe_test_")

    def tearDown(self):
        """Finishes the test."""
        # Deleting the output directory
        self.output_dir.cleanup()

    def test_file_sorter(self):
        """Tests the 'file_sorter' function."""
        filenames = [
            "chr1.1000_100000.impute2",
            "chr2.10000_1002300.impute2",
            "chr1.1_100.impute2",
            "/some/path/chr1.1_10.impute2",
            "chr25_1.3_40.impute2",
            "chr23.100_400.impute2",
            "some/path/to_file/chr25_2.1000_1500.impute2.gz",
            "chr25_2.1600_1650.impute2.gz",
            "chr1.100000_2000000.impute2.some_extension",
        ]
        expected_filenames = [
            "/some/path/chr1.1_10.impute2",
            "chr1.1_100.impute2",
            "chr1.1000_100000.impute2",
            "chr1.100000_2000000.impute2.some_extension",
            "chr2.10000_1002300.impute2",
            "chr23.100_400.impute2",
            "chr25_1.3_40.impute2",
            "some/path/to_file/chr25_2.1000_1500.impute2.gz",
            "chr25_2.1600_1650.impute2.gz",
        ]
        expected_results = [(1, 1000, 100000), (2, 10000, 1002300),
                            (1, 1, 100), (1, 1, 10), (25, 3, 40),
                            (23, 100, 400), (25, 1000, 1500),
                            (25, 1600, 1650),
                            (1, 100000, 2000000)]

        # Trying the function
        for filename, expected in zip(filenames, expected_results):
            self.assertEqual(expected, cli.file_sorter(filename))

        # Trying the sort function
        filenames.sort(key=cli.file_sorter)
        self.assertEqual(filenames, expected_filenames)

    def test_get_chromosome_length(self):
        """Tests the 'get_chromosome_length' function."""
        # The expected chromosome
        expected_chrom = {3, 6, 9, 23, 25}
        expected_length = {}

        # Writing some legend file for different chromosome
        legend_template = os.path.join(self.output_dir.name,
                                       "chr{chrom}.legend")
        legend_chr23 = os.path.join(self.output_dir.name,
                                    "chr23_nonPAR.legend")
        legend_par1 = os.path.join(self.output_dir.name, "chr23_PAR1.legend")
        legend_par2 = os.path.join(self.output_dir.name, "chr23_PAR2.legend")
        for chrom in expected_chrom:
            if chrom == 23:
                # Getting the positions and expected length
                positions = sorted([randint(5000, 100000) for i in range(100)])
                expected_length[chrom] = (min(positions), max(positions))

                # Writing the positions to file
                with open(legend_chr23, "w") as o_file:
                    print("id", "position", file=o_file)
                    for i, position in enumerate(positions):
                        print("marker_{}".format(i+1), position, file=o_file)

                # Continuing to next chromosome
                continue

            if chrom == 25:
                # Getting the positions and expected length
                positions = sorted([randint(1, 4999) for i in range(100)])
                expected_length[chrom] = [max(positions)]

                # Writing the positions to file
                with open(legend_par1, "w") as o_file:
                    print("id", "position", file=o_file)
                    for i, position in enumerate(positions):
                        print("marker_{}".format(i+1), position, file=o_file)

                # Getting the positions and expected length
                positions = sorted([
                    randint(100000, 200000) for i in range(100)
                ])
                expected_length[chrom].extend([min(positions), max(positions)])
                expected_length[chrom] = tuple(expected_length[chrom])

                # Writing the positions to file
                with open(legend_par2, "w") as o_file:
                    print("id", "position", file=o_file)
                    for i, position in enumerate(positions):
                        print("marker_{}".format(i+1), position, file=o_file)

                # Continuing to next chromosome
                continue

            # Getting the position and saving the expected length
            positions = sorted([randint(1, 2000000) for i in range(1000)])
            expected_length[chrom] = max(positions)

            # Saving the file
            with open(legend_template.format(chrom=chrom), "w") as o_file:
                print("id", "position", file=o_file)
                for i, position in enumerate(positions):
                    print("marker_{}".format(i+1), position, file=o_file)

        # Getting the chromosome length
        chrom_length = cli.get_chromosome_length(
            required_chrom=expected_chrom,
            legend=legend_template,
            legend_chr23=legend_chr23,
            legend_par1=legend_par1,
            legend_par2=legend_par2,
            out_dir=self.output_dir.name,
        )

        # Checking the expected length
        self.assertEqual(expected_length, chrom_length)

        # Checking the file was created
        self.assertTrue(os.path.isfile(os.path.join(self.output_dir.name,
                                                    "chromosome_lengths.txt")))

        # Tests that we correctly read the file
        expected_chrom = {6: expected_length[6], 9: expected_length[9],
                          23: expected_length[23], 25: expected_length[25]}

        # Writing the file
        chrom_filename = os.path.join(self.output_dir.name,
                                      "chromosome_lengths.txt")
        with open(chrom_filename, "w") as o_file:
            for chrom, length in expected_chrom.items():
                if (chrom == 23) or (chrom == 25):
                    print(chrom, *length, sep="\t", file=o_file)
                else:
                    print(chrom, length, sep="\t", file=o_file)

        # Comparing what we got
        chrom_length = cli.get_chromosome_length(
            required_chrom=expected_chrom.keys(),
            legend=legend_template,
            legend_chr23=legend_chr23,
            legend_par1=legend_par1,
            legend_par2=legend_par2,
            out_dir=self.output_dir.name,
        )
        self.assertEqual(expected_chrom, chrom_length)

        # Removing some autosomes from the file
        del expected_chrom[9]
        del expected_chrom[23]
        del expected_chrom[25]
        with open(chrom_filename, "w") as o_file:
            for chrom, length in expected_chrom.items():
                if (chrom == 23) or (chrom == 25):
                    print(chrom, *length, sep="\t", file=o_file)
                else:
                    print(chrom, length, sep="\t", file=o_file)
        expected_chrom[9] = expected_length[9]
        expected_chrom[23] = expected_length[23]
        expected_chrom[25] = expected_length[25]

        # Tests that a warning is logged if there is a missing chromosome
        with self.assertLogs(level="WARNING") as cm:
            chrom_length = cli.get_chromosome_length(
                required_chrom=sorted(expected_chrom.keys()),
                legend=legend_template,
                legend_chr23=legend_chr23,
                legend_par1=legend_par1,
                legend_par2=legend_par2,
                out_dir=self.output_dir.name,
            )
        log_m = [
            "WARNING:root:missing length for chromosome 9",
            "WARNING:root:missing length for chromosome 23",
            "WARNING:root:missing length for chromosome 25",
        ]
        self.assertEqual(log_m, cm.output)

        # Testing the content
        self.assertEqual(expected_chrom, chrom_length)

    @unittest.skipIf(not cli.HAS_PYFAIDX,
                     "optional requirement (pyfaidx) not satisfied")
    def test_get_chrom_encoding(self):
        """Tests the 'get_chrom_encoding' function."""
        # Creating the reference file (fasta file) and index (using samtools)
        fasta_content = [[">{}".format(i), "ACGT"] for i in range(1, 25)]
        fasta_content.append([">26", "ACGT"])
        fasta_content.append([">Unaligned", "ACGT"])
        index_content = [
            ["1", "4", "3", "4", "5"],
            ["2", "4", "11", "4", "5"],
            ["3", "4", "19", "4", "5"],
            ["4", "4", "27", "4", "5"],
            ["5", "4", "35", "4", "5"],
            ["6", "4", "43", "4", "5"],
            ["7", "4", "51", "4", "5"],
            ["8", "4", "59", "4", "5"],
            ["9", "4", "67", "4", "5"],
            ["10", "4", "76", "4", "5"],
            ["11", "4", "85", "4", "5"],
            ["12", "4", "94", "4", "5"],
            ["13", "4", "103", "4", "5"],
            ["14", "4", "112", "4", "5"],
            ["15", "4", "121", "4", "5"],
            ["16", "4", "130", "4", "5"],
            ["17", "4", "139", "4", "5"],
            ["18", "4", "148", "4", "5"],
            ["19", "4", "157", "4", "5"],
            ["20", "4", "166", "4", "5"],
            ["21", "4", "175", "4", "5"],
            ["22", "4", "184", "4", "5"],
            ["23", "4", "193", "4", "5"],
            ["24", "4", "202", "4", "5"],
            ["26", "4", "210", "4", "5"],
            ["Unaligned", "4", "226", "4", "5"],
        ]
        reference_filename = os.path.join(self.output_dir.name, "ref.fasta")
        with open(reference_filename, "w") as o_file:
            for chromosome_fasta in fasta_content:
                print(*chromosome_fasta, sep="\n", file=o_file)
        with open(reference_filename + ".fai", "w") as o_file:
            for line in index_content:
                print(*line, sep="\t", file=o_file)

        # Reading the reference using pyfaidx
        reference = pyfaidx.Fasta(reference_filename, as_raw=True)

        # The expected result
        expected = {str(i): str(i) for i in range(1, 25)}
        expected["26"] = "26"

        # The observed result
        observed = cli.get_chrom_encoding(reference)
        self.assertEqual(expected, observed)
        reference.close()

        # Replacing 23 and 24 to X and Y
        fasta_content[22][0] = ">X"
        fasta_content[23][0] = ">Y"
        index_content[22][0] = "X"
        index_content[23][0] = "Y"

        # Writing to file
        with open(reference_filename, "w") as o_file:
            for chromosome_fasta in fasta_content:
                print(*chromosome_fasta, sep="\n", file=o_file)
        with open(reference_filename + ".fai", "w") as o_file:
            for line in index_content:
                print(*line, sep="\t", file=o_file)

        # Reading the reference using pyfaidx
        reference = pyfaidx.Fasta(reference_filename, as_raw=True)

        # The expected result
        expected["23"] = "X"
        expected["24"] = "Y"

        # The observed result
        observed = cli.get_chrom_encoding(reference)
        self.assertEqual(expected, observed)
        reference.close()

        # Adding chr everywhere
        for i in range(len(fasta_content)):
            fasta_content[i][0] = ">chr" + fasta_content[i][0][1:]
            index_content[i][0] = "chr" + index_content[i][0]

        # Writing to file
        with open(reference_filename, "w") as o_file:
            for chromosome_fasta in fasta_content:
                print(*chromosome_fasta, sep="\n", file=o_file)
        with open(reference_filename + ".fai", "w") as o_file:
            for line in index_content:
                print(*line, sep="\t", file=o_file)

        # Reading the reference using pyfaidx
        reference = pyfaidx.Fasta(reference_filename, as_raw=True)

        # The expected result
        expected = {str(i): "chr{}".format(i) for i in range(1, 23)}
        expected["23"] = "chrX"
        expected["24"] = "chrY"
        expected["26"] = "chr26"

        # The observed result
        observed = cli.get_chrom_encoding(reference)
        self.assertEqual(expected, observed)
        reference.close()

        # Replacing 23 and 24 to X and Y
        fasta_content[22][0] = ">chr23"
        fasta_content[23][0] = ">chr24"
        index_content[22][0] = "chr23"
        index_content[23][0] = "chr24"

        # Writing to file
        with open(reference_filename, "w") as o_file:
            for chromosome_fasta in fasta_content:
                print(*chromosome_fasta, sep="\n", file=o_file)
        with open(reference_filename + ".fai", "w") as o_file:
            for line in index_content:
                print(*line, sep="\t", file=o_file)

        # Reading the reference using pyfaidx
        reference = pyfaidx.Fasta(reference_filename, as_raw=True)

        # The expected result
        expected["23"] = "chr23"
        expected["24"] = "chr24"

        # The observed result
        observed = cli.get_chrom_encoding(reference)
        self.assertEqual(expected, observed)
        reference.close()

        # Finally, removing from chromosome 18 to trigger warnings
        fasta_content = fasta_content[:18]
        index_content = index_content[:18]

        # Writing to file
        with open(reference_filename, "w") as o_file:
            for chromosome_fasta in fasta_content:
                print(*chromosome_fasta, sep="\n", file=o_file)
        with open(reference_filename + ".fai", "w") as o_file:
            for line in index_content:
                print(*line, sep="\t", file=o_file)

        # Reading the reference using pyfaidx
        reference = pyfaidx.Fasta(reference_filename, as_raw=True)

        # The observed result
        with self.assertLogs(level="WARNING") as cm:
            cli.get_chrom_encoding(reference)
        log_m = [
            "WARNING:root:{}: chromosome not in reference".format(i)
            for i in range(19, 27) if i != 25
        ]
        self.assertEqual(log_m, cm.output)
        reference.close()

    @unittest.skipIf(not cli.HAS_PYFAIDX,
                     "optional requirement (pyfaidx) not satisfied")
    def test_is_reversed(self):
        """Tests the 'is_reversed' function."""
        # Creating the reference file (fasta file) and index (using samtools)
        fasta_content = (
            ">1\n"
            "ACGT\n"
            ">2\n"
            "ACGT\n"
            ">3\n"
            "acgt\n"
        )
        index_content = (
            "1\t4\t3\t4\t5\n"
            "2\t4\t11\t4\t5\n"
            "3\t4\t19\t4\t5\n"
        )
        reference_filename = os.path.join(self.output_dir.name, "ref.fasta")
        with open(reference_filename, "w") as o_file:
            o_file.write(fasta_content)
        with open(reference_filename + ".fai", "w") as o_file:
            o_file.write(index_content)

        # Reading the reference using pyfaidx
        reference = pyfaidx.Fasta(reference_filename, as_raw=True)

        # The chromosome encoding
        encoding = {"1": "1", "2": "2", "3": "3"}

        # Testing invalid allele (should return False)
        self.assertFalse(cli.is_reversed(
            "1", 1, "I", "D", reference, encoding),
        )
        self.assertFalse(cli.is_reversed(
            "1", 1, "Z", "A", reference, encoding),
        )
        self.assertFalse(cli.is_reversed(
            "1", 1, "A", "K", reference, encoding),
        )

        # Testing invalid chromosome (should return False)
        self.assertFalse(cli.is_reversed(
            "23", 1, "A", "C", reference, encoding),
        )

        # Testing invalid position (should return False)
        self.assertFalse(cli.is_reversed(
            "1", 100, "A", "C", reference, encoding),
        )

        # Testing valid input, without strand problem (should return False)
        self.assertFalse(cli.is_reversed(
            "1", 3, "G", "T", reference, encoding),
        )
        self.assertFalse(cli.is_reversed(
            "2", 4, "G", "T", reference, encoding),
        )
        self.assertFalse(cli.is_reversed(
            "3", 2, "g", "c", reference, encoding),
        )

        # Testing valid input, but strand problem (should return True)
        self.assertTrue(cli.is_reversed("1", 1, "T", "G", reference, encoding))
        self.assertTrue(cli.is_reversed("2", 2, "t", "g", reference, encoding))
        self.assertTrue(cli.is_reversed("3", 3, "T", "C", reference, encoding))
        self.assertTrue(cli.is_reversed("1", 4, "A", "C", reference, encoding))

        # Closing the reference
        reference.close()

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
        frequency_barh = ""
        if cli.HAS_MATPLOTLIB:
            frequency_barh = os.path.join(self.output_dir.name,
                                          "frequency_barh.pdf")
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
            "frequency_barh":       frequency_barh,
        }

        # Creating the files for the test
        filename_template = os.path.join(self.output_dir.name, "chr{chrom}",
                                         "final_impute2",
                                         "chr{chrom}.imputed.{suffix}")
        for chrom in autosomes:
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
        observed = cli.gather_maf_stats(
            required_chrom=autosomes,
            o_dir=self.output_dir.name,
        )

        # Checking the observed results
        self.assertEqual(len(expected_results), len(observed))
        for expected_key, expected_value in expected_results.items():
            self.assertTrue(expected_key in observed)
            self.assertEqual(expected_value, observed[expected_key])

        # If matplotlib is installed, checking we have a figure (and not
        # otherwise)
        if cli.HAS_MATPLOTLIB:
            self.assertTrue(os.path.isfile(frequency_barh))
        else:
            self.assertFalse(os.path.isfile(frequency_barh))

        # Testing an invalid entry
        changed_filename = filename_template.format(chrom=1, suffix="maf")
        with open(changed_filename, "w") as o_file:
            print(*header, sep="\t", file=o_file)
            print("marker_1", "0.6", sep="\t", file=o_file)

        # This should raise an exception
        with self.assertRaises(GenipeError) as cm:
            cli.gather_maf_stats(
                required_chrom=autosomes,
                o_dir=self.output_dir.name,
            )
        self.assertEqual("{}: {}: invalid MAF".format("marker_1",
                                                      round(0.6, 3)),
                         str(cm.exception))

        # Testing an invalid entry
        changed_filename = filename_template.format(chrom=1, suffix="maf")
        with open(changed_filename, "w") as o_file:
            print(*header, sep="\t", file=o_file)
            print("marker_1", "-0.01", sep="\t", file=o_file)

        # This should raise an exception
        with self.assertRaises(GenipeError) as cm:
            cli.gather_maf_stats(
                required_chrom=autosomes,
                o_dir=self.output_dir.name,
            )
        self.assertEqual("{}: {}: invalid MAF".format("marker_1",
                                                      round(-0.01, 3)),
                         str(cm.exception))

        # Testing a good site with NA MAF
        changed_filename = filename_template.format(chrom=1, suffix="maf")
        with open(changed_filename, "w") as o_file:
            print(*header, sep="\t", file=o_file)
            print("marker_1", "NA", sep="\t", file=o_file)

        # This should issue a warning
        with self.assertLogs(level="WARNING") as cm:
            cli.gather_maf_stats(
                required_chrom=autosomes,
                o_dir=self.output_dir.name,
            )
        log_m = "WARNING:root:chr1: good sites with invalid MAF (NaN)"
        self.assertEqual(1, len(cm.output))
        self.assertEqual(log_m, cm.output[0])

        # Clearing the good sites file to see if we have a warning
        for chrom in autosomes:
            filename = filename_template.format(chrom=chrom,
                                                suffix="good_sites")
            with open(filename, "w") as o_file:
                pass

        # This should issue a warning
        with self.assertLogs(level="WARNING") as cm:
            cli.gather_maf_stats(
                required_chrom=autosomes,
                o_dir=self.output_dir.name,
            )
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
        with self.assertRaises(GenipeError) as cm:
            cli.gather_maf_stats(
                required_chrom=autosomes,
                o_dir=self.output_dir.name,
            )
        self.assertEqual("{}: no such file".format(removed_filename),
                         str(cm.exception))

        # Deleting a MAF sites file, and checking we have an error
        removed_filename = filename_template.format(chrom=1, suffix="maf")
        os.remove(removed_filename)
        self.assertFalse(os.path.isfile(removed_filename))

        # This should raise an exception
        with self.assertRaises(GenipeError) as cm:
            cli.gather_maf_stats(
                required_chrom=autosomes,
                o_dir=self.output_dir.name,
            )
        self.assertEqual("{}: no such file".format(removed_filename),
                         str(cm.exception))

    @unittest.skip("Test not implemented")
    def test_gather_execution_time(self):
        """Tests the 'gather_execution_time' function."""
        self.fail("Test not implemented")
