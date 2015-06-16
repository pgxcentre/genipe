
# This file is part of genipe.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import random
import platform
import unittest
from tempfile import TemporaryDirectory
from itertools import zip_longest as zip

import numpy as np
import pandas as pd
from pkg_resources import resource_filename

from ..tools.imputed_stats import *
from ..tools.imputed_stats import _get_result_from_linear_logistic_mixedlm

if HAS_STATSMODELS:
    # patsy is installed only if statsmodels is
    import patsy


__author__ = ["Louis-Philippe Lemieux Perreault", "Marc-Andre Legault"]
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["TestImputedStats", "TestImputedStatsCox", "TestImputedStatsLinear",
           "TestImputedStatsLogistic", "TestImputedStatsMixedLM",
           "TestImputedStatsSkat"]


def clean_logging_handlers():
    handlers = list(logging.root.handlers)
    for handler in handlers:
        logging.root.removeHandler(handler)


def reverse_dosage(dosage):
    """Finds values for d1, d2 and d3 from dosage."""
    d1, d2, d3 = 0, 0, 0
    if np.isnan(dosage):
        d1 = random.uniform(0.5, 0.8)
        d2 = random.uniform(0, 0.2)
        d3 = 1 - d1 - d2

    elif dosage == 1:
        d1 = 0
        d2 = dosage
        d3 = 0

    elif dosage == 2:
        d1 = 0
        d2 = 0
        d3 = 1

    elif dosage < 1:
        d2 = dosage
        d1 = 1 - d2
        d3 = 0

    elif dosage > 1.5:
        d1 = 0
        d3 = random.uniform(0.98, 1)
        d2 = ((dosage / 2) - d3) * 2

    else:
        d1 = 0
        d2 = random.uniform(0.98, 1)
        d3 = (dosage - d2) / 2

    return d1, d2, d3


def create_input_files(i_filename, output_dirname, analysis_type,
                       pheno_name="y", tte="y", event="y_d",
                       interaction=None, nb_process=None, sample_col=None,
                       use_ml=False):
    """Creates input files for the imputed_stats script."""
    # Reading the data
    data = pd.read_csv(i_filename, sep="\t", compression="bz2")

    # Adding a sample column
    if sample_col is None:
        data["sample_id"] = ["sample_{}".format(i+1) for i in range(len(data))]
        sample_col = "sample_id"
    else:
        assert sample_col in data.columns

    # Saving the phenotypes to file
    pheno_filename = os.path.join(output_dirname, "phenotypes.txt")
    data.to_csv(pheno_filename, sep="\t", index=False, na_rep="999999")

    # Creating the sample file (watch out for duplicated data)
    seen_samples = {}
    sample_filename = os.path.join(output_dirname, "samples.txt")
    with open(sample_filename, "w") as o_file:
        print("ID_1 ID_2 missing father mother sex plink_pheno", file=o_file)
        print("0 0 0 D D D B", file=o_file)
        for sample, gender in data[[sample_col, "gender"]].values:
            if sample in seen_samples:
                assert gender == seen_samples[sample]
                continue

            print(sample, sample, "0", "0", "0", gender, "-9", file=o_file)
            seen_samples[sample] = gender

    # Creating the IMPUTE2 file (watch out for duplicated data)
    impute2_filename = os.path.join(output_dirname, "impute2.txt")
    with open(impute2_filename, "w") as o_file:
        seen_samples = {}
        print("22 marker_1 1 T C", end="", file=o_file)
        for sample, dosage in data[[sample_col, "snp1"]].values:
            if sample in seen_samples:
                continue

            d1, d2, d3 = reverse_dosage(dosage)
            print("", d1, d2, d3, sep=" ", end="", file=o_file)
            seen_samples[sample] = dosage
        print(file=o_file)

        seen_samples = {}
        print("22 marker_2 2 G A", end="", file=o_file)
        for sample, dosage in data[[sample_col, "snp2"]].values:
            if sample in seen_samples:
                continue

            d1, d2, d3 = reverse_dosage(dosage)
            print("", d1, d2, d3, sep=" ", end="", file=o_file)
            seen_samples[sample] = dosage
        print(file=o_file)

        seen_samples = {}
        print("22 marker_3 3 AT A", end="", file=o_file)
        for sample, dosage in data[[sample_col, "snp3"]].values:
            if sample in seen_samples:
                continue

            d1, d2, d3 = reverse_dosage(dosage)
            print("", d1, d2, d3, sep=" ", end="", file=o_file)
            seen_samples[sample] = dosage
        print(file=o_file)

    # The prefix of the output files
    o_prefix = os.path.join(output_dirname, "test_imputed_stats")

    # The tool's options
    options = [
        analysis_type,
        "--impute2", impute2_filename,
        "--sample", sample_filename,
        "--pheno", pheno_filename,
        "--out", o_prefix,
        "--gender-column", "gender",
        "--covar", "C1,C2,C3,age,gender",
        "--missing-value", "999999",
        "--sample-column", sample_col,
    ]

    # Is this Cox or linear?
    if analysis_type == "cox":
        options.extend(["--time-to-event", tte, "--event", event])
    else:
        options.extend(["--pheno-name", pheno_name])

    # Is this mixedlm?
    if analysis_type == "mixedlm":
        if use_ml:
            options.append("--use-ml")

    # Is there interaction?
    if interaction is not None:
        options.extend(["--interaction", interaction])

    # Multi processes are required?
    if nb_process is not None:
        options.extend(["--nb-process", str(nb_process)])

    return o_prefix, options


class TestImputedStats(unittest.TestCase):

    def setUp(self):
        """Setup the tests."""
        # Creating the temporary directory
        self.output_dir = TemporaryDirectory(prefix="genipe_test_")

    def tearDown(self):
        """Finishes the test."""
        # Deleting the output directory
        self.output_dir.cleanup()

    def test_read_phenotype(self):
        """Tests the 'read_phenotype' function."""
        # A dummy object for options
        class Dummy(object):
            pass

        # The content of the phenotype file
        phenotype_content = (
            "sample_id\tTime_To_Event\tEvent\tC1\tC2\tC3\tC4\tGender\t"
            "Pheno_Lin\tPheno_Logit\tInter\n"
            "sample_1\t10\t0\t0.3\t0.2\t0.45\t0.01\t1\t0.01\t0\t0.00001\n"
            "sample_2\t2\t1\t0.9\t0.1\t0.42\t0.012\t2\t0.15\t1\t0.00332\n"
            "sample_3\t8\t1\t0.4\t0.67\t999\t0.001\t1\t0.0\t0\t0.000001\n"
        )
        filename = os.path.join(self.output_dir.name, "phenotypes.txt")
        with open(filename, "w") as o_file:
            o_file.write(phenotype_content)

        # Need an object for the function's option
        args = Dummy()
        args.missing_value = None
        args.sample_column = "sample_id"
        args.covar = ["C1", "C2", "C3", "Gender"]
        args.analysis_type = "cox"
        args.tte = "Time_To_Event"
        args.event = "Event"
        args.interaction = None
        args.chrx = False
        args.gender_column = "Gender"

        # The expected value
        expected_shape = (3, 6)
        expected_columns = {"C1", "C2", "C3", "Gender", "Time_To_Event",
                            "Event"}
        expected_index = ["sample_1", "sample_2", "sample_3"]
        expected_c1 = np.array([0.3, 0.9, 0.4], dtype=float)
        expected_c2 = np.array([0.2, 0.1, 0.67], dtype=float)
        expected_c3 = np.array([0.45, 0.42, 999], dtype=float)
        expected_gender = np.array([1, 2, 1], dtype=int)
        expected_tte = np.array([10, 2, 8], dtype=int)
        expected_event = np.array([0, 1, 1], dtype=int)
        expected_remove_g = False

        # The observed values
        observed_p, observed_remove_g = read_phenotype(filename, args)
        self.assertTrue(isinstance(observed_p, pd.DataFrame))
        self.assertEqual(expected_shape, observed_p.shape)
        self.assertEqual(expected_columns, set(observed_p.columns))
        self.assertEqual(expected_index, list(observed_p.index))
        self.assertTrue(np.allclose(expected_c1, observed_p.C1.values))
        self.assertTrue(np.allclose(expected_c2, observed_p.C2.values))
        self.assertTrue(np.allclose(expected_c3, observed_p.C3.values))
        self.assertTrue((expected_gender == observed_p.Gender.values).all())
        self.assertTrue(
            (expected_tte == observed_p.Time_To_Event.values).all()
        )
        self.assertTrue((expected_event == observed_p.Event.values).all())
        self.assertEqual(expected_remove_g, observed_remove_g)

        # Modifying the missing value to 999
        args.missing_value = "999"

        # The expected results
        expected_shape = (2, 6)

        # The observed values
        observed_p, observed_remove_g = read_phenotype(filename, args)
        self.assertTrue(isinstance(observed_p, pd.DataFrame))
        self.assertEqual(expected_shape, observed_p.shape)
        self.assertEqual(expected_columns, set(observed_p.columns))
        self.assertEqual(expected_index[:-1], list(observed_p.index))
        self.assertTrue(np.allclose(expected_c1[:-1], observed_p.C1.values))
        self.assertTrue(np.allclose(expected_c2[:-1], observed_p.C2.values))
        self.assertTrue(np.allclose(expected_c3[:-1], observed_p.C3.values))
        self.assertTrue(
            (expected_gender[:-1] == observed_p.Gender.values).all()
        )
        self.assertTrue(
            (expected_tte[:-1] == observed_p.Time_To_Event.values).all()
        )
        self.assertTrue(
            (expected_event[:-1] == observed_p.Event.values).all()
        )
        self.assertEqual(expected_remove_g, observed_remove_g)

        # Changing from Cox to linear should remove a column
        args.missing_value = None
        args.analysis_type = "linear"
        del args.tte
        del args.event
        args.pheno_name = "Pheno_Lin"

        # The expected results
        expected_shape = (3, 5)
        expected_columns = {"C1", "C2", "C3", "Gender", "Pheno_Lin"}
        expected_pheno = np.array([0.01, 0.15, 0], dtype=float)

        # The observed values
        observed_p, observed_remove_g = read_phenotype(filename, args)
        self.assertTrue(isinstance(observed_p, pd.DataFrame))
        self.assertEqual(expected_shape, observed_p.shape)
        self.assertEqual(expected_columns, set(observed_p.columns))
        self.assertEqual(expected_index, list(observed_p.index))
        self.assertTrue(np.allclose(expected_c1, observed_p.C1.values))
        self.assertTrue(np.allclose(expected_c2, observed_p.C2.values))
        self.assertTrue(np.allclose(expected_c3, observed_p.C3.values))
        self.assertTrue((expected_gender == observed_p.Gender.values).all())
        self.assertTrue((expected_pheno == observed_p.Pheno_Lin.values).all())
        self.assertEqual(expected_remove_g, observed_remove_g)

        # Changing from linear to logistic shouldn't change a thing
        args.analysis_type = "logistic"
        args.pheno_name = "Pheno_Logit"

        # The expected results
        expected_columns = {"C1", "C2", "C3", "Gender", "Pheno_Logit"}
        expected_pheno = np.array([0, 1, 0], dtype=int)

        # The observed values
        observed_p, observed_remove_g = read_phenotype(filename, args)
        self.assertTrue(isinstance(observed_p, pd.DataFrame))
        self.assertEqual(expected_shape, observed_p.shape)
        self.assertEqual(expected_columns, set(observed_p.columns))
        self.assertEqual(expected_index, list(observed_p.index))
        self.assertTrue(np.allclose(expected_c1, observed_p.C1.values))
        self.assertTrue(np.allclose(expected_c2, observed_p.C2.values))
        self.assertTrue(np.allclose(expected_c3, observed_p.C3.values))
        self.assertTrue((expected_gender == observed_p.Gender.values).all())
        self.assertTrue(
            (expected_pheno == observed_p.Pheno_Logit.values).all()
        )
        self.assertEqual(expected_remove_g, observed_remove_g)

        # Adding an interaction
        args.interaction = "Inter"

        # The expected results
        expected_shape = (3, 6)
        expected_columns = {"C1", "C2", "C3", "Gender", "Pheno_Logit", "Inter"}
        expected_inter = np.array([0.00001, 0.00332, 0.000001], dtype=float)

        # The observed values
        observed_p, observed_remove_g = read_phenotype(filename, args)
        self.assertTrue(isinstance(observed_p, pd.DataFrame))
        self.assertEqual(expected_shape, observed_p.shape)
        self.assertEqual(expected_columns, set(observed_p.columns))
        self.assertEqual(expected_index, list(observed_p.index))
        self.assertTrue(np.allclose(expected_c1, observed_p.C1.values))
        self.assertTrue(np.allclose(expected_c2, observed_p.C2.values))
        self.assertTrue(np.allclose(expected_c3, observed_p.C3.values))
        self.assertTrue((expected_gender == observed_p.Gender.values).all())
        self.assertTrue(
            (expected_pheno == observed_p.Pheno_Logit.values).all()
        )
        self.assertTrue(np.allclose(expected_inter, observed_p.Inter.values))
        self.assertEqual(expected_remove_g, observed_remove_g)

        # Removing the gender in the covars, but setting chrx to true
        args.covar = ["C1", "C2", "C3"]
        args.chrx = True

        # The expected results
        expected_remove_g = True

        # The observed values
        observed_p, observed_remove_g = read_phenotype(filename, args)
        self.assertTrue(isinstance(observed_p, pd.DataFrame))
        self.assertEqual(expected_shape, observed_p.shape)
        self.assertEqual(expected_columns, set(observed_p.columns))
        self.assertEqual(expected_index, list(observed_p.index))
        self.assertTrue(np.allclose(expected_c1, observed_p.C1.values))
        self.assertTrue(np.allclose(expected_c2, observed_p.C2.values))
        self.assertTrue(np.allclose(expected_c3, observed_p.C3.values))
        self.assertTrue((expected_gender == observed_p.Gender.values).all())
        self.assertTrue(
            (expected_pheno == observed_p.Pheno_Logit.values).all()
        )
        self.assertTrue(np.allclose(expected_inter, observed_p.Inter.values))
        self.assertEqual(expected_remove_g, observed_remove_g)

        # Adding a sample (with unknown gender) (should be included since not
        # in covar, even though we ask for chrX)
        with open(filename, "a") as o_file:
            o_file.write(
                "sample_4\t8\t1\t0.4\t0.67\t999\t0.001\t0\t0.0\t0\t0.000001\n"
            )

        # The expected values
        expected_shape = (4, 6)
        expected_index.append("sample_4")
        expected_c1 = np.array([0.3, 0.9, 0.4, 0.4], dtype=float)
        expected_c2 = np.array([0.2, 0.1, 0.67, 0.67], dtype=float)
        expected_c3 = np.array([0.45, 0.42, 999, 999], dtype=float)
        expected_gender = np.array([1, 2, 1, 0], dtype=int)
        expected_pheno = np.array([0, 1, 0, 0], dtype=int)
        expected_inter = np.array([0.00001, 0.00332, 0.000001, 0.000001],
                                  dtype=float)

        # The observed values
        observed_p, observed_remove_g = read_phenotype(filename, args)
        self.assertTrue(isinstance(observed_p, pd.DataFrame))
        self.assertEqual(expected_shape, observed_p.shape)
        self.assertEqual(expected_columns, set(observed_p.columns))
        self.assertEqual(expected_index, list(observed_p.index))
        self.assertTrue(np.allclose(expected_c1, observed_p.C1.values))
        self.assertTrue(np.allclose(expected_c2, observed_p.C2.values))
        self.assertTrue(np.allclose(expected_c3, observed_p.C3.values))
        self.assertTrue((expected_gender == observed_p.Gender.values).all())
        self.assertTrue(
            (expected_pheno == observed_p.Pheno_Logit.values).all()
        )
        self.assertTrue(np.allclose(expected_inter, observed_p.Inter.values))
        self.assertEqual(expected_remove_g, observed_remove_g)

        # Sample shouldn't be included if chrx is False, but Gender in covars
        args.covar.append("Gender")
        args.chrx = False

        # The expected values
        expected_shape = (3, 6)
        expected_remove_g = False

        # The observed values
        observed_p, observed_remove_g = read_phenotype(filename, args)
        self.assertTrue(isinstance(observed_p, pd.DataFrame))
        self.assertEqual(expected_shape, observed_p.shape)
        self.assertEqual(expected_columns, set(observed_p.columns))
        self.assertEqual(expected_index[:-1], list(observed_p.index))
        self.assertTrue(np.allclose(expected_c1[:-1], observed_p.C1.values))
        self.assertTrue(np.allclose(expected_c2[:-1], observed_p.C2.values))
        self.assertTrue(np.allclose(expected_c3[:-1], observed_p.C3.values))
        self.assertTrue(
            (expected_gender[:-1] == observed_p.Gender.values).all()
        )
        self.assertTrue(
            (expected_pheno[:-1] == observed_p.Pheno_Logit.values).all()
        )
        self.assertTrue(
            np.allclose(expected_inter[:-1], observed_p.Inter.values)
        )
        self.assertEqual(expected_remove_g, observed_remove_g)

        # Removing gender in the covar should add a sample
        args.covar = args.covar[:-1]

        # The expected values
        expected_shape = (4, 5)
        expected_columns.remove("Gender")

        # The observed values
        observed_p, observed_remove_g = read_phenotype(filename, args)
        self.assertTrue(isinstance(observed_p, pd.DataFrame))
        self.assertEqual(expected_shape, observed_p.shape)
        self.assertEqual(expected_columns, set(observed_p.columns))
        self.assertEqual(expected_index, list(observed_p.index))
        self.assertTrue(np.allclose(expected_c1, observed_p.C1.values))
        self.assertTrue(np.allclose(expected_c2, observed_p.C2.values))
        self.assertTrue(np.allclose(expected_c3, observed_p.C3.values))
        self.assertTrue(
            (expected_pheno == observed_p.Pheno_Logit.values).all()
        )
        self.assertTrue(np.allclose(expected_inter, observed_p.Inter.values))
        self.assertEqual(expected_remove_g, observed_remove_g)

        # Setting chromosome to X should also include the sample
        args.chrx = True

        # The expected values
        expected_shape = (4, 6)
        expected_columns.add("Gender")
        expected_remove_g = True

        # The observed values
        observed_p, observed_remove_g = read_phenotype(filename, args)
        self.assertTrue(isinstance(observed_p, pd.DataFrame))
        self.assertEqual(expected_shape, observed_p.shape)
        self.assertEqual(expected_columns, set(observed_p.columns))
        self.assertEqual(expected_index, list(observed_p.index))
        self.assertTrue(np.allclose(expected_c1, observed_p.C1.values))
        self.assertTrue(np.allclose(expected_c2, observed_p.C2.values))
        self.assertTrue(np.allclose(expected_c3, observed_p.C3.values))
        self.assertTrue(
            (expected_pheno == observed_p.Pheno_Logit.values).all()
        )
        self.assertTrue(np.allclose(expected_inter, observed_p.Inter.values))
        self.assertEqual(expected_remove_g, observed_remove_g)

    def test_read_samples(self):
        """Tests the 'read_samples' function."""
        # Creating the sample file
        sample_content = (
            "ID_1 ID_2 missing father mother sex plink_pheno\n"
            "0 0 0 D D D B\n"
            "fam_1 sample_1 0 0 0 2 -9\n"
            "fam_1 sample_2 0 0 0 1 -9\n"
            "fam_2 sample_3 0 0 0 2 -9\n"
        )
        sample_filename = os.path.join(self.output_dir.name, "test.sample")
        with open(sample_filename, "w") as o_file:
            o_file.write(sample_content)

        # The expected values
        expected_columns = ["ID_1"]
        expected_index = ["sample_1", "sample_2", "sample_3"]
        expected_fam = ["fam_1", "fam_1", "fam_2"]

        # The observed values
        observed = read_samples(sample_filename)

        # Checking
        self.assertTrue(isinstance(observed, pd.DataFrame))
        self.assertEqual((3, 1), observed.shape)
        self.assertEqual(expected_columns, observed.columns)
        self.assertEqual(expected_index, list(observed.index.values))
        self.assertEqual(expected_fam, list(observed.ID_1.values))

        # Having a duplicated samples will trigger an exception
        sample_content = (
            "ID_1 ID_2 missing father mother sex plink_pheno\n"
            "0 0 0 D D D B\n"
            "fam_1 sample_1 0 0 0 2 -9\n"
            "fam_1 sample_2 0 0 0 1 -9\n"
            "fam_2 sample_2 0 0 0 2 -9\n"
        )
        with open(sample_filename, "w") as o_file:
            o_file.write(sample_content)

        # Checking
        with self.assertRaises(ValueError) as cm:
            read_samples(sample_filename)
        self.assertEqual("Index has duplicate keys: ['sample_2']",
                         str(cm.exception))

    def test_read_sites_to_extract(self):
        """Tests the 'test_read_sites_to_extract' function."""
        file_content = ["marker_{}".format(i) for i in range(100)] * 2
        filename = os.path.join(self.output_dir.name, "markers.txt")
        with open(filename, "w") as o_file:
            o_file.write("\n".join(file_content) + "\n")

        # The expected values
        expected = {"marker_{}".format(i) for i in range(100)}

        # The observed values
        observed = read_sites_to_extract(filename)

        # Checking
        self.assertEqual(expected, observed)

    @unittest.skip("Test not implemented")
    def test_compute_statistics(self):
        """Tests the 'compute_statistics' function."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_process_impute2_site(self):
        """Tests the 'process_impute2_site' function."""
        self.fail("Test not implemented")

    def test_samples_with_hetero_calls(self):
        """Tests the 'samples_with_hetero_calls' function."""
        data = [
            ("sample_1", 1.0, 0.0, 0.0),
            ("sample_2", 0.0, 1.0, 0.0),
            ("sample_3", 0.0, 0.0, 1.0),
            ("sample_4", 0.9, 0.1, 0.0),
            ("sample_5", 0.1, 0.8, 0.1),
            ("sample_6", 0.0, 0.4, 0.6),
            ("sample_7", 0.2, 0.5, 0.3),
            ("sample_8", 0.9, 0.05, 0.05),
            ("sample_9", 0.0, 1.0, 0.0),
        ]
        data = pd.DataFrame(data, columns=["sample_id", "D1", "D2", "D3"])
        data = data.set_index("sample_id", verify_integrity=True)

        # The expected results
        expected = ["sample_2", "sample_5", "sample_7", "sample_9"]

        # The observed results
        observed = samples_with_hetero_calls(data, "D2")

        # Checking
        self.assertTrue(isinstance(observed, pd.Index))
        self.assertEqual(expected, list(observed))

    def test_get_formula(self):
        """Tests the 'get_formula' function."""
        # Testing with only one phenotype (no covars, no interaction)
        expected = "pheno ~ _GenoD"
        observed = get_formula("pheno", [], None, "None", set())
        self.assertEqual(expected, observed)

        # Testing with one covar, no interaction
        expected = "pheno ~ _GenoD + C1"
        observed = get_formula("pheno", ["C1"], None, "None", set())
        self.assertEqual(expected, observed)

        # Testing with more than one covar, no interaction
        expected = "pheno ~ _GenoD + C1 + C2 + C3"
        observed = get_formula("pheno", ["C1", "C2", "C3"], None, "None",
                               set())
        self.assertEqual(expected, observed)

        # Testing with more than one covar (with gender), no interaction
        expected = "pheno ~ _GenoD + C1 + C2 + C3 + C(Gender)"
        observed = get_formula("pheno", ["C1", "C2", "C3", "Gender"], None,
                               "Gender", set())
        self.assertEqual(expected, observed)

        # Testing with without covar, but with interaction
        expected = "pheno ~ _GenoD + _GenoD*inter"
        observed = get_formula("pheno", [], "inter", "None", set())
        self.assertEqual(expected, observed)

        # Testing with one covar and interaction
        expected = "pheno ~ _GenoD + inter + _GenoD*inter"
        observed = get_formula("pheno", ["inter"], "inter", "None", set())
        self.assertEqual(expected, observed)

        # Testing with more than one covar and interaction
        expected = "pheno ~ _GenoD + C1 + C2 + C3 + inter + _GenoD*inter"
        observed = get_formula("pheno", ["C1", "C2", "C3", "inter"], "inter",
                               "None", set())
        self.assertEqual(expected, observed)

        # Testing with more than one covar and interaction and Gender
        expected = ("pheno ~ _GenoD + C1 + C2 + C3 + C(Gender) + "
                    "_GenoD*C(Gender)")
        observed = get_formula("pheno", ["C1", "C2", "C3", "Gender"], "Gender",
                               "Gender", set())
        self.assertEqual(expected, observed)

        # Testing with categorical values
        expected = ("pheno ~ _GenoD + C1 + C2 + C3 + C(Treatment) + "
                    "C(Population) + C(Gender) + _GenoD*C(Treatment)")
        observed = get_formula(
            "pheno",
            ["C1", "C2", "C3", "Treatment", "Population", "Gender"],
            "Treatment",
            "Gender",
            {"Treatment", "Population"},
        )
        self.assertEqual(expected, observed)

    @unittest.skip("Test not implemented")
    def test_get_result_from_linear_logistic(self):
        """Tests the '_get_result_from_linear_logistic' function."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_check_args(self):
        """Tests the 'check_args' function."""
        self.fail("Test not implemented")


@unittest.skipIf(not HAS_LIFELINES,
                 "optional requirement (lifelines) not satisfied")
class TestImputedStatsCox(unittest.TestCase):

    def setUp(self):
        """Setup the tests."""
        # Creating the temporary directory
        self.output_dir = TemporaryDirectory(prefix="genipe_test_")

    def tearDown(self):
        """Finishes the test."""
        # Deleting the output directory
        self.output_dir.cleanup()

    def test_fit_cox(self):
        """Tests the 'fit_cox' function."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_sim.txt.bz2",
        )

        # This dataset contains 3 markers + 5 covariables
        data = pd.read_csv(data_filename, sep="\t", compression="bz2")
        columns_to_keep = ["y", "y_d", "snp1", "C1", "C2", "C3", "age",
                           "gender"]

        # The expected results for the first marker (according to R)
        # P value was given by 2*pnorm(z) according to this site
        # https://stat.ethz.ch/pipermail/r-help/2008-October/178439.html
        expected_coef = -0.9892699611323815256
        expected_se = 0.0645633075569013448
        expected_min_ci = -1.11581171866669093
        expected_max_ci = -0.86272820359807223
        expected_z = -15.32247957186071652
        expected_p = 5.41131727236088009e-53

        # The observed results for the first marker
        observed = fit_cox(
            data=data[columns_to_keep].dropna(axis=0),
            time_to_event="y",
            event="y_d",
            result_col="snp1",
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_z, observed_p = observed

        # Comparing the results (P is not that similar though...)
        self.assertAlmostEqual(expected_coef, observed_coef, places=10)
        self.assertAlmostEqual(expected_se, observed_se, places=6)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=4)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=4)
        self.assertAlmostEqual(expected_z, observed_z, places=4)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=3)

        # The second marker
        columns_to_keep = ["y", "y_d", "snp2", "C1", "C2", "C3", "age",
                           "gender"]

        # The expected results for the first marker (according to R)
        expected_coef = 0.0499084338021939383
        expected_se = 0.0401904025600694215
        expected_min_ci = -0.0288633077397084936
        expected_max_ci = 0.12868017534409637
        expected_z = 1.241799798536472599
        expected_p = 0.21431043709814412423

        # The observed results for the first marker
        observed = fit_cox(
            data=data[columns_to_keep].dropna(axis=0),
            time_to_event="y",
            event="y_d",
            result_col="snp2",
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_z, observed_p = observed

        # Comparing the results (P is not that similar though...)
        self.assertAlmostEqual(expected_coef, observed_coef, places=10)
        self.assertAlmostEqual(expected_se, observed_se, places=9)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=4)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=4)
        self.assertAlmostEqual(expected_z, observed_z, places=8)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=8)

        # The third marker
        columns_to_keep = ["y", "y_d", "snp3", "C1", "C2", "C3", "age",
                           "gender"]

        # The expected results for the first marker (according to R)
        expected_coef = 1.114732193640146418
        expected_se = 0.0631205046371715733
        expected_min_ci = 0.991018277865296726
        expected_max_ci = 1.23844610941499611
        expected_z = 17.660381520202268035
        expected_p = 8.46654634146707403e-70

        # The observed results for the first marker
        observed = fit_cox(
            data=data[columns_to_keep].dropna(axis=0),
            time_to_event="y",
            event="y_d",
            result_col="snp3",
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_z, observed_p = observed

        # Comparing the results (P is not that similar though...)
        self.assertAlmostEqual(expected_coef, observed_coef, places=10)
        self.assertAlmostEqual(expected_se, observed_se, places=7)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=4)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=4)
        self.assertAlmostEqual(expected_z, observed_z, places=5)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=4)

    def test_fit_cox_interaction(self):
        """Tests the 'fit_cox' function with interaction."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_sim_inter.txt.bz2",
        )

        # This dataset contains 3 markers + 5 covariables
        data = pd.read_csv(data_filename, sep="\t", compression="bz2")

        # Computing the interaction
        data["interaction"] = data.snp1 * data.gender
        columns_to_keep = ["y", "y_d", "snp1", "C1", "C2", "C3", "age",
                           "gender", "interaction"]

        # The expected results for the first marker (according to R)
        expected_coef = -0.6446297920053343233
        expected_se = 0.1430770156758568168
        expected_min_ci = -0.925055589745486406
        expected_max_ci = -0.364203994265182296
        expected_z = -4.5054741249688419202
        expected_p = 6.62249088800859198e-06

        # The observed results
        observed = fit_cox(
            data=data[columns_to_keep].dropna(axis=0),
            time_to_event="y",
            event="y_d",
            result_col="interaction",
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_z, observed_p = observed

        # Comparing the results
        self.assertAlmostEqual(expected_coef, observed_coef, places=10)
        self.assertAlmostEqual(expected_se, observed_se, places=7)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=3)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=3)
        self.assertAlmostEqual(expected_z, observed_z, places=5)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=5)

    def test_full_fit_cox(self):
        """Tests the full pipeline for Cox's regression."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_sim.txt.bz2",
        )

        # Creating the input files
        o_prefix, options = create_input_files(
            i_filename=data_filename,
            output_dirname=self.output_dir.name,
            analysis_type="cox",
            tte="y",
            event="y_d",
        )

        # Executing the tool
        try:
            main(args=options)
        finally:
            clean_logging_handlers()

        # Making sure the output file exists
        self.assertTrue(os.path.isfile(o_prefix + ".cox.dosage"))

        # Reading the data
        observed = pd.read_csv(o_prefix + ".cox.dosage", sep="\t")

        # Checking the shape
        self.assertEqual((3, 13), observed.shape)

        # Checking all columns are present
        self.assertEqual(["chr", "pos", "snp", "major", "minor", "maf", "n",
                          "coef", "se", "lower", "upper", "z", "p"],
                         list(observed.columns))

        # Chromosomes
        self.assertEqual([22], observed.chr.unique())

        # Positions
        self.assertEqual([1, 2, 3], list(observed.pos))

        # Marker names
        self.assertEqual(["marker_1", "marker_2", "marker_3"],
                         list(observed.snp))

        # Major alleles
        self.assertEqual(["T", "G", "AT"], list(observed.major))

        # Minor alleles
        self.assertEqual(["C", "A", "A"], list(observed.minor))

        # Minor allele frequency
        expected = [1724 / 11526, 4604 / 11526, 1379 / 11526]
        for expected_maf, observed_maf in zip(expected, observed.maf):
            self.assertAlmostEqual(expected_maf, observed_maf, places=10)

        # The number of samples
        expected = [5763, 5763, 5763]
        for expected_n, observed_n in zip(expected, observed.n):
            self.assertEqual(expected_n, observed_n)

        # The coefficients
        expected = [-0.9892699611323815256, 0.0499084338021939383,
                    1.114732193640146418]
        for expected_coef, observed_coef in zip(expected, observed.coef):
            self.assertAlmostEqual(expected_coef, observed_coef, places=10)

        # The standard error
        expected = [0.0645633075569013448, 0.0401904025600694215,
                    0.0631205046371715733]
        places = [6, 9, 7]
        zipped = zip(expected, observed.se, places)
        for expected_se, observed_se, place in zipped:
            self.assertAlmostEqual(expected_se, observed_se, places=place)

        # The lower CI
        expected = [-1.11581171866669093, -0.0288633077397084936,
                    0.991018277865296726]
        for expected_min_ci, observed_min_ci in zip(expected, observed.lower):
            self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=4)

        # The upper CI
        expected = [-0.86272820359807223, 0.12868017534409637,
                    1.23844610941499611]
        for expected_max_ci, observed_max_ci in zip(expected, observed.upper):
            self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=4)

        # The Z statistics
        expected = [-15.32247957186071652, 1.241799798536472599,
                    17.660381520202268035]
        places = [4, 8, 5]
        zipped = zip(expected, observed.z, places)
        for expected_z, observed_z, place in zipped:
            self.assertAlmostEqual(expected_z, observed_z, places=place)

        # The p values
        expected = [5.41131727236088009e-53, 0.21431043709814412423,
                    8.46654634146707403e-70]
        places = [3, 8, 4]
        zipped = zip(expected, observed.p, places)
        for expected_p, observed_p, place in zipped:
            self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                                   places=place)

    @unittest.skipIf(platform.system() == "Darwin",
                     "multiprocessing not supported with Mac OS")
    def test_full_fit_cox_multiprocess(self):
        """Tests the full pipeline for Cox's regression with >1 processes."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_sim.txt.bz2",
        )

        # Creating the input files
        o_prefix, options = create_input_files(
            i_filename=data_filename,
            output_dirname=self.output_dir.name,
            analysis_type="cox",
            tte="y",
            event="y_d",
            nb_process=2,
        )

        # Executing the tool
        try:
            main(args=options)
        finally:
            clean_logging_handlers()

        # Making sure the output file exists
        self.assertTrue(os.path.isfile(o_prefix + ".cox.dosage"))

        # Reading the data
        observed = pd.read_csv(o_prefix + ".cox.dosage", sep="\t")

        # Checking the shape
        self.assertEqual((3, 13), observed.shape)

        # Checking all columns are present
        self.assertEqual(["chr", "pos", "snp", "major", "minor", "maf", "n",
                          "coef", "se", "lower", "upper", "z", "p"],
                         list(observed.columns))

        # Chromosomes
        self.assertEqual([22], observed.chr.unique())

        # Positions
        self.assertEqual([1, 2, 3], list(observed.pos))

        # Marker names
        self.assertEqual(["marker_1", "marker_2", "marker_3"],
                         list(observed.snp))

        # Major alleles
        self.assertEqual(["T", "G", "AT"], list(observed.major))

        # Minor alleles
        self.assertEqual(["C", "A", "A"], list(observed.minor))

        # Minor allele frequency
        expected = [1724 / 11526, 4604 / 11526, 1379 / 11526]
        for expected_maf, observed_maf in zip(expected, observed.maf):
            self.assertAlmostEqual(expected_maf, observed_maf, places=10)

        # The number of samples
        expected = [5763, 5763, 5763]
        for expected_n, observed_n in zip(expected, observed.n):
            self.assertEqual(expected_n, observed_n)

        # The coefficients
        expected = [-0.9892699611323815256, 0.0499084338021939383,
                    1.114732193640146418]
        for expected_coef, observed_coef in zip(expected, observed.coef):
            self.assertAlmostEqual(expected_coef, observed_coef, places=10)

        # The standard error
        expected = [0.0645633075569013448, 0.0401904025600694215,
                    0.0631205046371715733]
        places = [6, 9, 7]
        zipped = zip(expected, observed.se, places)
        for expected_se, observed_se, place in zipped:
            self.assertAlmostEqual(expected_se, observed_se, places=place)

        # The lower CI
        expected = [-1.11581171866669093, -0.0288633077397084936,
                    0.991018277865296726]
        for expected_min_ci, observed_min_ci in zip(expected, observed.lower):
            self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=4)

        # The upper CI
        expected = [-0.86272820359807223, 0.12868017534409637,
                    1.23844610941499611]
        for expected_max_ci, observed_max_ci in zip(expected, observed.upper):
            self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=4)

        # The Z statistics
        expected = [-15.32247957186071652, 1.241799798536472599,
                    17.660381520202268035]
        places = [4, 8, 5]
        zipped = zip(expected, observed.z, places)
        for expected_z, observed_z, place in zipped:
            self.assertAlmostEqual(expected_z, observed_z, places=place)

        # The p values
        expected = [5.41131727236088009e-53, 0.21431043709814412423,
                    8.46654634146707403e-70]
        places = [3, 8, 4]
        zipped = zip(expected, observed.p, places)
        for expected_p, observed_p, place in zipped:
            self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                                   places=place)

    def test_full_fit_cox_interaction(self):
        """Tests the full pipeline for Cox's regression with interaction."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_sim_inter.txt.bz2",
        )

        # Creating the input files
        o_prefix, options = create_input_files(
            i_filename=data_filename,
            output_dirname=self.output_dir.name,
            analysis_type="cox",
            tte="y",
            event="y_d",
            interaction="gender",
        )

        # Executing the tool
        try:
            main(args=options)
        finally:
            clean_logging_handlers()

        # Making sure the output file exists
        self.assertTrue(os.path.isfile(o_prefix + ".cox.dosage"))

        # Reading the data
        observed = pd.read_csv(o_prefix + ".cox.dosage", sep="\t")

        # Checking the shape
        self.assertEqual((3, 13), observed.shape)

        # Checking all columns are present
        self.assertEqual(["chr", "pos", "snp", "major", "minor", "maf", "n",
                          "coef", "se", "lower", "upper", "z", "p"],
                         list(observed.columns))

        # Chromosomes
        self.assertEqual([22], observed.chr.unique())

        # Positions
        self.assertEqual([1, 2, 3], list(observed.pos))

        # Marker names
        self.assertEqual(["marker_1", "marker_2", "marker_3"],
                         list(observed.snp))

        # Major alleles
        self.assertEqual(["T", "G", "AT"], list(observed.major))

        # Minor alleles
        self.assertEqual(["C", "A", "A"], list(observed.minor))

        # Minor allele frequency
        expected = [1724 / 11526, 4604 / 11526, 1379 / 11526]
        for expected_maf, observed_maf in zip(expected, observed.maf):
            self.assertAlmostEqual(expected_maf, observed_maf, places=10)

        # The number of samples
        expected = [5763, 5763, 5763]
        for expected_n, observed_n in zip(expected, observed.n):
            self.assertEqual(expected_n, observed_n)

        # The coefficients
        expected = [-0.6446297920053343233, -0.00626679211298179859,
                    -0.468685693217335109]
        places = [10, 8, 8]
        zipped = zip(expected, observed.coef, places)
        for expected_coef, observed_coef, place in zipped:
            self.assertAlmostEqual(expected_coef, observed_coef, places=place)

        # The standard error
        expected = [0.1430770156758568168, 0.0834709381328159472,
                    0.1209229684191305970]
        places = [7, 10, 10]
        zipped = zip(expected, observed.se, places)
        for expected_se, observed_se, place in zipped:
            self.assertAlmostEqual(expected_se, observed_se, places=place)

        # The lower CI
        expected = [-0.925055589745486406, -0.16986682460907207,
                    -0.705690356222505422]
        places = [3, 4, 3]
        zipped = zip(expected, observed.lower, places)
        for expected_min_ci, observed_min_ci, place in zipped:
            self.assertAlmostEqual(expected_min_ci, observed_min_ci,
                                   places=place)

        # The upper CI
        expected = [-0.364203994265182296, 0.157333240383108447,
                    -0.231681030212164824]
        places = [3, 4, 3]
        zipped = zip(expected, observed.upper, places)
        for expected_max_ci, observed_max_ci, place in zipped:
            self.assertAlmostEqual(expected_max_ci, observed_max_ci,
                                   places=place)

        # The Z statistics
        expected = [-4.5054741249688419202, -0.0750775330092768867,
                    -3.875902976453783122]
        places = [5, 7, 7]
        zipped = zip(expected, observed.z, places)
        for expected_z, observed_z, place in zipped:
            self.assertAlmostEqual(expected_z, observed_z, places=place)

        # The p values
        expected = [6.62249088800859198e-06, 9.40153023426109735e-01,
                    1.06230010246344264e-04]
        places = [5, 8, 7]
        zipped = zip(expected, observed.p, places)
        for expected_p, observed_p, place in zipped:
            self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                                   places=place)


@unittest.skipIf(not HAS_STATSMODELS,
                 "optional requirement (statsmodels) not satisfied")
class TestImputedStatsLinear(unittest.TestCase):

    def setUp(self):
        """Setup the tests."""
        # Creating the temporary directory
        self.output_dir = TemporaryDirectory(prefix="genipe_test_")

    def tearDown(self):
        """Finishes the test."""
        # Deleting the output directory
        self.output_dir.cleanup()

    def test_fit_linear(self):
        """Tests the 'fit_linear' function."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_sim.txt.bz2",
        )

        # This dataset contains 3 markers + 5 covariables
        data = pd.read_csv(data_filename, sep="\t", compression="bz2")

        # The formula for the first marker
        formula = "y ~ snp1 + C1 + C2 + C3 + age + C(gender)"
        columns_to_keep = ["y", "snp1", "C1", "C2", "C3", "age", "gender"]

        # The expected results for the first marker (according to R)
        # The data was simulated so that snp1 had a coefficient of 0.1
        expected_coef = 0.09930262321654575
        expected_se = 0.00302135517743109
        expected_min_ci = 0.09337963040899197
        expected_max_ci = 0.10522561602409949
        expected_t = 32.866914806414108
        expected_p = 2.7965174627917724e-217

        # The observed results for the first marker
        observed = fit_linear(
            data=data[columns_to_keep].dropna(axis=0),
            formula=formula,
            result_col="snp1",
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_t, observed_p = observed

        # Comparing the results
        self.assertAlmostEqual(expected_coef, observed_coef, places=10)
        self.assertAlmostEqual(expected_se, observed_se, places=10)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=10)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=10)
        self.assertAlmostEqual(expected_t, observed_t, places=10)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=10)

        # The formula for the second marker
        formula = "y ~ snp2 + C1 + C2 + C3 + age + C(gender)"
        columns_to_keep = ["y", "snp2", "C1", "C2", "C3", "age", "gender"]

        # The expected results for the second marker (according to R)
        # The data was simulated so that snp1 had a coefficient of 0
        expected_coef = -0.00279702443754753
        expected_se = 0.00240385609310785
        expected_min_ci = -0.0075094867313642991
        expected_max_ci = 0.0019154378562692411
        expected_t = -1.1635573550209353
        expected_p = 0.24465167231462448

        # The observed results for the first marker
        observed = fit_linear(
            data=data[columns_to_keep].dropna(axis=0),
            formula=formula,
            result_col="snp2",
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_t, observed_p = observed

        # Comparing the results
        self.assertAlmostEqual(expected_coef, observed_coef, places=10)
        self.assertAlmostEqual(expected_se, observed_se, places=10)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=10)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=10)
        self.assertAlmostEqual(expected_t, observed_t, places=10)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=10)

        # The formula for the third (and last) marker
        formula = "y ~ snp3 + C1 + C2 + C3 + age + C(gender)"
        columns_to_keep = ["y", "snp3", "C1", "C2", "C3", "age", "gender"]

        # The expected results for the second marker (according to R)
        # The data was simulated so that snp1 had a coefficient of -0.12
        expected_coef = -0.11731595824657762
        expected_se = 0.00327175651867383
        expected_min_ci = -0.12372983188610413
        expected_max_ci = -0.1109020846070511
        expected_t = -35.857178728608552
        expected_p = 2.4882495142044017e-254

        # The observed results for the first marker
        observed = fit_linear(
            data=data[columns_to_keep].dropna(axis=0),
            formula=formula,
            result_col="snp3",
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_t, observed_p = observed

        # Comparing the results
        self.assertAlmostEqual(expected_coef, observed_coef, places=10)
        self.assertAlmostEqual(expected_se, observed_se, places=10)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=10)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=10)
        self.assertAlmostEqual(expected_t, observed_t, places=10)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=10)

        # Asking for an invalid column should raise a KeyError
        with self.assertRaises(KeyError) as cm:
            fit_linear(
                data=data[columns_to_keep].dropna(axis=0),
                formula=formula,
                result_col="unknown",
            )

        with self.assertRaises(patsy.PatsyError) as cm:
            fit_linear(
                data=data[columns_to_keep].dropna(axis=0),
                formula=formula + " + unknown",
                result_col="snp4",
            )

    def test_fit_linear_interaction(self):
        """Tests the 'fit_linear' function with interaction."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_sim_inter.txt.bz2",
        )

        # This dataset contains 3 markers + 5 covariables
        data = pd.read_csv(data_filename, sep="\t", compression="bz2")

        # The formula for the first marker
        formula = "y ~ snp1 + C1 + C2 + C3 + age + C(gender) + snp1*C(gender)"
        columns_to_keep = ["y", "snp1", "C1", "C2", "C3", "age", "gender"]

        # The expected results for the first marker (according to R)
        # The data was simulated so that snp1 had a coefficient of 0.1
        expected_coef = 0.101959570845217173
        expected_se = 0.00588764130382715949
        expected_min_ci = 0.090417578486636035
        expected_max_ci = 0.11350156320379831
        expected_t = 17.3175581839437314
        expected_p = 1.5527088106478929e-65

        # The observed results
        observed = fit_linear(
            data=data[columns_to_keep].dropna(axis=0),
            formula=formula,
            result_col="snp1:C(gender)[T.2]",
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_t, observed_p = observed

        # Comparing the results
        self.assertAlmostEqual(expected_coef, observed_coef, places=10)
        self.assertAlmostEqual(expected_se, observed_se, places=10)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=10)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=10)
        self.assertAlmostEqual(expected_t, observed_t, places=10)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=10)

    def test_full_fit_linear(self):
        """Tests the full pipeline for linear regression."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_sim.txt.bz2",
        )

        # Creating the input files
        o_prefix, options = create_input_files(
            i_filename=data_filename,
            output_dirname=self.output_dir.name,
            analysis_type="linear",
        )

        # Executing the tool
        try:
            main(args=options)
        finally:
            clean_logging_handlers()

        # Making sure the output file exists
        self.assertTrue(os.path.isfile(o_prefix + ".linear.dosage"))

        # Reading the data
        observed = pd.read_csv(o_prefix + ".linear.dosage", sep="\t")

        # Checking the shape
        self.assertEqual((3, 13), observed.shape)

        # Checking all columns are present
        self.assertEqual(["chr", "pos", "snp", "major", "minor", "maf", "n",
                          "coef", "se", "lower", "upper", "t", "p"],
                         list(observed.columns))

        # Chromosomes
        self.assertEqual([22], observed.chr.unique())

        # Positions
        self.assertEqual([1, 2, 3], list(observed.pos))

        # Marker names
        self.assertEqual(["marker_1", "marker_2", "marker_3"],
                         list(observed.snp))

        # Major alleles
        self.assertEqual(["T", "G", "AT"], list(observed.major))

        # Minor alleles
        self.assertEqual(["C", "A", "A"], list(observed.minor))

        # Minor allele frequency
        expected = [1724 / 11526, 4604 / 11526, 1379 / 11526]
        for expected_maf, observed_maf in zip(expected, observed.maf):
            self.assertAlmostEqual(expected_maf, observed_maf, places=10)

        # The number of samples
        expected = [5763, 5763, 5763]
        for expected_n, observed_n in zip(expected, observed.n):
            self.assertEqual(expected_n, observed_n)

        # The coefficients
        expected = [0.09930262321654575, -0.00279702443754753,
                    -0.11731595824657762]
        for expected_coef, observed_coef in zip(expected, observed.coef):
            self.assertAlmostEqual(expected_coef, observed_coef, places=10)

        # The standard error
        expected = [0.00302135517743109, 0.00240385609310785,
                    0.00327175651867383]
        for expected_se, observed_se in zip(expected, observed.se):
            self.assertAlmostEqual(expected_se, observed_se, places=10)

        # The lower CI
        expected = [0.09337963040899197, -0.0075094867313642991,
                    -0.12372983188610413]
        for expected_min_ci, observed_min_ci in zip(expected, observed.lower):
            self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=10)

        # The upper CI
        expected = [0.10522561602409949, 0.0019154378562692411,
                    -0.1109020846070511]
        for expected_max_ci, observed_max_ci in zip(expected, observed.upper):
            self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=10)

        # The T statistics
        expected = [32.866914806414108, -1.1635573550209353,
                    -35.857178728608552]
        for expected_t, observed_t in zip(expected, observed.t):
            self.assertAlmostEqual(expected_t, observed_t, places=10)

        # The p values
        expected = [2.7965174627917724e-217, 0.24465167231462448,
                    2.4882495142044017e-254]
        for expected_p, observed_p in zip(expected, observed.p):
            self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                                   places=10)

    @unittest.skipIf(platform.system() == "Darwin",
                     "multiprocessing not supported with Mac OS")
    def test_full_fit_linear_multiprocess(self):
        """Tests the full pipeline for linear regression with >1 processes."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_sim.txt.bz2",
        )

        # Creating the input files
        o_prefix, options = create_input_files(
            i_filename=data_filename,
            output_dirname=self.output_dir.name,
            analysis_type="linear",
            nb_process=2,
        )

        # Executing the tool
        try:
            main(args=options)
        finally:
            clean_logging_handlers()

        # Making sure the output file exists
        self.assertTrue(os.path.isfile(o_prefix + ".linear.dosage"))

        # Reading the data
        observed = pd.read_csv(o_prefix + ".linear.dosage", sep="\t")

        # Checking the shape
        self.assertEqual((3, 13), observed.shape)

        # Checking all columns are present
        self.assertEqual(["chr", "pos", "snp", "major", "minor", "maf", "n",
                          "coef", "se", "lower", "upper", "t", "p"],
                         list(observed.columns))

        # Chromosomes
        self.assertEqual([22], observed.chr.unique())

        # Positions
        self.assertEqual([1, 2, 3], list(observed.pos))

        # Marker names
        self.assertEqual(["marker_1", "marker_2", "marker_3"],
                         list(observed.snp))

        # Major alleles
        self.assertEqual(["T", "G", "AT"], list(observed.major))

        # Minor alleles
        self.assertEqual(["C", "A", "A"], list(observed.minor))

        # Minor allele frequency
        expected = [1724 / 11526, 4604 / 11526, 1379 / 11526]
        for expected_maf, observed_maf in zip(expected, observed.maf):
            self.assertAlmostEqual(expected_maf, observed_maf, places=10)

        # The number of samples
        expected = [5763, 5763, 5763]
        for expected_n, observed_n in zip(expected, observed.n):
            self.assertEqual(expected_n, observed_n)

        # The coefficients
        expected = [0.09930262321654575, -0.00279702443754753,
                    -0.11731595824657762]
        for expected_coef, observed_coef in zip(expected, observed.coef):
            self.assertAlmostEqual(expected_coef, observed_coef, places=10)

        # The standard error
        expected = [0.00302135517743109, 0.00240385609310785,
                    0.00327175651867383]
        for expected_se, observed_se in zip(expected, observed.se):
            self.assertAlmostEqual(expected_se, observed_se, places=10)

        # The lower CI
        expected = [0.09337963040899197, -0.0075094867313642991,
                    -0.12372983188610413]
        for expected_min_ci, observed_min_ci in zip(expected, observed.lower):
            self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=10)

        # The upper CI
        expected = [0.10522561602409949, 0.0019154378562692411,
                    -0.1109020846070511]
        for expected_max_ci, observed_max_ci in zip(expected, observed.upper):
            self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=10)

        # The T statistics
        expected = [32.866914806414108, -1.1635573550209353,
                    -35.857178728608552]
        for expected_t, observed_t in zip(expected, observed.t):
            self.assertAlmostEqual(expected_t, observed_t, places=10)

        # The p values
        expected = [2.7965174627917724e-217, 0.24465167231462448,
                    2.4882495142044017e-254]
        for expected_p, observed_p in zip(expected, observed.p):
            self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                                   places=10)

    def test_full_fit_linear_interaction(self):
        """Tests the full pipeline for linear regression with interaction."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_sim_inter.txt.bz2",
        )

        # Creating the input files
        o_prefix, options = create_input_files(
            i_filename=data_filename,
            output_dirname=self.output_dir.name,
            analysis_type="linear",
            interaction="gender",
        )

        # Executing the tool
        try:
            main(args=options)
        finally:
            clean_logging_handlers()

        # Making sure the output file exists
        self.assertTrue(os.path.isfile(o_prefix + ".linear.dosage"))

        # Reading the data
        observed = pd.read_csv(o_prefix + ".linear.dosage", sep="\t")

        # Checking the shape
        self.assertEqual((3, 13), observed.shape)

        # Checking all columns are present
        self.assertEqual(["chr", "pos", "snp", "major", "minor", "maf", "n",
                          "coef", "se", "lower", "upper", "t", "p"],
                         list(observed.columns))

        # Chromosomes
        self.assertEqual([22], observed.chr.unique())

        # Positions
        self.assertEqual([1, 2, 3], list(observed.pos))

        # Marker names
        self.assertEqual(["marker_1", "marker_2", "marker_3"],
                         list(observed.snp))

        # Major alleles
        self.assertEqual(["T", "G", "AT"], list(observed.major))

        # Minor alleles
        self.assertEqual(["C", "A", "A"], list(observed.minor))

        # Minor allele frequency
        expected = [1724 / 11526, 4604 / 11526, 1379 / 11526]
        for expected_maf, observed_maf in zip(expected, observed.maf):
            self.assertAlmostEqual(expected_maf, observed_maf, places=10)

        # The number of samples
        expected = [5763, 5763, 5763]
        for expected_n, observed_n in zip(expected, observed.n):
            self.assertEqual(expected_n, observed_n)

        # The coefficients
        expected = [0.1019595708452171734, -0.010917110807672320352,
                    0.02657634353683377762]
        for expected_coef, observed_coef in zip(expected, observed.coef):
            self.assertAlmostEqual(expected_coef, observed_coef, places=10)

        # The standard error
        expected = [0.005887641303827159493, 0.006578361146066042525,
                    0.009435607764482155033]
        for expected_se, observed_se in zip(expected, observed.se):
            self.assertAlmostEqual(expected_se, observed_se, places=10)

        # The lower CI
        expected = [0.0904175784866360355, -0.0238131739612693419,
                    0.00807900188569928013]
        for expected_min_ci, observed_min_ci in zip(expected, observed.lower):
            self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=10)

        # The upper CI
        expected = [0.113501563203798311, 0.00197895234592469944,
                    0.0450736851879682751]
        for expected_max_ci, observed_max_ci in zip(expected, observed.upper):
            self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=10)

        # The T statistics
        expected = [17.31755818394373136, -1.65954871817898208519,
                    2.816601134785761129]
        for expected_t, observed_t in zip(expected, observed.t):
            self.assertAlmostEqual(expected_t, observed_t, places=10)

        # The p values
        expected = [1.55270881064789290e-65, 9.70597574168038796e-02,
                    4.87000487450673421e-03]
        for expected_p, observed_p in zip(expected, observed.p):
            self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                                   places=10)


@unittest.skipIf(not HAS_STATSMODELS,
                 "optional requirement (statsmodels) not satisfied")
class TestImputedStatsLogistic(unittest.TestCase):

    def setUp(self):
        """Setup the tests."""
        # Creating the temporary directory
        self.output_dir = TemporaryDirectory(prefix="genipe_test_")

    def tearDown(self):
        """Finishes the test."""
        # Deleting the output directory
        self.output_dir.cleanup()

    def test_fit_logistic(self):
        """Tests the 'fit_logistic' function."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_sim.txt.bz2",
        )

        # This dataset contains 3 markers + 5 covariables
        data = pd.read_csv(data_filename, sep="\t", compression="bz2")

        # The formula for the first marker
        formula = "y_d ~ snp1 + C1 + C2 + C3 + age + C(gender)"
        columns_to_keep = ["y_d", "snp1", "C1", "C2", "C3", "age", "gender"]

        # The expected results for the first marker (according to R)
        expected_coef = -0.514309712761157334
        expected_se = 0.1148545370213162609
        expected_min_ci = -0.741288870474147710
        expected_max_ci = -0.290848443930784406
        expected_z = -4.477922475676383129
        expected_p = 7.53729612963228856e-06

        # The observed results for the first marker
        observed = fit_logistic(
            data=data[columns_to_keep].dropna(axis=0),
            formula=formula,
            result_col="snp1",
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_z, observed_p, = observed

        # Comparing the results
        self.assertAlmostEqual(expected_coef, observed_coef, places=10)
        self.assertAlmostEqual(expected_se, observed_se, places=7)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=2)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=2)
        self.assertAlmostEqual(expected_z, observed_z, places=6)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=5)

        # The formula for the second marker
        formula = "y_d ~ snp2 + C1 + C2 + C3 + age + C(gender)"
        columns_to_keep = ["y_d", "snp2", "C1", "C2", "C3", "age", "gender"]

        # The expected results for the second marker (according to R)
        expected_coef = -0.0409615621727592721
        expected_se = 0.0898086129043482589
        expected_min_ci = -0.217201661656268197
        expected_max_ci = 0.135039323830941221
        expected_z = -0.456098372395372320
        expected_p = 6.48319240531225471e-01

        # The observed results for the first marker
        observed = fit_logistic(
            data=data[columns_to_keep].dropna(axis=0),
            formula=formula,
            result_col="snp2",
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_z, observed_p, = observed

        # Comparing the results
        self.assertAlmostEqual(expected_coef, observed_coef, places=9)
        self.assertAlmostEqual(expected_se, observed_se, places=5)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=3)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=3)
        self.assertAlmostEqual(expected_z, observed_z, places=4)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=5)

        # The formula for the third (and last) marker
        formula = "y_d ~ snp3 + C1 + C2 + C3 + age + C(gender)"
        columns_to_keep = ["y_d", "snp3", "C1", "C2", "C3", "age", "gender"]

        # The expected results for the second marker (according to R)
        expected_coef = 0.6806154974808061864
        expected_se = 0.1216909125194588076
        expected_min_ci = 0.442504664626330035
        expected_max_ci = 0.919770526738653338
        expected_z = 5.5929854036715633825
        expected_p = 2.23198061422542684e-08

        # The observed results for the first marker
        observed = fit_logistic(
            data=data[columns_to_keep].dropna(axis=0),
            formula=formula,
            result_col="snp3",
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_z, observed_p, = observed

        # Comparing the results
        self.assertAlmostEqual(expected_coef, observed_coef, places=10)
        self.assertAlmostEqual(expected_se, observed_se, places=7)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=3)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=2)
        self.assertAlmostEqual(expected_z, observed_z, places=5)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=5)

        # Asking for an invalid column should raise a KeyError
        with self.assertRaises(KeyError) as cm:
            fit_logistic(
                data=data[columns_to_keep].dropna(axis=0),
                formula=formula,
                result_col="unknown",
            )

        with self.assertRaises(patsy.PatsyError) as cm:
            fit_logistic(
                data=data[columns_to_keep].dropna(axis=0),
                formula=formula + " + unknown",
                result_col="snp4",
            )

    def test_fit_logistic_interaction(self):
        """Tests the 'fit_logistic' function with interaction."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_sim_inter.txt.bz2",
        )

        # This dataset contains 3 markers + 5 covariables
        data = pd.read_csv(data_filename, sep="\t", compression="bz2")

        # The formula for the first marker
        formula = ("y_d ~ snp1 + C1 + C2 + C3 + age + C(gender) + "
                   "snp1*C(gender)")
        columns_to_keep = ["y_d", "snp1", "C1", "C2", "C3", "age", "gender"]

        # The expected results for the first marker (according to R)
        expected_coef = -0.4998229445983128905
        expected_se = 0.2425924038476814093
        expected_min_ci = -0.9779101205637230620
        expected_max_ci = -0.0263169240328645603
        expected_z = -2.0603404586078508665
        expected_p = 3.93660046768015970e-02

        # The observed results
        observed = fit_logistic(
            data=data[columns_to_keep],
            formula=formula,
            result_col="snp1:C(gender)[T.2]",
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_z, observed_p = observed

        # Comparing the results
        self.assertAlmostEqual(expected_coef, observed_coef, places=10)
        self.assertAlmostEqual(expected_se, observed_se, places=7)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=2)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=2)
        self.assertAlmostEqual(expected_z, observed_z, places=6)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=6)

    def test_full_fit_logistic(self):
        """Tests the full pipeline for logistic regression."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_sim.txt.bz2",
        )

        # Creating the input files
        o_prefix, options = create_input_files(
            i_filename=data_filename,
            output_dirname=self.output_dir.name,
            analysis_type="logistic",
            pheno_name="y_d",
        )

        # Executing the tool
        try:
            main(args=options)
        finally:
            clean_logging_handlers()

        # Making sure the output file exists
        self.assertTrue(os.path.isfile(o_prefix + ".logistic.dosage"))

        # Reading the data
        observed = pd.read_csv(o_prefix + ".logistic.dosage", sep="\t")

        # Checking the shape
        self.assertEqual((3, 13), observed.shape)

        # Checking all columns are present
        self.assertEqual(["chr", "pos", "snp", "major", "minor", "maf", "n",
                          "coef", "se", "lower", "upper", "z", "p"],
                         list(observed.columns))

        # Chromosomes
        self.assertEqual([22], observed.chr.unique())

        # Positions
        self.assertEqual([1, 2, 3], list(observed.pos))

        # Marker names
        self.assertEqual(["marker_1", "marker_2", "marker_3"],
                         list(observed.snp))

        # Major alleles
        self.assertEqual(["T", "G", "AT"], list(observed.major))

        # Minor alleles
        self.assertEqual(["C", "A", "A"], list(observed.minor))

        # Minor allele frequency
        expected = [1778 / 11880, 4703 / 11760, 1427 / 11880]
        for expected_maf, observed_maf in zip(expected, observed.maf):
            self.assertAlmostEqual(expected_maf, observed_maf, places=10)

        # The number of samples
        expected = [5940, 5880, 5940]
        for expected_n, observed_n in zip(expected, observed.n):
            self.assertEqual(expected_n, observed_n)

        # The coefficients
        expected = [-0.514309712761163662, -0.0409615621727604101,
                    0.6806154974808032998]
        places = [10, 9, 10]
        zipped = zip(expected, observed.coef, places)
        for expected_coef, observed_coef, place in zipped:
            self.assertAlmostEqual(expected_coef, observed_coef, places=place)

        # The standard error
        expected = [0.1148545370213169270, 0.0898086129043478426,
                    0.1216909125194589325]
        places = [7, 5, 7]
        zipped = zip(expected, observed.se, places)
        for expected_se, observed_se, place in zipped:
            self.assertAlmostEqual(expected_se, observed_se, places=place)

        # The lower CI
        expected = [-0.741288870474143935, -0.217201661656271972,
                    0.442504664626339528]
        places = [2, 3, 3]
        zipped = zip(expected, observed.lower, places)
        for expected_min_ci, observed_min_ci, place in zipped:
            self.assertAlmostEqual(expected_min_ci, observed_min_ci,
                                   places=place)

        # The upper CI
        expected = [-0.290848443930769640, 0.135039323830934560,
                    0.919770526738648897]
        places = [2, 3, 2]
        zipped = zip(expected, observed.upper, places)
        for expected_max_ci, observed_max_ci, place in zipped:
            self.assertAlmostEqual(expected_max_ci, observed_max_ci,
                                   places=place)

        # The Z statistics
        expected = [-4.477922475676412439, -0.456098372395387086,
                    5.592985403671533184]
        places = [6, 4, 5]
        zipped = zip(expected, observed.z, places)
        for expected_z, observed_z, place in zipped:
            self.assertAlmostEqual(expected_z, observed_z, places=place)

        # The p values
        expected = [7.53729612963125518e-06, 6.48319240531214813e-01,
                    2.23198061422581463e-08]
        for expected_p, observed_p in zip(expected, observed.p):
            self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                                   places=5)

    @unittest.skipIf(platform.system() == "Darwin",
                     "multiprocessing not supported with Mac OS")
    def test_full_fit_logistic_multiprocess(self):
        """Tests the full pipeline, logistic regression with >1 processes."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_sim.txt.bz2",
        )

        # Creating the input files
        o_prefix, options = create_input_files(
            i_filename=data_filename,
            output_dirname=self.output_dir.name,
            analysis_type="logistic",
            pheno_name="y_d",
            nb_process=2,
        )

        # Executing the tool
        try:
            main(args=options)
        finally:
            clean_logging_handlers()

        # Making sure the output file exists
        self.assertTrue(os.path.isfile(o_prefix + ".logistic.dosage"))

        # Reading the data
        observed = pd.read_csv(o_prefix + ".logistic.dosage", sep="\t")

        # Checking the shape
        self.assertEqual((3, 13), observed.shape)

        # Checking all columns are present
        self.assertEqual(["chr", "pos", "snp", "major", "minor", "maf", "n",
                          "coef", "se", "lower", "upper", "z", "p"],
                         list(observed.columns))

        # Chromosomes
        self.assertEqual([22], observed.chr.unique())

        # Positions
        self.assertEqual([1, 2, 3], list(observed.pos))

        # Marker names
        self.assertEqual(["marker_1", "marker_2", "marker_3"],
                         list(observed.snp))

        # Major alleles
        self.assertEqual(["T", "G", "AT"], list(observed.major))

        # Minor alleles
        self.assertEqual(["C", "A", "A"], list(observed.minor))

        # Minor allele frequency
        expected = [1778 / 11880, 4703 / 11760, 1427 / 11880]
        for expected_maf, observed_maf in zip(expected, observed.maf):
            self.assertAlmostEqual(expected_maf, observed_maf, places=10)

        # The number of samples
        expected = [5940, 5880, 5940]
        for expected_n, observed_n in zip(expected, observed.n):
            self.assertEqual(expected_n, observed_n)

        # The coefficients
        expected = [-0.514309712761163662, -0.0409615621727604101,
                    0.6806154974808032998]
        places = [10, 9, 10]
        zipped = zip(expected, observed.coef, places)
        for expected_coef, observed_coef, place in zipped:
            self.assertAlmostEqual(expected_coef, observed_coef, places=place)

        # The standard error
        expected = [0.1148545370213169270, 0.0898086129043478426,
                    0.1216909125194589325]
        places = [7, 5, 7]
        zipped = zip(expected, observed.se, places)
        for expected_se, observed_se, place in zipped:
            self.assertAlmostEqual(expected_se, observed_se, places=place)

        # The lower CI
        expected = [-0.741288870474143935, -0.217201661656271972,
                    0.442504664626339528]
        places = [2, 3, 3]
        zipped = zip(expected, observed.lower, places)
        for expected_min_ci, observed_min_ci, place in zipped:
            self.assertAlmostEqual(expected_min_ci, observed_min_ci,
                                   places=place)

        # The upper CI
        expected = [-0.290848443930769640, 0.135039323830934560,
                    0.919770526738648897]
        places = [2, 3, 2]
        zipped = zip(expected, observed.upper, places)
        for expected_max_ci, observed_max_ci, place in zipped:
            self.assertAlmostEqual(expected_max_ci, observed_max_ci,
                                   places=place)

        # The Z statistics
        expected = [-4.477922475676412439, -0.456098372395387086,
                    5.592985403671533184]
        places = [6, 4, 5]
        zipped = zip(expected, observed.z, places)
        for expected_z, observed_z, place in zipped:
            self.assertAlmostEqual(expected_z, observed_z, places=place)

        # The p values
        expected = [7.53729612963125518e-06, 6.48319240531214813e-01,
                    2.23198061422581463e-08]
        for expected_p, observed_p in zip(expected, observed.p):
            self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                                   places=5)

    def test_full_fit_logistic_interaction(self):
        """Tests the full pipeline for logistic regression with interaction."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_sim_inter.txt.bz2",
        )

        # Creating the input files
        o_prefix, options = create_input_files(
            i_filename=data_filename,
            output_dirname=self.output_dir.name,
            analysis_type="logistic",
            pheno_name="y_d",
            interaction="gender",
        )

        # Executing the tool
        try:
            main(args=options)
        finally:
            clean_logging_handlers()

        # Making sure the output file exists
        self.assertTrue(os.path.isfile(o_prefix + ".logistic.dosage"))

        # Reading the data
        observed = pd.read_csv(o_prefix + ".logistic.dosage", sep="\t")

        # Checking the shape
        self.assertEqual((3, 13), observed.shape)

        # Checking all columns are present
        self.assertEqual(["chr", "pos", "snp", "major", "minor", "maf", "n",
                          "coef", "se", "lower", "upper", "z", "p"],
                         list(observed.columns))

        # Chromosomes
        self.assertEqual([22], observed.chr.unique())

        # Positions
        self.assertEqual([1, 2, 3], list(observed.pos))

        # Marker names
        self.assertEqual(["marker_1", "marker_2", "marker_3"],
                         list(observed.snp))

        # Major alleles
        self.assertEqual(["T", "G", "AT"], list(observed.major))

        # Minor alleles
        self.assertEqual(["C", "A", "A"], list(observed.minor))

        # Minor allele frequency
        expected = [1778 / 11880, 4703 / 11760, 1427 / 11880]
        for expected_maf, observed_maf in zip(expected, observed.maf):
            self.assertAlmostEqual(expected_maf, observed_maf, places=10)

        # The number of samples
        expected = [5940, 5880, 5940]
        for expected_n, observed_n in zip(expected, observed.n):
            self.assertEqual(expected_n, observed_n)

        # The coefficients
        expected = [-0.4998229445983128905, 0.479196060192496665,
                    -0.081500925621293116]
        places = [10, 9, 10]
        zipped = zip(expected, observed.coef, places)
        for expected_coef, observed_coef, place in zipped:
            self.assertAlmostEqual(expected_coef, observed_coef, places=place)

        # The standard error
        expected = [0.2425924038476814093, 0.1732540442374335965,
                    0.236081283702105793]
        places = [7, 6, 7]
        zipped = zip(expected, observed.se, places)
        for expected_se, observed_se, place in zipped:
            self.assertAlmostEqual(expected_se, observed_se, places=place)

        # The lower CI
        expected = [-0.9779101205637230620, 0.140220059149237547,
                    -0.544942749175234886]
        for expected_min_ci, observed_min_ci in zip(expected, observed.lower):
            self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=2)

        # The upper CI
        expected = [-0.0263169240328645603, 0.819710161131259829,
                    0.380931176186822595]
        places = [2, 2, 3]
        zipped = zip(expected, observed.upper, places)
        for expected_max_ci, observed_max_ci, place in zipped:
            self.assertAlmostEqual(expected_max_ci, observed_max_ci,
                                   places=place)

        # The Z statistics
        expected = [-2.0603404586078508665, 2.765857860932753098,
                    -0.345224002272595865]
        places = [6, 4, 7]
        zipped = zip(expected, observed.z, places)
        for expected_z, observed_z, place in zipped:
            self.assertAlmostEqual(expected_z, observed_z, places=place)

        # The p values
        expected = [3.93660046768015970e-02, 5.67732747277691768e-03,
                    7.29925975515044678e-01]
        places = [6, 4, 7]
        zipped = zip(expected, observed.p, places)
        for expected_p, observed_p, place in zipped:
            self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                                   places=place)


@unittest.skipIf(not HAS_STATSMODELS,
                 "optional requirement (statsmodels) not satisfied")
class TestImputedStatsMixedLM(unittest.TestCase):

    def setUp(self):
        """Setup the tests."""
        # Creating the temporary directory
        self.output_dir = TemporaryDirectory(prefix="genipe_test_")

    def tearDown(self):
        """Finishes the test."""
        # Deleting the output directory
        self.output_dir.cleanup()

    def test_fit_mixedlm(self):
        """Tests the 'fit_mixedlm' function."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_mixedlm.txt.bz2",
        )

        # This dataset contains 3 markers + 5 covariables
        data = pd.read_csv(data_filename, sep="\t", compression="bz2")
        data = data[data.gender != 0].set_index("SampleID").dropna(axis=0)

        # The formula for the first marker
        formula = "y ~ snp1 + C1 + C2 + C3 + age + C(gender)"
        columns_to_keep = ["y", "snp1", "C1", "C2", "C3", "age", "gender"]

        # The expected results for the first marker (according to R)
        expected_coef = 0.12265168980987724
        expected_se = 0.041715401481444134
        expected_min_ci = 0.04090301809382547
        expected_max_ci = 0.20440036152516114
        expected_z = 2.9402015911182162
        expected_p = 3.279987709238652e-03

        # The observed results for the first marker
        observed = fit_mixedlm(
            data=data[columns_to_keep],
            formula=formula,
            groups=data.index.values,
            result_col="snp1",
            use_ml=False,
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_z, observed_p, = observed

        # Comparing the results
        self.assertAlmostEqual(expected_coef, observed_coef, places=10)
        self.assertAlmostEqual(expected_se, observed_se, places=10)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=4)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=4)
        self.assertAlmostEqual(expected_z, observed_z, places=8)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=8)

        # The formula for the second marker
        formula = "y ~ snp2 + C1 + C2 + C3 + age + C(gender)"
        columns_to_keep = ["y", "snp2", "C1", "C2", "C3", "age", "gender"]

        # The expected results for the second marker (according to R)
        expected_coef = 0.007945559606616614
        expected_se = 0.030465888551648834
        expected_min_ci = -0.05175771144653040
        expected_max_ci = 0.06764883066030558
        expected_z = 0.2608018339313000
        expected_p = 7.942453301062471e-01

        # The observed results for the first marker
        observed = fit_mixedlm(
            data=data[columns_to_keep],
            formula=formula,
            groups=data.index.values,
            result_col="snp2",
            use_ml=False,
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_z, observed_p, = observed

        # Comparing the results
        self.assertAlmostEqual(expected_coef, observed_coef, places=9)
        self.assertAlmostEqual(expected_se, observed_se, places=10)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=4)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=4)
        self.assertAlmostEqual(expected_z, observed_z, places=8)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=8)

        # The formula for the third (and last) marker
        formula = "y ~ snp3 + C1 + C2 + C3 + age + C(gender)"
        columns_to_keep = ["y", "snp3", "C1", "C2", "C3", "age", "gender"]

        # The expected results for the second marker (according to R)
        expected_coef = -0.15454491556451655
        expected_se = 0.045842987304696339
        expected_min_ci = -0.24438231821546189
        expected_max_ci = -0.06470751291475632
        expected_z = -3.3711789883440764
        expected_p = 7.484721142754225e-04

        # The observed results for the first marker
        observed = fit_mixedlm(
            data=data[columns_to_keep],
            formula=formula,
            groups=data.index.values,
            result_col="snp3",
            use_ml=False,
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_z, observed_p, = observed

        # Comparing the results
        self.assertAlmostEqual(expected_coef, observed_coef, places=9)
        self.assertAlmostEqual(expected_se, observed_se, places=10)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=4)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=4)
        self.assertAlmostEqual(expected_z, observed_z, places=8)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=8)

        # Asking for an invalid column should raise a KeyError
        with self.assertRaises(KeyError) as cm:
            fit_mixedlm(
                data=data[columns_to_keep].dropna(axis=0),
                formula=formula,
                groups=data.index.values,
                result_col="unknown",
                use_ml=False,
            )

        with self.assertRaises(patsy.PatsyError) as cm:
            fit_mixedlm(
                data=data[columns_to_keep].dropna(axis=0),
                formula=formula + " + unknown",
                groups=data.index.values,
                result_col="snp4",
                use_ml=False,
            )

    def test_fit_mixedlm_use_ml(self):
        """Tests the 'fit_mixedlm' function (using ML instead of REML)."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_mixedlm.txt.bz2",
        )

        # This dataset contains 3 markers + 5 covariables
        data = pd.read_csv(data_filename, sep="\t", compression="bz2")
        data = data[data.gender != 0].set_index("SampleID").dropna(axis=0)

        # The formula for the first marker
        formula = "y ~ snp1 + C1 + C2 + C3 + age + C(gender)"
        columns_to_keep = ["y", "snp1", "C1", "C2", "C3", "age", "gender"]

        # The expected results for the first marker (according to R)
        expected_coef = 0.1226516898098772429
        expected_se = 0.04170695572479218849
        expected_min_ci = 0.0409030180938254675
        expected_max_ci = 0.2044003615251611394
        expected_z = 2.940796988857387273
        expected_p = 3.27369009811673273e-03

        # The observed results for the first marker
        observed = fit_mixedlm(
            data=data[columns_to_keep],
            formula=formula,
            groups=data.index.values,
            result_col="snp1",
            use_ml=True,
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_z, observed_p, = observed

        # Comparing the results
        self.assertAlmostEqual(expected_coef, observed_coef, places=10)
        self.assertAlmostEqual(expected_se, observed_se, places=10)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=5)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=5)
        self.assertAlmostEqual(expected_z, observed_z, places=8)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=8)

        # The formula for the second marker
        formula = "y ~ snp2 + C1 + C2 + C3 + age + C(gender)"
        columns_to_keep = ["y", "snp2", "C1", "C2", "C3", "age", "gender"]

        # The expected results for the second marker (according to R)
        expected_coef = 0.00794555960661661419
        expected_se = 0.03045972038661254669
        expected_min_ci = -0.0517577114465304008
        expected_max_ci = 0.0676488306603055811
        expected_z = 0.260854646916220345
        expected_p = 7.94204600712016484e-01

        # The observed results for the first marker
        observed = fit_mixedlm(
            data=data[columns_to_keep],
            formula=formula,
            groups=data.index.values,
            result_col="snp2",
            use_ml=True,
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_z, observed_p, = observed

        # Comparing the results
        self.assertAlmostEqual(expected_coef, observed_coef, places=9)
        self.assertAlmostEqual(expected_se, observed_se, places=10)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=5)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=5)
        self.assertAlmostEqual(expected_z, observed_z, places=8)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=8)

        # The formula for the third (and last) marker
        formula = "y ~ snp3 + C1 + C2 + C3 + age + C(gender)"
        columns_to_keep = ["y", "snp3", "C1", "C2", "C3", "age", "gender"]

        # The expected results for the second marker (according to R)
        expected_coef = -0.1545449155645165495
        expected_se = 0.0458337058714313184
        expected_min_ci = -0.2443823182154618856
        expected_max_ci = -0.0647075129147563211
        expected_z = -3.371861660020081253
        expected_p = 7.46619406491788595e-04

        # The observed results for the first marker
        observed = fit_mixedlm(
            data=data[columns_to_keep],
            formula=formula,
            groups=data.index.values,
            result_col="snp3",
            use_ml=True,
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_z, observed_p, = observed

        # Comparing the results
        self.assertAlmostEqual(expected_coef, observed_coef, places=9)
        self.assertAlmostEqual(expected_se, observed_se, places=10)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=5)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=5)
        self.assertAlmostEqual(expected_z, observed_z, places=8)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=8)

        # Asking for an invalid column should raise a KeyError
        with self.assertRaises(KeyError) as cm:
            fit_mixedlm(
                data=data[columns_to_keep].dropna(axis=0),
                formula=formula,
                groups=data.index.values,
                result_col="unknown",
                use_ml=False,
            )

        with self.assertRaises(patsy.PatsyError) as cm:
            fit_mixedlm(
                data=data[columns_to_keep].dropna(axis=0),
                formula=formula + " + unknown",
                groups=data.index.values,
                result_col="snp4",
                use_ml=False,
            )

    def test_fit_mixedlm_interaction(self):
        """Tests the 'fit_mixedlm' function with interaction."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_mixedlm.txt.bz2",
        )

        # This dataset contains 3 markers + 5 covariables
        data = pd.read_csv(data_filename, sep="\t", compression="bz2")
        data = data[data.gender != 0].set_index("SampleID").dropna(axis=0)

        # The formula for the first marker
        formula = "y ~ snp1 + C1 + C2 + C3 + age + C(gender) + snp1*C(gender)"
        columns_to_keep = ["y", "snp1", "C1", "C2", "C3", "age", "gender"]

        # The expected results for the first marker (according to R)
        expected_coef = -0.04619683922092099
        expected_se = 0.08342819867401344
        expected_min_ci = -0.20968434904114919
        expected_max_ci = 0.11729067060054199
        expected_z = -0.5537317112818183
        expected_p = 0.5797624694538319190

        # The observed results
        observed = fit_mixedlm(
            data=data[columns_to_keep],
            formula=formula,
            groups=data.index.values,
            result_col="snp1:C(gender)[T.2]",
            use_ml=False,
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_z, observed_p = observed

        # Comparing the results
        self.assertAlmostEqual(expected_coef, observed_coef, places=10)
        self.assertAlmostEqual(expected_se, observed_se, places=10)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=4)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=4)
        self.assertAlmostEqual(expected_z, observed_z, places=9)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=9)

    def test_fit_mixedlm_interaction_use_ml(self):
        """Tests the 'fit_mixedlm' function with interaction (using ML)."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_mixedlm.txt.bz2",
        )

        # This dataset contains 3 markers + 5 covariables
        data = pd.read_csv(data_filename, sep="\t", compression="bz2")
        data = data[data.gender != 0].set_index("SampleID").dropna(axis=0)

        # The formula for the first marker
        formula = "y ~ snp1 + C1 + C2 + C3 + age + C(gender) + snp1*C(gender)"
        columns_to_keep = ["y", "snp1", "C1", "C2", "C3", "age", "gender"]

        # The expected results for the first marker (according to R)
        expected_coef = -0.0461968392209209863
        expected_se = 0.08340889441074257615
        expected_min_ci = -0.2096843490411491873
        expected_max_ci = 0.1172906706005419908
        expected_z = -0.553859867671032235
        expected_p = 5.79674752511348590e-01

        # The observed results
        observed = fit_mixedlm(
            data=data[columns_to_keep],
            formula=formula,
            groups=data.index.values,
            result_col="snp1:C(gender)[T.2]",
            use_ml=True,
        )
        self.assertEqual(6, len(observed))
        observed_coef, observed_se, observed_min_ci, observed_max_ci, \
            observed_z, observed_p = observed

        # Comparing the results
        self.assertAlmostEqual(expected_coef, observed_coef, places=10)
        self.assertAlmostEqual(expected_se, observed_se, places=10)
        self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=4)
        self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=4)
        self.assertAlmostEqual(expected_z, observed_z, places=9)
        self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                               places=9)

    def test_full_fit_mixedlm(self):
        """Tests the full pipeline for mixed effect model."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_mixedlm.txt.bz2",
        )

        # Creating the input files
        o_prefix, options = create_input_files(
            i_filename=data_filename,
            output_dirname=self.output_dir.name,
            analysis_type="mixedlm",
            pheno_name="y",
            sample_col="SampleID",
        )

        # Executing the tool
        try:
            main(args=options)
        finally:
            clean_logging_handlers()

        # Making sure the output file exists
        self.assertTrue(os.path.isfile(o_prefix + ".mixedlm.dosage"))

        # Reading the data
        observed = pd.read_csv(o_prefix + ".mixedlm.dosage", sep="\t")

        # Checking the shape
        self.assertEqual((3, 13), observed.shape)

        # Checking all columns are present
        self.assertEqual(["chr", "pos", "snp", "major", "minor", "maf", "n",
                          "coef", "se", "lower", "upper", "z", "p"],
                         list(observed.columns))

        # Chromosomes
        self.assertEqual([22], observed.chr.unique())

        # Positions
        self.assertEqual([1, 2, 3], list(observed.pos))

        # Marker names
        self.assertEqual(["marker_1", "marker_2", "marker_3"],
                         list(observed.snp))

        # Major alleles
        self.assertEqual(["T", "G", "AT"], list(observed.major))

        # Minor alleles
        self.assertEqual(["C", "A", "A"], list(observed.minor))

        # Minor allele frequency
        expected = [1724 / 11526, 4605 / 11528, 1379 / 11528]
        for expected_maf, observed_maf in zip(expected, observed.maf):
            self.assertAlmostEqual(expected_maf, observed_maf, places=10)

        # The number of samples
        expected = [5763, 5764, 5764]
        for expected_n, observed_n in zip(expected, observed.n):
            self.assertEqual(expected_n, observed_n)

        # The coefficients
        expected = [0.12265168980987724, 0.007931606465841023,
                    -0.15449981039313940]
        for expected_coef, observed_coef in zip(expected, observed.coef):
            self.assertAlmostEqual(expected_coef, observed_coef, places=9)

        # The standard error
        expected = [0.041715401481444134, 0.030464140147980882,
                    0.04583966116454336]
        zipped = zip(expected, observed.se)
        for expected_se, observed_se in zipped:
            self.assertAlmostEqual(expected_se, observed_se, places=10)

        # The lower CI
        expected = [0.04090301809382547, -0.05176823980480887,
                    -0.24433069717817782]
        for expected_min_ci, observed_min_ci in zip(expected, observed.lower):
            self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=4)

        # The upper CI
        expected = [0.20440036152516114, 0.06763145273676965,
                    -0.06466892360820460]
        for expected_max_ci, observed_max_ci in zip(expected, observed.upper):
            self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=4)

        # The Z statistics
        expected = [2.9402015911182162, 0.2603587833863980,
                    -3.3704396251655511]
        for expected_z, observed_z in zip(expected, observed.z):
            self.assertAlmostEqual(expected_z, observed_z, places=8)

        # The p values
        expected = [3.279987709238652e-03, 7.945870329397586e-01,
                    7.504834930376347e-04]
        for expected_p, observed_p in zip(expected, observed.p):
            self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                                   places=8)

    def test_full_fit_mixedlm_use_ml(self):
        """Tests the full pipeline for mixed effect model (using ML)."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_mixedlm.txt.bz2",
        )

        # Creating the input files
        o_prefix, options = create_input_files(
            i_filename=data_filename,
            output_dirname=self.output_dir.name,
            analysis_type="mixedlm",
            pheno_name="y",
            sample_col="SampleID",
            use_ml=True,
        )

        # Executing the tool
        try:
            main(args=options)
        finally:
            clean_logging_handlers()

        # Making sure the output file exists
        self.assertTrue(os.path.isfile(o_prefix + ".mixedlm.dosage"))

        # Reading the data
        observed = pd.read_csv(o_prefix + ".mixedlm.dosage", sep="\t")

        # Checking the shape
        self.assertEqual((3, 13), observed.shape)

        # Checking all columns are present
        self.assertEqual(["chr", "pos", "snp", "major", "minor", "maf", "n",
                          "coef", "se", "lower", "upper", "z", "p"],
                         list(observed.columns))

        # Chromosomes
        self.assertEqual([22], observed.chr.unique())

        # Positions
        self.assertEqual([1, 2, 3], list(observed.pos))

        # Marker names
        self.assertEqual(["marker_1", "marker_2", "marker_3"],
                         list(observed.snp))

        # Major alleles
        self.assertEqual(["T", "G", "AT"], list(observed.major))

        # Minor alleles
        self.assertEqual(["C", "A", "A"], list(observed.minor))

        # Minor allele frequency
        expected = [1724 / 11526, 4605 / 11528, 1379 / 11528]
        for expected_maf, observed_maf in zip(expected, observed.maf):
            self.assertAlmostEqual(expected_maf, observed_maf, places=10)

        # The number of samples
        expected = [5763, 5764, 5764]
        for expected_n, observed_n in zip(expected, observed.n):
            self.assertEqual(expected_n, observed_n)

        # The coefficients
        expected = [0.1226516898098772429, 0.00793160646584102277,
                    -0.1544998103931393973]
        for expected_coef, observed_coef in zip(expected, observed.coef):
            self.assertAlmostEqual(expected_coef, observed_coef, places=9)

        # The standard error
        expected = [0.04170695572479218849, 0.03045797340709452752,
                    0.04583038201498106784]
        zipped = zip(expected, observed.se)
        for expected_se, observed_se in zipped:
            self.assertAlmostEqual(expected_se, observed_se, places=10)

        # The lower CI
        expected = [0.0409030180938254675, -0.0517682398048088729,
                    -0.2443306971781778192]
        places = [4, 5, 5]
        zipped = zip(expected, observed.lower, places)
        for expected_min_ci, observed_min_ci, place in zipped:
            self.assertAlmostEqual(expected_min_ci, observed_min_ci,
                                   places=place)

        # The upper CI
        expected = [0.2044003615251611394, 0.0676314527367696539,
                    -0.0646689236082046010]
        places = [4, 5, 5]
        zipped = zip(expected, observed.upper, places)
        for expected_max_ci, observed_max_ci, place in zipped:
            self.assertAlmostEqual(expected_max_ci, observed_max_ci,
                                   places=place)

        # The Z statistics
        expected = [2.940796988857387273, 0.260411497502769707,
                    -3.371122028671643900]
        for expected_z, observed_z in zip(expected, observed.z):
            self.assertAlmostEqual(expected_z, observed_z, places=8)

        # The p values
        expected = [3.27369009811673273e-03, 7.94546375098304614e-01,
                    7.48626890377357412e-04]
        for expected_p, observed_p in zip(expected, observed.p):
            self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                                   places=8)

    @unittest.skipIf(platform.system() == "Darwin",
                     "multiprocessing not supported with Mac OS")
    def test_full_fit_mixedlm_multiprocess(self):
        """Tests the full pipeline, mixed linear model with >1 processes."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_mixedlm.txt.bz2",
        )

        # Creating the input files
        o_prefix, options = create_input_files(
            i_filename=data_filename,
            output_dirname=self.output_dir.name,
            analysis_type="mixedlm",
            pheno_name="y",
            sample_col="SampleID",
            nb_process=2,
        )

        # Executing the tool
        try:
            main(args=options)
        finally:
            clean_logging_handlers()

        # Making sure the output file exists
        self.assertTrue(os.path.isfile(o_prefix + ".mixedlm.dosage"))

        # Reading the data
        observed = pd.read_csv(o_prefix + ".mixedlm.dosage", sep="\t")

        # Checking the shape
        self.assertEqual((3, 13), observed.shape)

        # Checking all columns are present
        self.assertEqual(["chr", "pos", "snp", "major", "minor", "maf", "n",
                          "coef", "se", "lower", "upper", "z", "p"],
                         list(observed.columns))

        # Chromosomes
        self.assertEqual([22], observed.chr.unique())

        # Positions
        self.assertEqual([1, 2, 3], list(observed.pos))

        # Marker names
        self.assertEqual(["marker_1", "marker_2", "marker_3"],
                         list(observed.snp))

        # Major alleles
        self.assertEqual(["T", "G", "AT"], list(observed.major))

        # Minor alleles
        self.assertEqual(["C", "A", "A"], list(observed.minor))

        # Minor allele frequency
        expected = [1724 / 11526, 4605 / 11528, 1379 / 11528]
        for expected_maf, observed_maf in zip(expected, observed.maf):
            self.assertAlmostEqual(expected_maf, observed_maf, places=10)

        # The number of samples
        expected = [5763, 5764, 5764]
        for expected_n, observed_n in zip(expected, observed.n):
            self.assertEqual(expected_n, observed_n)

        # The coefficients
        expected = [0.12265168980987724, 0.007931606465841023,
                    -0.15449981039313940]
        for expected_coef, observed_coef in zip(expected, observed.coef):
            self.assertAlmostEqual(expected_coef, observed_coef, places=9)

        # The standard error
        expected = [0.041715401481444134, 0.030464140147980882,
                    0.04583966116454336]
        zipped = zip(expected, observed.se)
        for expected_se, observed_se in zipped:
            self.assertAlmostEqual(expected_se, observed_se, places=10)

        # The lower CI
        expected = [0.04090301809382547, -0.05176823980480887,
                    -0.24433069717817782]
        for expected_min_ci, observed_min_ci in zip(expected, observed.lower):
            self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=4)

        # The upper CI
        expected = [0.20440036152516114, 0.06763145273676965,
                    -0.06466892360820460]
        for expected_max_ci, observed_max_ci in zip(expected, observed.upper):
            self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=4)

        # The Z statistics
        expected = [2.9402015911182162, 0.2603587833863980,
                    -3.3704396251655511]
        for expected_z, observed_z in zip(expected, observed.z):
            self.assertAlmostEqual(expected_z, observed_z, places=8)

        # The p values
        expected = [3.279987709238652e-03, 7.945870329397586e-01,
                    7.504834930376347e-04]
        for expected_p, observed_p in zip(expected, observed.p):
            self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                                   places=8)

    @unittest.skipIf(platform.system() == "Darwin",
                     "multiprocessing not supported with Mac OS")
    def test_full_fit_mixedlm_multiprocess_use_ml(self):
        """Tests the full pipeline, mixed linear with >1 processes (ML)."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_mixedlm.txt.bz2",
        )

        # Creating the input files
        o_prefix, options = create_input_files(
            i_filename=data_filename,
            output_dirname=self.output_dir.name,
            analysis_type="mixedlm",
            pheno_name="y",
            sample_col="SampleID",
            nb_process=2,
            use_ml=True,
        )

        # Executing the tool
        try:
            main(args=options)
        finally:
            clean_logging_handlers()

        # Making sure the output file exists
        self.assertTrue(os.path.isfile(o_prefix + ".mixedlm.dosage"))

        # Reading the data
        observed = pd.read_csv(o_prefix + ".mixedlm.dosage", sep="\t")

        # Checking the shape
        self.assertEqual((3, 13), observed.shape)

        # Checking all columns are present
        self.assertEqual(["chr", "pos", "snp", "major", "minor", "maf", "n",
                          "coef", "se", "lower", "upper", "z", "p"],
                         list(observed.columns))

        # Chromosomes
        self.assertEqual([22], observed.chr.unique())

        # Positions
        self.assertEqual([1, 2, 3], list(observed.pos))

        # Marker names
        self.assertEqual(["marker_1", "marker_2", "marker_3"],
                         list(observed.snp))

        # Major alleles
        self.assertEqual(["T", "G", "AT"], list(observed.major))

        # Minor alleles
        self.assertEqual(["C", "A", "A"], list(observed.minor))

        # Minor allele frequency
        expected = [1724 / 11526, 4605 / 11528, 1379 / 11528]
        for expected_maf, observed_maf in zip(expected, observed.maf):
            self.assertAlmostEqual(expected_maf, observed_maf, places=10)

        # The number of samples
        expected = [5763, 5764, 5764]
        for expected_n, observed_n in zip(expected, observed.n):
            self.assertEqual(expected_n, observed_n)

        # The coefficients
        expected = [0.1226516898098772429, 0.00793160646584102277,
                    -0.1544998103931393973]
        for expected_coef, observed_coef in zip(expected, observed.coef):
            self.assertAlmostEqual(expected_coef, observed_coef, places=9)

        # The standard error
        expected = [0.04170695572479218849, 0.03045797340709452752,
                    0.04583038201498106784]
        zipped = zip(expected, observed.se)
        for expected_se, observed_se in zipped:
            self.assertAlmostEqual(expected_se, observed_se, places=10)

        # The lower CI
        expected = [0.0409030180938254675, -0.0517682398048088729,
                    -0.2443306971781778192]
        places = [4, 5, 5]
        zipped = zip(expected, observed.lower, places)
        for expected_min_ci, observed_min_ci, place in zipped:
            self.assertAlmostEqual(expected_min_ci, observed_min_ci,
                                   places=place)

        # The upper CI
        expected = [0.2044003615251611394, 0.0676314527367696539,
                    -0.0646689236082046010]
        places = [4, 5, 5]
        zipped = zip(expected, observed.upper, places)
        for expected_max_ci, observed_max_ci, place in zipped:
            self.assertAlmostEqual(expected_max_ci, observed_max_ci,
                                   places=place)

        # The Z statistics
        expected = [2.940796988857387273, 0.260411497502769707,
                    -3.371122028671643900]
        for expected_z, observed_z in zip(expected, observed.z):
            self.assertAlmostEqual(expected_z, observed_z, places=8)

        # The p values
        expected = [3.27369009811673273e-03, 7.94546375098304614e-01,
                    7.48626890377357412e-04]
        for expected_p, observed_p in zip(expected, observed.p):
            self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                                   places=8)

    def test_full_fit_mixedlm_interaction(self):
        """Tests the full pipeline for mixed linear model with interaction."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_mixedlm.txt.bz2",
        )

        # Creating the input files
        o_prefix, options = create_input_files(
            i_filename=data_filename,
            output_dirname=self.output_dir.name,
            analysis_type="mixedlm",
            pheno_name="y",
            sample_col="SampleID",
            interaction="gender",
        )

        # Executing the tool
        try:
            main(args=options)
        finally:
            clean_logging_handlers()

        # Making sure the output file exists
        self.assertTrue(os.path.isfile(o_prefix + ".mixedlm.dosage"))

        # Reading the data
        observed = pd.read_csv(o_prefix + ".mixedlm.dosage", sep="\t")

        # Checking the shape
        self.assertEqual((3, 13), observed.shape)

        # Checking all columns are present
        self.assertEqual(["chr", "pos", "snp", "major", "minor", "maf", "n",
                          "coef", "se", "lower", "upper", "z", "p"],
                         list(observed.columns))

        # Chromosomes
        self.assertEqual([22], observed.chr.unique())

        # Positions
        self.assertEqual([1, 2, 3], list(observed.pos))

        # Marker names
        self.assertEqual(["marker_1", "marker_2", "marker_3"],
                         list(observed.snp))

        # Major alleles
        self.assertEqual(["T", "G", "AT"], list(observed.major))

        # Minor alleles
        self.assertEqual(["C", "A", "A"], list(observed.minor))

        # Minor allele frequency
        expected = [1724 / 11526, 4605 / 11528, 1379 / 11528]
        for expected_maf, observed_maf in zip(expected, observed.maf):
            self.assertAlmostEqual(expected_maf, observed_maf, places=10)

        # The number of samples
        expected = [5763, 5764, 5764]
        for expected_n, observed_n in zip(expected, observed.n):
            self.assertEqual(expected_n, observed_n)

        # The coefficients
        expected = [-0.04619683922092099, 0.011409877542661442,
                    -0.003218416854476458]
        places = [10, 9, 10]
        zipped = zip(expected, observed.coef, places)
        for expected_coef, observed_coef, place in zipped:
            self.assertAlmostEqual(expected_coef, observed_coef, places=place)

        # The standard error
        expected = [0.08342819867401344, 0.060927215816902323,
                    0.091755870420978478]
        for expected_se, observed_se in zip(expected, observed.se):
            self.assertAlmostEqual(expected_se, observed_se, places=10)

        # The lower CI
        expected = [-0.20968434904114919, -0.10798427523465945,
                    -0.18302499859617510]
        for expected_min_ci, observed_min_ci in zip(expected, observed.lower):
            self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=4)

        # The upper CI
        expected = [0.11729067060054199, 0.13080403032042959,
                    0.17658816488845774]
        for expected_max_ci, observed_max_ci in zip(expected, observed.upper):
            self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=4)

        # The Z statistics
        expected = [-0.5537317112818183, 0.18727062101360839,
                    -0.03507586860339586]
        places = [9, 8, 9]
        zipped = zip(expected, observed.z, places)
        for expected_z, observed_z, place in zipped:
            self.assertAlmostEqual(expected_z, observed_z, places=place)

        # The p values
        expected = [0.5797624694538319190, 0.851448456503608009,
                    9.720192436335240e-01]
        places = [9, 8, 9]
        zipped = zip(expected, observed.p, places)
        for expected_p, observed_p, place in zipped:
            self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                                   places=place)

    def test_full_fit_mixedlm_interaction_use_ml(self):
        """Tests the full pipeline for mixed linear model (interaction, ML)."""
        # Reading the data
        data_filename = resource_filename(
            __name__,
            "data/regression_mixedlm.txt.bz2",
        )

        # Creating the input files
        o_prefix, options = create_input_files(
            i_filename=data_filename,
            output_dirname=self.output_dir.name,
            analysis_type="mixedlm",
            pheno_name="y",
            sample_col="SampleID",
            interaction="gender",
            use_ml=True,
        )

        # Executing the tool
        try:
            main(args=options)
        finally:
            clean_logging_handlers()

        # Making sure the output file exists
        self.assertTrue(os.path.isfile(o_prefix + ".mixedlm.dosage"))

        # Reading the data
        observed = pd.read_csv(o_prefix + ".mixedlm.dosage", sep="\t")

        # Checking the shape
        self.assertEqual((3, 13), observed.shape)

        # Checking all columns are present
        self.assertEqual(["chr", "pos", "snp", "major", "minor", "maf", "n",
                          "coef", "se", "lower", "upper", "z", "p"],
                         list(observed.columns))

        # Chromosomes
        self.assertEqual([22], observed.chr.unique())

        # Positions
        self.assertEqual([1, 2, 3], list(observed.pos))

        # Marker names
        self.assertEqual(["marker_1", "marker_2", "marker_3"],
                         list(observed.snp))

        # Major alleles
        self.assertEqual(["T", "G", "AT"], list(observed.major))

        # Minor alleles
        self.assertEqual(["C", "A", "A"], list(observed.minor))

        # Minor allele frequency
        expected = [1724 / 11526, 4605 / 11528, 1379 / 11528]
        for expected_maf, observed_maf in zip(expected, observed.maf):
            self.assertAlmostEqual(expected_maf, observed_maf, places=10)

        # The number of samples
        expected = [5763, 5764, 5764]
        for expected_n, observed_n in zip(expected, observed.n):
            self.assertEqual(expected_n, observed_n)

        # The coefficients
        expected = [-0.0461968392209209863, 0.01140987754266144238,
                    -0.00321841685447645789]
        places = [10, 9, 10]
        zipped = zip(expected, observed.coef, places)
        for expected_coef, observed_coef, place in zipped:
            self.assertAlmostEqual(expected_coef, observed_coef, places=place)

        # The standard error
        expected = [0.08340889441074257615, 0.06091312045157387667,
                    0.09173464292030783507]
        for expected_se, observed_se in zip(expected, observed.se):
            self.assertAlmostEqual(expected_se, observed_se, places=10)

        # The lower CI
        expected = [-0.2096843490411491873, -0.1079842752346594492,
                    -0.1830249985961751036]
        for expected_min_ci, observed_min_ci in zip(expected, observed.lower):
            self.assertAlmostEqual(expected_min_ci, observed_min_ci, places=4)

        # The upper CI
        expected = [0.1172906706005419908, 0.1308040303204295873,
                    0.1765881648884577437]
        for expected_max_ci, observed_max_ci in zip(expected, observed.upper):
            self.assertAlmostEqual(expected_max_ci, observed_max_ci, places=4)

        # The Z statistics
        expected = [-0.553859867671032235, 0.1873139556482306678,
                    -0.0350839851992706467]
        places = [9, 8, 9]
        zipped = zip(expected, observed.z, places)
        for expected_z, observed_z, place in zipped:
            self.assertAlmostEqual(expected_z, observed_z, places=place)

        # The p values
        expected = [5.79674752511348590e-01, 0.85141448161648325410,
                    9.72012771510516016e-01]
        places = [9, 8, 9]
        zipped = zip(expected, observed.p, places)
        for expected_p, observed_p, place in zipped:
            self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                                   places=place)


@unittest.skipIf(not HAS_SKAT, "SKAT is not installed")
class TestImputedStatsSkat(unittest.TestCase):

    tmp_dir = None
    args = None

    @classmethod
    def setUpClass(cls):
        cls.tmp_dir = TemporaryDirectory(prefix="genipe_test_")
        cls.args = cls.setup_skat_files(cls.tmp_dir.name)

    @classmethod
    def tearDownClass(cls):
        # Cleaning the temporary directory
        cls.tmp_dir.cleanup()

    def test_continuous(self):
        o_prefix = os.path.join(self.tmp_dir.name, "skat_test_continuous")
        args = self.args + [
            "--pheno-name", "outcome_continuous",
            "--outcome-type", "continuous",
            "--out", o_prefix,
        ]

        # Executing the tool
        try:
            main(args=args)
        finally:
            clean_logging_handlers()

        # The observed values
        results_filename = o_prefix + ".skat.dosage"
        observed = pd.read_csv(results_filename, header=0, sep="\t")
        self.assertEqual((3, 3), observed.shape)

        # The SNP set ID
        expected = ["set1", "set2", "set3"]
        for expected_id, observed_id in zip(expected, observed.snp_set_id):
            self.assertEqual(expected_id, observed_id)

        # The p values
        expected = [0.002877041, 0.09319144, 0.002877041]
        for expected_p, observed_p in zip(expected, observed.p_value):
            self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                                   places=10)

        # The Q statistics
        expected = [298041.5, 24870.14, 298041.5]
        for expected_q, observed_q in zip(expected, observed.q_value):
            self.assertAlmostEqual(expected_q, observed_q, places=10)

    def test_discrete(self):
        o_prefix = os.path.join(self.tmp_dir.name, "skat_test_discrete")
        args = self.args + [
            "--pheno-name", "outcome_discrete",
            "--outcome-type", "discrete",
            "--out", o_prefix,
        ]

        # Executing the tool
        try:
            main(args=args)
        finally:
            clean_logging_handlers()

        # The observed values
        results_filename = o_prefix + ".skat.dosage"
        observed = pd.read_csv(results_filename, header=0, sep="\t")
        self.assertEqual((3, 3), observed.shape)

        # The SNP set ID
        expected = ["set1", "set2", "set3"]
        for expected_id, observed_id in zip(expected, observed.snp_set_id):
            self.assertEqual(expected_id, observed_id)

        # The p values
        expected = [0.1401991, 0.5868036, 0.1401991]
        for expected_p, observed_p in zip(expected, observed.p_value):
            self.assertAlmostEqual(expected_p, observed_p, places=10)

        # The Q statistics
        expected = [34934.79, 2776.797, 34934.79]
        for expected_q, observed_q in zip(expected, observed.q_value):
            self.assertAlmostEqual(expected_q, observed_q, places=10)

    def test_continuous_multiprocess(self):
        o_prefix = os.path.join(self.tmp_dir.name, "skat_test_continuous_mp")
        args = self.args + [
            "--pheno-name", "outcome_continuous",
            "--outcome-type", "continuous",
            "--nb-process", "2",
            "--out", o_prefix,
        ]

        # Executing the tool
        try:
            main(args=args)
        finally:
            clean_logging_handlers()

        # The observed values
        results_filename = o_prefix + ".skat.dosage"
        observed = pd.read_csv(results_filename, header=0, sep="\t")
        self.assertEqual((3, 3), observed.shape)

        # The SNP set ID
        expected = ["set1", "set2", "set3"]
        for expected_id, observed_id in zip(expected, observed.snp_set_id):
            self.assertEqual(expected_id, observed_id)

        # The p values
        expected = [0.002877041, 0.09319144, 0.002877041]
        for expected_p, observed_p in zip(expected, observed.p_value):
            self.assertAlmostEqual(np.log10(expected_p), np.log10(observed_p),
                                   places=10)

        # The Q statistics
        expected = [298041.5, 24870.14, 298041.5]
        for expected_q, observed_q in zip(expected, observed.q_value):
            self.assertAlmostEqual(expected_q, observed_q, places=10)

    def test_discrete_multiprocess(self):
        o_prefix = os.path.join(self.tmp_dir.name, "skat_test_discrete_mp")
        args = self.args + [
            "--pheno-name", "outcome_discrete",
            "--outcome-type", "discrete",
            "--nb-process", "2",
            "--out", o_prefix,
        ]

        # Executing the tool
        try:
            main(args=args)
        finally:
            clean_logging_handlers()

        # The observed values
        results_filename = o_prefix + ".skat.dosage"
        observed = pd.read_csv(results_filename, header=0, sep="\t")
        self.assertEqual((3, 3), observed.shape)

        # The SNP set ID
        expected = ["set1", "set2", "set3"]
        for expected_id, observed_id in zip(expected, observed.snp_set_id):
            self.assertEqual(expected_id, observed_id)

        # The p values
        expected = [0.1401991, 0.5868036, 0.1401991]
        for expected_p, observed_p in zip(expected, observed.p_value):
            self.assertAlmostEqual(expected_p, observed_p, places=10)

        # The Q statistics
        expected = [34934.79, 2776.797, 34934.79]
        for expected_q, observed_q in zip(expected, observed.q_value):
            self.assertAlmostEqual(expected_q, observed_q, places=10)

    @staticmethod
    def setup_skat_files(out_directory):
        """Parses the SKAT example files into the format expected by genipe."""
        # Read the SKAT example files.
        skat_x = pd.read_csv(  # Covariates
            resource_filename(
                __name__,
                os.path.join("data", "skat_x_matrix.csv.bz2"),
            ),
            sep=",",
            compression="bz2",
            header=None,
            names=("covar_discrete", "covar_continuous"),
        )

        skat_z = pd.read_csv(  # Genotypes
            resource_filename(
                __name__,
                os.path.join("data", "skat_z_matrix.csv.bz2"),
            ),
            sep=",",
            compression="bz2",
            header=None,
        )

        skat_y_continuous = pd.read_csv(  # Outcome, continuous
            resource_filename(
                __name__,
                os.path.join("data", "skat_y_continuous_vector.csv.bz2"),
            ),
            sep=",",
            compression="bz2",
            header=None,
            names=("outcome_continuous", ),
        )

        skat_y_discrete = pd.read_csv(  # Outcome, discrete
            resource_filename(
                __name__,
                os.path.join("data", "skat_y_discrete_vector.csv.bz2"),
            ),
            sep=",",
            compression="bz2",
            header=None,
            names=("outcome_discrete", ),
        )

        # Write the Impute2 file.
        filename = os.path.join(out_directory, "impute2.txt")
        variants = []
        with open(filename, "w") as output_file:
            for col_idx in range(skat_z.shape[1]):

                chrom = random.randint(1, 22)
                pos = random.randint(1000, 90000000)
                a1, a2 = random.sample("ATGC", 2)

                variant_name = "{chrom}:{pos}:{a2}".format(
                    chrom=chrom, pos=pos, a2=a2
                )
                variants.append(variant_name)

                row = "{chrom} {name} {pos} {a1} {a2}".format(
                    chrom=chrom, name=variant_name, pos=pos, a1=a1, a2=a2,
                )

                dosage_list = []
                for sample_idx in range(skat_z.shape[0]):
                    dosage = reverse_dosage(
                        skat_z.iloc[sample_idx, col_idx]
                    )

                    dosage_list.extend(dosage)

                # Make sure all the probabilities are in the [0, 1] interval.
                for i, proba in enumerate(dosage_list):
                    if proba < 0:
                        dosage_list[i] = 0
                    elif proba > 1:
                        dosage_list[i] = 1

                print(row, *dosage_list, sep=" ", file=output_file)

        # Create the snp set file.
        filename = os.path.join(out_directory, "snp_sets.txt")
        with open(filename, "w") as output_file:
            print("variant", "snp_set", sep="\t", file=output_file)
            # The first SNP set
            for variant in variants:
                print(variant, "set1", sep="\t", file=output_file)

            # The second SNP set (which contains only the first 10 markers, and
            # the last 10)
            for variant in variants[:10] + [variants[30]] + variants[-10:]:
                print(variant, "set2", sep="\t", file=output_file)

            # The last SNP set (which is the same as the first one)
            for variant in variants:
                print(variant, "set3", sep="\t", file=output_file)

        # Create a list of samples.
        samples = ["sample_{}".format(i + 1) for i in range(skat_z.shape[0])]

        # Get the filename for the phenotype file.
        filename = os.path.join(out_directory, "phenotypes.txt")

        # Before we can write it, we need to add a samples column.
        skat_x["sample"] = samples
        skat_x = skat_x[["sample", "covar_discrete", "covar_continuous"]]

        # We also add the outcomes to this (genipe uses a single file for both
        # the outcome and the covariates).
        skat_x["outcome_continuous"] = skat_y_continuous.values
        skat_x["outcome_discrete"] = skat_y_discrete.values

        # Write the phenotype file.
        skat_x.to_csv(filename, sep="\t", index=False)

        # Get the filename for the samples file.
        filename = os.path.join(out_directory, "samples.txt")

        # Now we write the samples file by imitating the Impute2 format.
        with open(filename, "w") as f:
            print("ID_1 ID_2 missing father mother sex phenotype", file=f)
            print("0 0 0 D D D B", file=f)
            for sample in samples:
                print(sample, sample, "0", "0", "0", random.choice([1, 2]),
                      "-9", file=f, sep=" ")

        # Finally, prepare the arguments for the script.
        args = [
            "skat",
            "--impute2", os.path.join(out_directory, "impute2.txt"),
            "--sample", os.path.join(out_directory, "samples.txt"),
            "--pheno", os.path.join(out_directory, "phenotypes.txt"),
            "--snp-set", os.path.join(out_directory, "snp_sets.txt"),
            "--covar", "covar_discrete,covar_continuous",
            "--sample-column", "sample",
            "--gender-column", "None",
        ]

        return args
