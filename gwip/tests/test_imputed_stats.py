
# This file is part of gwip.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import unittest
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd

from ..tools.imputed_stats import *
from ..tools.imputed_stats import _get_result_from_linear_logistic


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


class TestImputedStats(unittest.TestCase):

    def setUp(self):
        """Setup the tests."""
        # Creating the temporary directory
        self.output_dir = TemporaryDirectory(prefix="gwip_test_")

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
            "sample_id\tTime_To_Event\tCensure\tC1\tC2\tC3\tC4\tGender\t"
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
        args.censure = "Censure"
        args.interaction = None
        args.chrx = False
        args.gender_column = "Gender"

        # The expected value
        expected_shape = (3, 6)
        expected_columns = {"C1", "C2", "C3", "Gender", "Time_To_Event",
                            "Censure"}
        expected_index = ["sample_1", "sample_2", "sample_3"]
        expected_c1 = np.array([0.3, 0.9, 0.4], dtype=float)
        expected_c2 = np.array([0.2, 0.1, 0.67], dtype=float)
        expected_c3 = np.array([0.45, 0.42, 999], dtype=float)
        expected_gender = np.array([1, 2, 1], dtype=int)
        expected_tte = np.array([10, 2, 8], dtype=int)
        expected_censure = np.array([0, 1, 1], dtype=int)
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
        self.assertTrue((expected_censure == observed_p.Censure.values).all())
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
            (expected_censure[:-1] == observed_p.Censure.values).all()
        )
        self.assertEqual(expected_remove_g, observed_remove_g)

        # Changing from Cox to linear should remove a column
        args.missing_value = None
        args.analysis_type = "linear"
        del args.tte
        del args.censure
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
        self.assertTrue((expected_pheno == observed_p.Pheno_Logit.values).all())
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
        self.assertTrue((expected_pheno == observed_p.Pheno_Logit.values).all())
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
        self.assertTrue((expected_pheno == observed_p.Pheno_Logit.values).all())
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
        self.assertTrue((expected_pheno == observed_p.Pheno_Logit.values).all())
        self.assertTrue(np.allclose(expected_inter, observed_p.Inter.values))
        self.assertEqual(expected_remove_g, observed_remove_g)

        # Sample shouldn't be included if chrx is False, but Gender in covars
        args.covar.append("Gender")
        args.chrx=False

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
        self.assertTrue((expected_gender[:-1] == observed_p.Gender.values).all())
        self.assertTrue((expected_pheno[:-1] == observed_p.Pheno_Logit.values).all())
        self.assertTrue(np.allclose(expected_inter[:-1], observed_p.Inter.values))
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
        self.assertTrue((expected_pheno == observed_p.Pheno_Logit.values).all())
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
        self.assertTrue((expected_pheno == observed_p.Pheno_Logit.values).all())
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
        observed = get_formula("pheno", [], None)
        self.assertEqual(expected, observed)

        # Testing with one covar, no interaction
        expected = "pheno ~ _GenoD + C1"
        observed = get_formula("pheno", ["C1"], None)
        self.assertEqual(expected, observed)

        # Testing with more than one covar, no interaction
        expected = "pheno ~ _GenoD + C1 + C2 + C3"
        observed = get_formula("pheno", ["C1", "C2", "C3"], None)
        self.assertEqual(expected, observed)

        # Testing with without covar, but with interaction
        expected = "pheno ~ _GenoD + _GenoD*inter"
        observed = get_formula("pheno", [], "inter")
        self.assertEqual(expected, observed)

        # Testing with one covar and interaction
        expected = "pheno ~ _GenoD + inter + _GenoD*inter"
        observed = get_formula("pheno", ["inter"], "inter")
        self.assertEqual(expected, observed)

        # Testing with more than one covar and interaction
        expected = "pheno ~ _GenoD + C1 + C2 + C3 + inter + _GenoD*inter"
        observed = get_formula("pheno", ["C1", "C2", "C3", "inter"], "inter")
        self.assertEqual(expected, observed)

    @unittest.skip("Test not implemented")
    def test_fit_cox(self):
        """Tests the 'fit_cox' function."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_fit_linear(self):
        """Tests the 'fit_linear' function."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_fit_logistic(self):
        """Tests the 'fit_logistic' function."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_get_result_from_linear_logistic(self):
        """Tests the '_get_result_from_linear_logistic' function."""
        self.fail("Test not implemented")

    @unittest.skip("Test not implemented")
    def test_check_args(self):
        """Tests the 'check_args' function."""
        self.fail("Test not implemented")
