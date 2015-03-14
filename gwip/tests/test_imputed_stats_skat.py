
# This file is part of gwip.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import random
import unittest
from tempfile import TemporaryDirectory

import numpy as np
import pandas as pd
from pkg_resources import resource_filename

from .test_imputed_stats import reverse_dosage
from ..tools.imputed_stats import *


__author__ = "Marc-Andre Legault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["TestImputedStatsSkat"]


temp_directories = []


class TestImputedStatsSkat(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        global temp_directories

        self.output_dir = TemporaryDirectory(prefix="gwip_test_")
        temp_directories.append(self.output_dir)
        self.args = self.setup_skat_files()
        super().__init__(*args, **kwargs)

    @classmethod
    def tearDownClass(cls):
        """Finishes the test."""
        # Deleting the output directory
        for directory in temp_directories:
            directory.cleanup()

    def test_continuous(self):
        args = self.args + [
            "--pheno-name", "outcome_continuous",
            "--outcome-type", "continuous"
        ]
        main(args=args)
        
        results_filename = os.path.join(
            self.output_dir.name,
            "skat_test.skat.dosage"
        )

        p = pd.read_csv(results_filename, header=0, sep="\t")["p_value"][0]
        self.assertAlmostEqual(np.log10(0.002877041), np.log10(p), places=10)

    def test_discrete(self):
        args = self.args + [
            "--pheno-name", "outcome_discrete",
            "--outcome-type", "discrete"
        ]
        main(args=args)
        
        results_filename = os.path.join(
            self.output_dir.name,
            "skat_test.skat.dosage"
        )

        p = pd.read_csv(results_filename, header=0, sep="\t")["p_value"][0]
        self.assertAlmostEqual(np.log10(0.1401991), np.log10(p), places=10)

    def setup_skat_files(self):
        """Parses the SKAT example files into the format expected by gwip."""
        out_directory = self.output_dir.name

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
                a1, a2 = [random.choice("ATGC") for i in range(2)]

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
        filename = os.path.join(out_directory, "snp_set.txt")
        with open(filename, "w") as output_file:
            print("variant", "snp_set", sep="\t", file=output_file)
            for variant in variants:
                print(variant, "set1", sep="\t", file=output_file)

        # Create a list of samples.
        samples = ["sample_{}".format(i + 1) for i in range(skat_z.shape[0])]

        # Get the filename for the phenotype file.
        filename = os.path.join(out_directory, "phenotypes.txt")

        # Before we can write it, we need to add a samples column.
        skat_x["sample"] = samples
        skat_x = skat_x[["sample", "covar_discrete", "covar_continuous"]]

        # We also add the outcomes to this (gwip uses a single file for both
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
            "--out", os.path.join(out_directory, "skat_test"),
            "--snp-set", os.path.join(out_directory, "snp_set.txt"),
            "--covar", "covar_discrete,covar_continuous",
            "--sample-column", "sample",
            "--gender-column", "None",
        ]

        return args
