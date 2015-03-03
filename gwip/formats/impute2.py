
# This file is part of gwip.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import logging

import numpy as np

from ..error import ProgramError


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["matrix_from_line", "get_good_probs", "maf_from_probs",
           "dosage_from_probs"]


def matrix_from_line(impute2_line):
    """Generates the probability matrix from an IMPUTE2 line."""
    # Creating the array
    probabilities = np.array(impute2_line[5:], dtype=float)

    # Changing it's shape
    try:
        probabilities.shape = (len(probabilities) // 3, 3)
    except ValueError:
        raise ProgramError("{}: invalid number of "
                           "entries".format(len(probabilities)))

    return impute2_line[:5], probabilities


def get_good_probs(prob_matrix, min_prob=0.9):
    """Gathers good imputed genotypes (>= probability threshold)."""
    return np.amax(prob_matrix, axis=1) >= min_prob


def maf_from_probs(prob_matrix, a1, a2, gender=None, site_name=None):
    """Computes MAF from a probability matrix (and gender if chromosome X)."""
    # By default, the MAF is NA, and a1=major, a2=minor
    maf = "NA"
    major, minor = a1, a2

    if gender is None:
        # Not checking gender (this isn't chromosome X)
        if prob_matrix.shape[0] != 0:
            nb_geno = np.bincount(np.argmax(prob_matrix, axis=1), minlength=3)
            maf = ((nb_geno[2] * 2) + nb_geno[1]) / (np.sum(nb_geno) * 2)

    else:
        # Getting the males
        males = (gender == 1)

        # Male counts
        males_nb_geno = np.bincount(np.argmax(prob_matrix[males], axis=1),
                                    minlength=3)

        # Female counts
        females_nb_geno = np.bincount(np.argmax(prob_matrix[~males], axis=1),
                                      minlength=3)

        # There shouldn't be heterozygous genotypes for males
        if males_nb_geno[1] > 0:
            raise ProgramError("{}: heterozygous male "
                               "present".format(site_name))

        # Computing the frequencies
        maf = males_nb_geno[2] + (females_nb_geno[2] * 2) + females_nb_geno[1]
        maf /= (np.sum(males_nb_geno) + (np.sum(females_nb_geno) * 2))

    # Is this the MAF?
    if maf != "NA" and maf > 0.5:
        minor, major = a1, a2
        maf = 1 - maf

    return maf, minor, major


def dosage_from_probs(homo_probs, hetero_probs, scale=2):
    """Computes dosage from probability matrix (for the minor allele)."""
    return (homo_probs + (hetero_probs / 2)) * scale
