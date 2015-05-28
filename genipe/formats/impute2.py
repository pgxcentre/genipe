
# This file is part of genipe.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import numpy as np

from ..error import ProgramError


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["matrix_from_line", "get_good_probs", "maf_from_probs",
           "dosage_from_probs", "hard_calls_from_probs"]


def matrix_from_line(impute2_line):
    """Generates the probability matrix from an IMPUTE2 line.

    Args:
        impute2_line (list): a single line from IMPUTE2's result (split by
                             space)

    Returns:
        tuple: a tuple containing the marker's information (first five values
               of the line) and the matrix probability (numpy array, float)

    The shape of the matrix is n x 3 where n is the number of samples.
    The columns represent the probability for AA, AB and BB.

    Note
    ----
        The ``impute2_line`` variable is a list of str, corresponding to a line
        from the IMPUTE2's result, split by space.

    """
    # Creating the array and changing it's shape
    probabilities = np.array(impute2_line[5:], dtype=float)
    probabilities.shape = (len(probabilities) // 3, 3)
    return impute2_line[:5], probabilities


def get_good_probs(prob_matrix, min_prob=0.9):
    """Gathers good imputed genotypes (>= probability threshold).

    Args:
        prob_matrix (numpy.array): the probability matrix
        min_prob (float): the probability threshold

    Returns:
        numpy.array: a mask array containing the positions where the
                     probabilities are equal or higher to the threshold

    """
    return np.amax(prob_matrix, axis=1) >= min_prob


def maf_from_probs(prob_matrix, a1, a2, gender=None, site_name=None):
    """Computes MAF from a probability matrix (and gender if chromosome X).

    Args:
        prob_matrix (numpy.array): the probability matrix
        a1 (str): the first allele
        a2 (str): the second allele
        gender (numpy.array): the gender of the samples
        site_name (str): the name for this site

    Returns:
        tuple: a tuple containing three values: the minor allele frequency, the
               minor and the major allele.

    When 'gender' is not None, we assume that the MAF on chromosome X is
    required (hence, males count as 1, and females as 2 alleles). There is also
    an Exception raised if there are any heterozygous males.

    """
    # By default, the MAF is NA, and a1=major, a2=minor
    maf = "NA"
    major, minor = a1, a2

    # If there are no data, we return default values
    if prob_matrix.shape[0] == 0:
        return maf, minor, major

    if gender is None:
        # Not checking gender (this isn't chromosome X)
        nb_geno = np.bincount(np.argmax(prob_matrix, axis=1), minlength=3)
        maf = ((nb_geno[2] * 2) + nb_geno[1]) / (np.sum(nb_geno) * 2)

    else:
        # Getting the males and females
        males = (gender == 1)
        females = (gender == 2)

        # Male counts
        males_nb_geno = np.bincount(np.argmax(prob_matrix[males], axis=1),
                                    minlength=3)

        # Female counts
        females_nb_geno = np.bincount(np.argmax(prob_matrix[females], axis=1),
                                      minlength=3)

        # The total number of genotypes
        total_geno_males = np.sum(males_nb_geno)
        total_geno_females = np.sum(females_nb_geno)

        # If there are no genotypes, we return default values
        if (total_geno_males + total_geno_females) == 0:
            return maf, minor, major

        # There shouldn't be heterozygous genotypes for males
        if males_nb_geno[1] > 0:
            raise ProgramError("{}: heterozygous male "
                               "present".format(site_name))

        # Computing the frequencies
        maf = males_nb_geno[2] + (females_nb_geno[2] * 2) + females_nb_geno[1]
        maf /= (total_geno_males + (total_geno_females * 2))

    # Is this the MAF?
    if maf != "NA" and maf > 0.5:
        minor, major = a1, a2
        maf = 1 - maf

    return maf, minor, major


def dosage_from_probs(homo_probs, hetero_probs, scale=2):
    """Computes dosage from probability matrix (for the minor allele).

    Args:
        homo_probs (numpy.array): the probabilities for the homozygous genotype
        hetero_probs (numpy.array): the probabilities for the heterozygous
                                    genotype
        scale (int): the scale value

    Returns:
        numpy.array: the dosage computed from the probabilities

    """
    return (homo_probs + (hetero_probs / 2)) * scale


def hard_calls_from_probs(a1, a2, probs):
    """Computes hard calls from probability matrix.

    Arg:
        a1 (str): the first allele
        a2 (str): the second allele
        probs (numpy.array): the probability matrix

    Returns:
        numpy.array: the hard calls computed from the probabilities

    """
    possible_geno = np.array([
        " ".join([a1] * 2),     # Homo A1
        " ".join([a1, a2]),     # Hetero
        " ".join([a2] * 2),     # Homo A2
    ])

    return possible_geno[np.argmax(probs, axis=1)]
