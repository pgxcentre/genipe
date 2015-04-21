
# This file is part of gwip.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import os
import re
import sys
import logging
import argparse
from collections import defaultdict

import numpy as np

from .. import __version__
from ..formats.impute2 import *
from ..error import ProgramError


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


def main(args=None):
    """The main function.

    Args:
        args (argparse.Namespace): the arguments to be parsed (if
                                   :py:func:`main` is called by another
                                   modulel)

    """
    # Creating the option parser
    desc = ("Concatenate IMPUTE2 output files and retrieve some "
            "statistics. This script is part of the 'gwip' package, "
            "version {}.".format(__version__))
    parser = argparse.ArgumentParser(description=desc)

    # Files that need closing
    logging_fh = None

    try:
        # Parsing the options
        args = parse_args(parser, args)

        # Getting the output directory (dirname of the output prefix
        out_dir = os.path.dirname(args.prefix)

        # Adding the logging capability
        log_file = args.prefix + ".log"
        logging_fh = logging.FileHandler(log_file, mode="w")
        logging.basicConfig(
            format="[%(asctime)s %(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            level=logging.DEBUG if args.debug else logging.INFO,
            handlers=[logging.StreamHandler(), logging_fh]
        )
        logging.info("Logging everything into '{}'".format(log_file))
        logging.info("Program arguments: '{}'".format(" ".join(sys.argv[1:])))

        # Checking the options
        check_args(args)

        # Concatenate and extract information
        concatenate_files(args.impute2, args.prefix, args.chrom, args)

    # Catching the Ctrl^C
    except KeyboardInterrupt:
        logging.info("Cancelled by user")
        sys.exit(0)

    # Catching the ProgramError
    except ProgramError as e:
        logging.error(e)
        parser.error(e.message)

    except Exception as e:
        logging.error(e)
        raise

    finally:
        if logging_fh is not None:
            logging_fh.close()


def concatenate_files(i_filenames, out_prefix, real_chrom, options):
    """Concatenates and extracts information from IMPUTE2 GEN file(s).

    Args:
        i_filenames (list): the list of input filenames (to concatenate)
        out_prefix (str): the output prefix for the output files
        real_chrom (str): the chromosome contained in all the input files
        options (argparse.Namespace): the options

    This function will create the following seven files:

    +-----------------------+-------------------------------------------------+
    | File name             | Description                                     |
    +=======================+=================================================+
    | ``.impute2``          | Imputation results (merged from all the input   |
    |                       | files).                                         |
    +-----------------------+-------------------------------------------------+
    | ``.alleles``          | Description of the reference and alternative    |
    |                       | allele at each sites.                           |
    +-----------------------+-------------------------------------------------+
    | ``.imputed_sites``    | List of imputed sites (excluding sites that     |
    |                       | were previously genotyped in the study cohort). |
    +-----------------------+-------------------------------------------------+
    | ``.completion_rates`` | Number of missing values and completion rate    |
    |                       | for all sites (using the probability threshold  |
    |                       | set by the user, where the default is higher    |
    |                       | and equal to 0.9).                              |
    +-----------------------+-------------------------------------------------+
    | ``.good_sites``       | List of sites which pass the completion rate    |
    |                       | threshold (set by the user, where the default   |
    |                       | is higher and equal to 0.98) using the          |
    |                       | probability threshold (set by the user, where   |
    |                       | the default is higher and equal to 0.9).        |
    +-----------------------+-------------------------------------------------+
    | ``.map``              | A map file describing the genomic location of   |
    |                       | all sites.                                      |
    +-----------------------+-------------------------------------------------+
    | ``.maf``              | File containing the minor allele frequency      |
    |                       | (along with minor allele identification) for    |
    |                       | all sites using the probabilitty threshold of   |
    |                       | 0.9. When no genotypes are available (because   |
    |                       | they are all below the threshold), the MAF is   |
    |                       | ``NA``.                                         |
    +-----------------------+-------------------------------------------------+

    """
    # Opening output files
    impute2_o_file = open(out_prefix + ".impute2", "w")
    alleles_o_file = open(out_prefix + ".alleles", "w")
    imputed_sites_o_file = open(out_prefix + ".imputed_sites", "w")
    completion_o_file = open(out_prefix + ".completion_rates", "w")
    good_sites_o_file = open(out_prefix + ".good_sites", "w")
    map_o_file = open(out_prefix + ".map", "w")
    maf_o_file = open(out_prefix + ".maf", "w")

    # Printing the headers
    print("name", "nb_missing", "completion_rate", sep="\t",
          file=completion_o_file)
    print("name", "a1", "a2", sep="\t", file=alleles_o_file)
    print("name", "major", "minor", "maf", sep="\t", file=maf_o_file)

    # The markers that were already seen
    already_seen = defaultdict(int)

    # Opening the input file(s) one by one
    chr23_par_already_warned = False
    for i_filename in i_filenames:
        logging.info("Working with {}".format(i_filename))

        # Getting the expected number of lines from summary file
        summary = None
        with open(i_filename + "_summary", "r") as i_file:
            summary = i_file.read()
        r = re.search(r"-Output file\n --\d+ type 0 SNPs\n --\d+ type 1 SNPs"
                      r"\n --\d+ type 2 SNPs\n --\d+ type 3 SNPs\n"
                      r" --(\d+) total SNPs", summary)
        if r is None:
            raise ProgramError("{}: unknown "
                               "format".format(i_filename + "_summary"))
        nb_expected = int(r.group(1))
        logging.info("  - expecting {:,d} lines".format(nb_expected))

        # The input file
        i_file = open(i_filename, "r")

        nb_line = 0
        for line in i_file:
            # The number of line
            nb_line += 1

            # Splitting the line
            row = line.rstrip("\r\n").split(" ")

            # Gathering genotypes
            (chrom, name, pos, a1, a2), geno = matrix_from_line(row)

            # Checking the name of the marker
            if name == ".":
                name = "{}:{}".format(real_chrom, pos)

            # Have we seen this marker?
            if name in already_seen:
                new_name = "{}_{}".format(name, already_seen[name])
                while new_name in already_seen:
                    already_seen[name] += 1
                    new_name = "{}_{}".format(name, already_seen[name])
                name = new_name
            already_seen[name] += 1

            # Checking the chromosome
            if chrom == "---":
                # This site is imputed
                print(name, file=imputed_sites_o_file)

            elif chrom != real_chrom:
                # The chromosome in the impute2 file is not the same as
                # the one required in the options
                if (chrom == "23") and (real_chrom == "25"):
                    # This might be that we are in the PAR
                    # (pseudo-autosomal) region, so we print a warning
                    # and continue
                    if not chr23_par_already_warned:
                        logging.warning("WARNING: asked for chromosome 25, "
                                        "but in chromosome 23: be sure to be "
                                        "in the pseudo-autosomal region")
                        chr23_par_already_warned = True

                else:
                    # This is a problem, so we quit
                    raise ProgramError("{} != {}: not same "
                                       "chromosome".format(chrom, real_chrom))

            # Adding the information to the coding file
            print(name, a1, a2, sep="\t", file=alleles_o_file)

            # Checking the completion rate and saving it
            good_calls = get_good_probs(geno, options.probability)
            nb = np.sum(good_calls)
            comp = 0
            if geno.shape[0] != 0:
                comp = nb / geno.shape[0]
            print(name, geno.shape[0] - nb, comp, sep="\t",
                  file=completion_o_file)

            # Computing the MAF and saving it
            maf, minor, major = maf_from_probs(geno[good_calls], a1, a2)
            print(name, major, minor, maf, sep="\t", file=maf_o_file)

            if comp >= options.completion:
                # The completion is over the threshold
                print(name, file=good_sites_o_file)

            # Saving the map file
            print(real_chrom, name, "0", pos, sep="\t", file=map_o_file)

            # Saving the data
            print(real_chrom, name, pos, a1, a2, *row[5:], sep=" ",
                  file=impute2_o_file)

        # Closing the input file
        i_file.close()

    if nb_line != nb_expected:
        logging.warning("  - number of lines ({:,d}) is not as expected "
                        "({:,d})".format(nb_line, nb_expected))

    # Closing output files
    impute2_o_file.close()
    alleles_o_file.close()
    imputed_sites_o_file.close()
    completion_o_file.close()
    good_sites_o_file.close()
    map_o_file.close()
    maf_o_file.close()


def check_args(args):
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the options to verify

    Note
    ----
        If there is a problem, a :py:class:`gwip.error.ProgramError` is raised.

    """
    # Checking the input files
    for filename in args.impute2:
        if not os.path.isfile(filename):
            raise ProgramError("{}: no such file".format(filename))

        summary_file = filename + "_summary"
        if not os.path.isfile(summary_file):
            raise ProgramError("{}: no such file".format(summary_file))

    # Checking the chromosome
    valid_chromosome = [str(i) for i in range(1, 24)]
    valid_chromosome.append("25")
    if args.chrom not in valid_chromosome:
        raise ProgramError("{}: invalid chromosome".format(args.chrom))
    if args.chrom == "23":
        logging.warning("MAF computation is wrong for chromosome 23 (males "
                        "have 2 alleles in the computation, instead of 1)...")

    return True


def parse_args(parser, args=None):
    """Parses the command line options and arguments.

    Args:
        parser (argparse.ArgumentParser): the argument parser
        args (list): the list of arguments (if not taken from ``sys.argv``)

    Returns:
        argparse.Namespace: the list of options and arguments

    Note
    ----
        The only check that is done here is by the parser itself. Values are
        verified later by the :py:func:`check_args` function.

    """
    # The parser object
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s (part of GWIP version {})".format(__version__),
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="set the logging level to debug",
    )

    # The input files
    group = parser.add_argument_group("Input Files")
    group.add_argument(
        "-i",
        "--impute2",
        type=str,
        metavar="GEN",
        required=True,
        nargs="+",
        help="IMPUTE2 file(s) to merge.",
    )

    # The options
    group = parser.add_argument_group("Options")
    group.add_argument(
        "--chr",
        type=str,
        metavar="CHR",
        required=True,
        dest="chrom",
        help="The chromosome on witch the imputation was made.",
    )
    group.add_argument(
        "--probability",
        type=float,
        metavar="FLOAT",
        default=0.9,
        help="The probability threshold for no calls. [<%(default).1f]",
    )
    group.add_argument(
        "--completion",
        type=float,
        metavar="FLOAT",
        default=0.98,
        help="The completion rate threshold for site exclusion. "
             "[<%(default).2f]",
    )

    # The output files
    group = parser.add_argument_group("Output Files")
    group.add_argument(
        "--prefix",
        type=str,
        metavar="FILE",
        default="imputed",
        help="The prefix for the output files. [%(default)s]",
    )

    if args is not None:
        return parser.parse_args(args)

    return parser.parse_args()


# Calling the main, if necessary
if __name__ == "__main__":
    main()
