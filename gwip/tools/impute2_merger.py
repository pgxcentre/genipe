
# This file is part of gwip.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


import os
import re
import sys
import logging
import argparse
from collections import defaultdict

import numpy as np

from .. import __version__
from ..error import ProgramError


def main():
    """The main function."""
    # Creating the option parser
    desc = ("Concatenate IMPUTE2 files and retrieve some statistics "
            "(gwip version {}).".format(__version__))
    parser = argparse.ArgumentParser(description=desc)

    try:
        # Parsing the options
        args = parse_args(parser)

        # Getting the output directory (dirname of the output prefix
        out_dir = os.path.dirname(args.prefix)

        # Adding the logging capability
        log_file = args.prefix + ".log"
        logging.basicConfig(
            format="[%(asctime)s %(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            level=logging.DEBUG if args.debug else logging.INFO,
            handlers=[logging.StreamHandler(),
                      logging.FileHandler(log_file, mode="w")]
        )
        logging.info("Logging everything into '{}'".format(log_file))

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


def concatenate_files(i_filenames, out_prefix, real_chrom, options):
    """Concatenates and extracts information from IMPUTE2 GEN file(s)."""
    # Opening output files
    impute2_o_file = open(out_prefix + ".impute2", "w")
    alleles_o_file = open(out_prefix + ".alleles", "w")
    imputed_sites_o_file = open(out_prefix + ".imputed_sites", "w")
    completion_o_file = open(out_prefix + ".completion_rates", "w")
    good_sites_o_file = open(out_prefix + ".good_sites", "w")
    map_o_file = open(out_prefix + ".map", "w")

    # Printing the headers
    print("name", "nb_missing", "completion_rate", sep="\t",
          file=completion_o_file)
    print("name", "a1", "a2", sep="\t", file=alleles_o_file)

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

        already_seen = defaultdict(int)
        nb_line = 0
        for line in i_file:
            # The number of line
            nb_line += 1

            # Splitting the line
            row = line.rstrip("\r\n").split(" ")

            # Gathering site information
            chrom, name, pos, a1, a2 = row[:5]

            # Gathering genotypes
            geno = np.array(row[5:], dtype=float)
            geno.shape = (len(geno) // 3, 3)

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
            max_probs = np.amax(geno, axis=1)
            nb = np.sum(max_probs >= options.probability)
            comp = nb / geno.shape[0]
            print(name, geno.shape[0] - nb, comp, sep="\t",
                  file=completion_o_file)

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


def check_args(args):
    """Checks the arguments and options."""
    # Checking the input files
    for filename in args.impute2:
        if not os.path.isfile(filename):
            raise ProgramError("{}: no such file".format(filename))

        summary_file = filename + "_summary"
        if not os.path.isfile(summary_file):
            raise ProgramError("{}: no such file".format(summary_file))

    return True


def parse_args(parser):
    """Parses the command line options and arguments."""
    # The parser object
    parser.add_argument("--version", action="version",
                        version="%(prog)s (part of GWIP "
                                "version {})".format(__version__))
    parser.add_argument("--debug", action="store_true",
                        help="Set the logging level to debug")

    # The input files
    group = parser.add_argument_group("Input Files")
    group.add_argument("-i", "--impute2", type=str, metavar="GEN",
                       required=True, nargs="+", help="IMPUTE2 file(s)")

    # The options
    group = parser.add_argument_group("Options")
    group.add_argument("--chr", type=str, metavar="CHR", required=True,
                       dest="chrom", help=("The chromosome on witch the "
                                           "imputation was made"))
    group.add_argument("--probability", type=float, metavar="FLOAT",
                       default=0.9, help=("The probability threshold for no "
                                          "calls [%(default).1f]"))
    group.add_argument("--completion", type=float, metavar="FLOAT",
                       default=0.98, help=("The site completion rate "
                                           "threshold [%(default).2f]"))

    # The output files
    group = parser.add_argument_group("Output Files")
    group.add_argument("--prefix", type=str, metavar="FILE", default="imputed",
                       help="The prefix for the output files [%(default)s]")

    return parser.parse_args()


# Calling the main, if necessary
if __name__ == "__main__":
    main()
