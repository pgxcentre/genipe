
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

from ..formats.index import *
from ..error import ProgramError
from .. import __version__, chromosomes


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


def main(args=None):
    """The main function."""
    # Creating the option parser
    desc = ("Extract imputed markers located in a specific genomic region. "
            "This script is part of the 'gwip' package, version "
            "{}).".format(__version__))
    parser = argparse.ArgumentParser(description=desc)

    # Files that need closing
    logging_fh = None

    try:
        # Parsing the options
        args = parse_args(parser, args)

        # Getting the output directory (dirname of the output prefix
        out_dir = os.path.dirname(args.out)

        # Adding the logging capability
        log_file = args.out + ".log"
        logging_fh = logging.FileHandler(log_file, mode="w")
        logging.basicConfig(
            format="[%(asctime)s %(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            level=logging.DEBUG if args.debug else logging.INFO,
            handlers=[logging.StreamHandler(), logging_fh]
        )
        logging.info("Logging everything into '{}'".format(log_file))

        # Checking the options
        check_args(args)

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


def check_args(args):
    """Checks the arguments and options."""
    # Checking that the impute2 files exists
    if not os.path.isfile(args.impute2):
        raise ProgramError("{}: no such file".format(args.impute2))

    # Is there something to extract?
    if not args.genomic and not args.maf and not args.rate:
        raise ProgramError("nothing to extreact: use '--genomic', '--maf' or "
                           "'--rate'")

    # If genomic, we check the format
    if args.genomic is not None:
        genomic_match = re.match(r"(.+):(\d+)-(\d+)$", args.genomic)
        if not genomic_match:
            raise ProgramError("{}: no a valid genomic "
                               "region".format(args.genomic))
        chrom = genomic_match.group(1).replace("chr", "")
        start = int(genomic_match.group(2))
        end = int(genomic_match.group(3))

        if chrom not in map(str, chromosomes):
            raise ProgramError("{}: invalid chromosome".format(chrom))

        if end < start:
            start, end = end, start

        args.genomic = (chrom, start, end)

    # If MAF, we check what's required
    if args.maf is not None:
        if args.maf < 0 or args.maf > 0.5:
            raise ProgramError("{}: invalid MAF".format(args.maf))
        if args.maf_file is None:
            raise ProgramError("needs '--maf-file' when using '--maf'")
        if args.map_file is None:
            raise ProgramError("needs '--map-file' when using '--maf'")

    # If completion rate, we check what's required
    if args.rate is not None:
        if args.rate < 0 or args.rate > 1:
            raise ProgramError("{}: invalid rate".format(args.rate))
        if args.rate_file is None:
            raise ProgramError("needs '--rate-file' when using '--rate'")
        if args.map_file is None:
            raise ProgramError("needs '--map-file' when using '--rate'")

    # Checking the other files
    for filename in [args.maf_file, args.rate_file, args.map_file]:
        if filename and not os.path.isfile(filename):
            raise ProgramError("{}: no such file".format(filename))

    return True


def parse_args(parser, args=None):
    """Parses the command line options and arguments."""
    # Adding the version option for the main parser
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s, part of gwip version {}".format(__version__),
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="set the logging level to debug",
    )

    # The input files
    group = parser.add_argument_group("Input Files")
    group.add_argument(
        "--impute2",
        type=str,
        metavar="FILE",
        required=True,
        help="The output from IMPUTE2.",
    )
    group.add_argument(
        "--maf-file",
        type=str,
        metavar="FILE",
        help="The file containing the MAF",
    )
    group.add_argument(
        "--rate-file",
        type=str,
        metavar="FILE",
        help="The file containing the completion rates.",
    )
    group.add_argument(
        "--map-file",
        type=str,
        metavar="FILE",
        help="The file containing the mapping information.",
    )

    # The output files
    group = parser.add_argument_group("Output Options")
    group.add_argument(
        "--out",
        metavar="PREFIX",
        default="impute2_extractor",
        help="The prefix of the output files. [%(default)s]",
    )
    group.add_argument(
        "--format",
        metavar="FORMAT",
        nargs="+",
        default="probs",
        help="The output format. Can specify either 'probs' for probabilities "
             "(same as impute2 format, i.e. 3 values per sample), 'dosage' "
             "for dosage values (one value between 0 and 2 by sample), or "
             "'calls' for hard calls. [%(default)s]",
    )

    # What to extract
    group = parser.add_argument_group("Extraction Options")
    group.add_argument(
        "--genomic",
        type=str,
        metavar="CHR:START-END",
        help="The range to extract (e.g. 22 1000000 1500000).",
    )
    group.add_argument(
        "--maf",
        type=float,
        metavar="FLOAT",
        help="Extract markers with a minor allele frequency equal or higher "
             "than the specified threshold. Requires the two options "
             "'--maf-file' and '--map-file'.",
    )
    group.add_argument(
        "--rate",
        type=float,
        metavar="FLOAT",
        help="Extract markers with a completion rate equal or higher to the "
             "specified threshold. Requires the two options '--rate-file' and "
             "'--map-file'.",
    )

    if args is not None:
        return parser.parse_args(args)

    return parser.parse_args()


# Calling the main, if necessary
if __name__ == "__main__":
    main()
