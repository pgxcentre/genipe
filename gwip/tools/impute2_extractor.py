
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

import pandas as pd

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

        # Retrieving the file index
        impute2_index = get_impute2_index(args.impute2)

        # Gathering what needs to be extracted
        pos_to_extract = gather_extraction_positions(
            maf=args.maf,
            rate=args.rate,
            maf_filename=args.maf_file,
            rate_filename=args.rate_file,
            map_filename=args.map_file,
            extract_filename=args.extract,
        )
        logging.info("Keeping a total of {:,d} "
                     "markers".format(len(pos_to_extract)))

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


def get_impute2_index(i_filename):
    """Gets the impute2 file index."""
    if not os.path.isfile(get_index_fn(i_filename)):
        # The index does not exists, so we create it
        logging.info("Creating index for '{}'".format(i_filename))
        return build_index(
            i_filename,
            chrom_col=0,
            pos_col=2,
            delimiter=' ',
        )
    else:
        logging.info("Retrieving index for '{}'".format(i_filename))
        return get_index(i_filename)


def gather_extraction_positions(maf, maf_filename, rate, rate_filename,
                                map_filename, extract_filename):
    "Gather positions that are required."""
    map_data = None
    markers_to_keep = None

    # Reading the map data
    if map_filename is not None:
        logging.info("Reading the MAP data")
        map_data = pd.read_csv(map_filename, sep="\t",
                               names=["chrom", "name", "cm", "pos"])
        map_data = map_data.set_index("name", verify_integrity=True)
        markers_to_keep = map_data.index

        logging.info("  - read {:,d} markers".format(len(markers_to_keep)))

    # Extracting using names?
    if extract_filename is not None:
        logging.info("Reading the markers to extract")
        with open(extract_filename, "r") as i_file:
            markers_to_keep = markers_to_keep.intersection({
                i for i in i_file.read().splitlines()
            })
        logging.info("  - kept {:,d} markers".format(len(markers_to_keep)))

    # Extracting using MAF?
    if maf is not None:
        logging.info("Reading the MAF data")

        maf_data = pd.read_csv(maf_filename, sep="\t")
        maf_data = maf_data.set_index("name", verify_integrity=True)
        maf_data = maf_data[maf_data.maf >= maf]

        markers_to_keep = markers_to_keep.intersection(maf_data.index)
        logging.info("  - kept {:,d} markers".format(len(markers_to_keep)))

    # Extracting using completion rate?
    if rate is not None:
        logging.info("Reading the completion rate data")

        rate_data = pd.read_csv(rate_filename, sep="\t")
        rate_data = rate_data.set_index("name", verify_integrity=True)
        rate_data = rate_data[rate_data.completion_rate >= rate]

        markers_to_keep = markers_to_keep.intersection(rate_data.index)
        logging.info("  - kept {:,d} markers".format(len(markers_to_keep)))

    # Needed to extract something?
    if markers_to_keep is None:
        return {}

    # Keeping only what is required
    map_data = map_data.loc[markers_to_keep, :].sort(["chrom", "pos"])

    # Extracting the positions
    return [
        tuple(i) for i in map_data[["chrom", "pos"]].values
    ]


def check_args(args):
    """Checks the arguments and options."""
    # Checking that the impute2 files exists
    if not os.path.isfile(args.impute2):
        raise ProgramError("{}: no such file".format(args.impute2))

    # Is there something to extract?
    if not args.genomic and not args.maf and not args.rate:
        if args.extract is None:
            raise ProgramError("nothing to extract: use '--extract', "
                               "'--genomic', '--maf' or '--rate'")

    elif args.extract is not None:
        raise ProgramError("'--extract' can only be used alone")

    # If extract, check the file
    if args.extract is not None:
        if args.map_file is None:
            raise ProgramError("needs '--map-file' when using '--extract'")
        if not os.path.isfile(args.extract):
            raise ProgramError("{}: no such file".format(args.extract))

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
        "--extract",
        type=str,
        metavar="FILE",
        help="File containing marker names to extract. Requires the option "
             "'--map-file'.",
    )
    group.add_argument(
        "--genomic",
        type=str,
        metavar="CHR:START-END",
        help="The range to extract (e.g. 22 1000000 1500000). Can be used in "
             "combination with '--rate' and '--maf'.",
    )
    group.add_argument(
        "--maf",
        type=float,
        metavar="FLOAT",
        help="Extract markers with a minor allele frequency equal or higher "
             "than the specified threshold. Requires the two options "
             "'--maf-file' and '--map-file'. Can be used in combination with "
             "'--rate' and '--genomic'.",
    )
    group.add_argument(
        "--rate",
        type=float,
        metavar="FLOAT",
        help="Extract markers with a completion rate equal or higher to the "
             "specified threshold. Requires the two options '--rate-file' and "
             "'--map-file'. Can be used in combination with '--maf' and "
             "'--genomic'.",
    )

    if args is not None:
        return parser.parse_args(args)

    return parser.parse_args()


# Calling the main, if necessary
if __name__ == "__main__":
    main()
