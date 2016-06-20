
# This file is part of genipe.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import os
import re
import sys
import shlex
import shutil
import logging
import argparse
from collections import namedtuple

import pandas as pd
from numpy import nan

from ..formats import index
from ..formats import impute2
from ..error import GenipeError
from .. import __version__, chromosomes

# Check if pyplink is installed
try:
    from pyplink import PyPlink
    HAS_PYPLINK = True
except ImportError:
    HAS_PYPLINK = False


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
    desc = ("Extract imputed markers located in a specific genomic region. "
            "This script is part of the 'genipe' package, version "
            "{}.".format(__version__))
    parser = argparse.ArgumentParser(description=desc)

    # Files that need closing
    logging_fh = None

    try:
        # Parsing the options
        args = parse_args(parser, args)

        # The logging handlers
        handlers = [logging.StreamHandler()]
        if not args.index_only:
            log_file = args.out + ".log"
            logging_fh = logging.FileHandler(log_file, mode="w")
            handlers.append(logging_fh)

        # Adding the logging capability
        logging.basicConfig(
            format="[%(asctime)s %(levelname)s] %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
            level=logging.DEBUG if args.debug else logging.INFO,
            handlers=handlers,
        )

        # First log
        if not args.index_only:
            logging.info("Logging everything into '{}'".format(log_file))
        logging.info("Program arguments: {}".format(
            " ".join(shlex.quote(part) for part in sys.argv[1:])
        ))

        # Checking the options
        check_args(args)

        if args.index_only:
            return index_file(args.impute2)

        # Gathering what needs to be extracted
        to_extract = gather_extraction(
            fn=args.impute2,
            maf=args.maf,
            rate=args.rate,
            info=args.info,
            extract_filename=args.extract,
            genomic_range=args.genomic,
        )

        # Extraction
        extract_markers(
            fn=args.impute2,
            to_extract=to_extract,
            out_prefix=args.out,
            out_format=args.out_format,
            prob_t=args.prob,
            is_long=args.long_format,
        )

    # Catching the Ctrl^C
    except KeyboardInterrupt:
        logging.info("Cancelled by user")
        sys.exit(0)

    # Catching the GenipeError
    except GenipeError as e:
        logging.error(e)
        parser.error(e.message)

    except Exception as e:
        logging.error(e)
        raise

    finally:
        if logging_fh is not None:
            logging_fh.close()


def index_file(fn):
    """Indexes the impute2 file.

    Args:
        fn (str): the name of the impute2 file

    This function uses the :py:func:`genipe.formats.index.get_index` to create
    the index file if it's missing.

    Note
    ----
        We won't catch the :py:class:`genipe.error.GenipeError` exception if
        it's raised, since the message will be relevant to the user.

    """
    # For each input file
    index.get_index(fn, cols=[0, 1, 2], names=["chrom", "name", "pos"],
                    sep=" ")


def extract_markers(fn, to_extract, out_prefix, out_format, prob_t, is_long):
    """Extracts according to names.

    Args:
        fn (str): the name of the input file
        to_extract (set): the list of markers to extract for each input file
        out_prefix (str): the output prefix
        out_format (list): the output format(s)
        prob_t (float): the probability threshold
        is_long (bool): True if format needs to be long

    """
    # The output files (probabilities)
    o_files = {
        suffix: open(out_prefix + "." + suffix, "w")
        for suffix in out_format if suffix not in {"bed"}
    }

    # If there is the 'bed' format, we actually need pyplink
    if "bed" in out_format:
        o_files["bed"] = (
            PyPlink(out_prefix, "w"),
            open(out_prefix + ".bim", "w"),
        )

    # Creating a fam (if bed)
    samples = get_samples(get_file_prefix(fn) + ".sample")
    sample_names = ["{}/{}".format(id_1, id_2) for id_1, id_2
                    in zip(samples.ID_1, samples.ID_2)]

    # Writing the header (for dosage)
    if "dosage" in o_files:
        if is_long:
            print("fid", "iid", "chrom", "pos", "name", "minor", "major",
                  "dosage", sep="\t", file=o_files["dosage"])
        else:
            print("chrom", "pos", "name", "minor", "major",
                  *sample_names, sep="\t", file=o_files["dosage"])

    # Writing the header (for calls)
    if "calls" in o_files:
        if is_long:
            print("fid", "iid", "chrom", "name", "cm", "pos", "call",
                  sep="\t", file=o_files["calls"])
        else:
            print("chrom", "name", "cm", "pos",
                  *sample_names, sep="\t", file=o_files["calls"])

    # Extracted positions
    all_extracted = set()

    # Reading the impute2 file
    extracted = set()

    # Finding the name of the file containing the index
    file_index = index.get_index(fn, cols=[0, 1, 2],
                                 names=["chrom", "name", "pos"], sep=" ")

    # Keeping only required values from the index
    file_index = file_index[file_index.name.isin(to_extract)]

    # Getting all the markers value
    logging.info("Extracting {:,d} markers".format(len(file_index)))
    with index.get_open_func(fn)(fn, "r") as i_file:
        for seek_value in file_index.seek.values:
            # Seeking
            i_file.seek(int(seek_value))

            # Reading the line
            line = i_file.readline()
            row = line.rstrip("\n").split(" ")

            # The marker name
            name = row[1]

            # Printing the data
            print_data(o_files, prob_t, samples.ID_1, samples.ID_2, line=line,
                       row=row, is_long=is_long)

            # Saving statistics
            extracted.add(name)

    logging.info("Extracted {:,d} markers".format(len(extracted)))
    if len(to_extract - extracted) > 0:
        logging.warning("Missing {:,d} "
                        "markers".format(len(to_extract - extracted)))

    # Keeping track of what has been extracted
    all_extracted |= extracted

    # Extracting the companion files (if impute2 and files are present)
    if "impute2" in o_files:
        extract_companion_files(
            i_prefix=get_file_prefix(fn),
            to_extract=to_extract,
            o_prefix=out_prefix,
        )

    # Writing the FAM file if bed
    if "bed" in o_files:
        cols = ["ID_1", "ID_2", "father", "mother", "sex", "plink_pheno"]
        samples[cols].to_csv(out_prefix + ".fam", sep=" ", index=False,
                             header=False)

    # Closing the files
    for o_format, o_file in o_files.items():
        if o_format == "bed":
            o_file[0].close()
            o_file[1].close()
        else:
            o_file.close()

    # Extraction complete
    logging.info("Extraction of {:,d} markers "
                 "completed".format(len(all_extracted)))


def get_samples(fn):
    """Reads the sample files, and extract the information.

    Args:
        fn (str): the name of the sample file

    Returns:
        pandas.DataFrame: the sample information

    """
    # Reading the file
    sample = pd.read_csv(fn, sep=" ")

    # Removing the first data row
    sample = sample.iloc[1:, ].reset_index(drop=True)

    return sample


def extract_companion_files(i_prefix, o_prefix, to_extract):
    """Extract markers from companion files (if they exists).

    Args:
        i_prefix (str): the prefix of the input file
        o_prefix (str): the prefix of the output file
        to_extract (set): the set of markers to extract

    """
    file_info = [
        dict(suffix=".alleles", header=True, name="name"),
        dict(suffix=".completion_rates", header=True, name="name"),
        dict(suffix=".good_sites", header=False, index=0),
        dict(suffix=".impute2_info", header=True, name="name"),
        dict(suffix=".imputed_sites", header=False, index=0),
        dict(suffix=".maf", header=True, name="name"),
        dict(suffix=".map", header=False, index=1),
    ]

    for info in file_info:
        # The name of the input file
        i_fn = i_prefix + info["suffix"]

        if not os.path.isfile(i_fn):
            # The file doesn't exist, so we continue
            continue

        # The name of the output file
        o_fn = o_prefix + info["suffix"]

        # If the file doesn't have a header, we just read line per line
        header = None
        with open(i_fn, "r") as i_file, open(o_fn, "w") as o_file:
            for i, line in enumerate(i_file):
                row = line.rstrip("\r\n").split(info.get("sep", "\t"))
                if info["header"] and i == 0:
                    header = {name: i for i, name in enumerate(row)}
                    if info["name"] not in header:
                        raise GenipeError("{}: missing column {}".format(
                            i_fn,
                            info["name"],
                        ))

                    info["index"] = header[info["name"]]
                    o_file.write(line)
                    continue

                if row[info["index"]] in to_extract:
                    o_file.write(line)

    # We need to copy the sample file
    sample_fn = i_prefix + ".sample"
    if os.path.isfile(sample_fn):
        o_fn = o_prefix + ".sample"
        shutil.copyfile(sample_fn, o_fn)


def print_data(o_files, prob_t, fid, iid, is_long, *, line=None, row=None):
    """Prints an impute2 line.

    Args:
        o_files (dict): the output files
        prob_t (float): the probability threshold
        fid (list): the list of family IDs
        iid (list): the list of sample IDs
        is_long (bool): True if the format is long (dosage, calls)
        line (str): the impute2 line
        row (list): the impute2 line, split by spaces

    """
    # Probabilities?
    if "impute2" in o_files:
        o_files["impute2"].write(line)

    # Require more?
    a1 = None
    a2 = None
    pos = None
    name = None
    chrom = None
    good_calls = None
    probabilities = None
    if ("dosage" in o_files) or ("calls" in o_files) or ("bed" in o_files):
        # Getting the informations
        marker_info, probabilities = impute2.matrix_from_line(row)
        chrom, name, pos, a1, a2 = marker_info

        # Getting the good calls
        good_calls = impute2.get_good_probs(probabilities, min_prob=prob_t)

    # Dosage?
    if "dosage" in o_files:
        # Getting the maf
        maf, minor, major = impute2.maf_from_probs(
            prob_matrix=probabilities[good_calls, :],
            a1=0,
            a2=2,
        )
        dosage = impute2.dosage_from_probs(
            homo_probs=probabilities[:, minor],
            hetero_probs=probabilities[:, 1],
            scale=2,
        )
        dosage[~good_calls] = nan

        alleles = [a1, nan, a2]

        if is_long:
            for sample_f, sample_i, sample_d in zip(fid, iid, dosage):
                print(sample_f, sample_i, chrom, pos, name, alleles[minor],
                      alleles[major], sample_d, sep="\t",
                      file=o_files["dosage"])
        else:
            print(chrom, pos, name, alleles[minor], alleles[major], *dosage,
                  sep="\t", file=o_files["dosage"])

    # Bed?
    if "bed" in o_files:
        geno, minor, major = impute2.additive_from_probs(a1, a2, probabilities)
        geno[~good_calls] = -1
        o_files["bed"][0].write_genotypes(geno)
        print(chrom, name, "0", pos, minor, major, sep="\t",
              file=o_files["bed"][1])

    # Hard calls?
    if "calls" in o_files:
        calls = impute2.hard_calls_from_probs(a1, a2, probabilities)
        calls[~good_calls] = "0 0"

        if is_long:
            for sample_f, sample_i, sample_c in zip(fid, iid, calls):
                print(sample_f, sample_i, chrom, name, "0", pos, sample_c,
                      sep="\t", file=o_files["calls"])
        else:
            print(chrom, name, "0", pos, *calls, sep="\t",
                  file=o_files["calls"])


def gather_extraction(fn, maf, rate, info, extract_filename, genomic_range):
    """Gather positions that are required.

    Args:
        fn (str): the impute2 filename
        maf (float): the minor allele frequency threshold (might be ``None``)
        rate (float): the call rate threshold (might be ``None``)
        info (float): the marker information value threshold (might be
                      ``None``)
        extract_filename (str): the name of the file containing marker names to
                                extract (might be ``None``)
        genomic_range (str): the genomic range for extraction

    Returns:
        set: the set of markers to extract

    If extraction by marker name is required, only those markers will be
    extracted. Otherwise, ``maf``, ``rate``, ``info`` or ``genomic_range`` can
    be specified (alone or together) to extract markers according to minor
    allele frequency, call rate and genomic location.

    """
    logging.info("Gathering information about {}".format(fn))

    # The prefix of all the input files
    prefix = get_file_prefix(fn)

    # Reading the map file
    logging.info("Reading MAP data")
    map_data = pd.read_csv(prefix + ".map", sep="\t", usecols=[0, 1, 3],
                           names=["chrom", "name", "pos"])
    map_data = map_data.set_index("name", verify_integrity=True)
    logging.info("MAP data contained {:,d} markers".format(len(map_data)))

    # If extraction, we only require a list of marker names
    if extract_filename is not None:
        available_markers = map_data.index
        marker_list = None
        with open(extract_filename, "r") as i_file:
            marker_list = set(i_file.read().splitlines())

        return set(available_markers.intersection(marker_list))

    # Do we require a genomic location?
    if genomic_range is not None:
        logging.info("Keeping markers in required genomic region")
        map_data = map_data[(
            (map_data.chrom == genomic_range.chrom) &
            (map_data.pos >= genomic_range.start) &
            (map_data.pos <= genomic_range.end)
        )]
        logging.info("Required genomic region contained {:,d} "
                     "markers".format(len(map_data)))

    # Do we require a certain MAF?
    if maf is not None:
        logging.info("Reading MAF data")
        maf_data = pd.read_csv(prefix + ".maf", sep="\t")
        maf_data = maf_data.set_index("name", verify_integrity=True)

        # Merging
        map_data = pd.merge(
            map_data,
            maf_data[maf_data.maf >= maf],
            how="inner",
            left_index=True,
            right_index=True,
        )
        logging.info("{:,d} markers with maf >= {}".format(len(map_data), maf))

    # Do we required a certain completion rate?
    if rate is not None:
        logging.info("Reading completion rates")
        rate_data = pd.read_csv(prefix + ".completion_rates", sep="\t",
                                usecols=[0, 2])
        rate_data = rate_data.set_index("name", verify_integrity=True)
        map_data = pd.merge(
            map_data,
            rate_data[rate_data.completion_rate >= rate],
            how="inner",
            left_index=True,
            right_index=True,
        )
        logging.info("{:,d} markers with completion rate >= "
                     "{}".format(len(map_data), rate))

    # Do we required a certain information value?
    if info is not None:
        logging.info("Reading information values")
        info_data = pd.read_csv(prefix + ".impute2_info", sep="\t")
        info_data = info_data.set_index("name", verify_integrity=True)
        map_data = pd.merge(
            map_data,
            info_data[info_data["info"] >= info],
            how="inner",
            left_index=True,
            right_index=True,
        )
        logging.info("{:,d} markers with information value >= "
                     "{}".format(len(map_data), info))

    # Extracting the names
    to_extract = set(map_data.index)

    if len(to_extract) == 0:
        logging.warning("No marker left for analysis")
        sys.exit(0)

    return to_extract


def get_file_prefix(fn):
    """Gets the filename prefix.

    Args:
        fn (str): the name of the file from which the prefix is required

    Returns:
        str: the prefix of the file

    This function removes the extension from the file name, and return its
    prefix (*e.g.* ``test.impute2`` returns ``test``, and
    ``../test.impute2.gz`` returns ``../test``).

    """
    prefix = os.path.splitext(fn)[0]
    if prefix.endswith("impute2"):
        prefix = os.path.splitext(prefix)[0]
    return prefix


def check_args(args):
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the options to verify

    Note
    ----
        If there is a problem, a :py:class:`genipe.error.GenipeError` is
        raised.

    Note
    ----
        Noting is checked (apart from the impute2 files) if indexation is asked
        (``--index`` option).

    """
    # Checking that the impute2 file exists
    if not os.path.isfile(args.impute2):
        raise GenipeError("{}: no such file".format(args.impute2))

    if args.index_only:
        return True

    # Is there something to extract?
    if not args.genomic and not args.maf and not args.rate and not args.info:
        if args.extract is None:
            raise GenipeError("nothing to extract: use '--extract', "
                              "'--genomic', '--maf' or '--rate'")

    elif args.extract is not None:
        raise GenipeError("'--extract' can only be use alone")

    # What extensions to look for after checking arguments
    extensions = set()

    # If extract, check the file
    if args.extract is not None:
        if not os.path.isfile(args.extract):
            raise GenipeError("{}: no such file".format(args.extract))

    # If genomic, we check the format
    if args.genomic is not None:
        genomic_match = re.match(r"(.+):(\d+)-(\d+)$", args.genomic)
        if not genomic_match:
            raise GenipeError("{}: no a valid genomic "
                              "region".format(args.genomic))
        chrom = int(genomic_match.group(1).replace("chr", ""))
        start = int(genomic_match.group(2))
        end = int(genomic_match.group(3))

        if chrom not in chromosomes:
            raise GenipeError("{}: invalid chromosome".format(chrom))

        if end < start:
            start, end = end, start

        GenomicRange = namedtuple("GenomicRange", ["chrom", "start", "end"])
        args.genomic = GenomicRange(chrom, start, end)

    # If MAF, we check what's required
    if args.maf is not None:
        extensions.add("map")
        extensions.add("maf")
        if args.maf < 0 or args.maf > 0.5:
            raise GenipeError("{}: invalid MAF".format(args.maf))

    # If completion rate, we check what's required
    if args.rate is not None:
        extensions.add("map")
        extensions.add("completion_rates")
        if args.rate < 0 or args.rate > 1:
            raise GenipeError("{}: invalid rate".format(args.rate))

    # If info, we check what's required
    if args.info is not None:
        extensions.add("map")
        extensions.add("impute2_info")
        if args.info < 0 or args.info > 1:
            raise GenipeError("{}: invalid information "
                              "value".format(args.info))

    # Checking the probability threshold
    if args.prob < 0 or args.prob > 1:
        raise GenipeError("{}: invalid probability "
                          "threshold".format(args.prob))

    # Checking the companion files (for impute2)
    f_prefix = get_file_prefix(args.impute2)
    for f_extension in extensions:
        fn = f_prefix + "." + f_extension
        if not os.path.isfile(fn):
            raise GenipeError("{}: no such file".format(fn))

    # Checking the output format
    for out_format in args.out_format:
        if out_format not in {"impute2", "dosage", "calls", "bed"}:
            raise GenipeError("{}: invalid output format".format(out_format))

        if out_format == "bed":
            if not HAS_PYPLINK:
                raise GenipeError("missing optional module: pyplink")

        if out_format in {"bed", "calls", "dosage"}:
            if not f_prefix + ".sample":
                raise GenipeError("{}: sample file missing".format(
                    f_prefix + ".sample"),
                )

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
    # Adding the version option for the main parser
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s, part of genipe version {}".format(__version__),
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

    # Indexation options
    group = parser.add_argument_group("Indexation Options")
    group.add_argument(
        "--index",
        dest="index_only",
        action="store_true",
        help="Only perform the indexation.",
    )

    # The output files
    group = parser.add_argument_group("Output Options")
    group.add_argument(
        "--out",
        type=str,
        metavar="PREFIX",
        default="impute2_extractor",
        help="The prefix of the output files. [%(default)s]",
    )
    group.add_argument(
        "--format",
        type=str,
        metavar="FORMAT",
        nargs="+",
        default=["impute2"],
        dest="out_format",
        help="The output format. Can specify either 'impute2' for "
             "probabilities (same as impute2 format, i.e. 3 values per "
             "sample), 'dosage' for dosage values (one value between 0 and 2 "
             "by sample), 'calls' for hard calls, or 'bed' for Plink binary "
             "format (with hard calls). %(default)s",
    )
    group.add_argument(
        "--long",
        action="store_true",
        dest="long_format",
        help="Write the output file in the long format (one line per sample "
             "per marker). This option is only compatible with the 'calls' "
             "and 'dosage' format (option '--format').",
    )
    group.add_argument(
        "--prob",
        type=float,
        metavar="FLOAT",
        default=0.9,
        help="The probability threshold used when creating a file in the "
             "dosage or call format. [%(default).1f]",
    )

    # What to extract
    group = parser.add_argument_group("Extraction Options")
    group.add_argument(
        "--extract",
        type=str,
        metavar="FILE",
        help="File containing marker names to extract.",
    )
    group.add_argument(
        "--genomic",
        type=str,
        metavar="CHR:START-END",
        help="The range to extract (e.g. 22 1000000 1500000). Can be use in "
             "combination with '--rate', '--maf' and '--info'.",
    )
    group.add_argument(
        "--maf",
        type=float,
        metavar="FLOAT",
        help="Extract markers with a minor allele frequency equal or higher "
             "than the specified threshold. Can be use in combination with "
             "'--rate', '--info' and '--genomic'.",
    )
    group.add_argument(
        "--rate",
        type=float,
        metavar="FLOAT",
        help="Extract markers with a completion rate equal or higher to the "
             "specified threshold. Can be use in combination with '--maf', "
             "'--info' and '--genomic'.",
    )
    group.add_argument(
        "--info",
        type=float,
        metavar="FLOAT",
        help="Extract markers with an information equal or higher to the "
             "specified threshold. Can be use in combination with '--maf', "
             "'--rate' and '--genomic'.",
    )

    if args is not None:
        return parser.parse_args(args)

    return parser.parse_args()


# Calling the main, if necessary
if __name__ == "__main__":
    main()
