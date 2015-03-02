
# This file is part of gwip.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import os
import logging
import argparse

import pandas as pd

from .. import __version__
from ..error import ProgramError


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


def main(args=None):
    """The main function."""
    # Creating the option parser
    desc = ("Performs statistical analysis on imputed data (linear, logistic "
            "or Cox's regressions) (gwip version {})".format(__version__))
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

        # Reading the phenotype file
        logging.info("Reading phenotype file")
        phenotypes, remove_gender = read_phenotype(args.pheno, args)

        # Reading the sample file
        logging.info("Reading the sample file")
        samples = read_samples(args.sample)

        # Reading the sites to extract (if required)
        sites_to_extract = None
        if args.extract_sites is not None:
            sites_to_extract = read_sites_to_extract(args.extract_sites)

        # Computing the statistics
        compute_statistics(
            impute2_filename=args.impute2,
            samples=samples,
            markers_to_extract=sites_to_extract,
            phenotypes=phenotypes,
            remove_gender=remove_gender,
            out_prefix=args.out,
            other_opts=args,
        )

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


def read_phenotype(i_filename, opts):
    """Reads the phenotype file."""
    # Reading the data (and setting the index)
    pheno = pd.read_csv(i_filename, sep="\t", na_values=opts.missing_value)
    pheno = pheno.set_index(opts.sample_column, verify_integrity=True)

    # Finding the required column
    required_columns = opts.covar
    if opts.analysis_type == "cox":
        required_columns.extend([opts.tte, opts.censure])
    else:
        required_columns.append(opts.pheno_name)

    # Interaction?
    if opts.interaction is not None:
        if opts.interaction not in required_columns:
            required_columns.append(opts.interaction)

    # Do we need the gender column?
    remove_gender_column = False
    if opts.chrx and (opts.gender_column not in required_columns):
        # It's not already in the gender column, so we keep it, but we need to
        # remove it later on...
        required_columns.append(opts.gender_column)
        remove_gender_column = True

    # Extracting the required column
    pheno = pheno.loc[:, required_columns]

    # We need to exclude unknown gender if we are on chrX
    if opts.chrx:
        pheno = pheno[(pheno[opts.gender_column] != 1) |
                      (pheno[opts.gender_column] != 2)]

    # Returning the phenotypes
    return pheno.dropna(), remove_gender_column


def read_samples(i_filename):
    """Reads the sample file (produced by SHAPEIT)."""
    samples = pd.read_csv(i_filename, sep=" ", usecols=[0, 1])
    samples = samples.drop(samples.index[0])
    return samples.set_index("ID_2", verify_integrity=True)


def read_sites_to_extract(i_filename):
    """Reads the list of sites to extract."""
    markers_to_extract = None
    with open(i_filename, "r") as i_file:
        markers_to_extract = set(i_file.read().splitlines())
    return markers_to_extract


def compute_statistics(impute2_filename, samples, markers_to_extract,
                       phenotypes, remove_gender, out_prefix, other_opts):
    """Parses IMPUTE2 file while computing statistics."""
    # The name of the output file
    o_name = "{}.{}.dosage".format(out_prefix, other_opts.analysis_type)

    # Reading the IMPUTE2 file one line (site) at a time, creating a subprocess
    # if required
    proc = None
    i_file = None
    o_file = open(o_name, "w")
    sites_to_process = []
    pool = None

    # Multiprocessing?
    if other_opts.nb_process > 1:
        pool = Pool(processes=other_opts.nb_process)

    try:
        if impute2_filename.endswith(".gz"):
            proc = Popen(["gzip", "-d", "-c", impute2_filename], stdout=PIPE)
            i_file = proc.stdout

        else:
            i_file = open(impute2_filename, "rb")

        # Printing the header of the output file
        print("chr", "pos", "snp", "major", "minor", "maf", "n", "coef", "se",
              "lower", "upper",
              "z" if other_opts.analysis_type == "cox" else "t", "p", sep="\t",
              file=o_file)

    except Exception as e:
        if pool is not None:
            pool.terminate()
        raise

    finally:
        # Closing the input file
        i_file.close()

        # Finishing the rows if required
        if other_opts.nb_process > 1:
            if len(sites_to_process) > 0:
                pass
            pool.close()

        # Closing the output file
        o_file.close()

        # Closing the proc
        if proc is not None:
            if proc.wait() != 0:
                raise ProgramError("{}: problem while reading the GZ "
                                   "file".format(impute2_filename))

def check_args(args):
    """Checks the arguments and options."""
    # Checking the required input files
    for filename in [args.impute2, args.sample, args.pheno]:
        if not os.path.isfile(filename):
            raise ProgramError("{}: no such file".format(filename))

    # Checking the optional input files
    for filename in [args.extract_sites]:
        if filename is not None:
            if not os.path.isfile(filename):
                raise ProgramError("{}: no such file".format(filename))

    # Checking the number of process
    if args.nb_process < 1:
        raise ProgramError("{}: invalid number of "
                           "processes".format(args.nb_process))

    # Checking the number of lines to read
    if args.nb_lines < 1:
        raise ProgramError("{}: invalid number of lines to "
                           "read".format(args.nb_lines))

    # Checking the MAF threshold
    if args.maf < 0 or args.maf > 1:
        raise ProgramError("{}: invalid MAF".format(args.maf))

    # Checking the probability threshold
    if args.prob < 0 or args.prob > 1:
        raise ProgramError("{}: invalid probability "
                           "threshold".format(args.prob))

    # Reading all the variables in the phenotype file
    header = None
    with open(args.pheno, "r") as i_file:
        header = {name for name in i_file.readline().rstrip("\n").split("\t")}

    # Checking the restricted columns
    restricted_columns = {"_D1", "_D2", "_D3", "_MaxD", "_GenoD", "_Inter"}
    if len(restricted_columns & header) != 0:
        raise ProgramError("{}: {}: restricted variables".format(
            args.pheno,
            restricted_columns & header,
        ))

    # Checking the required columns
    variables_to_check = None
    if args.analysis_type == "cox":
        variables_to_check = {args.tte, args.censure}
    else:
        variables_to_check = {args.pheno_name}
    for variable in variables_to_check:
        if variable not in header:
            raise ProgramError("{}: {}: missing variable for {}".format(
                args.pheno,
                variable,
                args.analysis_type,
            ))

    # Checking the co-variables
    covar_list = []
    if args.covar != "":
        covar_list = args.covar.split(",")
    for covar in covar_list:
        if covar not in header:
            raise ProgramError("{}: {}: missing co-variable".format(
                args.pheno,
                covar,
            ))
    args.covar = covar_list

    # Checking the sample column
    if args.sample_column not in header:
        raise ProgramError("{}: {}: no such column (--sample-column)".format(
            args.pheno,
            args.sample_column,
        ))

    # Checking the gender column (only if required)
    if args.chrx:
        if args.gender_column not in header:
            raise ProgramError(
                "{}: {}: no such column (--gender-column)".format(
                    args.pheno,
                    args.gender_column,
                )
            )

    # Checking the interaction column (if required)
    if args.interaction is not None:
        if args.interaction not in header:
            raise ProgramError(
                "{}: {}: no such column (--interaction)".format(
                    args.pheno,
                    args.interaction,
                )
            )

    return True


def parse_args(parser, args=None):
    """Parses the command line options and arguments."""
    # The parser object
    parser.add_argument("--version", action="version",
                        version="%(prog)s (part of GWIP "
                                "version {})".format(__version__))
    parser.add_argument("--debug", action="store_true",
                        help="Set the logging level to debug.")

    # Sub parsers
    subparsers = parser.add_subparsers(help="Analysis type",
                                       dest="analysis_type")
    cox_parser = subparsers.add_parser("cox", help="Cox (survival regression)")
    lin_parser = subparsers.add_parser("linear", help="Linear regression")
    logit_parser = subparsers.add_parser("logistic", help="Logistic regrssion")

    # The input files
    group = parser.add_argument_group("Input Files")
    group.add_argument("--impute2", type=str, metavar="FILE", required=True,
                       help="The output from IMPUTE2.")
    group.add_argument("--sample", type=str, metavar="FILE", required=True,
                       help=("The sample file (the order should be the same "
                             "as in the IMPUTE2 files)."))
    group.add_argument("--pheno", type=str, metavar="FILE", required=True,
                       help="The file containing phenotypes and co variables.")
    group.add_argument("--extract-sites", type=str, metavar="FILE",
                       help="The 'good' sites to extract for analysis.")

    # The dosage options
    group = parser.add_argument_group("Dosage Options")
    group.add_argument("--scale", type=int, metavar="INT", default=2,
                       choices=[1, 2],
                       help=("Scale dosage so that values are in [0, n] "
                             "(possible values are 1 (no scaling) or 2). "
                             "[%(default)d]"))
    group.add_argument("--prob", type=float, metavar="FLOAT", default=0.9,
                       help=("The minimal probability for which a genotype "
                             "should be considered. [>=%(default).1f]"))
    group.add_argument("--maf", type=float, metavar="FLOAT", default=0.01,
                       help=("Minor allele frequency threshold for which "
                             "marker will be skipped. [<%(default).2f]"))

    # The phenotype options for cox
    group = cox_parser.add_argument_group("Phenotype Options")
    group.add_argument("--time-to-event", type=str, metavar="NAME",
                       required=True, dest="tte",
                       help="The time to event variable (in the pheno file).")
    group.add_argument("--censure", type=str, metavar="NAME", required=True,
                       help=("The censure value (1 if observed, 0 if "
                             "censored)"))

    # The phenotype options for linear regression
    group = lin_parser.add_argument_group("Phenotype Options")
    group.add_argument("--pheno-name", type=str, metavar="NAME", required=True,
                       help="The phenotype.")

    # The phenotype options for logistic regression
    group = logit_parser.add_argument_group("Phenotype Options")
    group.add_argument("--pheno-name", type=str, metavar="NAME", required=True,
                       help="The phenotype.")

    # The general phenotype options
    group = parser.add_argument_group("Phenotype Options")
    group.add_argument("--covar", type=str, metavar="NAME", default="",
                       help=("The co variable names (in the pheno file), "
                             "separated by coma."))
    group.add_argument("--missing-value", type=str, metavar="NAME",
                       help="The missing value.")
    group.add_argument("--sample-column", type=str, metavar="NAME",
                       default="sample_id",
                       help=("The name of the sample ID column (in the pheno "
                             "file). [%(default)s]"))
    group.add_argument("--interaction", type=str, metavar="NAME",
                       help=("Add an interaction between the genotype and "
                             "this co-variable"))

    # General options
    group = parser.add_argument_group("General Options")
    group.add_argument("--nb-process", type=int, metavar="INT", default=1,
                       help="The number of process to use. [%(default)d]")
    group.add_argument("--nb-lines", type=int, metavar="INT", default=1000,
                       help=("The number of line to read at a time."
                             "[%(default)d]"))
    group.add_argument("--chrx", action="store_true",
                       help=("The analysis is performed for the non pseudo-"
                             "autosomal region of the chromosome X (male "
                             "dosage will be divided by 2 to get values "
                             "[0, 0.5] instead of [0, 1]) (males are coded "
                             "as 1 and option '--gender-column' should be "
                             "used)."))
    group.add_argument("--gender-column", type=str, metavar="NAME",
                       default="Gender",
                       help=("The name of the gender column (use in "
                             "conjunction with the '--chrx' option) "
                             "[%(default)s]"))

    # The output files
    group = parser.add_argument_group("Output Options")
    group.add_argument("--out", metavar="FILE", default="imputed_stats",
                       help="The prefix for the output files. [%(default)s]")

    if args is not None:
        return parser.parse_args(args)

    return parser.parse_args()


# Calling the main, if necessary
if __name__ == "__main__":
    main()
