
# This file is part of gwip.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.

import os
import sys
import datetime
import logging
import argparse
import platform
import traceback
from multiprocessing import Pool
from subprocess import Popen, PIPE
from collections import namedtuple

import numpy as np
import pandas as pd
from lifelines import CoxPHFitter

import statsmodels.api as sm
import statsmodels.formula.api as smf

from .. import __version__
from ..formats.impute2 import *
from ..error import ProgramError


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


# An IMPUTE2 row to process
_Row = namedtuple("_Row", ("row", "samples", "pheno", "pheno_name", "formula",
                           "time_to_event", "censure", "inter_c", "is_chrx",
                           "gender_c", "del_g", "scale", "maf_t", "prob_t",
                           "analysis_type", "number_to_print"))


def main(args=None):
    """The main function."""
    # Creating the option parser
    desc = ("Performs statistical analysis on imputed data (either linear, "
            "logistic, SKAT or Cox's regressions). This script is part of the "
            "'gwip' package, version {}).".format(__version__))
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
        if args.analysis_type != "skat":
            compute_statistics(
                impute2_filename=args.impute2,
                samples=samples,
                markers_to_extract=sites_to_extract,
                phenotypes=phenotypes,
                remove_gender=remove_gender,
                out_prefix=args.out,
                options=args,
            )
        else:
            skat_parse_impute2(
                impute2_filename=args.impute2,
                samples=samples,
                markers_to_extract=sites_to_extract,
                phenotypes=phenotypes,
                remove_gender=remove_gender,
                out_prefix=args.out,
                args=args,
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
    """Reads the phenotype file.

    This file is expected to be a tab separated file of phenotypes and
    covariates. The columns to use will be determined by the `--sample-column`
    and the `--covar` options.

    For analysis including the X chromosome, the gender is automatically
    added as a covariate. The results are not shown to the user unless asked
    for.

    """
    # Reading the data (and setting the index)
    pheno = pd.read_csv(i_filename, sep="\t", na_values=opts.missing_value)
    pheno = pheno.set_index(opts.sample_column, verify_integrity=True)

    # Finding the required column
    required_columns = opts.covar.copy()
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

    # We need to exclude unknown gender only if gender was required
    if (opts.gender_column in required_columns) and (not remove_gender_column):
        pheno = pheno[(pheno[opts.gender_column] == 1) |
                      (pheno[opts.gender_column] == 2)]

    # Extracting the required column
    pheno = pheno.loc[:, required_columns]

    # Returning the phenotypes (dropping the nan values)
    return pheno.dropna(axis=0), remove_gender_column


def read_samples(i_filename):
    """Reads the sample file (produced by SHAPEIT).

    The expected format for this file is a tab separated file with a first
    row containing the following columns: ::

        ID_1	ID_2	missing	father	mother	sex	plink_pheno

    The subsequent row will be discarded and should contain: ::

        0	0	0 D	D	D	B

    We are mostly interested in the sample IDs corresponding to the `ID_2`
    column.

    Their uniqueness is verified by pandas.

    """
    samples = pd.read_csv(i_filename, sep=" ", usecols=[0, 1])
    samples = samples.drop(samples.index[0], axis=0)
    return samples.set_index("ID_2", verify_integrity=True)


def skat_read_snp_set(i_filename):
    """Reads the SKAT SNP set file.

    This file has to be supplied by the user. The recognized columns are:
    `variant, snp_set, weight`. The `weight` column is optional and can
    be used to specify a custom weighting scheme for SKAT. If nothing is
    specified, the default Beta weights are used.

    The file has to be tab delimited.

    """
    skat_info = pd.read_csv(i_filename, sep="\t", header=0)
    if "variant" not in skat_info.columns:
        raise ProgramError("The SKAT SNP set file needs to have a 'variant' "
                           "column containing the ID of every variant of "
                           "interest.")

    if "snp_set" not in skat_info.columns:
        raise ProgramError("The SKAT SNP set file needs to have a 'snp_set' "
                           "column containing the SNP set ID for every "
                           "variant. The user is free to choose the SNP ID")

    return skat_info


def read_sites_to_extract(i_filename):
    """Reads the list of sites to extract.

    :param i_filename: The input filename containing the IDs of the variants
                       to consider for the analysis.
    :type i_filename: str

    :returns: A set containing the variants.
    :rtype: set

    The expected file format is simply a list of variants. Every row should
    correspond to a single variant identifier.

    e.g. ::

        3:60069:t
        rs123456:A
        3:60322:A

    Typically, this is used to analyze only variants that passed some QC
    threshold. The _gwip_ pipeline generates this file at the 'merge_impute2'
    step.

    """
    markers_to_extract = None
    with open(i_filename, "r") as i_file:
        markers_to_extract = set(i_file.read().splitlines())
    return markers_to_extract


def skat_parse_impute2(impute2_filename, samples, markers_to_extract,
                       phenotypes, remove_gender, out_prefix, args):
    """Read the impute2 file and generate the files required by SKAT."""

    # Read the SNP set and create the output files.
    snp_set = skat_read_snp_set(args.snp_sets)

    # We will use a temporary directory for our analysis.
    dir_name = os.path.abspath(
        datetime.datetime.today().strftime("%Y-%m-%d_%H.%M.%S.skat")
    )

    if not os.path.isdir(dir_name):
        os.makedirs(dir_name)
    else:
        # This should not happen because the folder name contains a timestamp
        raise ProgramError("A folder named '{}' already exists.".format(
            dir_name
        ))

    # Open genotype CSV files for every SNP set. Those CSV files will be
    # read in R and used by SKAT.
    genotype_files = {}
    for set_id in snp_set["snp_set"].unique():
        filename = os.path.join(dir_name, "{}.genotypes.csv".format(set_id))
        genotype_files[set_id] = open(filename, "w")
        # We write the column headers for every file.
        print("", *samples.index, sep=",", file=genotype_files[set_id])

    # The markers of interest are the markers we want to include in the
    # analysis. Concretly they are the markers that were included in the snp
    # sets. If a list of markers to extract is provided by the user, we also
    # look at this to filter ou undesired markers.
    #
    # We also fill a set of written_markers to stop whenever all the markers
    # were written to file.
    markers_of_interest = set(snp_set["variant"])
    written_markers = set()
    if markers_to_extract is not None:
        markers_of_interest = markers_of_interest & markers_to_extract

    # Open the file. We use subprocess if it's gunzipped because it's faster.
    # We use gzip -d -c instead of zcat because the default Mac OS zcat has
    # weird behavior.
    if impute2_filename.endswith(".gz"):
        proc = Popen(["gzip", "-d", "-c", impute2_filename], stdout=PIPE)
        i_file = proc.stdout
    else:
        i_file = open(impute2_filename, "rb")

    for line in i_file:
        # If we already found everything, we will stop here.
        if len(written_markers) == len(markers_of_interest):
            if written_markers == markers_of_interest:
                logging.info("Found all the necessary markers, skipping the "
                             "rest of the file.")
                break

        line = line.decode("ascii")

        # TODO: Add the gender and do QC (men with hetero on chrX).
        line = _skat_parse_line(line, markers_of_interest, samples)
        if line is not None:
            name, dosage = line
            written_markers.add(name)
            _skat_write_marker(name, dosage, snp_set, genotype_files)

    # Close the genotype files.
    for file_handle in genotype_files.values():
        file_handle.close()

    # Write the covariate file.
    phenotype_df = samples.join(phenotypes)
    phenotype_df.index.name = "sample"

    if args.covar:
        # Make sure the samples are consistent by merging phenotype and
        # samples.
        phenotype_df[args.covar].to_csv(
            os.path.join(dir_name, "covariates.csv"),
            sep=",",
        )

    # Write the phenotype file.
    phenotype_df[[args.pheno_name]].to_csv(
        os.path.join(
            dir_name,
            "{}.{}.csv".format(args.pheno_name, args.outcome_type)
        ),
        sep=",",
    )


def _skat_parse_line(line, markers_of_interest, samples, gender=None):
    """Parses a single line of the Impute2 file.

    :param line: A line from the Impute2 file.
    :type line: str

    :param markers_of_interest: A set of markers that are required for the
                                analysis.
    :type markers_of_interest: set

    :param samples: A DataFrame containing the samples IDs. This is useful
                    to make sure we return a dosage vector with the appropriate
                    data.
    :type samples: :py:class:`pandas.DataFrame`

    :returns: Either None if the marker is not of interest or a tuple of `(name
              , dosage_vector)` where `name` is a `str` representing the
              variant ID and `dosage_vector` is a numpy array containing the
              dosage values for every sample in the `samples` dataframe.
    :rtype: tuple or None

    """
    line = line.split(" ")
    # info_tuple contains: chrom, name, pos, a1, a2
    # proba_matrix is a matrix of sample x (aa, ab, bb)
    info_tuple, proba_matrix = matrix_from_line(line)

    chrom, name, pos, a1, a2 = info_tuple

    # If this marker is not of interest, we don't bother continuing.
    if name not in markers_of_interest:
        return None

    # We need to compute the dosage vector.
    # This is given wrt to the minor and major alleles, so we need to get the
    # MAF to identify those.
    maf, minor, major = maf_from_probs(proba_matrix, a1, a2, gender, name)

    # The minor allele columns could be either the first or the third.
    # We identify it here.
    minor_allele_col = 0 if a1 == minor else 2

    # We don't pass a scale parameter, because we want additive coding.
    dosage = dosage_from_probs(
        homo_probs=proba_matrix[:, minor_allele_col],
        hetero_probs=proba_matrix[:, 1],
    )

    return (name, dosage)


def _skat_write_marker(name, dosage, snp_set, genotype_files):
    """Write the dosage information to the appropriate genotype file.

    :param name: The name of the marker.
    :type name: str

    :param dosage: The dosage vector.
    :type dosage: :py:class:`numpy.ndarray`

    :param snp_set: The dataframe that allows us to identify the correct SNP
                    set for the specified variant.
    :type snp_set: :py:class:`pandas.DataFrame`

    :param genotype_files: The dict containing the opened CSV files for the
                           genotypes.
    :type genotype_files: dict

    """
    # Identify the correct SNP set.
    this_snp_set = snp_set.loc[snp_set["variant"] == name, "snp_set"].unique()
    for set_id in this_snp_set:
        file_object = genotype_files[set_id]
        print(name, *dosage, sep=",", file=file_object)

def compute_statistics(impute2_filename, samples, markers_to_extract,
                       phenotypes, remove_gender, out_prefix, options):
    """Parses IMPUTE2 file while computing statistics.

    This function takes care of parallelism. It reads the Impute2 file and
    fills a queue that will trigger the analysis when full.

    If the number of process to launch is 1, the rows are analyzed as they
    come.

    """
    # The name of the output file
    o_name = "{}.{}.dosage".format(out_prefix, options.analysis_type)

    # Do we need to create a formula?
    formula = None
    if options.analysis_type != "cox":
        formula = get_formula(
            phenotype=options.pheno_name,
            covars=options.covar,
            interaction=options.interaction,
        )
        logging.info("Using: {}".format(formula))

    # Reading the IMPUTE2 file one line (site) at a time, creating a subprocess
    # if required
    proc = None
    i_file = None
    o_file = open(o_name, "w")
    pool = None

    # Multiprocessing?
    if options.nb_process > 1:
        pool = Pool(processes=options.nb_process)

    try:
        if impute2_filename.endswith(".gz"):
            proc = Popen(["gzip", "-d", "-c", impute2_filename], stdout=PIPE)
            i_file = proc.stdout

        else:
            i_file = open(impute2_filename, "rb")

        # Printing the header of the output file
        header = ("chr", "pos", "snp", "major", "minor", "maf", "n", "coef",
                  "se", "lower", "upper",
                  "t" if options.analysis_type == "linear" else "z", "p")
        print(*header, sep="\t", file=o_file)

        # The sites to process (if multiprocessing)
        sites_to_process = []

        # Reading the file
        for line in i_file:
            row = line.decode().rstrip("\r\n").split(" ")

            # Is this site required?
            if markers_to_extract and (row[1] not in markers_to_extract):
                continue

            # Constructing the row object
            site = _Row(
                row=row,
                samples=samples,
                pheno=phenotypes,
                pheno_name=vars(options).get("pheno_name", None),
                formula=formula,
                time_to_event=vars(options).get("tte", None),
                censure=vars(options).get("censure", None),
                inter_c=options.interaction,
                is_chrx=options.chrx,
                gender_c=options.gender_column,
                del_g=remove_gender,
                scale=options.scale,
                maf_t=options.maf,
                prob_t=options.prob,
                analysis_type=options.analysis_type,
                number_to_print=len(header),
            )

            # Is there more than one process
            if options.nb_process > 1:
                # Saving this site to process later
                sites_to_process.append(site)

                # Is there enough sites to process?
                if len(sites_to_process) >= options.nb_lines:
                    for result in pool.map(process_impute2_site,
                                           sites_to_process):
                        print(*result, sep="\t", file=o_file)
                    sites_to_process = []

            else:
                # Processing this row
                print(*process_impute2_site(site), sep="\t", file=o_file)

        if len(sites_to_process) > 0:
            for result in pool.map(process_impute2_site, sites_to_process):
                print(*result, sep="\t", file=o_file)

    except Exception as e:
        if pool is not None:
            pool.terminate()
        logging.error(traceback.format_exc())
        raise

    finally:
        # Closing the input file
        i_file.close()

        # Finishing the rows if required
        if (options.nb_process > 1) and (pool is not None):
            pool.close()

        # Closing the output file
        o_file.close()

        # Closing the proc
        if proc is not None:
            if proc.wait() != 0:
                raise ProgramError("{}: problem while reading the GZ "
                                   "file".format(impute2_filename))


def process_impute2_site(site_info):
    """Process an IMPUTE2 site (a line in an IMPUTE2 file)."""
    # Getting the probability matrix and site information
    (chrom, name, pos, a1, a2), geno = matrix_from_line(site_info.row)

    # The name of the dosage column
    dosage_columns = ["_D1", "_D2", "_D3"]

    # Allele encoding
    allele_encoding = {dosage_columns[0]: a1, dosage_columns[-1]: a2}

    # Creating the sample data frame
    samples = site_info.samples
    for i, col_name in enumerate(dosage_columns):
        samples[col_name] = geno[:, i]

    # Merging with phenotypes
    data = pd.merge(
        site_info.pheno,
        samples,
        how="inner",
        left_index=True,
        right_index=True
    ).dropna(axis=0)[list(site_info.pheno.columns) + dosage_columns]

    # Keeping only good quality markers
    data = data[
        get_good_probs(data[dosage_columns].values, site_info.prob_t)
    ]

    # Checking gender if required
    gender = None
    if site_info.is_chrx:
        # We want to exclude males with heterozygous calls for the rest of the
        # analysis
        invalid_rows = samples_with_hetero_calls(
            data.loc[data[site_info.gender_c] == 1, dosage_columns],
            dosage_columns[1]
        )
        if len(invalid_rows) > 0:
            logging.warning("There were {:,d} males with heterozygous "
                            "calls for {}".format(len(invalid_rows), name))
            logging.debug(data.shape)
            data = data.drop(invalid_rows, axis=0)
            logging.debug(data.shape)
            logging.debug(invalid_rows.isin(data.index))

        # Getting the genders
        gender = data[site_info.gender_c].values

    # Computing the frequency
    maf, minor, major = maf_from_probs(data[dosage_columns].values,
                                       dosage_columns[0], dosage_columns[-1],
                                       gender, name)

    # What we want to print
    to_return = [chrom, pos, name, allele_encoding[major],
                 allele_encoding[minor], maf, data.shape[0]]

    # If the marker is too rare, we continue with the rest
    if (maf == "NA") or (maf < site_info.maf_t):
        to_return.extend(["NA"] * (site_info.number_to_print - len(to_return)))
        return to_return

    # Computing the dosage on the minor allele
    data["_GenoD"] = dosage_from_probs(
        homo_probs=data[minor],
        hetero_probs=data[dosage_columns[1]],
        scale=site_info.scale,
    )

    # We might need to specifically add the interaction column (if Cox and if
    # interaction)
    if (site_info.inter_c is not None) and (site_info.analysis_type == "cox"):
        data["_Inter"] = data._GenoD * data[site_info.inter_c]

    # Removing the unwanted columns
    unwanted_columns = dosage_columns
    if site_info.del_g:
        unwanted_columns.append(site_info.gender_c)
    data = data.drop(unwanted_columns, axis=1)

    # The column to get the result from
    result_from_column = "_GenoD"
    if site_info.inter_c is not None:
        result_from_column = "_GenoD:" + site_info.inter_c
        if site_info.analysis_type == "cox":
            result_from_column = "_Inter"

    # Fitting
    results = _fit_map[site_info.analysis_type](
        data=data,
        time_to_event=site_info.time_to_event,
        censure=site_info.censure,
        formula=site_info.formula,
        result_col=result_from_column,
    )

    # Extending the list to return
    if len(results) == 0:
        results = ["NA"] * (site_info.number_to_print - len(to_return))
    to_return.extend(results)

    return to_return


def samples_with_hetero_calls(data, hetero_c):
    """Gets male and heterozygous calls."""
    return data[data.idxmax(axis=1) == hetero_c].index


def get_formula(phenotype, covars, interaction):
    """Creates the linear/logistic regression formula (for statsmodel)."""
    # The phenotype and genetics
    formula = "{} ~ _GenoD".format(phenotype)

    # Are there any covars?
    if len(covars) > 0:
        formula += " + " + " + ".join(covars)

    # Is there an interaction term?
    if interaction is not None:
        formula += " + _GenoD*{}".format(interaction)

    return formula


def fit_cox(data, time_to_event, censure, result_col, **kwargs):
    """Fit a Cox' regression to the data."""
    cf = CoxPHFitter(alpha=0.95, tie_method="Efron", normalize=False)

    failed = False
    try:
        cf.fit(data, duration_col=time_to_event, event_col=censure)
    except np.linalg.linalg.LinAlgError:
        failed = True

    # The results
    result = []
    if not failed:
        result.extend([
            float(cf.hazards_[result_col]),
            float(cf._compute_standard_errors()[result_col]),
            cf._compute_confidence_intervals().loc["lower-bound", result_col],
            cf._compute_confidence_intervals().loc["upper-bound", result_col],
            cf._compute_z_values()[result_col],
            cf._compute_p_values()[cf.hazards_.columns == result_col][0]
        ])

    return result


def fit_linear(data, formula, result_col, **kwargs):
    """Fit a linear regression to the data."""
    return _get_result_from_linear_logistic(
        smf.ols(formula=formula, data=data).fit(),
        result_col=result_col,
    )


def fit_logistic(data, formula, result_col, **kwargs):
    """Fit a logistic regression to the data."""
    return _get_result_from_linear_logistic(
        smf.glm(formula=formula, data=data,
                family=sm.families.Binomial()).fit(),
        result_col=result_col,
    )


def _get_result_from_linear_logistic(fit_result, result_col):
    """Gets results from either a linear or logistic regression."""
    conf_int = fit_result.conf_int().loc[result_col, :].values
    assert len(conf_int) == 2
    return [
        fit_result.params[result_col],
        fit_result.bse[result_col],
        conf_int[0],
        conf_int[1],
        fit_result.tvalues[result_col],
        fit_result.pvalues[result_col],
    ]


_fit_map = {"cox": fit_cox, "linear": fit_linear, "logistic": fit_logistic}


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
    if args.nb_process > 1 and platform.system() == "Darwin":
        raise ProgramError("multiprocessing is not supported on Mac OS")

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
    if args.gender_column != "None":
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
    # Adding the version option for the main parser
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s, part of gwip version {}".format(__version__),
    )

    # Creating a parent parser (for common options between analysis type)
    p_parser = argparse.ArgumentParser(add_help=False)

    # The parser object
    p_parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s (part of gwip version {})".format(__version__),
    )
    p_parser.add_argument(
        "--debug",
        action="store_true",
        help="set the logging level to debug",
    )

    # The input files
    group = p_parser.add_argument_group("Input Files")
    group.add_argument(
        "--impute2",
        type=str,
        metavar="FILE",
        required=True,
        help="The output from IMPUTE2.",
    )
    group.add_argument(
        "--sample",
        type=str,
        metavar="FILE",
        required=True,
        help="The sample file (the order should be the same as in the IMPUTE2 "
             "files).",
    )
    group.add_argument(
        "--pheno",
        type=str,
        metavar="FILE",
        required=True,
        help="The file containing phenotypes and co variables.",
    )
    group.add_argument(
        "--extract-sites",
        type=str,
        metavar="FILE",
        help="A list of sites to extract for analysis (optional).",
    )

    # The output files
    group = p_parser.add_argument_group("Output Options")
    group.add_argument(
        "--out",
        metavar="FILE",
        default="imputed_stats",
        help="The prefix for the output files. [%(default)s]",
    )

    # General options
    group = p_parser.add_argument_group("General Options")
    group.add_argument(
        "--nb-process",
        type=int,
        metavar="INT",
        default=1,
        help="The number of process to use. [%(default)d]",
    )
    group.add_argument(
        "--nb-lines",
        type=int,
        metavar="INT",
        default=1000,
        help="The number of line to read at a time. [%(default)d]",
    )
    group.add_argument(
        "--chrx",
        action="store_true",
        help="The analysis is performed for the non pseudo-autosomal region "
             "of the chromosome X (male dosage will be divided by 2 to get "
             "values [0, 0.5] instead of [0, 1]) (males are coded as 1 and "
             "option '--gender-column' should be used).",
    )
    group.add_argument(
        "--gender-column",
        type=str,
        metavar="NAME",
        default="Gender",
        help="The name of the gender column (use to exclude samples with "
             "unknown gender (i.e. not 1, male, or 2, female). If gender not "
             "available, use 'None'. [%(default)s]",
    )

    # The dosage options
    group = p_parser.add_argument_group("Dosage Options")
    group.add_argument(
        "--scale",
        type=int,
        metavar="INT",
        default=2,
        choices=[1, 2],
        help="Scale dosage so that values are in [0, n] (possible values are "
             "1 (no scaling) or 2). [%(default)d]",
    )
    group.add_argument(
        "--prob",
        type=float,
        metavar="FLOAT",
        default=0.9,
        help="The minimal probability for which a genotype should be "
             "considered. [>=%(default).1f]",
    )
    group.add_argument(
        "--maf",
        type=float,
        metavar="FLOAT",
        default=0.01,
        help="Minor allele frequency threshold for which marker will be "
             "skipped. [<%(default).2f]",
    )

    # The general phenotype options
    group = p_parser.add_argument_group("Phenotype Options")
    group.add_argument(
        "--covar",
        type=str,
        metavar="NAME",
        default="",
        help="The co variable names (in the phenotype file), separated by "
             "coma.",
    )
    group.add_argument(
        "--missing-value",
        type=str,
        metavar="NAME",
        help="The missing value in the phenotype file.",
    )
    group.add_argument(
        "--sample-column",
        type=str,
        metavar="NAME",
        default="sample_id",
        help="The name of the sample ID column (in the phenotype file). "
             "[%(default)s]",
    )
    group.add_argument(
        "--interaction",
        type=str,
        metavar="NAME",
        help="Add an interaction between the genotype and this variable.",
    )

    # Sub parsers
    subparsers = parser.add_subparsers(
        title="Statistical Analysis Type",
        description="The type of statistical analysis to be performed on the "
                    "imputed data.",
        dest="analysis_type",
    )

    subparsers.required = True

    # The Cox parser
    cox_parser = subparsers.add_parser(
        "cox",
        help="Cox's regression (survival analysis).",
        description="Performs a Cox's regression (survival analysis) on "
                    "imputed data. This script is part of the 'gwip' "
                    "package, version {}).".format(__version__),
        parents=[p_parser],
    )

    group = cox_parser.add_argument_group("Cox's Regression Options")
    group.add_argument(
        "--time-to-event",
        type=str,
        metavar="NAME",
        required=True,
        dest="tte",
        help="The time to event variable (in the pheno file).",
    )
    group.add_argument(
        "--censure",
        type=str,
        metavar="NAME",
        required=True,
        help="The censure value (1 if observed, 0 if censored)",
    )

    # The linear parser
    lin_parser = subparsers.add_parser(
        "linear",
        help="Linear regression (ordinary least squares).",
        description="Performs a linear regression (ordinary least squares) on "
                    "imputed data. This script is part of the 'gwip' "
                    "package, version {}).".format(__version__),
        parents=[p_parser],
    )

    # The phenotype options for linear regression
    group = lin_parser.add_argument_group("Linear Regression Options")
    group.add_argument(
        "--pheno-name",
        type=str,
        metavar="NAME",
        required=True,
        help="The phenotype.",
    )

    # The logistic parser
    logit_parser = subparsers.add_parser(
        "logistic",
        help="Logistic regression.",
        description="Performs a logistic regression on imputed data. "
                    "This script is part of the 'gwip' package, "
                    "version {}).".format(__version__),
        parents=[p_parser],
    )

    # The phenotype options for logistic regression
    group = logit_parser.add_argument_group("Logistic Regression Options")
    group.add_argument(
        "--pheno-name",
        type=str,
        metavar="NAME",
        required=True,
        help="The phenotype.",
    )

    # The SKAT parser.
    skat_parser = subparsers.add_parser(
        "skat",
        help="SKAT analysis.",
        description="Uses the SKAT R package to analyze user defined gene "
                    "sets. This script is part of the 'gwip' package, "
                    "version {}).".format(__version__),
        parents=[p_parser],
    )

    # Additional options for SKAT analysis.
    group = skat_parser.add_argument_group("SKAT Options")

    # The gene set file.
    group.add_argument(
        "--snp-sets",
        type=str,
        metavar="FILE",
        required=True,
        help="A file indicating a snp_set and an optional weight for every "
             "variant."
    )

    # The outcome type: discrete or continuous.
    group.add_argument(
        "--outcome-type",
        type=str,
        choices=("continuous", "discrete"),
        required=True,
        help="The variable type for the outcome. This will be passed to SKAT."
    )

    # SKAT-O flag.
    group.add_argument(
        "--skat-o",
        action="store_true",
        help="By default, the regular SKAT is used. Setting this flag will "
             "use the SKAT-O algorithm instead."
    )

    group.add_argument(
        "--pheno-name",
        type=str,
        metavar="NAME",
        required=True,
        help="The phenotype.",
    )

    if args is not None:
        return parser.parse_args(args)

    return parser.parse_args()


# Calling the main, if necessary
if __name__ == "__main__":
    main()
