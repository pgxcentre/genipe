
# This file is part of genipe.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.

import re
import os
import sys
import stat
import shlex
import datetime
import logging
import argparse
import platform
import traceback
from shutil import which
from multiprocessing import Pool
from subprocess import Popen, PIPE
from collections import namedtuple

import jinja2
import pandas as pd
from numpy.linalg.linalg import LinAlgError

from .. import __version__
from ..formats import impute2
from ..error import GenipeError


__author__ = ["Louis-Philippe Lemieux Perreault", "Marc-Andre Legault"]
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


# Check if lifelines is installed
try:
    from lifelines import CoxPHFitter
    HAS_LIFELINES = True
except ImportError:
    HAS_LIFELINES = False

# Check if statsmodels is installed
try:
    import statsmodels.api as sm
    import statsmodels.formula.api as smf
    from patsy import dmatrices
    HAS_STATSMODELS = True
except ImportError:
    HAS_STATSMODELS = False


def _has_skat():
    """Checks if the SKAT R library is installed.

    Returns:
        bool: True if SKAT is installed, False otherwise.

    """
    proc = Popen(
        ["Rscript", "-e", 'is.element("SKAT", installed.packages()[,1])'],
        stdout=PIPE,
    )
    out = proc.communicate()[0].decode().strip()

    return out.endswith("TRUE")


# Check if R and SKAT is installed
HAS_R = which("Rscript") is not None
HAS_SKAT = _has_skat() if HAS_R else False


# An IMPUTE2 row to process
_Row = namedtuple("_Row", ("row", "samples", "pheno", "pheno_name", "use_ml",
                           "categorical", "formula", "time_to_event", "event",
                           "inter_c", "is_chrx", "gender_c", "del_g", "scale",
                           "maf_t", "prob_t", "analysis_type",
                           "number_to_print", "random_effects", "mixedlm_p"))

# The Cox's regression required values
_COX_REQ_COLS = ["coef", "se(coef)", "lower 0.95", "upper 0.95", "z", "p"]


def main(args=None):
    """The main function.

    Args:
        args (argparse.Namespace): the arguments to be parsed (if
                                   :py:func:`main` is called by another
                                   modulel)

    """
    # Creating the option parser
    desc = ("Performs statistical analysis on imputed data (either SKAT "
            "analysis, or linear, logistic or survival regression). This "
            "script is part of the 'genipe' package, version "
            "{}.".format(__version__))
    parser = argparse.ArgumentParser(description=desc)

    # Files that need closing
    logging_fh = None

    try:
        # Parsing the options
        args = parse_args(parser, args)

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
        logging.info("Program arguments: {}".format(
            " ".join(shlex.quote(part) for part in sys.argv[1:])
        ))

        # Checking the options
        check_args(args)

        # Reading the phenotype file
        logging.info("Reading phenotype file")
        phenotypes, remove_gender = read_phenotype(
            args.pheno,
            args,
            check_duplicated=args.analysis_type != "mixedlm",
        )
        logging.info("  - {:,d} samples with phenotype".format(
            len(phenotypes.index.unique()),
        ))

        # Reading the sample file
        logging.info("Reading the sample file")
        samples = read_samples(args.sample)
        logging.info("  - {:,d} samples with imputation data".format(
            len(samples),
        ))

        # Reading the sites to extract (if required)
        sites_to_extract = None
        if args.extract_sites is not None:
            logging.info("Reading the sites to extract for analysis")
            sites_to_extract = read_sites_to_extract(args.extract_sites)
            logging.info("  - {:,d} sites to extract".format(
                len(sites_to_extract),
            ))

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

        logging.info("Analysis completed")

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


def read_phenotype(i_filename, opts, check_duplicated=True):
    """Reads the phenotype file.

    Args:
        i_filename (str): the name of the input file
        opts (argparse.Namespace): the options
        check_duplicated (bool): whether or not to check for duplicated samples

    Returns:
        pandas.DataFrame: the phenotypes

    This file is expected to be a tab separated file of phenotypes and
    covariates. The columns to use will be determined by the
    ``--sample-column`` and the ``--covar`` options.

    For analysis including the X chromosome, the gender is automatically
    added as a covariate. The results are not shown to the user unless asked
    for.

    """
    # Reading the data (and setting the index)
    pheno = pd.read_csv(i_filename, sep="\t", na_values=opts.missing_value)
    pheno[opts.sample_column] = pheno[opts.sample_column].astype(str)
    pheno = pheno.set_index(opts.sample_column,
                            verify_integrity=check_duplicated)

    # Finding the required column
    required_columns = opts.covar.copy()
    if opts.analysis_type == "cox":
        required_columns.extend([opts.tte, opts.event])
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
        # Log the number of male/female/missing.
        sex_counts = pheno.groupby(opts.gender_column).size()

        # Computing the number of unknown sex
        unknown = sex_counts.loc[
            ~((sex_counts.index == 1) | (sex_counts.index == 2))
        ]
        nb_unknown = 0
        if unknown.shape[0] > 0:
            nb_unknown = unknown.sum()

        logging.info("  - {:,d} males / {:,d} females ({:,d} unknown)".format(
            sex_counts.get(1, 0), sex_counts.get(2, 0), nb_unknown,
        ))

        pheno = pheno[(pheno[opts.gender_column] == 1) |
                      (pheno[opts.gender_column] == 2)]

    # Extracting the required column
    pheno = pheno.loc[:, required_columns]

    # Returning the phenotypes (dropping the nan values)
    return pheno.dropna(axis=0), remove_gender_column


def read_samples(i_filename):
    """Reads the sample file (produced by SHAPEIT).

    Args:
        i_filename (str): the name of the input file

    Returns:
        pandas.DataFrame: the list of samples

    This file contains the list of samples that are contained in the
    ``impute2`` file (with same order). The expected format for this file is a
    tab separated file with a first row containing the following columns: ::

        ID_1	ID_2	missing	father	mother	sex	plink_pheno

    The subsequent row will be discarded and should contain: ::

        0	0	0 D	D	D	B

    Notes
    -----
        We are mostly interested in the sample IDs corresponding to the
        ``ID_2`` column. Their uniqueness is verified by pandas.

    """
    samples = pd.read_csv(i_filename, sep=" ", usecols=[0, 1])
    samples = samples.drop(samples.index[0], axis=0)
    samples["ID_2"] = samples["ID_2"].astype(str)
    return samples.set_index("ID_2", verify_integrity=True)


def skat_read_snp_set(i_filename):
    """Reads the SKAT SNP set file.

    Args:
        i_filename (str): the name of the input file

    Returns:
        pandas.DataFrame: the SNP set for the SKAT analysis

    This file has to be supplied by the user. The recognized columns are:
    ``variant``, ``snp_set`` and ``weight``. The ``weight`` column is optional
    and can be used to specify a custom weighting scheme for SKAT. If nothing
    is specified, the default Beta weights are used.

    The file has to be tab delimited.

    """
    skat_info = pd.read_csv(i_filename, sep="\t", header=0)
    if "variant" not in skat_info.columns:
        raise GenipeError("The SKAT SNP set file needs to have a 'variant' "
                          "column containing the ID of every variant of "
                          "interest.")

    if "snp_set" not in skat_info.columns:
        raise GenipeError("The SKAT SNP set file needs to have a 'snp_set' "
                          "column containing the SNP set ID for every "
                          "variant. The user is free to choose the SNP ID")

    return skat_info


def read_sites_to_extract(i_filename):
    """Reads the list of sites to extract.

    Args:
        i_filename (str): The input filename containing the IDs of the variants
                          to consider for the analysis.

    Returns:
        set: A set containing the variants.

    The expected file format is simply a list of variants. Every row should
    correspond to a single variant identifier. ::

        3:60069:t
        rs123456:A
        3:60322:A

    Typically, this is used to analyze only variants that passed some QC
    threshold. The :py:mod:`genipe` pipeline generates this file at the
    'merge_impute2' step.

    """
    markers_to_extract = None
    with open(i_filename, "r") as i_file:
        markers_to_extract = set(i_file.read().splitlines())
    return markers_to_extract


def skat_parse_impute2(impute2_filename, samples, markers_to_extract,
                       phenotypes, remove_gender, out_prefix, args):
    """Read the impute2 file and run the SKAT analysis.

    Args:
        impute2_filename (str): the name of the input file
        samples (pandas.DataFrame): the samples
        markers_to_extract (set): the set of markers to analyse
        phenotypes (pandas.DataFrame): the phenotypes
        remove_gender (bool): whether or not to remove the gender column
        out_prefix (str): the output prefix
        args (argparse.Namespace): the options


    This function does most of the "dispatching" to run SKAT. It writes the
    input files to the disk, runs the generated R scripts to do the actual
    analysis and then writes the results to disk.

    """

    # We keep track of the files that are generated because we will need the
    # paths to generate the R script correctly.
    r_files = {"snp_sets": [], "covariates": None, "outcome": None,
               "weights": None}

    # Read the SNP set and create the output files.
    snp_set = skat_read_snp_set(args.snp_sets)

    # We will use a temporary directory for our analysis.
    dir_name = "{}.{}".format(
        args.out,
        datetime.datetime.today().strftime("skat.%Y.%m.%d")
    )

    dir_name = os.path.abspath(dir_name)

    if not os.path.isdir(dir_name):
        os.makedirs(dir_name)
    else:
        # This should not happen because the folder name contains a timestamp
        raise GenipeError("A folder named '{}' already exists.".format(
            dir_name
        ))

    # If weights were provided, we will write a weight vector to disk for R.
    if "weight" in snp_set.columns:
        logging.info("SKAT will use the provided weights (from the SNP set "
                     "file).")
        weight_filename = os.path.join(dir_name, "weights.csv")
        snp_set[["weight"]].to_csv(weight_filename, index=False, header=False)
        r_files["weights"] = weight_filename

    # Open genotype CSV files for every SNP set. Those CSV files will be
    # read in R and used by SKAT.
    genotype_files = {}
    snp_sets = snp_set["snp_set"].unique()
    for set_id in snp_sets:
        filename = os.path.join(dir_name, "{}.genotypes.csv".format(set_id))
        r_files["snp_sets"].append(filename)  # Track the filenames.
        genotype_files[set_id] = open(filename, "w")
        # We write the column headers for every file.
        print("", *samples.index, sep=",", file=genotype_files[set_id])

    # The markers of interest are the markers we want to include in the
    # analysis. Concretely, they are the markers that were included in the snp
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

    # We will warn the user if some variants were not found.
    number_missing = len(markers_of_interest) - len(written_markers)

    if number_missing > 0:
        logging.warning("{} markers of interest were not found in the Impute2 "
                        "file.".format(number_missing))

    i_file.close()

    # Close the genotype files.
    for file_handle in genotype_files.values():
        file_handle.close()

    # Write the covariate file.
    phenotype_df = samples.join(phenotypes)
    phenotype_df.index.name = "sample"

    if args.covar:
        # Make sure the samples are consistent by merging phenotype and
        # samples.
        filename = os.path.join(dir_name, "covariates.csv")
        phenotype_df[args.covar].to_csv(
            filename,
            sep=",",
        )
        r_files["covariates"] = filename

    # Write the phenotype file.
    filename = os.path.join(
        dir_name,
        "{}.{}.csv".format(args.pheno_name, args.outcome_type)
    )
    phenotype_df[[args.pheno_name]].to_csv(
        filename,
        sep=",",
    )
    r_files["outcome"] = filename

    r_scripts = _skat_generate_r_script(dir_name, r_files, args)

    # Run the SKAT analysis by calling Rscript either in different subprocesses
    # or linearly.
    logging.info("Launching SKAT using {} processes on {} SNP sets.".format(
        args.nb_process, len(snp_sets)
    ))

    results = []
    if args.nb_process > 1:
        with Pool(processes=args.nb_process) as pool:
            results = pool.map(_skat_run_job, r_scripts)
    else:
        for script in r_scripts:
            results.append(_skat_run_job(script))

    # Finally, write the SKAT output to disk.
    output_filename = args.out + ".skat.dosage"

    logging.info("SKAT completed, writing the results to disk ({}).".format(
        output_filename
    ))

    with open(output_filename, "w") as f:
        # Write the header.
        if args.skat_o:
            print("snp_set_id", "p_value", sep="\t", file=f)
        else:
            print("snp_set_id", "p_value", "q_value", sep="\t", file=f)

        # The order in the snp_sets list should have been consistent
        # across the whole analysis. We just want to make sure that we
        # actually have the right number of results.
        assert len(snp_sets) == len(results)

        # Write a tab separated file containing the set and p-values.
        for i, (p_value, q_value) in enumerate(results):
            set_id = snp_sets[i]
            if args.skat_o:
                print(set_id, p_value, sep="\t", file=f)
            else:
                print(set_id, p_value, q_value, sep="\t", file=f)


def _skat_run_job(script_filename):
    """Calls Rscript with the generated script and parses the results.

    Args:
        script_filename (str): the name of the script

    Returns:
        tuple: two values: the *p-value* and the *q-value* (for SKAT-O, the
               *q-value* is set to None)

    The results should be somewhere in the standard output. The expected
    format is: ::

        _PYTHON_HOOK_QVAL:[0.123]
        _PYTHON_HOOK_PVAL:[0.123]

    If the template script is modified, this format should still be respected.

    It is also noteworthy that this function uses ``Rscript`` to run the
    analysis. Hence, it should be in the path when using the imputed_stats
    skat mode.

    """
    # Parse the p-value.
    proc = Popen(
        ["Rscript", script_filename],
        stdout=PIPE,
        stderr=PIPE,
    )
    out, err = proc.communicate()
    if err:
        logging.info("SKAT Warning: " + err.decode("utf-8"))

    # Decoding
    out = out.decode("utf-8")

    # Getting the p value
    p_match = re.search(r"_PYTHON_HOOK_PVAL:\[(.+)\]", out)
    if p_match is None:
        raise GenipeError("SKAT did not return properly. See script "
                          "'{}' for details.".format(script_filename))

    # Getting the Q value
    q_match = re.search(r"_PYTHON_HOOK_QVAL:\[(.+)\]", out)
    if q_match is None:
        raise GenipeError("SKAT did not return properly. See script "
                          "'{}' for details.".format(script_filename))

    if q_match.group(1) == "NA":
        # For SKAT-O, the Q will be set to NA.
        return float(p_match.group(1)), None
    else:
        return float(p_match.group(1)), float(q_match.group(1))


def _skat_generate_r_script(dir_name, r_files, args):
    """Uses jinja2 to generate an R script to do the SKAT analysis.

    Args:
        dir_name (str): the output directory name to write the scripts in
        r_files (dict): contains the different input files required by the R
                        script
        args (argparse.Namespace): the parsed arguments

    Returns:
        list: the list of all script names

    """
    jinja_env = jinja2.Environment(
        loader=jinja2.PackageLoader("genipe", "script_templates")
    )
    template = jinja_env.get_template("run_skat.R")

    scripts = []

    # Create one file per SNP set for parallelism.
    for snp_set_file in r_files["snp_sets"]:
        rendered_script = template.render(
            version=__version__,

            snp_set_file=snp_set_file,
            covariate_file=r_files["covariates"],
            outcome_file=r_files["outcome"],
            weights=r_files["weights"],

            outcome_type="C" if args.outcome_type == "continuous" else "D",
            skat_o=args.skat_o,
        )
        # Write the rendered script to disk.
        script_filename = "run_skat_{}R".format(
            # We parse the set id name.
            os.path.basename(snp_set_file)[:-len(".genotype.csv")]
        )
        script_filename = os.path.join(dir_name, script_filename)

        with open(script_filename, "w") as f:
            f.write(rendered_script)

        scripts.append(script_filename)

    return scripts


def _skat_parse_line(line, markers_of_interest, samples, gender=None):
    """Parses a single line of the Impute2 file.

    Args:
        line (str): a line from the Impute2 file
        markers_of_interest (set): a set of markers that are required for the
                                   analysis
        samples (pandas.DataFrame): contains the samples IDs (this is useful to
                                    make sure we return a dosage vector with
                                    the appropriate data)

    Returns:
        tuple: Either None if the marker is not of interest or a tuple of
               ``(name , dosage_vector)`` where ``name`` is a ``str``
               representing the variant ID and ``dosage_vector`` is a numpy
               array containing the dosage values for every sample in the
               ``samples`` dataframe.

    """
    line = line.split(" ")
    # info_tuple contains: chrom, name, pos, a1, a2
    # proba_matrix is a matrix of sample x (aa, ab, bb)
    info_tuple, proba_matrix = impute2.matrix_from_line(line)

    chrom, name, pos, a1, a2 = info_tuple

    # If this marker is not of interest, we don't bother continuing.
    if name not in markers_of_interest:
        return None

    # We need to compute the dosage vector.
    # This is given wrt to the minor and major alleles, so we need to get the
    # MAF to identify those.
    maf, minor_i, major_i = impute2.maf_from_probs(
        prob_matrix=proba_matrix,
        a1=0,
        a2=2,
        gender=gender,
        site_name=name,
    )

    # We don't pass a scale parameter, because we want additive coding.
    dosage = impute2.dosage_from_probs(
        homo_probs=proba_matrix[:, minor_i],
        hetero_probs=proba_matrix[:, 1],
    )

    return (name, dosage)


def _skat_write_marker(name, dosage, snp_set, genotype_files):
    """Write the dosage information to the appropriate genotype file.

    Args:
        name (str): the name of the marker
        dosage (numpy.array): the dosage vector
        snp_set (pandas.DataFrame: the dataframe that allows us to identify the
                                   correct SNP set for the specified variant
        genotype_files (dict): a dictionary containing the opened CSV files for
                               the genotypes

    """
    # Identify the correct SNP set.
    this_snp_set = snp_set.loc[snp_set["variant"] == name, "snp_set"].unique()
    for set_id in this_snp_set:
        file_object = genotype_files[set_id]
        print(name, *dosage, sep=",", file=file_object)


def _extract_mixedlm_random_effect(fitted):
    """Extracts the random effects from a MixedLM fit object.

    Args:
        fitted (MixedLMResultsWrapper): The fitted object.

    Returns:
        pandas.DataFrame: The random effects as a DataFrame (with a column
                          named "RE").

    Note
    ====
        Depending of the version of StatsModels, the object might be a pandas
        DataFrame or a dictionary...

    """
    # Getting the random effects
    random_effects = fitted.random_effects

    # If it's a dictionary, we need to create a DataFrame
    if isinstance(random_effects, dict):
        return pd.DataFrame(random_effects).T.rename(
            columns={"groups": "RE", "Group": "RE"},
        )

    return random_effects.rename(columns={"Intercept": "RE"})


def compute_statistics(impute2_filename, samples, markers_to_extract,
                       phenotypes, remove_gender, out_prefix, options):
    """Parses IMPUTE2 file while computing statistics.

    Args:
        impute2_filename (str): the name of the input file
        samples (pandas.DataFrame): the list of samples
        markers_to_extract (set): the set of markers to extract
        phenotypes (pandas.DataFrame): the phenotypes
        remove_gender (bool): whether or not to remove the gender column
        out_prefix (str): the output prefix
        options (argparse.Namespace): the options

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
            gender_c=options.gender_column,
            categorical=options.categorical,
        )
        logging.info("{}: '{}'".format(options.analysis_type, formula))
        if options.analysis_type == "mixedlm" and options.use_ml:
            logging.info("  - using ML")

    else:
        # This is Cox
        formula = get_formula(
            phenotype=options.tte + " + " + options.event,
            covars=options.covar,
            interaction=options.interaction,
            gender_c=options.gender_column,
            categorical=options.categorical,
        )

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
        if options.analysis_type == "linear":
            header = header + ("adj.r-squared", )
        if options.analysis_type == "mixedlm":
            header = header + ("type", )

        print(*header, sep="\t", file=o_file)

        # The sites to process (if multiprocessing)
        sites_to_process = []

        # We need to compute the random effects if it's a MixedLM analysis
        random_effects = None
        if options.analysis_type == "mixedlm" and options.interaction is None:
            random_effects = _extract_mixedlm_random_effect(smf.mixedlm(
                formula=formula.replace("_GenoD + ", ""),
                data=phenotypes,
                groups=phenotypes.index,
            ).fit(reml=not options.use_ml))

        # Reading the file
        nb_processed = 0
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
                use_ml=vars(options).get("use_ml", None),
                pheno_name=vars(options).get("pheno_name", None),
                formula=formula,
                time_to_event=vars(options).get("tte", None),
                event=vars(options).get("event", None),
                inter_c=options.interaction,
                is_chrx=options.chrx,
                gender_c=options.gender_column,
                del_g=remove_gender,
                scale=options.scale,
                maf_t=options.maf,
                prob_t=options.prob,
                analysis_type=options.analysis_type,
                number_to_print=len(header),
                categorical=options.categorical,
                random_effects=random_effects,
                mixedlm_p=vars(options).get("p_threshold", None),
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

                    # Logging
                    nb_processed += options.nb_lines
                    logging.info("Processed {:,d} lines".format(nb_processed))

                    # Resetting the sites to process
                    sites_to_process = []

            else:
                # Processing this row
                print(*process_impute2_site(site), sep="\t", file=o_file)

        if len(sites_to_process) > 0:
            for result in pool.map(process_impute2_site, sites_to_process):
                print(*result, sep="\t", file=o_file)

            # Logging
            nb_processed += len(sites_to_process)
            logging.info("Processed {:,d} lines".format(nb_processed))

    except Exception:
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
                raise GenipeError("{}: problem while reading the GZ "
                                  "file".format(impute2_filename))


def process_impute2_site(site_info):
    """Process an IMPUTE2 site (a line in an IMPUTE2 file).

    Args:
        site_info (list): the impute2 line (split by space)

    Returns:
        list: the results of the analysis

    """
    # Getting the probability matrix and site information
    (chrom, name, pos, a1, a2), geno = impute2.matrix_from_line(site_info.row)

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
        impute2.get_good_probs(data[dosage_columns].values, site_info.prob_t)
    ]

    # FIXME: Quick and dirty fix for mixedlm...
    # If the analysis type is not MixedLM, then t_data is just a pointer to the
    # real data. Otherwise, t_data is the grouped values (first occurrence,
    # which shouldn't cause a problem, since gender and probabilities are the
    # same for each measurement).
    t_data = data
    if site_info.analysis_type == "mixedlm":
        t_data = data.groupby(level=0).first()

    # Checking gender if required
    gender = None
    if site_info.is_chrx:
        # We want to exclude males with heterozygous calls for the rest of the
        # analysis
        invalid_rows = samples_with_hetero_calls(
            t_data.loc[t_data[site_info.gender_c] == 1, dosage_columns],
            dosage_columns[1]
        )
        if len(invalid_rows) > 0:
            logging.warning("There were {:,d} males with heterozygous "
                            "calls for {}".format(len(invalid_rows), name))
            logging.debug(t_data.shape)
            data = data.drop(invalid_rows, axis=0)
            t_data = t_data.drop(invalid_rows, axis=0)
            logging.debug(t_data.shape)

        # Getting the genders
        gender = t_data[site_info.gender_c].values

    # Computing the frequency
    maf, minor, major = impute2.maf_from_probs(
        prob_matrix=t_data[dosage_columns].values,
        a1=dosage_columns[0],
        a2=dosage_columns[-1],
        gender=gender,
        site_name=name,
    )

    # What we want to print
    to_return = [chrom, pos, name, allele_encoding[major],
                 allele_encoding[minor], maf, t_data.shape[0]]

    # If the marker is too rare, we continue with the rest
    if (maf == "NA") or (maf < site_info.maf_t):
        to_return.extend(["NA"] * (site_info.number_to_print - len(to_return)))
        return to_return

    # Computing the dosage on the minor allele
    data["_GenoD"] = impute2.dosage_from_probs(
        homo_probs=data[minor],
        hetero_probs=data[dosage_columns[1]],
        scale=site_info.scale,
    )

    # Removing the unwanted columns
    unwanted_columns = dosage_columns
    if site_info.del_g:
        unwanted_columns.append(site_info.gender_c)
    data = data.drop(unwanted_columns, axis=1)

    # The column to get the result from
    result_from_column = "_GenoD"
    if site_info.inter_c is not None:
        if ((site_info.inter_c == site_info.gender_c) or
                (site_info.inter_c in site_info.categorical)):
            result_from_column = "_GenoD:C({})[T.{}]".format(
                site_info.inter_c,
                sorted(data[site_info.inter_c].unique())[-1],
            )
        else:
            result_from_column = "_GenoD:" + site_info.inter_c

    # Fitting
    results = []
    try:
        results = _fit_map[site_info.analysis_type](
            data=data,
            groups=data.index.values,
            time_to_event=site_info.time_to_event,
            event=site_info.event,
            formula=site_info.formula,
            result_col=result_from_column,
            use_ml=site_info.use_ml,
            random_effects=site_info.random_effects,
            mixedlm_p=site_info.mixedlm_p,
            interaction=site_info.inter_c is not None,
        )
    except LinAlgError as e:
        # Something strange happened...
        logging.warning("{}: numpy LinAlgError: {}".format(name, str(e)))

    # Extending the list to return
    if len(results) == 0:
        results = ["NA"] * (site_info.number_to_print - len(to_return))
    to_return.extend(results)

    return to_return


def samples_with_hetero_calls(data, hetero_c):
    """Gets male and heterozygous calls.

    Args:
        data (pandas.DataFrame): the probability matrix
        hetero_c (str): the name of the heterozygous column

    Returns:
        pandas.Index: samples where call is heterozygous

    Note
    ----
        If there are no data (i.e. no males), an empty list is returned.

    """
    if data.shape[0] == 0:
        return []
    return data[data.idxmax(axis=1) == hetero_c].index


def get_formula(phenotype, covars, interaction, gender_c, categorical):
    """Creates the linear/logistic regression formula (for statsmodel).

    Args:
        phenotype (str): the phenotype column
        covars (list): the list of co variable columns
        interaction (str): the interaction column

    Returns:
        str: the formula for the statistical analysis

    Note
    ----
        The phenotype column needs to be specified. The list of co variables
        might be empty (if no co variables are necessary). The interaction
        column can be set to ``None`` if there is no interaction.

    Note
    ----
        The gender column should be categorical (hence, the formula requires
        the gender to be included into ``C()``, *e.g.* ``C(Gender)``).

    """
    # The phenotype and genetics
    formula = "{} ~ _GenoD".format(phenotype)

    # Are there any covars?
    for covar in covars:
        if (covar == gender_c) or (covar in categorical):
            covar = "C({})".format(covar)
        formula += " + " + covar

    # Is there an interaction term?
    if interaction is not None:
        if (interaction == gender_c) or (interaction in categorical):
            interaction = "C({})".format(interaction)
        formula += " + _GenoD*{}".format(interaction)

    return formula


def fit_cox(data, time_to_event, event, formula, result_col, **kwargs):
    """Fit a Cox' proportional hazard to the data.

    Args:
        data (pandas.DataFrame): the data to analyse
        time_to_event (str): the time to event column for the survival analysis
        event (str): the event column for the survival analysis
        formula (str): the formula for the data preparation
        result_col (str): the column that will contain the results

    Returns:
        numpy.array: the results from the survival analysis

    Note
    ----
        The tie method used is ``Efron``. Normalization is set to ``False``.

    """
    # Preparing the data using Patsy
    y, X = dmatrices(formula, data=data, return_type="dataframe")
    data = pd.merge(y, X.drop("Intercept", axis=1), left_index=True,
                    right_index=True)

    # Fitting
    cf = CoxPHFitter(alpha=0.95, tie_method="Efron")
    cf.fit(data, duration_col=time_to_event, event_col=event)
    return cf.summary.loc[result_col, _COX_REQ_COLS].values


def fit_linear(data, formula, result_col, **kwargs):
    """Fit a linear regression to the data.

    Args:
        data (pandas.DataFrame): the data to analyse
        formula (str): the formula for the linear regression
        result_col (str): the column that will contain the results

    Returns:
        list: the results from the linear regression

    """
    return _get_result_from_linear(
        smf.ols(formula=formula, data=data).fit(),
        result_col=result_col,
    )


def fit_logistic(data, formula, result_col, **kwargs):
    """Fit a logistic regression to the data.

    Args:
        data (pandas.DataFrame): the data to analyse
        formula (str): the formula for the logistic regression
        result_col (str): the column that will contain the results

    Returns:
        list: the results from the logistic regression

    """
    return _get_result_from_logistic_mixedlm(
        smf.glm(formula=formula, data=data,
                family=sm.families.Binomial()).fit(),
        result_col=result_col,
    )


def fit_mixedlm(data, formula, use_ml, groups, result_col, random_effects,
                mixedlm_p, interaction, **kwargs):
    """Fit a linear mixed effects model to the data.

    Args:
        data (pandas.DataFrame): the data to analyse
        formula (str): the formula for the linear mixed effects model
        use_ml (bool): whether to use ML instead of REML
        groups (str): the column containing the groups
        result_col (str): the column that will contain the results
        random_effects (pandas.Series): the random effects
        mixedlm_p (float): the p-value threshold for which loci will be
                           computed with the real MixedLM analysis
        interaction (bool): Whether there is an interaction or not

    Returns:
        list: the results from the linear mixed effects model

    """
    # We perform the optimization if there is no interaction
    if not interaction:
        # Getting the single copy of each genotypes
        geno = data.reset_index()[["index", "_GenoD"]]
        indexes = geno[["index"]].drop_duplicates().index
        geno = geno.loc[indexes, :]

        # Merging with the random effects
        t_data = pd.merge(random_effects, geno.set_index("index"),
                          left_index=True, right_index=True)

        # Approximating the results
        approximate_r = _get_result_from_linear(
            smf.ols(formula="RE ~ _GenoD", data=t_data).fit(),
            result_col="_GenoD",
        )

        # If the approximated p-value is higher or equal to the threshold, we
        # return the approximation
        if approximate_r[5] >= mixedlm_p:
            result = ["NA"] * (len(approximate_r) - 2)
            result.append(approximate_r[5])
            result.append("TS-MixedLM")
            return result

    # If we get here, it's because the p-value was low enough so we compute the
    # real statistics
    result = _get_result_from_logistic_mixedlm(
        smf.mixedlm(formula=formula, data=data,
                    groups=groups).fit(reml=not use_ml),
        result_col=result_col,
    )
    result.append("MixedLM")
    return result


_fit_map = {
    "cox": fit_cox,
    "linear": fit_linear,
    "logistic": fit_logistic,
    "mixedlm": fit_mixedlm,
}


def _get_result_from_linear(fit_result, result_col):
    """Gets results from either a linear, a logistic or a mixedlm regression.

    Args:
        fit_result (RegressionResults): the results from the regression
        result_col (str): the name of the result column

    Returns:
        list: the regression results

    """
    conf_int = fit_result.conf_int().loc[result_col, :].values
    assert len(conf_int) == 2
    return [
        fit_result.params[result_col],
        fit_result.bse[result_col],
        conf_int[0],
        conf_int[1],
        fit_result.tvalues[result_col],
        fit_result.pvalues[result_col],
        fit_result.rsquared_adj,
    ]


def _get_result_from_logistic_mixedlm(fit_result, result_col):
    """Gets results from either a linear, a logistic or a mixedlm regression.

    Args:
        fit_result (RegressionResults): the results from the regression
        result_col (str): the name of the result column

    Returns:
        list: the regression results

    """
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


def is_file_like(fn):
    """Checks if the path is like a file (it might be a named pipe).

    Args:
        fn (str): the path to check

    Returns:
        bool: True if path is like a file, False otherwise.

    """
    return os.path.isfile(fn) or stat.S_ISFIFO(os.stat(fn).st_mode)


def check_args(args):
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the options to verify

    Note
    ----
        If there is a problem, a :py:class:`genipe.error.GenipeError` is
        raised.

    """
    # Checking if we have the requirements
    if args.analysis_type == "cox":
        if not HAS_LIFELINES:
            raise GenipeError("missing optional module: lifelines")

    elif args.analysis_type in {"linear", "logistic", "mixedlm"}:
        if not HAS_STATSMODELS:
            raise GenipeError("missing optional module: statsmodels")

    elif args.analysis_type == "skat":
        if not HAS_R:
            raise GenipeError("R is not installed")
        if not HAS_SKAT:
            raise GenipeError("R library missing: SKAT")

    # Checking the required input files
    for filename in [args.impute2, args.sample, args.pheno]:
        if not os.path.isfile(filename):
            raise GenipeError("{}: no such file".format(filename))

    # Checking the optional input files
    for filename in [args.extract_sites]:
        if filename is not None:
            if not is_file_like(filename):
                raise GenipeError("{}: no such file".format(filename))

    # Checking the number of process
    on_mac_os = platform.system() == "Darwin"
    if args.nb_process < 1:
        raise GenipeError("{}: invalid number of "
                          "processes".format(args.nb_process))
    if (args.nb_process > 1 and on_mac_os and args.analysis_type != "skat"):
        raise GenipeError("multiprocessing is not supported on Mac OS when "
                          "using linear regression, logistic regression or "
                          "Cox.")

    # Checking the number of lines to read
    if args.nb_lines < 1:
        raise GenipeError("{}: invalid number of lines to "
                          "read".format(args.nb_lines))

    # Checking the MAF threshold
    if args.maf < 0 or args.maf > 1:
        raise GenipeError("{}: invalid MAF".format(args.maf))

    # Checking the probability threshold
    if args.prob < 0 or args.prob > 1:
        raise GenipeError("{}: invalid probability "
                          "threshold".format(args.prob))

    # Reading all the variables in the phenotype file
    header = None
    with open(args.pheno, "r") as i_file:
        header = {name for name in i_file.readline().rstrip("\n").split("\t")}

    # Checking the restricted columns
    restricted_columns = {"_D1", "_D2", "_D3", "_MaxD", "_GenoD", "_Inter"}
    if len(restricted_columns & header) != 0:
        raise GenipeError("{}: {}: restricted variables".format(
            args.pheno,
            restricted_columns & header,
        ))

    # Checking the required columns
    variables_to_check = None
    if args.analysis_type == "cox":
        variables_to_check = {args.tte, args.event}
    else:
        variables_to_check = {args.pheno_name}
    for variable in variables_to_check:
        if variable not in header:
            raise GenipeError("{}: {}: missing variable for {}".format(
                args.pheno,
                variable,
                args.analysis_type,
            ))

    # The categorical variables
    categorical_set = set()
    if args.categorical != "":
        categorical_set = set(args.categorical.split(","))
    for name in categorical_set:
        if name not in header:
            raise GenipeError("{}: {}: missing categorical value".format(
                args.pheno,
                name,
            ))
        if args.analysis_type != "cox":
            if args.pheno_name in categorical_set:
                raise GenipeError("{}: {}: should not be in categorical "
                                  "list".format(args.pheno, args.pheno_name))
    args.categorical = categorical_set

    # Checking the co-variables
    covar_list = []
    if args.covar != "":
        covar_list = args.covar.split(",")
    for covar in covar_list:
        if covar not in header:
            raise GenipeError("{}: {}: missing co-variable".format(
                args.pheno,
                covar,
            ))
    args.covar = covar_list

    # Checking the sample column
    if args.sample_column not in header:
        raise GenipeError("{}: {}: no such column (--sample-column)".format(
            args.pheno,
            args.sample_column,
        ))

    # Checking the gender column (only if required)
    if args.gender_column != "None":
        if args.gender_column not in header:
            raise GenipeError(
                "{}: {}: no such column (--gender-column)".format(
                    args.pheno,
                    args.gender_column,
                )
            )

    # Checking the interaction column (if required)
    if args.interaction is not None:
        if args.interaction not in header:
            raise GenipeError(
                "{}: {}: no such column (--interaction)".format(
                    args.pheno,
                    args.interaction,
                )
            )
        if args.interaction in args.categorical:
            logging.warning("interaction term is categorical: the last "
                            "category will be used in the results")

        if args.analysis_type == "mixedlm":
            logging.warning("when using interaction, mixedlm optimization "
                            "cannot be performed, analysis will be slow")

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

    # Creating a parent parser (for common options between analysis type)
    p_parser = argparse.ArgumentParser(add_help=False)

    # The parser object
    p_parser.add_argument(
        "-v", "--version", action="version",
        version="%(prog)s (part of genipe version {})".format(__version__),
    )
    p_parser.add_argument(
        "--debug", action="store_true",
        help="set the logging level to debug",
    )

    # The input files
    group = p_parser.add_argument_group("Input Files")
    group.add_argument(
        "--impute2", type=str, metavar="FILE", required=True,
        help="The output from IMPUTE2.",
    )
    group.add_argument(
        "--sample", type=str, metavar="FILE", required=True,
        help="The sample file (the order should be the same as in the IMPUTE2 "
             "files).",
    )
    group.add_argument(
        "--pheno", type=str, metavar="FILE", required=True,
        help="The file containing phenotypes and co variables.",
    )
    group.add_argument(
        "--extract-sites", type=str, metavar="FILE",
        help="A list of sites to extract for analysis (optional).",
    )

    # The output files
    group = p_parser.add_argument_group("Output Options")
    group.add_argument(
        "--out", metavar="FILE", default="imputed_stats",
        help="The prefix for the output files. [%(default)s]",
    )

    # General options
    group = p_parser.add_argument_group("General Options")
    group.add_argument(
        "--nb-process", type=int, metavar="INT", default=1,
        help="The number of process to use. [%(default)d]",
    )
    group.add_argument(
        "--nb-lines", type=int, metavar="INT", default=1000,
        help="The number of line to read at a time. [%(default)d]",
    )
    group.add_argument(
        "--chrx", action="store_true",
        help="The analysis is performed for the non pseudo-autosomal region "
             "of the chromosome X (male dosage will be divided by 2 to get "
             "values [0, 0.5] instead of [0, 1]) (males are coded as 1 and "
             "option '--gender-column' should be used).",
    )
    group.add_argument(
        "--gender-column", type=str, metavar="NAME", default="Gender",
        help="The name of the gender column (use to exclude samples with "
             "unknown gender (i.e. not 1, male, or 2, female). If gender not "
             "available, use 'None'. [%(default)s]",
    )

    # The dosage options
    group = p_parser.add_argument_group("Dosage Options")
    group.add_argument(
        "--scale", type=int, metavar="INT", default=2, choices=[1, 2],
        help="Scale dosage so that values are in [0, n] (possible values are "
             "1 (no scaling) or 2). [%(default)d]",
    )
    group.add_argument(
        "--prob", type=float, metavar="FLOAT", default=0.9,
        help="The minimal probability for which a genotype should be "
             "considered. [>=%(default).1f]",
    )
    group.add_argument(
        "--maf", type=float, metavar="FLOAT", default=0.01,
        help="Minor allele frequency threshold for which marker will be "
             "skipped. [<%(default).2f]",
    )

    # The general phenotype options
    group = p_parser.add_argument_group("Phenotype Options")
    group.add_argument(
        "--covar", type=str, metavar="NAME", default="",
        help="The co variable names (in the phenotype file), separated by "
             "coma.",
    )
    group.add_argument(
        "--categorical", type=str, metavar="NAME", default="",
        help="The name of the variables that are categorical (note that the "
             "gender is always categorical). The variables are separated by "
             "coma.",
    )
    group.add_argument(
        "--missing-value", type=str, metavar="NAME",
        help="The missing value in the phenotype file.",
    )
    group.add_argument(
        "--sample-column", type=str, metavar="NAME", default="sample_id",
        help="The name of the sample ID column (in the phenotype file). "
             "[%(default)s]",
    )
    group.add_argument(
        "--interaction", type=str, metavar="NAME",
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
        help="Cox's proportional hazard model (survival regression).",
        description="Performs a survival regression on imputed data using "
                    "Cox's proportional hazard model. This script is part of "
                    "the 'genipe' package, version {}.".format(__version__),
        parents=[p_parser],
    )

    group = cox_parser.add_argument_group("Cox's Proportional Hazard Model "
                                          "Options")
    group.add_argument(
        "--time-to-event", type=str, metavar="NAME", required=True, dest="tte",
        help="The time to event variable (in the pheno file).",
    )
    group.add_argument(
        "--event", type=str, metavar="NAME", required=True,
        help="The event variable (1 if observed, 0 if not observed)",
    )

    # The linear parser
    lin_parser = subparsers.add_parser(
        "linear",
        help="Linear regression (ordinary least squares).",
        description="Performs a linear regression (ordinary least squares) on "
                    "imputed data. This script is part of the 'genipe' "
                    "package, version {}.".format(__version__),
        parents=[p_parser],
    )

    # The phenotype options for linear regression
    group = lin_parser.add_argument_group("Linear Regression Options")
    group.add_argument(
        "--pheno-name", type=str, metavar="NAME", required=True,
        help="The phenotype.",
    )

    # The logistic parser
    logit_parser = subparsers.add_parser(
        "logistic",
        help="Logistic regression (GLM with binomial distribution).",
        description="Performs a logistic regression on imputed data using a "
                    "GLM with a binomial distribution. This script is part of "
                    "the 'genipe' package, version {}.".format(__version__),
        parents=[p_parser],
    )

    # The phenotype options for logistic regression
    group = logit_parser.add_argument_group("Logistic Regression Options")
    group.add_argument(
        "--pheno-name", type=str, metavar="NAME", required=True,
        help="The phenotype.",
    )

    # The mixed effect model parser
    mixedlm_parser = subparsers.add_parser(
        "mixedlm",
        help="Linear mixed effect model (random intercept).",
        description="Performs a linear mixed effects regression on imputed "
                    "data using a random intercept for each group. A p-value "
                    "approximation is performed so that computation time is "
                    "acceptable for imputed data. This script is part of the "
                    "'genipe' package, version {}.".format(__version__),
        parents=[p_parser],
    )

    # The phenotype options for the mixed effect model
    group = mixedlm_parser.add_argument_group("Linear Mixed Effects Options")
    group.add_argument(
        "--pheno-name", type=str, metavar="NAME", required=True,
        help="The phenotype.",
    )

    group.add_argument(
        "--use-ml", action="store_true",
        help="Fit the standard likelihood using maximum likelihood (ML) "
             "estimation instead of REML (default is REML).",
    )

    group.add_argument(
        "--p-threshold", type=float, metavar="FLOAT", default=1e-4,
        help="The p-value threshold for which the real MixedLM analysis will "
             "be performed. [<%(default).4f]",
    )

    # The SKAT parser.
    skat_parser = subparsers.add_parser(
        "skat",
        help="SKAT analysis.",
        description="Uses the SKAT R package to analyze user defined gene "
                    "sets. This script is part of the 'genipe' package, "
                    "version {}.".format(__version__),
        parents=[p_parser],
    )

    # Additional options for SKAT analysis.
    group = skat_parser.add_argument_group("SKAT Options")

    # The SNP set file.
    group.add_argument(
        "--snp-sets", type=str, metavar="FILE", required=True,
        help="A file indicating a snp_set and an optional weight for every "
             "variant."
    )

    # The outcome type: discrete or continuous.
    group.add_argument(
        "--outcome-type", type=str, choices=("continuous", "discrete"),
        required=True,
        help="The variable type for the outcome. This will be passed to SKAT."
    )

    # SKAT-O flag.
    group.add_argument(
        "--skat-o", action="store_true",
        help="By default, the regular SKAT is used. Setting this flag will "
             "use the SKAT-O algorithm instead."
    )

    group.add_argument(
        "--pheno-name", type=str, metavar="NAME", required=True,
        help="The phenotype.",
    )

    if args is not None:
        return parser.parse_args(args)

    return parser.parse_args()


# Calling the main, if necessary
if __name__ == "__main__":
    main()
