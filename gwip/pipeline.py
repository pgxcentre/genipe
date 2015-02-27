
# This file is part of gwip.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import os
import re
import sys
import json
import logging
import argparse
from glob import glob
from math import floor
from shutil import which
from urllib.request import urlopen
from subprocess import Popen, PIPE
from collections import defaultdict

import pandas as pd

from .db import *
from . import __version__
from . import chromosomes
from .task import launcher
from .error import ProgramError
from .config import parse_drmaa_config
from .reporting import generate_report


try:
    import matplotlib as mpl
    mpl.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


def main():
    """The main function.

    This is what the pipeline should do:
        1- Exclude markers that have ambiguous alleles [A/T or C/G] (Plink)
        2- Exclude duplicated markers (only keep one) (Plink)
        3- Split the dataset by chromosome (Plink)
        4- Find markers that need to be flip (strand problem) (SHAPEIT)
        5- Flip those markers (Plink)
        6- Find markers with strand problem (SHAPEIT)
        7- Exclude markers with strand problem (Plink)
        8- Phase using SHAPEIT
        9- Impute using IMPUTE2
        10- Merge IMPUTE2 files

    """
    # Creating the option parser
    desc = ("Execute the genome-wide imputation pipeline "
            "(gwip version {}).".format(__version__))
    parser = argparse.ArgumentParser(description=desc)

    # We run the script
    try:
        # Parsing the options
        args = parse_args(parser)

        # Creating the output directory if it doesn't exist
        if not os.path.isdir(args.out_dir):
            os.mkdir(args.out_dir)

        # Adding the logging capability
        log_file = os.path.join(args.out_dir, "gwip.log")
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

        # Getting the task options
        args.task_options = None
        if args.use_drmaa:
            logging.info("Task will be launched with DRMAA")
            args.task_options = parse_drmaa_config(args.drmaa_config)
            args.preamble = read_preamble(args.preamble)

        # Creating the database
        db_name = create_task_db(args.out_dir)

        # Creating the output directories
        for chrom in chromosomes:
            chr_dir = os.path.join(args.out_dir, "chr{}".format(chrom))
            if not os.path.isdir(chr_dir):
                os.mkdir(chr_dir)

        # Creating a data structure to gather run information
        run_information = {}

        # Getting shapeit version
        run_information["shapeit_version"] = get_shapeit_version(
            "shapeit" if args.shapeit_bin is None else args.shapeit_bin
        )

        # Getting impute2 version
        run_information["impute2_version"] = get_impute2_version(
            "impute2" if args.impute2_bin is None else args.impute2_bin
        )

        # Getting Plink version
        run_information["plink_version"] = get_plink_version(
            "plink" if args.plink_bin is None else args.plink_bin
        )

        # Excluding markers prior to phasing (ambiguous markers [A/T and [G/C]
        # and duplicated markers
        numbers = exclude_markers_before_phasing(args.bfile, db_name, args)
        run_information.update(numbers)

        # Computing the marker missing rate
        missing_rate = compute_marker_missing_rate(args.bfile, db_name, args)

        # Checking the strand
        numbers = check_strand(
            os.path.join(args.out_dir, "chr{chrom}", "chr{chrom}"),
            "_1",
            db_name,
            args,
        )
        run_information.update(numbers)

        # Flipping the markers
        flip_markers(
            os.path.join(args.out_dir, "chr{chrom}", "chr{chrom}"),
            os.path.join(args.out_dir, "chr{chrom}", "chr{chrom}.to_flip"),
            db_name,
            args,
        )

        # Checking the strand
        numbers = check_strand(
            os.path.join(args.out_dir, "chr{chrom}", "chr{chrom}.flipped"),
            "_2",
            db_name,
            args,
            exclude=True,
        )
        run_information.update(numbers)

        # The final marker exclusion
        numbers = final_exclusion(
            os.path.join(args.out_dir, "chr{chrom}", "chr{chrom}.flipped"),
            os.path.join(args.out_dir, "chr{chrom}", "chr{chrom}.to_exclude"),
            db_name,
            args,
        )
        run_information.update(numbers)

        # Phasing the data
        samples = phase_markers(
            os.path.join(args.out_dir, "chr{chrom}", "chr{chrom}.final"),
            os.path.join(args.out_dir, "chr{chrom}",
                         "chr{chrom}.final.phased"),
            db_name,
            args,
        )

        # Gathering the chromosome length from Ensembl REST API
        chromosome_length = get_chromosome_length(args.out_dir)

        # Performs the imputation
        impute_markers(os.path.join(args.out_dir, "chr{chrom}",
                                    "chr{chrom}.final.phased.haps"),
                       os.path.join(args.out_dir, "chr{chrom}",
                                    "chr{chrom}.{start}_{end}.impute2"),
                       chromosome_length, db_name, args)

        # Getting the weighed average for cross-validation
        numbers = get_cross_validation_results(
            os.path.join(args.out_dir, "chr{chrom}",
                         "chr{chrom}.*.impute2_summary"),
        )
        run_information.update(numbers)

        # Merging the impute2 files
        merge_impute2_files(os.path.join(args.out_dir, "chr{chrom}",
                                         "chr{chrom}.*.impute2"),
                            os.path.join(args.out_dir, "chr{chrom}",
                                         "final_impute2",
                                         "chr{chrom}.imputed"),
                            args.probability, args.completion, db_name, args)

        # Gathering the imputation statistics
        numbers = gather_imputation_stats(args.probability, args.completion,
                                          len(samples), missing_rate,
                                          args.out_dir)
        run_information.update(numbers)

        # Gathering the MAF statistics
        numbers = gather_maf_stats(args.out_dir)
        run_information.update(numbers)

        # Checking that the number of good sites equals the number of sites
        # with MAF (for internal validation)
        nb_good_sites = run_information["nb_good_sites"]
        nb_sites_with_maf = run_information["nb_marker_with_maf"]
        if nb_good_sites != nb_sites_with_maf:
            logging.warning("All {} (good) imputed sites should have a MAF, "
                            "but only {} with MAF".format(nb_good_sites,
                                                          nb_sites_with_maf))

        # Gathering the execution time
        exec_time = gather_execution_time(db_name)
        run_information.update(exec_time)

        # Creating the output directory (if it doesn't exits)
        report_dir = os.path.join(args.out_dir, "report")
        if not os.path.isdir(report_dir):
            os.mkdir(report_dir)

        # Generating the report
        generate_report(report_dir, args, run_information)

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


def phase_markers(prefix, o_prefix, db_name, options):
    """Phase markers using shapeit."""
    commands_info = []
    base_command = [
        "shapeit" if options.shapeit_bin is None else options.shapeit_bin,
        "-phase",
        "--thread", str(options.shapeit_thread),
    ]

    for chrom in chromosomes:
        # The current output prefix
        c_prefix = o_prefix.format(chrom=chrom)

        remaining_command = [
            "-B", prefix.format(chrom=chrom),
            "-M", options.map_template.format(chrom=chrom),
            "-O", c_prefix,
            "-L", c_prefix + ".log",
        ]
        commands_info.append({
            "task_id": "shapeit_phase_chr{}".format(chrom),
            "name": "SHAPEIT phase chr{}".format(chrom),
            "command": base_command + remaining_command,
            "task_db": db_name,
            "o_files": [c_prefix + ext for ext in (".haps", ".sample")],
        })

    # Executing command
    logging.info("Phasing markers")
    launcher.launch_tasks(commands_info, options.thread, hpc=options.use_drmaa,
                          hpc_options=options.task_options,
                          out_dir=options.out_dir, preamble=options.preamble)
    logging.info("Done phasing markers")

    # Checking that all the sample files are the same
    compare_with = None
    for chrom in chromosomes:
        filename = o_prefix.format(chrom=chrom) + ".sample"
        compare_to = None
        with open(filename, "r") as i_file:
            compare_to = i_file.read()
            if chrom == 1:
                compare_with = compare_to

        if compare_with != compare_to:
            raise ProgramError("phased sample files are different...")

    # Returning the samples
    return [i.split(" ")[0] for i in compare_with.splitlines()[2:]]


def impute_markers(phased_haplotypes, out_prefix, chrom_length, db_name,
                   options):
    """Imputes the markers using IMPUTE2."""
    commands_info = []
    base_command = [
        "impute2" if options.impute2_bin is None else options.impute2_bin,
        "-use_prephased_g",
        "-Ne", "20000",
    ]

    # Are there any filtering rules?
    if options.filtering_rules is not None:
        base_command.append("-filt_rules_l")
        for rule in options.filtering_rules:
            base_command.append(rule)

    # Each chromosome have multiple segments
    for chrom in chromosomes:
        assert str(chrom) in chrom_length

        length = chrom_length[str(chrom)]
        start = 1
        while start < length:
            end = start + floor(options.segment_length) - 1

            # The current output prefix
            c_prefix = out_prefix.format(chrom=chrom, start=start, end=end)

            # The task ID
            task_id = "impute2_chr{}_{}_{}".format(chrom, start, end)

            # The command for this segment
            remaining_command = [
                "-known_haps_g", phased_haplotypes.format(chrom=chrom),
                "-h", options.hap_template.format(chrom=chrom),
                "-l", options.legend_template.format(chrom=chrom),
                "-m", options.map_template.format(chrom=chrom),
                "-int", str(start), str(end),
                "-o", c_prefix,
            ]
            commands_info.append({
                "task_id": task_id,
                "name": "IMPUTE2 chr{} from {} to {}".format(chrom, start,
                                                             end),
                "command": base_command + remaining_command,
                "task_db": db_name,
                "o_files": [c_prefix + "_summary", ],
            })

            # The new starting position
            start = end + 1

            # Adding the walltime for this particular task_id
            if options.use_drmaa:
                if task_id not in options.task_options:
                    # Sending the chromosome specific instead
                    value = options.task_options["impute2_chr{}".format(chrom)]
                    options.task_options[task_id] = value

    # Executing the commands
    logging.info("Imputing markers")
    launcher.launch_tasks(commands_info, options.thread, hpc=options.use_drmaa,
                          hpc_options=options.task_options,
                          out_dir=options.out_dir, preamble=options.preamble)
    logging.info("Done imputing markers")


def merge_impute2_files(in_glob, o_prefix, probability_t, completion_t,
                        db_name, options):
    """Merges impute2 files."""
    commands_info = []
    base_command = [
        "impute2-merger",
        "--probability", str(probability_t),
        "--completion", str(completion_t),
    ]

    for chrom in chromosomes:
        # The current output prefix
        c_prefix = o_prefix.format(chrom=chrom)

        # Checking that the output directory exists
        if not os.path.isdir(os.path.dirname(c_prefix)):
            os.mkdir(os.path.dirname(c_prefix))

        remaining_command = [
            "--prefix", c_prefix,
            "--chr", str(chrom),
            "-i",
        ]
        filenames = sorted(glob(in_glob.format(chrom=chrom)), key=file_sorter)
        remaining_command.extend(filenames)
        commands_info.append({
            "task_id": "merge_impute2_chr{}".format(chrom),
            "name": "Merge imputed chr{}".format(chrom),
            "command": base_command + remaining_command,
            "task_db": db_name,
            "o_files": [c_prefix + ext for ext in (".alleles",
                                                   ".completion_rates",
                                                   ".good_sites",
                                                   ".impute2",
                                                   ".imputed_sites",
                                                   ".map",
                                                   ".maf")],
        })

    # Executing command
    logging.info("Merging impute2 files")
    launcher.launch_tasks(commands_info, options.thread, hpc=options.use_drmaa,
                          hpc_options=options.task_options,
                          out_dir=options.out_dir, preamble=options.preamble)
    logging.info("Done merging reports")


def file_sorter(filename):
    """Helps in filename sorting."""
    r = re.search(r"chr(\d+)\.(\d+)_(\d+)\.impute2", filename)
    return (int(r.group(1)), int(r.group(2)), int(r.group(3)))


def get_chromosome_length(out_dir):
    """Gets the chromosome length from Ensembl REST API."""
    chrom_length = None
    filename = os.path.join(out_dir, "chromosome_lengths.txt")
    if not os.path.isfile(filename):
        logging.info("Gathering chromosome length (Ensembl, GRCh37)")

        # The URL
        url = ("http://grch37.rest.ensembl.org/info/assembly/homo_sapiens"
               "?content-type=application/json")
        result = json.loads(urlopen(url).read().decode())

        # Checking the build
        if not result["assembly_name"].startswith("GRCh37") or \
           result["default_coord_system_version"] != "GRCh37":
            raise ProgramError("{}: wrong "
                               "build".format(result["assembly_name"]))

        # Gathering the chromosome length
        chrom_length = {}
        req_chrom = {str(i) for i in range(23)} | {"X"}
        for region in result["top_level_region"]:
            if region["name"] in req_chrom:
                chrom_length[region["name"]] = region["length"]

        # Saving to file
        with open(filename, "w") as o_file:
            for chrom in sorted(chrom_length.keys()):
                print(chrom, chrom_length[chrom], sep="\t", file=o_file)

    else:
        # Gathering from file
        logging.info("Gathering chromosome length ({})".format(filename))
        with open(filename, "r") as i_file:
            chrom_length = {}
            for line in i_file:
                row = line.rstrip("\n").split("\t")
                chrom_length[row[0]] = int(row[1])

    # Checking we have all the required data
    required_chrom = {str(i) for i in chromosomes}
    if (set(chrom_length) & required_chrom) != required_chrom:
        missing = ", ".join(sorted(required_chrom - set(chrom_length)))
        raise ProgramError("missing chromosomes: {}".format(missing))

    return chrom_length


def check_strand(prefix, id_suffix, db_name, options, exclude=False):
    """Checks the strand using SHAPEIT2."""
    # Creating the command to launch
    base_command = [
        "shapeit" if options.shapeit_bin is None else options.shapeit_bin,
        "-check"
    ]

    # The output suffix
    suffix = "alignments"
    commands_info = []
    if exclude:
        # This is for exclusion
        suffix = "to_exclude.alignments"

    # The output file prefix
    o_prefix = os.path.join(options.out_dir, "chr{chrom}",
                            "chr{chrom}." + suffix)

    for chrom in chromosomes:
        # The current output prefix
        c_prefix = o_prefix.format(chrom=chrom)

        remaining_command = [
            "-B", prefix.format(chrom=chrom),
            "-M", options.map_template.format(chrom=chrom),
            "--input-ref",
            options.hap_template.format(chrom=chrom),
            options.legend_template.format(chrom=chrom),
            options.sample_file,
            "--output-log", c_prefix,
        ]
        commands_info.append({
            "task_id": "shapeit_check_chr{}{}".format(chrom, id_suffix),
            "name": "SHAPEIT check strand chr{}".format(chrom),
            "command": base_command + remaining_command,
            "task_db": db_name,
            "o_files": [c_prefix + ".snp.strand", ],
        })

    # Executing command
    logging.info("Checking strand of markers")
    launcher.launch_tasks(commands_info, options.thread, check_rc=False,
                          hpc=options.use_drmaa,
                          hpc_options=options.task_options,
                          out_dir=options.out_dir, preamble=options.preamble)
    logging.info("Done checking strand of markers")

    # The output suffix
    o_suffix = "to_flip"
    what = "flip"
    if exclude:
        o_suffix = "to_exclude"
        what = "exclude"

    # The name of the files
    filename = o_prefix + ".snp.strand"
    o_filename = os.path.join(options.out_dir, "chr{chrom}",
                              "chr{chrom}.{o_suffix}")

    # For each chromosome, we find markers to change strand
    nb_total = 0
    for chrom in chromosomes:
        # The SNP to print in the output file
        to_write = set()

        chrom_filename = filename.format(chrom=chrom)
        chrom_o_filename = o_filename.format(chrom=chrom, o_suffix=o_suffix)

        # Checking the input file exists
        if not os.path.isfile(chrom_filename):
            raise ProgramError("{}: no such file".format(chrom_filename))

        # Markers to flip
        with open(chrom_filename, "r") as i_file:
            # Reading the header
            header = i_file.readline().rstrip("\r\n").split("\t")
            header = {name: i + 1 for i, name in enumerate(header)}

            # Checking header
            for name in ("type", "main_id"):
                if name not in header:
                    raise ProgramError("{}: no column named "
                                       "{}".format(chrom_filename, name))

            # Reading the file
            for line in i_file:
                row = line.rstrip("\r\n").split("\t")
                if row[header["type"]] == "Strand":
                    to_write.add(row[header["main_id"]])

        # The number of markers to flip
        to_flip = len(to_write)

        with open(chrom_o_filename, "w") as o_file:
            print(*to_write, sep="\n", file=o_file)

        nb_total += to_flip
        logging.info("chr{}: {:,d} markers to {}".format(chrom, to_flip, what))

    # Logging the last one
    logging.info("After strand check: {:,d} markers "
                 "to {}".format(nb_total, what))

    return {"nb_{}".format(what): nb_total}


def flip_markers(prefix, to_flip, db_name, options):
    """Flip markers."""
    # The commands to run
    commands_info = []
    base_command = [
        "plink" if options.plink_bin is None else options.plink_bin,
        "--noweb",
        "--make-bed",
    ]

    # The output prefix
    o_prefix = os.path.join(options.out_dir, "chr{chrom}",
                            "chr{chrom}.flipped")

    for chrom in chromosomes:
        # The current output prefix
        c_prefix = o_prefix.format(chrom=chrom)

        remaining_command = [
            "--bfile", prefix.format(chrom=chrom),
            "--flip", to_flip.format(chrom=chrom),
            "--out", c_prefix,
        ]
        commands_info.append({
            "task_id": "plink_flip_chr{}".format(chrom),
            "name": "plink flip chr{}".format(chrom),
            "command": base_command + remaining_command,
            "task_db": db_name,
            "o_files": [c_prefix + ext for ext in (".bed", ".bim", ".fam")],
        })

    # Executing command
    logging.info("Flipping markers")
    launcher.launch_tasks(commands_info, options.thread, hpc=options.use_drmaa,
                          hpc_options=options.task_options,
                          out_dir=options.out_dir, preamble=options.preamble)
    logging.info("Done flipping markers")


def final_exclusion(prefix, to_exclude, db_name, options):
    """Flip markers."""
    # The commands to run
    commands_info = []
    base_command = [
        "plink" if options.plink_bin is None else options.plink_bin,
        "--noweb",
        "--make-bed",
    ]

    # The output prefix
    o_prefix = os.path.join(options.out_dir, "chr{chrom}", "chr{chrom}.final")

    # The output files (for statistics)
    bims = []

    for chrom in chromosomes:
        # The current output prefix
        c_prefix = o_prefix.format(chrom=chrom)

        remaining_command = [
            "--bfile", prefix.format(chrom=chrom),
            "--exclude", to_exclude.format(chrom=chrom),
            "--out", c_prefix,
        ]
        commands_info.append({
            "task_id": "plink_final_exclude_chr{}".format(chrom),
            "name": "plink final exclude chr{}".format(chrom),
            "command": base_command + remaining_command,
            "task_db": db_name,
            "o_files": [c_prefix + ext for ext in (".bed", ".bim", ".fam")],
        })

        bims.append(c_prefix + ".bim")

    # Executing command
    logging.info("Final marker exclusion")
    launcher.launch_tasks(commands_info, options.thread, hpc=options.use_drmaa,
                          hpc_options=options.task_options,
                          out_dir=options.out_dir, preamble=options.preamble)
    logging.info("Done final marker exclusion")

    nb_markers = 0
    for bim in bims:
        with open(bim, "r") as i_file:
            for line in i_file:
                nb_markers += 1

    return {"nb_phasing_markers": "{:,d}".format(nb_markers)}


def compute_marker_missing_rate(prefix, db_name, options):
    """Compute (using Plink) marker missing rate."""
    # The output prefix
    o_prefix = os.path.join(options.out_dir, "missing")
    if not os.path.isdir(o_prefix):
        os.mkdir(o_prefix)
    o_prefix = os.path.join(o_prefix, "missing")

    # The command to run
    commands_info = []
    command = [
        "plink" if options.plink_bin is None else options.plink_bin,
        "--noweb",
        "--bfile", prefix,
        "--missing",
        "--out", o_prefix,
    ]

    commands_info.append({
        "task_id": "plink_missing_rate",
        "name": "plink missing rate",
        "command": command,
        "task_db": db_name,
        "o_files": [o_prefix + ext for ext in (".lmiss", ".imiss")],
    })

    # Executing command
    logging.info("Computing missing rate")
    launcher.launch_tasks(commands_info, options.thread, hpc=options.use_drmaa,
                          hpc_options=options.task_options,
                          out_dir=options.out_dir, preamble=options.preamble)
    logging.info("Done computing missing rate")

    # Getting the missing rate
    logging.info("Reading the missing rate")
    return pd.read_csv(o_prefix + ".lmiss", delim_whitespace=True)


def exclude_markers_before_phasing(prefix, db_name, options):
    """Finds and excludes ambiguous markers (A/T and G/C) or duplicated."""
    # The ambiguous genotypes
    ambiguous_genotypes = {"AT", "TA", "GC", "CG"}

    # The positions that have been written to file
    kept_positions = set()
    nb_ambiguous = 0
    nb_kept = 0
    nb_dup = 0

    # Logging
    logging.info("Finding markers to exclude")

    # Counting the total number of markers and samples
    nb_markers = 0
    nb_samples = 0
    nb_special_markers = 0

    with open(prefix + ".fam", "r") as i_file:
        for line in i_file:
            nb_samples += 1

    o_filename = os.path.join(options.out_dir, "markers_to_exclude.txt")
    with open(o_filename, "w") as o_file, \
            open(prefix + ".bim", "r") as i_file:
        for line in i_file:
            nb_markers += 1
            row = line.rstrip("\r\n").split("\t")

            is_special = False
            if row[0] in {"23", "24", "25", "26"}:
                nb_special_markers += 1
                is_special = True

            # Checking the alleles
            if row[4] + row[5] in ambiguous_genotypes:
                if not is_special:
                    nb_ambiguous += 1
                print(row[1], file=o_file)
                logging.debug("  - {}: {}: "
                              "ambiguous".format(row[1], row[4] + row[5]))
                continue

            # Checking if we already have this marker
            if (row[0], row[3]) in kept_positions:
                if not is_special:
                    nb_dup += 1
                print(row[1], file=o_file)
                logging.debug("  - {}: duplicated".format(row[1]))
                continue

            # We keep this marker
            kept_positions.add((row[0], row[3]))
            if not is_special:
                nb_kept += 1

    # Logging
    logging.info("  - {:,d} special markers".format(nb_special_markers))
    logging.info("  - {:,d} ambiguous markers removed".format(nb_ambiguous))
    logging.info("  - {:,d} duplicated markers removed".format(nb_dup))
    logging.info("  - {:,d} markers kept".format(nb_kept))

    # The commands to run
    commands_info = []
    base_command = [
        "plink" if options.plink_bin is None else options.plink_bin,
        "--noweb",
        "--exclude", os.path.join(options.out_dir, "markers_to_exclude.txt"),
        "--make-bed",
    ]

    # The output prefix
    o_prefix = os.path.join(options.out_dir, "chr{chrom}", "chr{chrom}")

    for chrom in chromosomes:
        # The current output prefix
        c_prefix = o_prefix.format(chrom=chrom)

        remaining_command = [
            "--bfile", prefix,
            "--chr", str(chrom),
            "--out", c_prefix,
        ]
        commands_info.append({
            "task_id": "plink_exclude_chr{}".format(chrom),
            "name": "plink exclude chr{}".format(chrom),
            "command": base_command + remaining_command,
            "task_db": db_name,
            "o_files": [c_prefix + ext for ext in (".bed", ".bim", ".fam")],
        })

    # Executing command
    logging.info("Excluding and splitting markers")
    launcher.launch_tasks(commands_info, options.thread, hpc=options.use_drmaa,
                          hpc_options=options.task_options,
                          out_dir=options.out_dir, preamble=options.preamble)
    logging.info("Done excluding and splitting markers")

    return {
        "initial_nb_markers": "{:,d}".format(nb_markers),
        "initial_nb_samples": "{:,d}".format(nb_samples),
        "nb_ambiguous":       "{:,d}".format(nb_ambiguous),
        "nb_duplicates":      "{:,d}".format(nb_dup),
        "nb_special_markers": "{:,d}".format(nb_special_markers),
    }


def get_cross_validation_results(glob_pattern):
    """Creates a weighed mean for each chromosome for cross-validation."""
    logging.info("Gathering cross-validation statistics")

    # The regular expressions used here
    nb_genotypes_re = re.compile(r"^In the current analysis, IMPUTE2 masked, "
                                 r"imputed, and evaluated (\d+) genotypes")
    table_header_re = re.compile(r"Interval\s+#Genotypes\s+%Concordance\s+"
                                 r"Interval\s+%Called\s+%Concordance")
    split_re = re.compile(r"\s+")

    # The intervals
    table_1_intervals = ["[0.0-0.1]", "[0.1-0.2]", "[0.2-0.3]", "[0.3-0.4]",
                         "[0.4-0.5]", "[0.5-0.6]", "[0.6-0.7]", "[0.7-0.8]",
                         "[0.8-0.9]", "[0.9-1.0]"]
    table_2_intervals = ["[>=0.0]", "[>=0.1]", "[>=0.2]", "[>=0.3]", "[>=0.4]",
                         "[>=0.5]", "[>=0.6]", "[>=0.7]", "[>=0.8]", "[>=0.9]"]

    # The final data per chromosome
    per_chrom_table_1 = {}
    per_chrom_table_2 = {}
    chrom_nb_geno = {}

    # The final data for all the dataset
    tot_concordance = defaultdict(float)
    tot_nb_geno = defaultdict(int)
    tot_weight = defaultdict(int)
    tot_cumm_called = defaultdict(int)
    tot_cumm_concordance = defaultdict(float)
    tot_cumm_weight = defaultdict(int)
    final_nb_genotypes = 0

    # For each chromosome
    for chrom in chromosomes:
        filenames = glob(glob_pattern.format(chrom=chrom))

        # The total number of genotypes for this chromosome
        tot_chrom_nb_genotypes = 0

        # The numbers for table 1
        tot_chrom_concordance = defaultdict(float)
        tot_chrom_nb_geno = defaultdict(int)
        tot_chrom_weight = defaultdict(int)

        # The numbers for table 2
        tot_chrom_cumm_called = defaultdict(int)
        tot_chrom_cumm_concordance = defaultdict(float)
        tot_chrom_cumm_weight = defaultdict(int)

        for filename in filenames:
            with open(filename, "r") as i_file:
                line = i_file.readline()
                nb_genotypes = nb_genotypes_re.search(line)
                while (line != "") and (nb_genotypes is None):
                    line = i_file.readline()
                    nb_genotypes = nb_genotypes_re.search(line)

                # Check if we have data
                if nb_genotypes is None:
                    continue

                # We have the number of genotypes
                nb_genotypes = int(nb_genotypes.group(1))
                if nb_genotypes == 0:
                    continue

                # Increasing the total number of genotypes for this chromosome
                tot_chrom_nb_genotypes += nb_genotypes
                final_nb_genotypes += nb_genotypes

                # Looking for the headers of the two tables
                line = i_file.readline()
                table_header = table_header_re.search(line)
                while (line != "") and (table_header is None):
                    line = i_file.readline()
                    table_header = table_header_re.search(line)

                # We now should have data
                if table_header is None:
                    raise ProgramError("Problem with {}".format(filename))

                # Now reading the data inside the tables
                for line in i_file:
                    # Getting the value for the tables
                    table = split_re.split(line.strip())

                    # The first table
                    interval, nb_geno, concordance = table[:3]

                    # Casting to values
                    nb_geno = int(nb_geno)
                    concordance = float(concordance)

                    # The total number of genotypes for this interval
                    tot_nb_geno[interval] += nb_geno
                    tot_chrom_nb_geno[interval] += nb_geno

                    # The concordance for this interval
                    tot_concordance[interval] += (concordance * nb_geno)
                    tot_chrom_concordance[interval] += (concordance * nb_geno)

                    # The weight for this interval
                    tot_weight[interval] += nb_geno
                    tot_chrom_weight[interval] += nb_geno

                    # The second table
                    table = table[3:]
                    interval = "".join(table[:3])
                    pc_called = float(table[-2])
                    pc_concordance = float(table[-1])

                    # The interval number of genotypes
                    nb_called = pc_called * nb_genotypes / 100
                    tot_chrom_cumm_called[interval] += nb_called
                    tot_cumm_called[interval] += nb_called

                    # The concordance for this interval
                    weighted_value = nb_called * pc_concordance
                    tot_chrom_cumm_concordance[interval] += weighted_value
                    tot_cumm_concordance[interval] += weighted_value

                    # The weight for this interval
                    tot_chrom_cumm_weight[interval] += nb_called
                    tot_cumm_weight[interval] += nb_called

        # Computing the weighted average for the first table
        table_1_data = []
        for interval in table_1_intervals:
            nb_geno = tot_chrom_nb_geno[interval]
            concordance = tot_chrom_concordance[interval]
            weight = tot_chrom_weight[interval]

            # Computing the weighted concordance
            weighted_concordance = 0
            if weight != 0:
                weighted_concordance = concordance / weight

            # Saving the data for the current chromosome
            table_1_data.append([
                interval,
                "{:,d}".format(nb_geno),
                "{:.1f}".format(weighted_concordance),
            ])

        # Computing the weighed average for the second table
        table_2_data = []
        for interval in table_2_intervals:
            nb_called = tot_chrom_cumm_called[interval]
            concordance = tot_chrom_cumm_concordance[interval]
            weight = tot_chrom_cumm_weight[interval]

            # Computing the weighted concordance
            weighted_concordance = 0
            if weight != 0:
                weighted_concordance = concordance / weight

            # Saving the data for the current chromosome
            table_2_data.append([
                interval,
                "{:.1f}".format(nb_called / tot_chrom_nb_genotypes * 100),
                "{:.1f}".format(weighted_concordance),
            ])

        # Saving the tables
        per_chrom_table_1[chrom] = table_1_data
        per_chrom_table_2[chrom] = table_2_data
        chrom_nb_geno[chrom] = tot_chrom_nb_genotypes

    # The final table 1
    table_1_data = []
    for interval in table_1_intervals:
        nb_geno = tot_nb_geno[interval]
        concordance = tot_concordance[interval]
        weight = tot_weight[interval]

        # Computing the weighted concordance
        weighted_concordance = 0
        if weight != 0:
            weighted_concordance = concordance / weight

        table_1_data.append([
            interval,
            "{:,d}".format(nb_geno),
            "{:.1f}".format(weighted_concordance),
        ])

    # The final table 2
    table_2_data = []
    for interval in table_2_intervals:
        nb_called = tot_cumm_called[interval]
        concordance = tot_cumm_concordance[interval]
        weight = tot_cumm_weight[interval]

        # Computing the weighted concordance
        weighted_concordance = 0
        if weight != 0:
            weighted_concordance = concordance / weight

        # Saving the data for all the chromosome
        table_2_data.append([
            interval,
            "{:.1f}".format(nb_called / final_nb_genotypes * 100),
            "{:.1f}".format(weighted_concordance),
        ])

    # Returning the data
    return {
        "cross_validation_final_nb_genotypes": final_nb_genotypes,
        "cross_validation_nb_genotypes_chrom": chrom_nb_geno,
        "cross_validation_table_1":            table_1_data,
        "cross_validation_table_2":            table_2_data,
        "cross_validation_table_1_chrom":      per_chrom_table_1,
        "cross_validation_table_2_chrom":      per_chrom_table_2,
    }


def gather_imputation_stats(prob_t, completion_t, nb_samples, missing, o_dir):
    """Gathers imputation statistics from the merged dataset."""
    logging.info("Gathering imputation statistics")

    # The stats we want to gather
    # For all the sites
    tot_nb_sites = 0
    sum_rates = 0

    # For the best sites
    tot_good_sites = 0
    sum_good_rates = 0

    # All imputed sites
    genotyped_sites = []

    # For each chromosome, get the statistics
    filename_template = os.path.join(o_dir, "chr{chrom}", "final_impute2",
                                     "chr{chrom}.imputed.{suffix}")
    for chrom in chromosomes:
        logging.info("  - chromosome {}".format(chrom))

        # First, we read the imputed sites
        filename = filename_template.format(chrom=chrom,
                                            suffix="imputed_sites")
        imputed_sites = None
        with open(filename, "r") as i_file:
            imputed_sites = {i for i in i_file.read().splitlines()}

        # Then, we read the completion using a DataFrame
        filename = filename_template.format(chrom=chrom,
                                            suffix="completion_rates")
        completion_data = pd.read_csv(filename, sep="\t")

        # Saving the stats for all sites
        tot_nb_sites += completion_data.shape[0]
        sum_rates += completion_data.completion_rate.sum()

        # Saving the stats for good sites
        good_sites = completion_data.completion_rate >= completion_t
        tot_good_sites += completion_data[good_sites].shape[0]
        sum_good_rates += completion_data[good_sites].completion_rate.sum()

        # Saving the imputed sites call rates
        imputed_sites = completion_data.name.isin(imputed_sites)
        genotyped_sites.append(completion_data[~imputed_sites])

    # Concatenating the completion rates for the genotyped sites
    genotyped_sites = pd.concat(genotyped_sites)
    now_incomplete_sites = genotyped_sites.nb_missing > 0

    # Getting the original missing rates for the genotyped sites
    sites = missing.SNP.isin(genotyped_sites.name)
    incomplete_sites = missing.N_MISS > 0

    # Gathering genotyped sites statistics
    nb_genotyped_sites = genotyped_sites.shape[0]
    nb_incomplete_sites = missing[sites & incomplete_sites].shape[0]
    pct_incomplete_sites = nb_incomplete_sites / nb_genotyped_sites * 100
    nb_sites_now_complete = genotyped_sites[~now_incomplete_sites].shape[0]
    nb_missing_geno = missing[sites & incomplete_sites].N_MISS.sum()
    now_nb_missing_geno = genotyped_sites[now_incomplete_sites]
    now_nb_missing_geno = now_nb_missing_geno.nb_missing.sum()
    nb_geno_now_complete = nb_missing_geno - now_nb_missing_geno
    pct_geno_now_complete = nb_geno_now_complete / nb_missing_geno * 100

    # Computing the rates
    comp_rate = sum_rates / tot_nb_sites
    good_comp_rate = sum_good_rates / tot_good_sites

    # Other statistics
    pct_good_sites = tot_good_sites / tot_nb_sites * 100
    mean_missing = nb_samples * (1 - good_comp_rate)

    return {
        "prob_threshold":             "{:.1f}".format(prob_t * 100),
        "nb_imputed":                 "{:,d}".format(tot_nb_sites),
        "average_comp_rate":          "{:.1f}".format(comp_rate * 100),
        "rate_threshold":             "{:.1f}".format(completion_t * 100),
        "nb_good_sites":              "{:,d}".format(tot_good_sites),
        "pct_good_sites":             "{:.1f}".format(pct_good_sites),
        "average_comp_rate_cleaned":  "{:.1f}".format(good_comp_rate * 100),
        "mean_missing":               "{:.1f}".format(mean_missing),
        "nb_samples":                 "{:,d}".format(nb_samples),
        "nb_genotyped":               "{:,d}".format(nb_genotyped_sites),
        "nb_genotyped_not_complete":  "{:,d}".format(nb_incomplete_sites),
        "pct_genotyped_not_complete": "{:.1f}".format(pct_incomplete_sites),
        "nb_geno_now_complete":       "{:,d}".format(nb_geno_now_complete),
        "nb_missing_geno":            "{:,d}".format(nb_missing_geno),
        "pct_geno_now_complete":      "{:.1f}".format(pct_geno_now_complete),
        "nb_site_now_complete":       "{:,d}".format(nb_sites_now_complete)
    }


def gather_maf_stats(o_dir):
    """Gather minor allele frequencies from imputation."""
    logging.info("Gathering imputation statistics")

    # The statistics we want to gather
    nb_marker_with_maf = 0   # The number of markers for which we have a MAF
    nb_maf_geq_01 = 0        # The number of markers with MAF >= 0.01
    nb_maf_geq_05 = 0        # The number of markers with MAF >= 0.05
    nb_maf_lt_05 = 0         # The number of markers with MAF < 0.05
    nb_maf_lt_01 = 0         # The number of markers with MAF < 0.01
    nb_maf_geq_01_lt_05 = 0  # The number of markers with 0.01 <= MAF < 0.05
    nb_maf_nan = 0           # The number of markers without MAF

    # For each chromosome, get the MAF statistics
    filename_template = os.path.join(o_dir, "chr{chrom}", "final_impute2",
                                     "chr{chrom}.imputed.{suffix}")
    for chrom in chromosomes:
        logging.info("  - chromosome {}".format(chrom))

        # The name of the file
        maf_filename = filename_template.format(chrom=chrom, suffix="maf")
        good_sites_filename = filename_template.format(chrom=chrom,
                                                       suffix="good_sites")

        # Checking the file exists
        for filename in [maf_filename, good_sites_filename]:
            if not os.path.isfile(filename):
                raise ProgramError("{}: no such file".format(filename))

        # Reading the list of good sites
        good_sites = None
        with open(good_sites_filename, "r") as i_file:
            good_sites = set(i_file.read().splitlines())

        # Reading the file using pandas
        maf = pd.read_csv(maf_filename, sep="\t")

        # Keeping only the good sites
        maf = maf[maf.name.isin(good_sites)]

        # Excluding sites with no MAF (NaN values)
        null_maf = maf.maf.isnull()
        nb_nan = null_maf.sum()
        maf = maf[~null_maf]

        # There should not be any NaN sites...
        if nb_nan > 0:
            logging.warning("chr{}: good sites with invalid MAF "
                            "(NaN)".format(chrom))

        # Checking we have MAF (and not just frequencies)
        maf_description = maf.maf.describe()
        if maf_description["max"] > 0.5:
            bad = maf.loc[maf.maf.idxmax(), ["name", "maf"]]
            raise ProgramError("{}: {}: invalid MAF".format(str(bad["name"]),
                                                            round(bad.maf, 3)))
        if maf_description["max"] < 0:
            bad = maf.loc[maf.maf.idxmin(), ["name", "maf"]]
            raise ProgramError("{}: {}: invalid MAF".format(str(bad["name"]),
                                                            round(bad.maf, 3)))

        # Some of the true/false we need to keep (to not compute multiple time)
        maf_geq_01 = maf.maf >= 0.01
        maf_lt_05 = maf.maf < 0.05

        # Updating the statistics
        nb_marker_with_maf += maf.shape[0]
        nb_maf_nan += nb_nan
        nb_maf_geq_01 += maf_geq_01.sum()
        nb_maf_geq_05 += (maf.maf >= 0.05).sum()
        nb_maf_lt_05 += maf_lt_05.sum()
        nb_maf_lt_01 += (maf.maf < 0.01).sum()
        nb_maf_geq_01_lt_05 += (maf_geq_01 & maf_lt_05).sum()

    # Checking
    nb_total = nb_maf_lt_01 + nb_maf_geq_01_lt_05 + nb_maf_geq_05
    if nb_total != nb_marker_with_maf:
        raise ProgramError("something went wrong")

    # Computing the percentages
    pct_maf_geq_01 = 0
    pct_maf_geq_05 = 0
    pct_maf_lt_05 = 0
    pct_maf_lt_01 = 0
    pct_maf_geq_01_lt_05 = 0

    if nb_marker_with_maf > 0:
        pct_maf_geq_01 = nb_maf_geq_01 / nb_marker_with_maf * 100
        pct_maf_geq_05 = nb_maf_geq_05 / nb_marker_with_maf * 100
        pct_maf_lt_05 = nb_maf_lt_05 / nb_marker_with_maf * 100
        pct_maf_lt_01 = nb_maf_lt_01 / nb_marker_with_maf * 100
        pct_maf_geq_01_lt_05 = nb_maf_geq_01_lt_05 / nb_marker_with_maf * 100

    else:
        logging.warning("There were no marker with MAF (something went wrong)")

    # Generating a pie chart if matplotlib is installed
    frequency_pie = ""
    if (nb_marker_with_maf > 0) and HAS_MATPLOTLIB:
        # Creating a figure and axe
        figure, axe = plt.subplots(1, 1, figsize=(6, 9))

        # The colors
        colors = ["#0099CC", "#669900", "#FF8800"]

        # The data for the pie chart
        labels = [
            "{:.1f}%".format(nb_maf_lt_01 / nb_marker_with_maf * 100),
            "{:.1f}%".format(nb_maf_geq_01_lt_05 / nb_marker_with_maf * 100),
            "{:.1f}%".format(nb_maf_geq_05 / nb_marker_with_maf * 100),
        ]
        sizes = [nb_maf_lt_01, nb_maf_geq_01_lt_05, nb_maf_geq_05]
        explode = (0.05, 0.05, 0.05)
        wedges, texts = axe.pie(sizes, explode=explode, labels=labels,
                                colors=colors, startangle=90,
                                wedgeprops={"linewidth": 0},
                                textprops={"fontsize": 12, "weight": "bold"})

        # Changing the label parameters
        for text in texts:
            text.set_bbox({"boxstyle": "round", "fc": "#C0C0C0",
                           "ec": "#C0C0C0"})

        # Shrink current axis by 50%
        bbox = axe.get_position()
        axe.set_position([bbox.x0, bbox.y0, bbox.width * 0.5,
                         bbox.height * 0.5])

        # If no restriction was performed while imputing, there will be a high
        # majority of ultra rare variants... We need to move the text box if
        # they touch
        # Getting the renderer and the bboxes
        renderer = figure.canvas.get_renderer()
        bbox_1 = texts[1].get_window_extent(renderer=renderer)
        bbox_2 = texts[2].get_window_extent(renderer=renderer)
        while bbox_1.overlaps(bbox_2):
            # Moving the first label a bit
            x, y = texts[1].get_position()
            texts[1].set_x(x + 0.006)
            bbox_1 = texts[1].get_window_extent(renderer=renderer)

            # Moving the second label (MAF >= 0.05)
            x, y = texts[2].get_position()
            texts[2].set_x(x - 0.036)
            bbox_2 = texts[2].get_window_extent(renderer=renderer)

        # Adding a legend in the margin with custom patches
        ultra_rare = mpatches.Patch(color=colors[0], linewidth=3)
        rare = mpatches.Patch(color=colors[1], linewidth=3)
        common = mpatches.Patch(color=colors[2], linewidth=3)
        axe.legend(
            [ultra_rare, rare, common],
            [r"$MAF < 1\%$", r"$1\% \leq MAF < 5\%$", r"$MAF \geq 5\%$"],
            bbox_to_anchor=(1.3, 1),
            loc="upper left",
            ncol=1,
            frameon=False,
            fontsize=21,
        )

        # Setting the axis to equal size
        axe.axis("equal")

        # Saving and closing the figure
        frequency_pie = os.path.join(o_dir, "frequency_pie.pdf")
        plt.savefig(frequency_pie, bbox_inches="tight", figure=figure)
        plt.close(figure)

    return {
        "nb_maf_nan":           "{:,d}".format(nb_maf_nan),
        "nb_marker_with_maf":   "{:,d}".format(nb_marker_with_maf),
        "nb_maf_geq_01":        "{:,d}".format(nb_maf_geq_01),
        "nb_maf_geq_05":        "{:,d}".format(nb_maf_geq_05),
        "nb_maf_lt_05":         "{:,d}".format(nb_maf_lt_05),
        "nb_maf_lt_01":         "{:,d}".format(nb_maf_lt_01),
        "nb_maf_geq_01_lt_05":  "{:,d}".format(nb_maf_geq_01_lt_05),
        "pct_maf_geq_01":       "{:.1f}".format(pct_maf_geq_01),
        "pct_maf_geq_05":       "{:.1f}".format(pct_maf_geq_05),
        "pct_maf_lt_05":        "{:.1f}".format(pct_maf_lt_05),
        "pct_maf_lt_01":        "{:.1f}".format(pct_maf_lt_01),
        "pct_maf_geq_01_lt_05": "{:.1f}".format(pct_maf_geq_01_lt_05),
        "frequency_pie":        frequency_pie,
    }


def gather_execution_time(db_name):
    """Gather all the execution times."""
    # Getting all the execution time from the DB
    exec_time = get_all_runtimes(db_name)

    # Getting the execution time for the steps
    plink_exclude_exec_time = []
    shapeit_check_1_exec_time = []
    plink_flip_exec_time = []
    shapeit_check_2_exec_time = []
    plink_final_exec_time = []
    shapeit_phase_exec_time = []
    impute2_exec_time = []
    merge_impute2_exec_time = []
    for chrom in chromosomes:
        # Getting the time for 'plink_exclude'
        seconds = exec_time["plink_exclude_chr{}".format(chrom)]
        plink_exclude_exec_time.append([chrom, seconds])

        # Getting the time for 'shapeit_check_1'
        seconds = exec_time["shapeit_check_chr{}_1".format(chrom)]
        shapeit_check_1_exec_time.append([chrom, seconds])

        # Getting the time for 'plink_flip'
        seconds = exec_time["plink_flip_chr{}".format(chrom)]
        plink_flip_exec_time.append([chrom, seconds])

        # Getting the time for 'shapeit_check_2'
        seconds = exec_time["shapeit_check_chr{}_2".format(chrom)]
        shapeit_check_2_exec_time.append([chrom, seconds])

        # Getting the time for 'plink_final_exclude'
        seconds = exec_time["plink_final_exclude_chr{}".format(chrom)]
        plink_final_exec_time.append([chrom, seconds])

        # Getting the time for 'shapeit_phase'
        seconds = exec_time["shapeit_phase_chr{}".format(chrom)]
        shapeit_phase_exec_time.append([chrom, seconds])

        # Getting the execution times for the imputation step
        chr_imputation_tasks = [
            i for i in exec_time.keys()
            if i.startswith("impute2_chr{}_".format(chrom))
        ]
        seconds = [
            exec_time[task_name] for task_name in chr_imputation_tasks
        ]
        impute2_exec_time.append([
            chrom,
            len(chr_imputation_tasks),
            int(round(sum(seconds) / len(seconds), 0)),
            max(seconds),
        ])

        # Getting the time for 'merge_impute2'
        seconds = exec_time["merge_impute2_chr{}".format(chrom)]
        merge_impute2_exec_time.append([chrom, seconds])

    # Getting the execution time for the second step (plink missing)
    plink_missing_exec_time = exec_time["plink_missing_rate"]

    # Returning the time
    return {
        "plink_exclude_exec_time":   plink_exclude_exec_time,
        "shapeit_check_1_exec_time": shapeit_check_1_exec_time,
        "plink_missing_exec_time":   plink_missing_exec_time,
        "plink_flip_exec_time":      plink_flip_exec_time,
        "shapeit_check_2_exec_time": shapeit_check_2_exec_time,
        "plink_final_exec_time":     plink_final_exec_time,
        "shapeit_phase_exec_time":   shapeit_phase_exec_time,
        "merge_impute2_exec_time":   merge_impute2_exec_time,
        "impute2_exec_time":         impute2_exec_time,
    }


def read_preamble(filename):
    """Reads the preamble file."""
    if filename is None:
        return ""

    preamble = None
    with open(filename, "r") as i_file:
        preamble = i_file.read()

    while not preamble.endswith("\n\n"):
        preamble += "\n"

    if not preamble.startswith("\n"):
        preamble = "\n" + preamble

    return preamble


def get_shapeit_version(binary):
    """Gets the SHAPEIT version from the binary."""
    # Running the command
    command = [binary, "--version"]
    proc = Popen(command, stdout=PIPE)
    output = proc.communicate()[0].decode()

    # Finding the version
    version = re.search(r"Version : ([\S]+)", output)
    if version is None:
        version = "unknown"
    else:
        version = version.group(1)

    logging.info("Will be using SHAPEIT version {}".format(version))

    return version


def get_impute2_version(binary):
    """Gets the IMPUTE2 version from the binary."""
    # Running the command
    command = [binary]
    proc = Popen(command, stdout=PIPE)
    output = proc.communicate()[0].decode()

    # Deleting the output files automatically created by IMPUTE2
    for filename in ["test.impute2_summary", "test.impute2_warnings"]:
        if os.path.isfile(filename):
            os.remove(filename)

    # Finding the version
    version = re.search(r"IMPUTE version ([\S]+)", output)
    if version is None:
        version = "unknown"
    else:
        version = version.group(1)

    logging.info("Will be using IMPUTE2 version {}".format(version))

    return version


def get_plink_version(binary):
    """Gets the Plink version from the binary."""
    # Running the command
    command = [binary, "--noweb"]
    proc = Popen(command, stdout=PIPE, stderr=PIPE)
    output = proc.communicate()[0].decode()

    # Deleting the output file automatically created by Plink
    if os.path.isfile("plink.log"):
        os.remove("plink.log")

    # Finding the version
    version = re.search(r"\|\s+PLINK!\s+\|\s+(\S+)\s+\|", output)
    if version is None:
        version = "unknown"
    else:
        version = version.group(1)

    logging.info("Will be using Plink version {}".format(version))

    return version


def check_args(args):
    """Checks the arguments and options."""
    # Checking the presence of the BED, BIM and BAM files
    for suffix in (".bed", ".bim", ".fam"):
        if not os.path.isfile(args.bfile + suffix):
            raise ProgramError("{}: no such file".format(args.bfile + suffix))

    # Checking the thread
    if args.thread < 1:
        raise ProgramError("thread should be one or more")
    if args.shapeit_thread < 1:
        raise ProgramError("thread should be one or more")

    # Checking IMPUTE2's files
    for template in (args.hap_template, args.legend_template,
                     args.map_template):
        for chrom in chromosomes:
            # Checking the haplotype file
            filename = template.format(chrom=chrom)
            if not os.path.isfile(filename):
                raise ProgramError("{}: no such file".format(filename))
    if not os.path.isfile(args.sample_file):
        raise ProgramError("{}: no such file".format(args.sample_file))

    # Checking the SHAPEIT binary if required
    if args.shapeit_bin is not None:
        if not os.path.isfile(args.shapeit_bin):
            raise ProgramError("{}: no such file".format(args.shapeit_bin))
    else:
        if which("shapeit") is None:
            raise ProgramError("shapeit: not in the path (use --shapeit-bin)")

    # Checking the IMPUTE2 binary if required
    if args.impute2_bin is not None:
        if not os.path.isfile(args.impute2_bin):
            raise ProgramError("{}: no such file".format(args.impute2_bin))
    else:
        if which("impute2") is None:
            raise ProgramError("impute2: not in the path (use --impute2-bin)")

    # Checking that Plink is in the path
    if args.plink_bin is not None:
        if not os.path.isfile(args.plink_bin):
            raise ProgramError("{}: no such file".format(args.plink_bin))
    else:
        if which("plink") is None:
            raise ProgramError("plink: not in the path (use --plink-bin)")

    # Checking the segment length
    if args.segment_length <= 0:
        raise ProgramError("{}: invalid segment "
                           "length".format(args.segment_length))
    if args.segment_length < 1e3:
        # This is too small.. We continue with a warning
        logging.warning("segment length ({:g} bp) is too "
                        "small".format(args.segment_length))
    if args.segment_length > 5e6:
        # This is too big... We continue with a warning
        logging.warning("segment length ({:g} bp) is more than "
                        "5Mb".format(args.segment_length))

    # Checking the preamble file (if required)
    if args.preamble is not None:
        if not os.path.isfile(args.preamble):
            raise ProgramError("{}: no such file".format(args.preamble))

    # Checking the DRMAA configuration file
    if args.use_drmaa:
        if args.drmaa_config is not None:
            if not os.path.isfile(args.drmaa_config):
                raise ProgramError("{}: no such "
                                   "file".format(args.drmaa_config))

    return True


def parse_args(parser):
    """Parses the command line options and arguments."""
    parser.add_argument("-v", "--version", action="version",
                        version="%(prog)s {}".format(__version__))
    parser.add_argument("--debug", action="store_true",
                        help="Set the logging level to debug")
    parser.add_argument("--thread", type=int, default=1,
                        help="The number of thread [%(default)d]")

    # The input files
    group = parser.add_argument_group("Input Options")
    group.add_argument("--bfile", type=str, metavar="PREFIX", required=True,
                       help="The prefix of the binary pedfiles (input data)")

    # The output options
    group = parser.add_argument_group("Output Options")
    group.add_argument("--output-dir", type=str, metavar="DIR", default="gwip",
                       dest="out_dir", help=("The name of the output "
                                             "directory [%(default)s]"))

    # The HPC options
    group = parser.add_argument_group("HPC Options")
    group.add_argument("--use-drmaa", action="store_true",
                       help="Launch tasks using DRMAA",)
    group.add_argument("--drmaa-config", type=str, metavar="FILE",
                       help=("The configuration file for tasks (use this "
                             "option when launching tasks using DRMAA). This "
                             "file should describe the walltime and the "
                             "number of nodes/processors to use for each "
                             "task."))
    group.add_argument("--preamble", type=str, metavar="FILE",
                       help=("This option should be used when using DRMAA on "
                             "a HPC to load required module and set "
                             "environment variables. The content of the file "
                             "will be added between the 'shebang' line and "
                             "the tool command."))

    # The SHAPEIT software options
    group = parser.add_argument_group("SHAPEIT Options")
    group.add_argument("--shapeit-bin", type=str, metavar="BINARY",
                       help="The SHAPEIT binary if it's not in the path")
    group.add_argument("--shapeit-thread", type=int, default=1,
                       help="The number of thread for phasing [%(default)d]")

    # The Plink option
    group = parser.add_argument_group("Plink Options")
    group.add_argument("--plink-bin", type=str, metavar="BINARY",
                       help="The Plink binary if it's not in the path")

    # The IMPUTE2 software options
    group = parser.add_argument_group("IMPUTE2 Options")
    group.add_argument("--impute2-bin", type=str, metavar="BINARY",
                       help="The IMPUTE2 binary if it's not in the path")
    group.add_argument("--segment-length", type=float, metavar="BP",
                       default=5e6, help=("The length of a single segment for "
                                          "imputation [%(default).1g)]"))
    group.add_argument("--hap-template", type=str, metavar="TEMPLATE",
                       required=True,
                       help=("The template for IMPUTE2's haplotype files "
                             "(replace the chromosome number by '{chrom}', "
                             "e.g. '1000GP_Phase3_chr{chrom}.hap.gz')"))
    group.add_argument("--legend-template", type=str, metavar="TEMPLATE",
                       required=True,
                       help=("The template for IMPUTE2's legend files "
                             "(replace the chromosome number by '{chrom}', "
                             "e.g. '1000GP_Phase3_chr{chrom}.legend.gz')"))
    group.add_argument("--map-template", type=str, metavar="TEMPLATE",
                       required=True,
                       help=("The template for IMPUTE2's map files (replace "
                             "the chromosome number by '{chrom}', e.g. "
                             "'genetic_map_chr{chrom}_combined_b37.txt')"))
    group.add_argument("--sample-file", type=str, metavar="FILE",
                       required=True, help="The name of IMPUTE2's sample file")
    group.add_argument("--filtering-rules", type=str, metavar="RULE",
                       nargs="+", help="IMPUTE2 filtering rules (if required)")

    # The impute2 file merger options
    group = parser.add_argument_group("IMPUTE2 Merger Options")
    group.add_argument("--probability", type=float, metavar="FLOAT",
                       default=0.9, help=("The probability threshold for no "
                                          "calls [%(default).1f]"))
    group.add_argument("--completion", type=float, metavar="FLOAT",
                       default=0.98, help=("The site completion rate "
                                           "threshold [%(default).2f]"))

    # The automatic report options
    group = parser.add_argument_group("Automatic Report Options")
    group.add_argument("--report-number", type=str, metavar="NB",
                       default="GWIP automatic report",
                       help="The report number")
    group.add_argument("--report-title", type=str, metavar="TITLE",
                       default="GWIP: Automatic genome-wide imputation",
                       help="The report title")
    group.add_argument("--report-author", type=str, metavar="AUTHOR",
                       default="Automatically generated by GWIP",
                       help="The report author")

    return parser.parse_args()


if __name__ == "__main__":
    main()
