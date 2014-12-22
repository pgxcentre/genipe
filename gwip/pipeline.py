"""The main script that do everything."""

import os
import sys
import json
import logging
import argparse
from math import floor
from shutil import which
from urllib.request import urlopen

from .db import *
from . import __version__
from .task import launcher
from .error import ProgramError


def main():
    """The main function.
    
    This is what the pipeline should do:
        1- Exclude markers that have ambiguous alleles [A/T or C/G] (Plink)
        2- Exclude duplicated markers (only keep one) (Plink)
        2- Split the dataset by chromosome (Plink)
        3- Find markers that need to be flip (strand problem) (SHAPEIT)
        4- Flip those markers (Plink)
        5- Find markers with strand problem (SHAPEIT)
        6- Exclude markers with strand problem (Plink)
        7- Phase using SHAPEIT
        8- Impute using IMPUTE2
    
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
            args.task_options = get_task_options()

        # Creating the database
        db_name = create_task_db(args.out_dir)

        # Creating the output directories
        for chromosome in range(1, 24):
            chr_dir = os.path.join(args.out_dir, "chr{}".format(chromosome))
            if not os.path.isdir(chr_dir):
                os.mkdir(chr_dir)

        # Excluding markers prior to phasing (ambiguous markers [A/T and [G/C]
        # and duplicated markers
        exclude_markers_before_phasing(args.bfile, db_name, args)

        # Checking the strand
        check_strand(os.path.join(args.out_dir, "chr{chrom}", "chr{chrom}"),
                     "_1", db_name, args)

        # Flipping the markers
        flip_markers(
            os.path.join(args.out_dir, "chr{chrom}", "chr{chrom}"),
            os.path.join(args.out_dir, "chr{chrom}", "chr{chrom}.to_flip"),
            db_name,
            args,
        )

        # Checking the strand
        check_strand(
            os.path.join(args.out_dir, "chr{chrom}", "chr{chrom}.flipped"),
            "_2",
            db_name,
            args,
            exclude=True,
        )

        # The final marker exclusion
        final_exclusion(
            os.path.join(args.out_dir, "chr{chrom}", "chr{chrom}.flipped"),
            os.path.join(args.out_dir, "chr{chrom}", "chr{chrom}.to_exclude"),
            db_name,
            args,
        )

        # Phasing the data
        phase_markers(
            os.path.join(args.out_dir, "chr{chrom}", "chr{chrom}.final"),
            os.path.join(args.out_dir, "chr{chrom}", "chr{chrom}.final.phased"),
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

    # Catching the Ctrl^C
    except KeyboardInterrupt:
        print("Cancelled by user", file=sys.stderr)
        sys.exit(0)

    # Catching the ProgramError
    except ProgramError as e:
        logging.error(e)
        parser.error(e.message)


def get_task_options():
    """Gets the task options."""
    task_options = {}
    task_names = ["plink_exclude_chr1", "plink_exclude_chr10",
                  "plink_exclude_chr11", "plink_exclude_chr12",
                  "plink_exclude_chr13", "plink_exclude_chr14",
                  "plink_exclude_chr15", "plink_exclude_chr16",
                  "plink_exclude_chr17", "plink_exclude_chr18",
                  "plink_exclude_chr19", "plink_exclude_chr2",
                  "plink_exclude_chr20", "plink_exclude_chr21",
                  "plink_exclude_chr22", "plink_exclude_chr3",
                  "plink_exclude_chr4", "plink_exclude_chr5",
                  "plink_exclude_chr6", "plink_exclude_chr7",
                  "plink_exclude_chr8", "plink_exclude_chr9",
                  "plink_final_exclude_chr1", "plink_final_exclude_chr10",
                  "plink_final_exclude_chr11", "plink_final_exclude_chr12",
                  "plink_final_exclude_chr13", "plink_final_exclude_chr14",
                  "plink_final_exclude_chr15", "plink_final_exclude_chr16",
                  "plink_final_exclude_chr17", "plink_final_exclude_chr18",
                  "plink_final_exclude_chr19", "plink_final_exclude_chr2",
                  "plink_final_exclude_chr20", "plink_final_exclude_chr21",
                  "plink_final_exclude_chr22", "plink_final_exclude_chr3",
                  "plink_final_exclude_chr4", "plink_final_exclude_chr5",
                  "plink_final_exclude_chr6", "plink_final_exclude_chr7",
                  "plink_final_exclude_chr8", "plink_final_exclude_chr9",
                  "plink_flip_chr1", "plink_flip_chr10", "plink_flip_chr11",
                  "plink_flip_chr12", "plink_flip_chr13", "plink_flip_chr14",
                  "plink_flip_chr15", "plink_flip_chr16", "plink_flip_chr17",
                  "plink_flip_chr18", "plink_flip_chr19", "plink_flip_chr2",
                  "plink_flip_chr20", "plink_flip_chr21", "plink_flip_chr22",
                  "plink_flip_chr3", "plink_flip_chr4", "plink_flip_chr5",
                  "plink_flip_chr6", "plink_flip_chr7", "plink_flip_chr8",
                  "plink_flip_chr9", "shapeit_check_chr10_1",
                  "shapeit_check_chr10_2", "shapeit_check_chr11_1",
                  "shapeit_check_chr11_2", "shapeit_check_chr12_1",
                  "shapeit_check_chr12_2", "shapeit_check_chr13_1",
                  "shapeit_check_chr13_2", "shapeit_check_chr14_1",
                  "shapeit_check_chr14_2", "shapeit_check_chr15_1",
                  "shapeit_check_chr15_2", "shapeit_check_chr16_1",
                  "shapeit_check_chr16_2", "shapeit_check_chr17_1",
                  "shapeit_check_chr17_2", "shapeit_check_chr18_1",
                  "shapeit_check_chr18_2", "shapeit_check_chr19_1",
                  "shapeit_check_chr19_2", "shapeit_check_chr1_1",
                  "shapeit_check_chr1_2", "shapeit_check_chr20_1",
                  "shapeit_check_chr20_2", "shapeit_check_chr21_1",
                  "shapeit_check_chr21_2", "shapeit_check_chr22_1",
                  "shapeit_check_chr22_2", "shapeit_check_chr2_1",
                  "shapeit_check_chr2_2", "shapeit_check_chr3_1",
                  "shapeit_check_chr3_2", "shapeit_check_chr4_1",
                  "shapeit_check_chr4_2", "shapeit_check_chr5_1",
                  "shapeit_check_chr5_2", "shapeit_check_chr6_1",
                  "shapeit_check_chr6_2", "shapeit_check_chr7_1",
                  "shapeit_check_chr7_2", "shapeit_check_chr8_1",
                  "shapeit_check_chr8_2", "shapeit_check_chr9_1",
                  "shapeit_check_chr9_2"]

    for task_name in task_names:
        # Creating the task options
        task_options[task_name] = {
            "walltime": bytes("00:15:00", encoding="ascii"),
            "nodes": bytes("-l nodes=1:ppn=1", encoding="ascii"),
        }

    walltimes = {1: "24:00:00", 2: "24:00:00", 3: "20:00:00", 4: "20:00:00",
                 5: "18:00:00", 6: "18:00:00", 7: "15:00:00", 8: "15:00:00",
                 9: "12:00:00", 10: "15:00:00", 11: "15:00:00", 12: "12:00:00",
                 13: "10:00:00", 14: "09:00:00", 15: "08:00:00", 16: "08:00:00",
                 17: "08:00:00", 18: "08:00:00", 19: "06:00:00", 20: "07:00:00",
                 21: "05:00:00", 22: "05:00:00"}
    for chrom in range(1, 23):
        task_name = "shapeit_phase_chr{}".format(chrom)
        # Creating the task options
        task_options[task_name] = {
            "walltime": bytes(walltimes[chrom], encoding="ascii"),
            "nodes": bytes("-l nodes=1:ppn=1", encoding="ascii"),
        }

    # The time for an impute2 segment
    task_options["impute"] = {
        "walltime": bytes("120:00:00", encoding="ascii"),
        "nodes": bytes("-l nodes=3:ppn=1", encoding="ascii"),
    }

    return task_options


def phase_markers(prefix, o_prefix, db_name, options):
    """Phase markers using shapeit."""
    commands_info = []
    base_command = [
        "shapeit" if options.shapeit_bin is None else options.shapeit_bin,
        "-phase",
        "--thread", str(options.shapeit_thread),
    ]

    for chrom in range(1, 23):
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
                          out_dir=options.out_dir)
    logging.info("Done phasing markers")


def impute_markers(phased_haplotypes, out_prefix, chrom_length, db_name,
                   options):
    """Imputes the markers using IMPUTE2."""
    commands_info = []
    base_command = [
        "impute2" if options.impute2_bin is None else options.impute2_bin,
        "-use_prephased_g",
        "-Ne", "20000",
    ]

    # Each chromosome have multiple segments
    for chrom in range(1, 23):
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
                "name": "IMPUTE2 chr{} from {} to {}".format(chrom, start, end),
                "command": base_command + remaining_command,
                "task_db": db_name,
                "o_files": [c_prefix + "_summary", ],
            })

            # The new starting position
            start = end + 1

            # Adding the walltime for this particular task_id
            if options.use_drmaa:
                options.task_options[task_id] = options.task_options["impute"]

    # Executing the commands
    logging.info("Imputing markers")
    launcher.launch_tasks(commands_info, options.thread, hpc=options.use_drmaa,
                          hpc_options=options.task_options,
                          out_dir=options.out_dir)
    logging.info("Done imputing markers")


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
        chromosomes = {str(i) for i in range(23)} | {"X"}
        for region in result["top_level_region"]:
            if region["name"] in chromosomes:
                chrom_length[region["name"]] = region["length"]

        # Checking we have all the required data
        if len(chrom_length) != 23:
            raise ProgramError("missing chromosomes")

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
                chrom_length[row[0]] = row[1]

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

    for chrom in range(1, 23):
        # The current output prefix
        c_prefix = o_prefix.format(chrom=chrom)

        remaining_command = [
            "-B", prefix.format(chrom=chrom),
            "-M", options.map_template.format(chrom=chrom),
            "--input-ref", options.hap_template.format(chrom=chrom),
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
                          out_dir=options.out_dir)
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
    for chrom in range(1, 23):
        chrom_filename = filename.format(chrom=chrom)
        chrom_o_filename = o_filename.format(chrom=chrom, o_suffix=o_suffix)

        # Checking the input file exists
        if not os.path.isfile(chrom_filename):
            raise ProgramError("{}: no such file".format(chrom_filename))

        # Markers to flip
        to_flip = 0
        with open(chrom_filename, "r") as i_file,\
             open(chrom_o_filename, "w") as o_file:
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
                    print(row[header["main_id"]], file=o_file)
                    to_flip += 1

        nb_total += to_flip
        logging.info("chr{}: {:,d} markers to {}".format(chrom, to_flip, what))

    # Logging the last one
    logging.info("After strand check: {:,d} markers "
                 "to {}".format(nb_total, what))


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
    o_prefix = os.path.join(options.out_dir, "chr{chrom}", "chr{chrom}.flipped")

    for chrom in range(1, 23):
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
                          out_dir=options.out_dir)
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

    for chrom in range(1, 23):
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

    # Executing command
    logging.info("Final marker exclusion")
    launcher.launch_tasks(commands_info, options.thread, hpc=options.use_drmaa,
                          hpc_options=options.task_options,
                          out_dir=options.out_dir)
    logging.info("Done final marker exclusion")


def exclude_markers_before_phasing(prefix, db_name, options):
    """Finds and excludes ambiguous markers (A/T and G/C) or duplicated ones."""
    # The ambiguous genotypes
    ambiguous_genotypes = {"AT", "TA", "GC", "CG"}

    # The positions that have been written to file
    kept_positions = set()
    nb_ambiguous = 0
    nb_dup = 0

    # Logging
    logging.info("Finding markers to exclude")

    o_filename = os.path.join(options.out_dir, "markers_to_exclude.txt")
    with open(o_filename, "w") as o_file, \
         open(prefix + ".bim", "r") as i_file:
        for line in i_file:
            row = line.rstrip("\r\n").split("\t")

            # Checking the alleles
            if row[4] + row[5] in ambiguous_genotypes:
                nb_ambiguous += 1
                print(row[1], file=o_file)
                logging.debug("  - {}: {}: "
                             "ambiguous".format(row[1], row[4] + row[5]))
                continue

            # Checking if we already have this marker
            if (row[0], row[3]) in kept_positions:
                nb_dup += 1
                print(row[1], file=o_file)
                logging.debug("  - {}: duplicated".format(row[1]))
                continue

            # We keep this marker
            kept_positions.add((row[0], row[3]))

    # Logging
    logging.info("  - {:,d} markers kept".format(len(kept_positions)))
    logging.info("  - {:,d} ambiguous markers removed".format(nb_ambiguous))
    logging.info("  - {:,d} duplicated markers removed".format(nb_dup))

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

    for chrom in range(1, 23):
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
                          out_dir=options.out_dir)
    logging.info("Done excluding and splitting markers")


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
    for template in (args.hap_template, args.legend_template, args.map_template):
        for chrom in range(1, 23):
            # Checking the haplotype file
            filename = template.format(chrom=chrom)
            if not os.path.isfile(filename):
                raise ProgramError("{}: no such file".format(filename))
    if not os.path.isfile(args.sample_file):
        raise ProgramError("{}: no such file".format(filename))

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
            raise ProgramError("plink: not in the path")

    # Checking the segment length
    if args.segment_length < 0:
        raise ProgramError("{}: invalid segment "
                           "length".format(args.segment_length))
    if args.segment_length > 5e6:
        # This is too big... We continue with a warning
        logging.warning("segment length ({:g} bp) is more than "
                        "5Mb".format(args.segment_length))

    return True


def parse_args(parser):
    """Parses the command line options and arguments."""
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {}".format(__version__),
    )
    parser.add_argument(
        "--debug",
        action="store_true",
        help="Set the logging level to debug",
    )
    parser.add_argument(
        "--thread",
        type=int,
        default=1,
        help="The number of thread [%(default)d]",
    )

    # The input files
    group = parser.add_argument_group("Input Options")
    group.add_argument(
        "--bfile",
        type=str,
        metavar="PREFIX",
        required=True,
        help="The prefix of the binary pedfiles (input data)",
    )

    # The output options
    group = parser.add_argument_group("Output Options")
    group.add_argument(
        "--output-dir",
        type=str,
        metavar="DIR",
        default="gwip",
        dest="out_dir",
        help="The name of the output directory [%(default)s]",
    )

    # The HPC options
    group = parser.add_argument_group("HPC Options")
    group.add_argument(
        "--use-drmaa",
        action="store_true",
        help="Launch tasks using DRMAA",
    )

    # The SHAPEIT software options
    group = parser.add_argument_group("SHAPEIT Options")
    group.add_argument(
        "--shapeit-bin",
        type=str,
        metavar="BINARY",
        help="The SHAPEIT binary if it's not in the path",
    )
    group.add_argument(
        "--shapeit-thread",
        type=int,
        default=1,
        help="The number of thread for phasing [%(default)d]",
    )

    # The Plink option
    group = parser.add_argument_group("Plink Options")
    group.add_argument(
        "--plink-bin",
        type=str,
        metavar="BINARY",
        help="The Plink binary if it's not in the path",
    )

    # The IMPUTE2 software options
    group = parser.add_argument_group("IMPUTE2 Options")
    group.add_argument(
        "--impute2-bin",
        type=str,
        metavar="BINARY",
        help="The IMPUTE2 binary if it's not in the path",
    )
    group.add_argument(
        "--segment-length",
        type=float,
        metavar="BP",
        default=5e6,
        help=("The length of a single segment for imputation "
              "[%(default).1g)]"),
    )
    group.add_argument(
        "--hap-template",
        type=str,
        metavar="TEMPLATE",
        required=True,
        help=("The template for IMPUTE2's haplotype files (replace the "
              "chromosome number by '{chrom}', e.g. "
              "'1000GP_Phase3_chr{chrom}.hap.gz')"),
    )
    group.add_argument(
        "--legend-template",
        type=str,
        metavar="TEMPLATE",
        required=True,
        help=("The template for IMPUTE2's legend files (replace the chromosome "
              "number by '{chrom}', e.g. "
              "'1000GP_Phase3_chr{chrom}.legend.gz')"),
    )
    group.add_argument(
        "--map-template",
        type=str,
        metavar="TEMPLATE",
        required=True,
        help=("The template for IMPUTE2's map files (replace the chromosome "
              "number by '{chrom}', e.g. "
              "'genetic_map_chr{chrom}_combined_b37.txt')"),
    )
    group.add_argument(
        "--sample-file",
        type=str,
        metavar="FILE",
        required=True,
        help="The name of IMPUTE2's sample file",
    )

    return parser.parse_args()


if __name__ == "__main__":
    main()

