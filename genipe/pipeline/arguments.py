
# This file is part of genipe.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import os
import shlex
import logging
from shutil import which

from ..error import GenipeError
from .. import __version__, autosomes, chromosomes, HAS_PYFAIDX, HAS_DRMAA


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["parse_args", "check_args"]


def parse_args(parser):
    """Parses the command line options and arguments.

    Args:
        parser (argparse.ArgumentParser): the parser object

    Returns:
        argparse.Namespace: the parsed options and arguments

    """
    parser.add_argument(
        "-v", "--version", action="version",
        version="%(prog)s {}".format(__version__),
    )
    parser.add_argument(
        "--debug", action="store_true",
        help="set the logging level to debug",
    )
    parser.add_argument(
        "--thread", type=int, default=1,
        help="number of threads [%(default)d]",
    )

    # The input files
    group = parser.add_argument_group("Input Options")
    group.add_argument(
        "--bfile", type=str, metavar="PREFIX", required=True,
        help="The prefix of the binary pedfiles (input data).",
    )
    group.add_argument(
        "--reference", type=str, metavar="FILE",
        help="The human reference to perform an initial strand check (useful "
             "for genotyped markers not in the IMPUTE2 reference files) "
             "(optional).",
    )

    # The output options
    group = parser.add_argument_group("Output Options")
    possible_chromosomes = chromosomes + ("autosomes", )
    group.add_argument(
        "--chrom", type=str, nargs="+", metavar="CHROM", dest="required_chrom",
        choices=[str(c) for c in possible_chromosomes], default=chromosomes,
        help="The chromosomes to process. It is possible to write 'autosomes' "
             "to process all the autosomes (from chromosome 1 to  22, "
             "inclusively).",
    )
    group.add_argument(
        "--output-dir", type=str, metavar="DIR", default="genipe",
        dest="out_dir",
        help="The name of the output directory. [%(default)s]",
    )
    group.add_argument(
        "--bgzip", action="store_true",
        help="Use bgzip to compress the impute2 files.",
    )

    # The HPC options
    group = parser.add_argument_group("HPC Options")
    group.add_argument(
        "--use-drmaa", action="store_true",
        help="Launch tasks using DRMAA.",
    )
    group.add_argument(
        "--drmaa-config", type=str, metavar="FILE",
        help="The configuration file for tasks (use this option when "
             "launching tasks using DRMAA). This file should describe the "
             "walltime and the number of nodes/processors to use for each "
             "task.",
    )
    group.add_argument(
        "--preamble", type=str, metavar="FILE",
        help="This option should be used when using DRMAA on a HPC to load "
             "required module and set environment variables. The content of "
             "the file will be added between the 'shebang' line and the tool "
             "command.",
    )

    # The SHAPEIT software options
    group = parser.add_argument_group("SHAPEIT Options")
    group.add_argument(
        "--shapeit-bin", type=str, metavar="BINARY",
        help="The SHAPEIT binary if it's not in the path.",
    )
    group.add_argument(
        "--shapeit-thread", type=int, metavar="INT", default=1,
        help="The number of thread for phasing. [%(default)d]",
    )
    group.add_argument(
        "--shapeit-extra", type=str, metavar="OPTIONS",
        help="SHAPEIT extra parameters. Put extra parameters between single "
             "or normal quotes (e.g. --shapeit-extra '--states 100 "
             "--window 2').",
    )

    # The Plink option
    group = parser.add_argument_group("Plink Options")
    group.add_argument(
        "--plink-bin", type=str, metavar="BINARY",
        help="The Plink binary if it's not in the path.",
    )

    # The IMPUTE2 file options
    group = parser.add_argument_group("IMPUTE2 Autosomal Reference")
    group.add_argument(
        "--hap-template", type=str, metavar="TEMPLATE",
        help="The template for IMPUTE2's haplotype files (replace the "
             "chromosome number by '{chrom}', e.g. "
             "'1000GP_Phase3_chr{chrom}.hap.gz').",
    )
    group.add_argument(
        "--legend-template", type=str, metavar="TEMPLATE",
        help="The template for IMPUTE2's legend files (replace the chromosome "
             "number by '{chrom}', e.g. "
             "'1000GP_Phase3_chr{chrom}.legend.gz').",
    )
    group.add_argument(
        "--map-template", type=str, metavar="TEMPLATE",
        help="The template for IMPUTE2's map files (replace the chromosome "
             "number by '{chrom}', e.g. "
             "'genetic_map_chr{chrom}_combined_b37.txt').",
    )
    group.add_argument(
        "--sample-file", type=str, metavar="FILE", required=True,
        help="The name of IMPUTE2's sample file.",
    )

    # The IMPUTE2 sexual chromosome file options
    group = parser.add_argument_group("IMPUTE2 Chromosome X Reference")
    group.add_argument(
        "--hap-nonPAR", type=str, metavar="FILE", dest="hap_chr23",
        help="The IMPUTE2's haplotype file for the non-pseudoautosomal region "
             "of chromosome 23.",
    )
    group.add_argument(
        "--hap-PAR1", type=str, metavar="FILE", dest="hap_par1",
        help="The IMPUTE2's haplotype file for the first pseudoautosomal "
             "region of chromosome 23.",
    )
    group.add_argument(
        "--hap-PAR2", type=str, metavar="FILE", dest="hap_par2",
        help="The IMPUTE2's haplotype file for the second pseudoautosomal "
             "region of chromosome 23.",
    )
    group.add_argument(
        "--legend-nonPAR", type=str, metavar="FILE", dest="legend_chr23",
        help="The IMPUTE2's legend file for the non-pseudoautosomal region "
             "of chromosome 23.",
    )
    group.add_argument(
        "--legend-PAR1", type=str, metavar="FILE", dest="legend_par1",
        help="The IMPUTE2's legend file for the first pseudoautosomal "
             "region of chromosome 23.",
    )
    group.add_argument(
        "--legend-PAR2", type=str, metavar="FILE", dest="legend_par2",
        help="The IMPUTE2's legend file for the second pseudoautosomal "
             "region of chromosome 23.",
    )
    group.add_argument(
        "--map-nonPAR", type=str, metavar="FILE", dest="map_chr23",
        help="The IMPUTE2's map file for the non-pseudoautosomal region "
             "of chromosome 23.",
    )
    group.add_argument(
        "--map-PAR1", type=str, metavar="FILE", dest="map_par1",
        help="The IMPUTE2's map file for the first pseudoautosomal "
             "region of chromosome 23.",
    )
    group.add_argument(
        "--map-PAR2", type=str, metavar="FILE", dest="map_par2",
        help="The IMPUTE2's map file for the second pseudoautosomal "
             "region of chromosome 23.",
    )

    # The IMPUTE2 software options
    group = parser.add_argument_group("IMPUTE2 Options")
    group.add_argument(
        "--impute2-bin", type=str, metavar="BINARY",
        help="The IMPUTE2 binary if it's not in the path.",
    )
    group.add_argument(
        "--segment-length", type=float, metavar="BP", default=5e6,
        help="The length of a single segment for imputation. [%(default).1g]",
    )
    group.add_argument(
        "--filtering-rules", type=str, metavar="RULE", nargs="+",
        help="IMPUTE2 filtering rules (optional).",
    )
    group.add_argument(
        "--impute2-extra", type=str, metavar="OPTIONS",
        help="IMPUTE2 extra parameters. Put the extra parameters between "
             "single or normal quotes (e.g. --impute2-extra '-buffer 250 "
             "-Ne 20000').",
    )

    # The impute2 file merger options
    group = parser.add_argument_group("IMPUTE2 Merger Options")
    group.add_argument(
        "--probability", type=float, metavar="FLOAT", default=0.9,
        help="The probability threshold for no calls. [<%(default).1f]",
    )
    group.add_argument(
        "--completion", type=float, metavar="FLOAT", default=0.98,
        help="The completion rate threshold for site exclusion. "
             "[<%(default).2f]",
    )
    group.add_argument(
        "--info", type=float, metavar="FLOAT", default=0,
        help="The measure of the observed statistical information associated "
             "with the allele frequency estimate threshold for site "
             "exclusion. [<%(default).2f]",
    )

    # The automatic report options
    group = parser.add_argument_group("Automatic Report Options")
    group.add_argument(
        "--report-number", type=str, metavar="NB",
        default="genipe automatic report",
        help="The report number. [%(default)s]",
    )
    group.add_argument(
        "--report-title", type=str, metavar="TITLE",
        default="genipe: Automatic genome-wide imputation",
        help="The report title. [%(default)s]",
    )
    group.add_argument(
        "--report-author", type=str, metavar="AUTHOR",
        default="Automatically generated by genipe",
        help="The report author. [%(default)s]",
    )
    group.add_argument(
        "--report-background", type=str, metavar="BACKGROUND",
        default="The aim of this project is to perform genome-wide imputation "
                "using the study cohort.",
        help="The report background section (can either be a string or a file "
             "containing the background. [General background]",
    )

    return parser.parse_args()


def check_args(args):
    """Checks the arguments and options.

    Args:
        args (argparse.Namespace): the arguments and options

    Returns:
        bool: `True` if everything is OK

    If an option is invalid, a :py:class:`genipe.error.GenipeError` is raised.

    """
    # Checking the presence of the BED, BIM and BAM files
    for suffix in (".bed", ".bim", ".fam"):
        if not os.path.isfile(args.bfile + suffix):
            raise GenipeError("{}: no such file".format(args.bfile + suffix))

    # Checking the thread
    if args.thread < 1:
        raise GenipeError("thread should be one or more")
    if args.shapeit_thread < 1:
        raise GenipeError("thread should be one or more")

    # Checking the chromosome (if autosomes)
    if args.required_chrom == ["autosomes"]:
        args.required_chrom = tuple(autosomes)
    else:
        try:
            args.required_chrom = [int(c) for c in args.required_chrom]
        except ValueError:
            raise GenipeError(
                "{}: invalid chromosome(s) (if all autosomes are required, "
                "write only 'autosomes')".format(args.required_chrom)
            )

    # Checking IMPUTE2's templates (if required)
    for chrom in args.required_chrom:
        if chrom in autosomes:
            if args.hap_template is None:
                raise GenipeError(
                    "chr{} requires '--hap-template'".format(chrom)
                )

            if args.legend_template is None:
                raise GenipeError(
                    "chr{} requires '--legend-template'".format(chrom)
                )

            if args.map_template is None:
                raise GenipeError(
                    "chr{} requires '--map-template'".format(chrom)
                )

            for template in (args.hap_template, args.legend_template,
                             args.map_template):
                filename = template.format(chrom=chrom)
                if not os.path.isfile(filename):
                    raise GenipeError("{}: no such file".format(filename))

    # Checking the non pseudo-autosomal region of chromosome 23
    if 23 in args.required_chrom:
        if args.hap_chr23 is None:
            raise GenipeError("chr23 requires '--hap-nonPAR'")
        if not os.path.isfile(args.hap_chr23):
            raise GenipeError("{}: no such file".format(args.hap_chr23))

        if args.legend_chr23 is None:
            raise GenipeError("chr23 requires '--legend-nonPAR'")
        if not os.path.isfile(args.legend_chr23):
            raise GenipeError("{}: no such file".format(args.legend_chr23))

        if args.map_chr23 is None:
            raise GenipeError("chr23 requires '--map-nonPAR'")
        if not os.path.isfile(args.map_chr23):
            raise GenipeError("{}: no such file".format(args.map_chr23))

    # Checking the pseudo-autosomal region of chromosome 23
    if 25 in args.required_chrom:
        for i in ("1", "2"):
            if vars(args)["hap_par" + i] is None:
                raise GenipeError("chr25 requires '--hap-PAR" + i + "'")
            if not os.path.isfile(vars(args)["hap_par" + i]):
                raise GenipeError(
                    "{}: no such file".format(vars(args)["hap_par" + i])
                )

            if vars(args)["legend_par" + i] is None:
                raise GenipeError("chr25 requires '--legend-PAR" + i + "'")
            if not os.path.isfile(vars(args)["legend_par" + i]):
                raise GenipeError(
                    "{}: no such file".format(vars(args)["legend_par" + i])
                )

            if vars(args)["map_par" + i] is None:
                raise GenipeError("chr25 requires '--map-PAR" + i + "'")
            if not os.path.isfile(vars(args)["map_par" + i]):
                raise GenipeError(
                    "{}: no such file".format(vars(args)["map_par" + i])
                )

    # The final chromosomal requirement
    chrom_names = []
    args.required_chrom = tuple(sorted(args.required_chrom))
    for chrom in args.required_chrom:
        if chrom == 25:
            chrom_names.extend(["25_1", "25_2"])
            continue
        chrom_names.append(chrom)
    args.required_chrom_names = tuple(chrom_names)

    # Checking IMPUTE2's sample file
    if not os.path.isfile(args.sample_file):
        raise GenipeError("{}: no such file".format(args.sample_file))

    # Checking if bgzip is installed, if asking for compression
    if args.bgzip:
        if which("bgzip") is None:
            raise GenipeError("bgzip: no installed")

    # Checking the SHAPEIT binary if required
    if args.shapeit_bin is not None:
        if not os.path.isfile(args.shapeit_bin):
            raise GenipeError("{}: no such file".format(args.shapeit_bin))
    else:
        if which("shapeit") is None:
            raise GenipeError("shapeit: not in the path (use --shapeit-bin)")

    # Checking the IMPUTE2 binary if required
    if args.impute2_bin is not None:
        if not os.path.isfile(args.impute2_bin):
            raise GenipeError("{}: no such file".format(args.impute2_bin))
    else:
        if which("impute2") is None:
            raise GenipeError("impute2: not in the path (use --impute2-bin)")

    # Checking that Plink is in the path
    if args.plink_bin is not None:
        if not os.path.isfile(args.plink_bin):
            raise GenipeError("{}: no such file".format(args.plink_bin))
    else:
        if which("plink") is None:
            raise GenipeError("plink: not in the path (use --plink-bin)")

    # Checking the segment length
    if args.segment_length <= 0:
        raise GenipeError("{}: invalid segment "
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
            raise GenipeError("{}: no such file".format(args.preamble))

    # Checking the DRMAA configuration file
    if args.use_drmaa:
        # Checking the DRMAA module
        if not HAS_DRMAA:
            raise GenipeError("The --use-drmaa option was used, but the drmaa "
                              "module is not installed")

        # Checking the DRMAA_LIBRARY_PATH environment variable
        if "DRMAA_LIBRARY_PATH" not in os.environ:
            raise GenipeError("The DRMAA_LIBRARY_PATH environment variable is "
                              "not set (required by the drmaa module)")

        # Checking the DRMAA configuration file
        if args.drmaa_config is None:
            raise GenipeError("DRMAA configuration file was not provided "
                              "(--drmaa-config), but DRMAA is used "
                              "(--use-drmaa)")
        if not os.path.isfile(args.drmaa_config):
            raise GenipeError("{}: no such file".format(args.drmaa_config))

    # Checking the reference file (if required)
    if args.reference is not None:
        if not HAS_PYFAIDX:
            logging.warning("pyfaidx is not installed, can not perform "
                            "initial strand check")
            args.reference = None

        else:
            if not os.path.isfile(args.reference):
                raise GenipeError("{}: no such file".format(args.reference))

            if not os.path.isfile(args.reference + ".fai"):
                raise GenipeError("{}: should be indexed using "
                                  "FAIDX".format(args.reference))

    # The shapeit extra parameters (if required)
    if args.shapeit_extra is not None:
        # Proofing the command
        args.shapeit_extra = [
            shlex.quote(s) for s in args.shapeit_extra.split(" ")
        ]

        # Checking for protected options
        protected_args = {"-B", "--input-bed", "-M", "--input-map",
                          "-O", "--output-max", "-L", "--output-log",
                          "-phase", "--thread"}
        if len(protected_args & set(args.shapeit_extra)) != 0:
            raise GenipeError(
                "The following SHAPEIT options are hidden from the user: "
                "{}".format(", ".join(sorted(protected_args))),
            )

    # The impute2 extra parameters (if required)
    if args.impute2_extra is not None:
        # Proofing the command
        args.impute2_extra = [
            shlex.quote(s) for s in args.impute2_extra.split(" ")
        ]

        # Checking for protected options
        protected_args = {"-use_prephased_g", "-known_haps_g", "-h", "-l",
                          "-m", "-int", "-o"}
        if len(protected_args & set(args.impute2_extra)) != 0:
            raise GenipeError(
                "The following IMPUTE2 options are hidden from the user: "
                "{}".format(", ".join(sorted(protected_args))),
            )

    return True
