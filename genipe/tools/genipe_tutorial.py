
# This file is part of genipe.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import os
import sys
import stat
import shutil
import logging
import zipfile
import argparse
import platform
from glob import glob
from urllib.request import urlretrieve
from tempfile import TemporaryDirectory
from distutils.spawn import find_executable
from urllib.error import HTTPError, URLError
from subprocess import check_call, CalledProcessError

from .. import __version__
from ..error import GenipeError


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


# Logging
logging.basicConfig(
    level=logging.INFO,
    format="[%(asctime)s %(name)s %(levelname)s] %(message)s",
    datefmt="%Y-%m-%d %H:%M:%S",
)
logger = logging.getLogger("genipe-tutorial")


_SCRIPT = r"""#!/usr/bin/env bash
# Changing directory
cd {path}

# Launching the imputation with genipe
genipe-launcher \
    --chrom autosomes \
    --bfile {genotypes_prefix} \
    --shapeit-bin {shapeit_bin} \
    --impute2-bin {impute2_bin} \
    --plink-bin {plink_bin} \
    --reference {hg19_fasta} \
    --hap-template {hap_template} \
    --legend-template {legend_template} \
    --map-template {map_template} \
    --sample-file {sample_file} \
    --filtering-rules 'ALL<0.01' 'ALL>0.99' \
    --thread 4 \
    --report-title "Tutorial" \
    --report-number "Test Report"
"""


def main(args=None):
    """The main function.

    Args:
        args (argparse.Namespace): the arguments to be parsed (if
                                   :py:func:`main` is called by another
                                   modulel)

    """
    # Creating the option parser
    desc = ("Prepares the input files for the main pipeline tutorial. This "
            "script is part of the 'genipe' package, "
            "version {}.".format(__version__))
    parser = argparse.ArgumentParser(description=desc)

    try:
        # Parsing the options
        args = parse_args(parser, args)

        # Let's start
        logger.info("Preparing '{}' for genipe tutorial".format(args.path))
        is_ok = input("Proceed (y/n)? ").upper()
        if is_ok != "Y":
            raise KeyboardInterrupt

        # Getting the system type and architecture
        logger.info("Inferring operating system and architecture")
        os_name, architecture = get_os_info()
        logger.info("  - " + os_name)
        logger.info("  - " + architecture + " bits")

        # Creating the directories (if required)
        if not os.path.isdir(args.path):
            os.mkdir(args.path)
        for dirname in ("bin", "data", "hg19"):
            if not os.path.isdir(os.path.join(args.path, dirname)):
                os.mkdir(os.path.join(args.path, dirname))

        # Downloading Plink (if required)
        has_plink = check_files(
            os.path.join(args.path, "bin", "plink"),
        )
        if not has_plink:
            logger.info("Getting Plink")
            get_plink(
                os_name=os_name,
                arch=architecture,
                path=os.path.join(args.path, "bin"),
            )
        else:
            logger.info("Plink already present")

        # Downloading impute2
        has_impute2 = check_files(
            os.path.join(args.path, "bin", "impute2"),
        )
        if not has_impute2:
            logger.info("Getting impute2")
            get_impute2(
                os_name=os_name,
                arch=architecture,
                path=os.path.join(args.path, "bin"),
            )
        else:
            logger.info("Impute2 already present")

        # Downloading shapeit
        has_shapeit = check_files(
            os.path.join(args.path, "bin", "shapeit"),
        )
        if not has_shapeit:
            logger.info("Getting shapeit")
            get_shapeit(
                os_name=os_name,
                arch=architecture,
                path=os.path.join(args.path, "bin"),
            )
        else:
            logger.info("Shapeit already present")

        # Downloading the reference
        has_hg19 = check_files(
            os.path.join(args.path, "hg19", "hg19.fasta"),
            os.path.join(args.path, "hg19", "hg19.fasta.fai"),
        )
        if not has_hg19:
            logger.info("Downloading hg19 reference")
            get_hg19(path=os.path.join(args.path, "hg19"))
        else:
            logger.info("hg19 already downloaded")

        # Downloading the genotypes
        has_geno = check_files(
            os.path.join(args.path, "data", "hapmap_CEU_r23a_hg19.bed"),
            os.path.join(args.path, "data", "hapmap_CEU_r23a_hg19.bim"),
            os.path.join(args.path, "data", "hapmap_CEU_r23a_hg19.fam"),
        )
        if not has_geno:
            logger.info("Downloading genotypes")
            get_genotypes(path=os.path.join(args.path, "data"))
        else:
            logger.info("Genotypes already downloaded")

        # Downloading IMPUTE2 reference files
        has_impute2_ref = check_files(
            os.path.join(args.path, "1000GP_Phase3", "genipe_tut_done"),
        )
        if not has_impute2_ref:
            logger.info("Downloading IMPUTE2's reference files")
            get_impute2_ref(args.path)
        else:
            logger.info("Impute2 reference files already downloaded")

        # Generating the bash script
        generate_bash(args.path)

    # Catching the Ctrl^C
    except KeyboardInterrupt:
        logger.info("Cancelled by user")
        sys.exit(0)

    # Catching the GenipeError
    except GenipeError as e:
        logger.error(e)
        parser.error(e.message)

    except Exception as e:
        logger.error(e)
        raise


def generate_bash(path):
    """Generates a bash script to launch the imputation pipeline.

    Args:
        path (str): the path to write the bash script

    """
    fn = os.path.join(path, "execute.sh")
    path = os.path.abspath(path)
    with open(fn, "w") as f:
        f.write(_SCRIPT.format(
            path=path,
            genotypes_prefix=os.path.join(path, "data",
                                          "hapmap_CEU_r23a_hg19"),
            shapeit_bin=os.path.join(path, "bin", "shapeit"),
            impute2_bin=os.path.join(path, "bin", "impute2"),
            plink_bin=os.path.join(path, "bin", "plink"),
            hg19_fasta=os.path.join(path, "hg19", "hg19.fasta"),
            hap_template=os.path.join(path, "1000GP_Phase3",
                                      "1000GP_Phase3_chr{chrom}.hap.gz"),
            legend_template=os.path.join(path, "1000GP_Phase3",
                                         "1000GP_Phase3_chr{chrom}.legend.gz"),
            map_template=os.path.join(path, "1000GP_Phase3",
                                      "genetic_map_chr{chrom}_combined_b37"
                                      ".txt"),
            sample_file=os.path.join(path, "1000GP_Phase3",
                                     "1000GP_Phase3.sample"),
        ))

    # Making the script executable
    os.chmod(fn, stat.S_IRWXU)


def check_files(*filenames):
    """Checks that all files exists.

    Args:
        filenames (list): the list of file to check

    Returns:
        bool: True if all files exist, False otherwise

    """
    return all(os.path.isfile(fn) for fn in filenames)


def get_os_info():
    """Getting the OS information.

    Returns:
        tuple: first element is the name of the os, and the second is the
               system's architecture

    Note
    ----
        The tutorial does not work on the Windows operating system. The script
        will quit unless the operating system is Linux or Darwin (MacOSX).

    """
    # Getting the OS name
    os_name = platform.system()
    if os_name == "Windows":
        raise GenipeError("Windows OS it not compatible with tutorial")

    # Getting the system's architecture
    architecture = platform.architecture()[0][:2]
    if architecture != "64":
        raise GenipeError("{}: unknown architecture".format(architecture))

    return os_name, architecture


def get_impute2_ref(path):
    """Gets the impute2's reference files.

    Args:
        path (str): the path where to put the reference files

    """
    # The url for each platform
    url = "https://mathgen.stats.ox.ac.uk/impute/{filename}"

    # Getting the name of the file to download
    filename = "1000GP_Phase3.tgz"

    # Downloading Impute2 in a temporary directory
    logger.info("  - " + filename)
    tar_path = os.path.join(path, filename)
    download_file(url.format(filename=filename), tar_path)

    # Extracting genotypes
    logger.info("  - Extracting file")
    untar_file(path, tar_path)

    # Deleting the archive
    os.remove(tar_path)

    # Checking the directory exists
    if not os.path.isdir(os.path.join(path, "1000GP_Phase3")):
        raise GenipeError("Problem extracting the impute2 reference files")

    # Creating an empty file to say that the download and extraction was
    # completed
    done_fn = os.path.join(path, "1000GP_Phase3", "genipe_tut_done")
    with open(done_fn, "w"):
        pass


def get_genotypes(path):
    """Gets the genotypes files.

    Args:
        path (str): the path where to put the genotypes

    """
    # The url for each platform
    url = "http://pgxcentre.github.io/genipe/_static/tutorial/{filename}"

    # Getting the name of the file to download
    filename = "hapmap_CEU_r23a_hg19.tar.bz2"

    # Downloading genotypes in a temporary directory
    logger.info("  - " + filename)
    with TemporaryDirectory() as tmpdir:
        tar_path = os.path.join(tmpdir, filename)
        download_file(url.format(filename=filename), tar_path)

        # Extracting genotypes
        logger.info("  - Extracting file")
        untar_file(tmpdir, tar_path)

        # Finding the genotypes files
        os.remove(os.path.join(tmpdir, filename))
        genotypes_files = glob(os.path.join(tmpdir, "hapmap_CEU_r23a_hg19.*"))
        if len(genotypes_files) != 3:
            raise GenipeError("Unable to locate genotypes")
        for filename in genotypes_files:
            if not os.path.isfile(filename):
                raise GenipeError("Unable to locate genotypes")

        # Moving the files
        for filename in genotypes_files:
            shutil.move(filename, path)


def get_hg19(path):
    """Gets the hg19 reference file.

    Args:
        path (str): the path where to put the reference

    """
    # The url for each platform
    url = ("http://statgen.org/wp-content/uploads/Softwares/genipe/supp_files/"
           "{filename}")

    # Getting the name of the file to download
    filename = "hg19.tar.bz2"

    # Downloading hg19 in a temporary directory
    logger.info("  - " + filename)
    with TemporaryDirectory() as tmpdir:
        tar_path = os.path.join(tmpdir, filename)
        download_file(url.format(filename=filename), tar_path)

        # Extracting the reference
        logger.info("  - Extracting file")
        untar_file(tmpdir, tar_path)

        # Finding the hg19 file
        hg19_files = glob(os.path.join(tmpdir, "hg19.fasta*"))
        if len(hg19_files) != 2:
            raise GenipeError("Unable to locate hg19")
        for filename in hg19_files:
            if not os.path.isfile(filename):
                raise GenipeError("Unable to locate hg19")

        # Moving the files
        for filename in hg19_files:
            shutil.move(filename, path)


def get_plink(os_name, arch, path):
    """Gets Plink depending of the system, and puts it in 'path'.

    Args:
        os_name (str): the name of the OS
        arch (str): the architecture of the system
        path (str): the path where to put Plink

    Note
    ====
        If the binary is in the system path, it is copied to the destination
        path. Otherwise, we download it.

    """
    system_plink = find_executable("plink")
    if system_plink is not None:
        logger.info("  - Copying Plink from {}".format(system_plink))
        shutil.copy(system_plink, path, follow_symlinks=True)

    else:
        # The url for each platform
        url = "http://zzz.bwh.harvard.edu/plink/dist/{filename}"

        # Getting the name of the file to download
        filename = ""
        if os_name == "Darwin":
            filename = "plink-1.07-mac-intel.zip"
        elif os_name == "Linux":
            filename = "plink-1.07-x86_64.zip"
        if filename == "":
            raise GenipeError("Problem choosing a file to download for "
                              "{} {} bits".format(os_name, arch))

        # Downloading Plink in a temporary directory
        logger.info("  - " + filename)
        with TemporaryDirectory() as tmpdir:
            zip_path = os.path.join(tmpdir, filename)
            download_file(url.format(filename=filename), zip_path)

            # Unzipping Plink
            logger.info("  - Extracting file")
            with zipfile.ZipFile(zip_path, "r") as z:
                z.extractall(tmpdir)

            # Finding the plink file
            plink_file = glob(os.path.join(tmpdir, "*", "plink"))
            if len(plink_file) != 1 or not os.path.isfile(plink_file[0]):
                raise GenipeError("Unable to locate Plink")
            plink_file = plink_file[0]

            # Moving the file
            shutil.move(plink_file, path)
            plink_path = os.path.join(path, "plink")

        # Making the script executable
        os.chmod(plink_path, stat.S_IRWXU)


def get_impute2(os_name, arch, path):
    """Gets impute2 depending of the system, and puts it in 'path'.

    Args:
        os_name (str): the name of the OS
        arch (str): the architecture of the system
        path (str): the path where to put impute2

    Note
    ====
        If the binary is in the system path, it is copied to the destination
        path. Otherwise, we download it.

    """
    system_impute2 = find_executable("impute2")
    if system_impute2 is not None:
        logger.info("  - Copying impute2 from {}".format(system_impute2))
        shutil.copy(system_impute2, path, follow_symlinks=True)

    else:
        # The url for each platform
        url = "https://mathgen.stats.ox.ac.uk/impute/{filename}"

        # Getting the name of the file to download
        filename = ""
        if os_name == "Darwin":
            filename = "impute_v2.3.2_MacOSX_Intel.tgz"
        elif os_name == "Linux":
            filename = "impute_v2.3.2_x86_64_static.tgz"
        if filename == "":
            raise GenipeError("Problem choosing a file to download for "
                              "{} {} bits".format(os_name, arch))

        # Downloading Impute2 in a temporary directory
        logger.info("  - " + filename)
        with TemporaryDirectory() as tmpdir:
            tar_path = os.path.join(tmpdir, filename)
            download_file(url.format(filename=filename), tar_path)

            # Extracting impute2
            logger.info("  - Extracting file")
            untar_file(tmpdir, tar_path)

            # Finding the impute2 file
            impute2_file = glob(os.path.join(tmpdir, "*", "impute2"))
            if len(impute2_file) != 1 or not os.path.isfile(impute2_file[0]):
                raise GenipeError("Unable to locate impute2")
            impute2_file = impute2_file[0]

            # Moving the file
            shutil.move(impute2_file, path)
            impute2_path = os.path.join(path, "impute2")

        # Making the script executable
        os.chmod(impute2_path, stat.S_IRWXU)


def get_shapeit(os_name, arch, path):
    """Gets shapeit depending of the system, and puts it in 'path'.

    Args:
        os_name (str): the name of the OS
        arch (str): the architecture of the system
        path (str): the path where to put shapeit

    Note
    ====
        If the binary is in the system path, it is copied to the destination
        path. Otherwise, we download it.

    """
    system_shapeit = find_executable("shapeit")
    if system_shapeit is not None:
        logger.info("  - Copying shapeit from {}".format(system_shapeit))
        shutil.copy(system_shapeit, path, follow_symlinks=True)

    else:
        # The url for each platform
        url = ("https://mathgen.stats.ox.ac.uk/genetics_software/"
               "shapeit/{filename}")

        # Getting the name of the file to download
        filename = ""
        if os_name == "Darwin":
            filename = "shapeit.v2.r837.MacOSX.tgz"
        elif os_name == "Linux":
            filename = "shapeit.v2.r837.GLIBCv2.12.Linux.static.tgz"
        if filename == "":
            raise GenipeError("Problem choosing a file to download for "
                              "{} {} bits".format(os_name, arch))

        # Downloading shapeit in a temporary directory
        logger.info("  - " + filename)
        with TemporaryDirectory() as tmpdir:
            tar_path = os.path.join(tmpdir, filename)
            download_file(url.format(filename=filename), tar_path)

            # Extracting shapeit
            logger.info("  - Extracting file")
            untar_file(tmpdir, tar_path)

            # Finding the shapeit file
            shapeit_file = glob(os.path.join(tmpdir, "*", "shapeit"))
            if len(shapeit_file) != 1 or not os.path.isfile(shapeit_file[0]):
                raise GenipeError("Unable to locate shapeit")
            shapeit_file = shapeit_file[0]

            # Moving the file
            shutil.move(shapeit_file, path)
            shapeit_path = os.path.join(path, "shapeit")

        # Making the script executable
        os.chmod(shapeit_path, stat.S_IRWXU)


def download_file(url, path):
    """Downloads a file from a URL to a path.

    Args:
        url (str): the url to download
        path (str): the path where to save the file

    """
    try:
        urlretrieve(url, path)
    except (HTTPError, URLError):
        raise GenipeError("URL not available: " + url)


def untar_file(path, fn):
    """Extracts a tar archive.

    Args:
        path (str): the path to where the file will be extracted
        fn (str): the name of the tar archive

    """
    try:
        check_call(["tar", "-C", path, "-xf", fn])
    except CalledProcessError:
        raise GenipeError("Could not extract {}".format(fn))


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
    # The options
    parser.add_argument(
        "--version",
        action="version",
        version="%(prog)s, part of genipe version {}".format(__version__),
    )

    parser.add_argument(
        "--tutorial-path",
        type=str,
        metavar="PATH",
        dest="path",
        default=os.path.join(os.environ["HOME"], "genipe_tutorial"),
        help="The path where the tutorial will be run. [%(default)s]",
    )

    if args is not None:
        return parser.parse_args(args)

    return parser.parse_args()


# Calling the main, if necessary
if __name__ == "__main__":
    main()
