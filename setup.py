#!/usr/bin/env python

# How to build source distribution
#   - python setup.py sdist --format bztar
#   - python setup.py sdist --format gztar
#   - python setup.py sdist --format zip
#   - python setup.py bdist_wheel

# How to build for conda
#   - bash conda_build.sh


import os
import sys
from setuptools import setup


MAJOR = 1
MINOR = 4
MICRO = 1
VERSION = "{0}.{1}.{2}".format(MAJOR, MINOR, MICRO)


def check_python_version():
    """Checks the python version, exits if < 3.4."""
    python_major, python_minor = sys.version_info[:2]

    if python_major != 3 or python_minor < 4:
        sys.stderr.write("genipe requires python 3 (version 3.4 or higher)\n")
        sys.exit(1)


def write_version_file(fn=None):
    if fn is None:
        fn = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            os.path.join("genipe", "version.py"),
        )

    content = ("\n# THIS FILE WAS GENERATED AUTOMATICALLY BY GENIPE SETUP.PY\n"
               'genipe_version = "{version}"\n')

    a = open(fn, "w")
    try:
        a.write(content.format(version=VERSION))
    finally:
        a.close()


def setup_package():
    # Checking the python version prior to installation
    check_python_version()

    # Saving the version into a file
    write_version_file()

    setup(
        name="genipe",
        version=VERSION,
        description="An automatic genome-wide imputation pipeline.",
        long_description=("This package provides tools to automatically "
                          "perform a genome-wide imputation analysis, "
                          "including the different imputation steps using "
                          "well known softwares, as well as downstream "
                          "statistical analysis. It also provides an "
                          "automatic report (using LaTeX), showing different "
                          "quality metrics about the imputation process."),
        author="Louis-Philippe Lemieux Perreault",
        author_email="louis-philippe.lemieux.perreault@statgen.org",
        url="https://github.com/pgxcentre/genipe",
        license="CC BY-NC 4.0",
        entry_points={
            "console_scripts": [
                "genipe-launcher=genipe.pipeline.cli:main",
                "impute2-merger=genipe.tools.impute2_merger:main",
                "impute2-extractor=genipe.tools.impute2_extractor:main",
                "imputed-stats=genipe.tools.imputed_stats:main",
                "genipe-tutorial=genipe.tools.genipe_tutorial:main",
            ],
        },
        install_requires=["numpy >= 1.11", "Jinja2 >= 2.9",
                          "pandas >= 0.19", "setuptools >= 12.0.5"],
        packages=["genipe", "genipe.pipeline", "genipe.task", "genipe.db",
                  "genipe.tools", "genipe.formats", "genipe.reporting",
                  "genipe.config", "genipe.tests"],
        package_data={"genipe.reporting": ["templates/*.tex",
                                           "templates/biblio/*",
                                           "templates/utils/*",
                                           "templates/parts/*.tex"],
                      "genipe.tests": ["data/*"],
                      "genipe": ["script_templates/*"]},
        test_suite="genipe.tests.test_suite",
        zip_safe=False,
        classifiers=["Development Status :: 5 - Production/Stable",
                     "Intended Audience :: Science/Research",
                     "License :: Free for non-commercial use",
                     "Operating System :: Unix",
                     "Operating System :: POSIX :: Linux",
                     "Operating System :: MacOS :: MacOS X",
                     "Programming Language :: Python",
                     "Programming Language :: Python :: 3.4",
                     "Programming Language :: Python :: 3.5",
                     "Topic :: Scientific/Engineering :: Bio-Informatics"],
        keywords="bioinformatics imputation pipeline analysis",
    )

    return


if __name__ == "__main__":
    setup_package()
