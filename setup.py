#!/usr/bin/env python

# How to build source distribution
#   - python setup.py sdist --format bztar
#   - python setup.py sdist --format gztar
#   - python setup.py sdist --format zip
#   - python setup.py bdist_wheel

# How to build for conda
#   - python setup.py bdist_conda
#   - conda convert -p all .../gwip-0.1-py34_0.tar.bz2 -o dist
#   - cd dist && conda index win-* osx-64 linux-*


import os
from setuptools import setup


MAJOR = 1
MINOR = 0
MICRO = 0
VERSION = "{}.{}.{}".format(MAJOR, MINOR, MICRO)


def write_version_file(fn=None):
    if fn is None:
        fn = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            os.path.join("gwip", "version.py"),
        )

    content = ("\n# THIS FILE WAS GENERATED AUTOMATICALLY BY GWIP SETUP.PY\n"
               'gwip_version = "{version}"\n')

    a = open(fn, "w")
    try:
        a.write(content.format(version=VERSION))
    finally:
        a.close()


def setup_package():
    # Saving the version into a file
    write_version_file()

    setup(
        name="gwip",
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
        url="https://github.com/pgxcentre/gwip",
        license="CC BY-NC 4.0",
        entry_points={
            "console_scripts": [
                "gwip-launcher=gwip.pipeline:main",
                "impute2-merger=gwip.tools.impute2_merger:main",
                "impute2-extractor=gwip.tools.impute2_extractor:main",
                "imputed-stats=gwip.tools.imputed_stats:main",
            ],
        },
        install_requires=["numpy >= 1.8.2", "jinja2 >= 2.7.3",
                          "pandas >= 0.15.2", "setuptools >= 12.0.5"],
        packages=["gwip", "gwip.task", "gwip.db", "gwip.tools", "gwip.formats",
                  "gwip.reporting", "gwip.config", "gwip.tests"],
        package_data={"gwip.reporting": ["templates/*.tex",
                                         "templates/biblio/*",
                                         "templates/utils/*",
                                         "templates/parts/*.tex"],
                      "gwip.tests": ["data/*"],
                      "gwip": ["script_templates/*"]},
        test_suite="gwip.tests.test_suite",
        classifiers=["Development Status :: 5 - Production/Stable",
                     "Intended Audience :: Science/Research",
                     "License :: Free for non-commercial use",
                     "Operating System :: Unix",
                     "Operating System :: POSIX :: Linux",
                     "Operating System :: MacOS :: MacOS X",
                     "Programming Language :: Python",
                     "Programming Language :: Python :: 3.3",
                     "Programming Language :: Python :: 3.4",
                     "Topic :: Scientific/Engineering :: Bio-Informatics"],
        keywords="bioinformatics imputation pipeline analysis",
    )

    return


if __name__ == "__main__":
    setup_package()
