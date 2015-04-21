#!/usr/bin/env python

# How to build source distribution
# python setup.py sdist --format bztar
# python setup.py sdist --format gztar
# python setup.py sdist --format zip


import os
from setuptools import setup


MAJOR = 0
MINOR = 1
VERSION = "{}.{}".format(MAJOR, MINOR)


def write_version_file(fn=os.path.join("gwip", "version.py")):
    content = """
# THIS FILE WAS GENERATED AUTOMATICALLY BY GWIP SETUP.PY
gwip_version = "{version}"
"""
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
        description="PGx Genome-Wide Imputation Pipeline.",
        author="Louis-Philippe Lemieux Perreault",
        author_email="louis-philippe.lemieux.perreault@statgen.org",
        url="http://www.statgen.org",
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
        classifiers=["Operating System :: Linux",
                     "Programming Language :: Python",
                     "Programming Language :: Python :: 3", ],
    )

    return


if __name__ == "__main__":
    setup_package()
