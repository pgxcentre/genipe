#!/usr/bin/env python

# How to build source distribution
# python setup.py sdist --format bztar
# python setup.py sdist --format gztar
# python setup.py sdist --format zip


import os
from setuptools import setup

from gwip import __version__


MAJOR = 0
MINOR = 1
VERSION = "{}.{}".format(MAJOR, MINOR)


def write_version_file(fn=os.path.join("gwip", "version.py")):
    content = """
# THIS FILE WAS GENERATED AUTOMATICALLY BY GWIP SETUP.PY
gwip_version = {version}
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
        license="MIT",
        entry_points={
            "console_scripts": [
                "gwip-launcher=gwip.pipeline:main",
                "impute2-merger=gwip.tool.impute2_merger:main",
            ],
        },
        install_requires=["numpy >= 1.8.2", "jinja2 >= 2.7.3", ],
        packages=["gwip", "gwip.task", "gwip.db", "gwip.tool",
                  "gwip.reporting", ],
        package_data={"gwip.reporting": ["template/*.tex", ], },
        classifiers=["Operating System :: Linux",
                     "Programming Language :: Python",
                     "Programming Language :: Python :: 3", ],
    )

    return


if __name__ == "__main__":
    setup_package()
