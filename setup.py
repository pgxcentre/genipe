#!/usr/bin/env python

# How to build source distribution
# python setup.py sdist --format bztar
# python setup.py sdist --format gztar
# python setup.py sdist --format zip

import os
from setuptools import setup

from gwip import __version__

setup(
    name="gwip",
    version=__version__,
    description="Genome-Wide Imputation Pipeline.",
    author="Louis-Philippe Lemieux Perreault",
    author_email="louis-philippe.lemieux.perreault@statgen.org",
    url="http://www.statgen.org",
    license="GPL",
    entry_points = {
        "console_scripts": ["gwip-launcher=gwip.pipeline:main",
                            "impute2-merger=gwip.tool.impute2_merger:main",],
    },
    install_requires=["numpy >= 1.8.2", ],
    packages=["gwip", "gwip.task", "gwip.db", "gwip.tool", ],
    classifiers=['Operating System :: Linux',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 3.4'],
)
