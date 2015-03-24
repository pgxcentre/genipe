
# This file is part of gwip, but came from gepyto.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import os
import gzip
import pickle
import logging

import numpy as np
import pandas as pd

from ..error import ProgramError


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["generate_index", "get_index", "goto", "get_index_fn", "has_index",
           "get_open_func"]


_CHECK_STRING = b"GWIP INDEX FILE"

try:
    from Bio.bgzf import BgzfReader
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False


def _seek_generator(f):
    """Yields seek position for each line."""
    yield 0
    for line in f:
        yield f.tell()


def generate_index(fn, cols=None, names=None, sep=" "):
    """Build a index for the given file."""
    logging.info("Generating index for '{}'".format(fn))

    # Some assertions
    assert cols is not None, "'cols' was not set"
    assert names is not None, "'names' was not set"
    assert len(cols) == len(names)

    # Getting the open function
    bgzip, open_func = get_open_func(fn)

    # Reading the required columns
    data = pd.read_csv(fn, sep=sep, engine="c", usecols=cols, names=names,
                       compression="gzip" if bgzip else None)

    # Getting the seek information
    f = open_func(fn, "rb")
    data["seek"] = np.fromiter(_seek_generator(f), dtype=np.uint)[:-1]
    f.close()

    # Saving the index to file
    with gzip.open(get_index_fn(fn), "wb") as o_file:
        o_file.write(_CHECK_STRING)
        pickle.dump(data, o_file)

    return data


def get_open_func(fn):
    """Get the opening function."""
    # The file might be compressed using bgzip
    bgzip = None
    with open(fn, "rb") as i_file:
        bgzip = i_file.read(3) == b"\x1f\x8b\x08"

    if bgzip and not HAS_BIOPYTHON:
        raise ProgramError("needs BioPython to index a bgzip file")

    open_func = open
    if bgzip:
        open_func = BgzfReader

    # Trying to read
    try:
        with open_func(fn, "r") as i_file:
            if bgzip:
                if not i_file.seekable():
                    raise ValueError
            pass

    except ValueError:
        raise ProgramError("{}: use bgzip for compression...".format(fn))

    return bgzip, open_func


def get_index(fn, cols, names, sep):
    """Restores the index for a given file."""
    if not has_index(fn):
        # The index doesn't exists, generate it
        return generate_index(fn, cols, names, sep)

    # Retrieving the index
    logging.info("Retrieving the index for '{}'".format(fn))
    file_index = None
    with gzip.open(get_index_fn(fn), "rb") as i_file:
        if i_file.read(len(_CHECK_STRING)) != _CHECK_STRING:
            raise ProgramError("{}: not a valid index file: "
                               "reindex".format(get_index_fn(fn)))
        file_index = pickle.load(i_file)

    # Checking the names are there
    if len(set(names) - (set(file_index.columns) - {'seek'})) != 0:
        raise ProgramError("{}: missing index columns: reindex".format(fn))

    if "seek" not in file_index.columns:
        raise ProgramError("{}: invalid index: reindex".format(fn))

    return file_index


def goto(f, seek):
    """Given a file and a seek position, go to the required line."""
    pass


def get_index_fn(fn):
    """Generates the index filename from the path to the indexed file."""
    return os.path.abspath("{}.idx".format(fn))


def has_index(fn):
    """Checks if the index exists, if not, create it."""
    return os.path.isfile(get_index_fn(fn))
