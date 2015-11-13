
# This file is part of genipe, but came from gepyto.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import io
import os
import zlib
import logging

import numpy as np
import pandas as pd

from ..error import GenipeError


__author__ = "Louis-Philippe Lemieux Perreault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


__all__ = ["get_index", "get_open_func"]


_CHECK_STRING = b"GENIPE INDEX FILE"

try:
    from Bio.bgzf import BgzfReader
    HAS_BIOPYTHON = True
except ImportError:
    HAS_BIOPYTHON = False


def _seek_generator(f):
    """Yields seek position for each line.

    Args:
        f (file): the file object

    """
    yield 0
    for line in f:
        yield f.tell()


def generate_index(fn, cols=None, names=None, sep=" "):
    """Build a index for the given file.

    Args:
        fn (str): the name of the file
        cols (list): a list containing column to keep (as int)
        names (list): the name corresponding to the column to keep (as str)
        sep (str): the field separator

    Returns:
        pandas.DataFrame: the index

    """
    logging.info("Generating index for '{}'".format(fn))

    # Some assertions
    assert cols is not None, "'cols' was not set"
    assert names is not None, "'names' was not set"
    assert len(cols) == len(names)

    # Getting the open function
    bgzip, open_func = get_open_func(fn, return_fmt=True)

    # Reading the required columns
    data = pd.read_csv(fn, sep=sep, engine="c", usecols=cols, names=names,
                       compression="gzip" if bgzip else None)

    # Getting the seek information
    f = open_func(fn, "rb")
    data["seek"] = np.fromiter(_seek_generator(f), dtype=np.uint)[:-1]
    f.close()

    # Saving the index to file
    write_index(get_index_fn(fn), data)

    return data


def get_open_func(fn, return_fmt=False):
    """Get the opening function.

    Args:
        fn (str): the name of the file
        return_fmt (bool): if the file format needs to be returned

    Returns:
        tuple: either a tuple containing two elements: a boolean telling if the
               format is bgzip, and the opening function.

    """
    # The file might be compressed using bgzip
    bgzip = None
    with open(fn, "rb") as i_file:
        bgzip = i_file.read(3) == b"\x1f\x8b\x08"

    if bgzip and not HAS_BIOPYTHON:
        raise GenipeError("needs BioPython to index a bgzip file")

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
        raise GenipeError("{}: use bgzip for compression...".format(fn))

    if return_fmt:
        return bgzip, open_func

    return open_func


def get_index(fn, cols, names, sep):
    """Restores the index for a given file.

    Args:
        fn (str): the name of the file
        cols (list): a list containing column to keep (as int)
        names (list): the name corresponding to the column to keep (as str)
        sep (str): the field separator

    Returns:
        pandas.DataFrame: the index

    If the index doesn't exist for the file, it is first created.

    """
    if not has_index(fn):
        # The index doesn't exists, generate it
        return generate_index(fn, cols, names, sep)

    # Retrieving the index
    logging.info("Retrieving the index for '{}'".format(fn))
    file_index = read_index(get_index_fn(fn))

    # Checking the names are there
    if len(set(names) - (set(file_index.columns) - {'seek'})) != 0:
        raise GenipeError("{}: missing index columns: reindex".format(fn))

    if "seek" not in file_index.columns:
        raise GenipeError("{}: invalid index: reindex".format(fn))

    return file_index


def write_index(fn, index):
    """Writes the index to file.

    Args:
        fn (str): the name of the file that will contain the index
        index (pandas.DataFrame): the index

    """
    with open(fn, "wb") as o_file:
        o_file.write(_CHECK_STRING)
        o_file.write(zlib.compress(bytes(
            index.to_csv(None, index=False, encoding="utf-8"),
            encoding="utf-8",
        )))


def read_index(fn):
    """Reads index from file.

    Args:
        fn (str): the name of the file containing the index

    Returns:
        pandas.DataFrame: the index of the file

    Before reading the index, we check the first couple of bytes to see if it
    is a valid index file.

    """
    index = None
    with open(fn, "rb") as i_file:
        if i_file.read(len(_CHECK_STRING)) != _CHECK_STRING:
            raise GenipeError("{}: not a valid index file".format(fn))

        index = pd.read_csv(io.StringIO(
            zlib.decompress(i_file.read()).decode(encoding="utf-8"),
        ))

    return index


def get_index_fn(fn):
    """Generates the index filename from the path to the indexed file.

    Args:
        fn (str): the name of the file for which we want an index

    Returns:
        str: the name of the file containing the index

    """
    return os.path.abspath("{}.idx".format(fn))


def has_index(fn):
    """Checks if the index exists, if not, create it.

    Args:
        fn (str): the name of the file for which we want the index

    Returns:
        bool: ``True`` if the file contains an index, ``False`` otherwise

    """
    return os.path.isfile(get_index_fn(fn))
