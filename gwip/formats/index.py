
# This file is part of gwip, but came from gepyto.
#
# This work is licensed under the Creative Commons Attribution-NonCommercial
# 4.0 International License. To view a copy of this license, visit
# http://creativecommons.org/licenses/by-nc/4.0/ or send a letter to Creative
# Commons, PO Box 1866, Mountain View, CA 94042, USA.


import os
import re
import bisect
import pickle
import logging

import numpy as np


__author__ = "Marc-Andre Legault"
__copyright__ = "Copyright 2014, Beaulieu-Saucier Pharmacogenomics Centre"
__license__ = "Attribution-NonCommercial 4.0 International (CC BY-NC 4.0)"


def build_index(fn, chrom_col, pos_col, delimiter='\t', skip_lines=0,
                index_rate=0.2, ignore_startswith=None):
    """Build a index for the given file.

    :param fn: The filename
    :type fn: str

    :param chrom_col: The column representing the chromosome (0 based).
    :type chrom_col: int

    :param pos_col: The column for the position on the chromosome (0 based).
    :type pos_col: int

    :param delimiter: The delimiter for the columns (default tab).
    :type delimiter: str

    :param skip_lines: Number of header lines to skip.
    :type skip_lines: int

    :param index_rate: The approximate rate of line indexing. As an example,
                       a file with 1000 lines and the default index_rate of
                       0.2 will have an index with ~200 entries.
    :type index_rate: float

    :param ignore_startswith: Ignore lines that start with a given string.
                              This can be used to skip headers, but will not
                              be used to parse the rest of the file.
    :type ignore_startswith: str

    :returns: The index dict.
    :rtype: dict

    Basically, the index will be a pickle in the same directory.
    The underlying data structure is:

    {chrom: [(pos, index), ...], ...  }

    Note that we make sure the indexed positions are increasing (file sorted)
    but we do not look at all the positions so the check is very partial.

    """

    idx_fn = get_index_fn(fn)
    idx = {}

    size = os.path.getsize(fn)  # Total filesize
    start = 0  # Byte position of the meat of the file (data).
    with open(fn, "r") as f:
        if skip_lines > 0:
            for i in range(skip_lines):
                # Skip header lines if needed.
                _ = f.readline()

        cur = f.tell()
        if ignore_startswith is not None:
            line = f.readline()
            while line.startswith(ignore_startswith):
                cur = f.tell()
                line = f.readline()
            f.seek(cur)

        start = f.tell()

        size -= start  # Adjust file size to remove header.

        # Estimate the line length using first 100 lines
        line_length = np.empty((100))
        prev = start
        for i in range(100):
            _ = f.readline()
            line_length[i] = f.tell() - prev
            prev = f.tell()

        line_length = np.mean(line_length)
        approx_num_lines = size / line_length

        # Go back to start
        f.seek(start)

        # Add the first line to the index.
        l1 = f.readline().rstrip().split(delimiter)
        chrom1, pos1 = (l1[chrom_col], l1[pos_col])
        chrom1 = chrom1.lstrip("chr")
        pos1 = int(pos1)
        idx[chrom1] = [(pos1, start), ]
        f.seek(start)

        # Compute the seek jump size.
        target_num_lines = index_rate * approx_num_lines
        seek_jump = size / target_num_lines

        # Now we will jump of seek_jump, get to a fresh line (because we will
        # probably land in the middle of a line).
        cur = f.tell() + seek_jump

        # We need to remember the previous position to make sure the file is
        # sorted.
        prev_chrom = None
        prev_pos = None

        while cur + seek_jump < size:
            # Jump in the file.
            f.seek(cur)
            # Throw away partial line.
            _ = f.readline()

            # Make sure we did not end at the end of a sequence.
            l1 = f.readline().rstrip().split(delimiter)
            l2_tell = f.tell()
            l2 = f.readline().rstrip().split(delimiter)

            if len(l1) == 1 or len(l2) == 1:
                # We're at the end of the file.
                assert f.readline() == ""
                break

            chrom1, pos1 = (l1[chrom_col], l1[pos_col])
            chrom2, pos2 = (l2[chrom_col], l2[pos_col])

            chrom1 = chrom1.lstrip("chr")
            chrom2 = chrom2.lstrip("chr")

            pos1 = int(pos1)
            pos2 = int(pos2)

            if chrom1 == chrom2 and pos1 == pos2:
                # We were in a repeated region so we will move further.
                cur += 2 * line_length
            else:
                # We know that l2 is the first occurence of this position
                # so we will index it.
                if chrom2 not in idx:
                    idx[chrom2] = []
                idx[chrom2].append((pos2, l2_tell))
                cur = l2_tell + seek_jump

                if prev_chrom is None:
                    prev_chrom = chrom2
                    prev_pos = pos2
                else:
                    # Make sure file is sorted.
                    unsrt = False
                    if prev_chrom == chrom2:
                        if prev_pos > pos2:
                            unsrt = True
                    else:
                        if prev_chrom > chrom2:
                            unsrt = True
                    if unsrt:
                        raise Exception("File is not sorted.")

    # Dump the index to disk.
    idx = {
        "idx": idx,
        "delim": delimiter,
        "chrom_col": chrom_col,
        "pos_col": pos_col,
    }
    with open(idx_fn, "wb") as f:
        pickle.dump(idx, f)

    return idx


def get_index(fn):
    """Restores the index for a given file.

    :param fn: The filname of the file to index.
    :type fn: str

    :returns: The index dictionary corresponding to the input file.
    :rtype: dict

    """

    with open(get_index_fn(fn), "rb") as f:
        idx = pickle.load(f)
    return idx


def goto(f, idx, chrom, pos):
    """Given a file and an index, go to the genomic coordinates.

    :param f: An open file.
    :type f: file

    :param idx: The index data structure (see build_index).
    :rtype idx: dict

    :param chrom: The queried chromosome.
    :param pos: The queried position on the chromosome.

    :returns: True if the position was found and the cursor moved, False if
              the queried chromosome, position wasn't found.
    :rtype: bool

    See https://docs.python.org/2/library/bisect.html

    """

    # Those are the types used by the indexing structure.
    chrom = str(chrom)
    chrom = chrom.lstrip("chr")
    num_chrom = None
    if re.match(r"^[0-9]+$", chrom):
        num_chrom = int(chrom)
    pos = int(pos)

    # Extract important information on index structure.
    delim = idx["delim"]
    chrom_col = idx["chrom_col"]
    pos_col = idx["pos_col"]
    idx = idx["idx"]

    # Look if we can use chromosome indexing.
    if chrom in idx:
        li = [i[0] for i in idx[chrom]]  # Sorted list of pos
        i = bisect.bisect_right(li, pos)
    else:
        i = None

    if not i:
        # We didn't find something before in the index for this chromosome.
        # OR this chromosome is not in the index.
        # We will look at the end of the previous chromosome if it is a
        # numeric chromosome.
        if num_chrom is not None:
            num_prev_chrom = num_chrom - 1
            while str(num_prev_chrom) not in idx and num_prev_chrom > 1:
                num_prev_chrom -= 1
            prev_chrom = str(num_prev_chrom)

            # A previous chromosome is in the index, start from there.
            if prev_chrom in idx:
                tell = idx[prev_chrom][-1][1]
            # This chromosome and the previous are not in the index, start
            # at the beginning.
            elif num_chrom != 1 and chrom not in idx:
                logging.warning("Sparse coverage index, try increasing the "
                                "index rate. Neither chromosome {} or "
                                "{} were sampled in the index.".format(
                                    num_chrom, prev_chrom
                                ))
                tell = 0  # This could be bad.
            # This is the first chromosome, start from the beginning, but don't
            # print a warning because it will probably be up top.
            elif num_chrom == 1:
                tell = 0
            else:
                raise Exception("The observed index structure was not "
                                "anticipated, please report it (query: "
                                "{}, {}).".format(chrom, pos))
        else:
            # This is a non numeric chromosome.
            # Find the last numeric chromosome and look after.
            last_num_chrom = None
            pat = r"^[0-9]+$"
            keys = sorted([int(i) for i in idx.keys() if re.match(pat, i)])
            last_num_chrom = str(keys[-1])
            # Set the tell to the last index of this chromosome.
            tell = idx[last_num_chrom][-1][1]
    else:
        tell = idx[chrom][i - 1][1]

    f.seek(tell)

    while True:
        # Iterate through the file.
        here = f.tell()
        line = f.readline().rstrip().split(delim)
        if len(line) == 1:
            # End of file.
            return False

        this_chrom = line[chrom_col]
        this_chrom = this_chrom.lstrip("chr")

        if re.match(r"^[0-9]+$", this_chrom):
            num_this_chrom = int(this_chrom)
        else:
            num_this_chrom = None

        this_pos = int(line[pos_col])

        # We're too far, didn't find it.
        if num_this_chrom is not None and num_chrom is not None:
            if num_this_chrom > num_chrom:
                return False
            elif num_this_chrom == num_chrom and this_pos > pos:
                return False

        # We got it!
        if this_chrom == chrom and this_pos == pos:
            f.seek(here)
            return True


def get_index_fn(fn):
    """Generates the index filename from the path to the indexed file.

    :param fn: The name of the file to index.
    :type fn: str

    """
    return os.path.abspath("{}.gtidx".format(fn))
