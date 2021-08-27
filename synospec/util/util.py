"""
Miscellaneous package utilities.

.. include:: ../include/links.rst
"""

from itertools import chain, combinations

from IPython import embed 

import numpy


def all_subclasses(cls):
    """
    Collect all the subclasses of the provided class.

    The search follows the inheritance to the highest-level class.  Intermediate
    base classes are included in the returned set, but not the base class itself.

    Thanks to:
    https://stackoverflow.com/questions/3862310/how-to-find-all-the-subclasses-of-a-class-given-its-name

    Args:
        cls (object):
            The base class

    Returns:
        :obj:`set`: The unique set of derived classes, including any
        intermediate base classes in the inheritance thread.
    """
    return set(cls.__subclasses__()).union(
            [s for c in cls.__subclasses__() for s in all_subclasses(c)])


def string_table(tbl, delimeter='print', has_header=True):
    """
    Provided the array of data, format it with equally spaced columns
    and add a header (first row) and contents delimeter.

    Args:
        tbl (`numpy.ndarray`_):
            Array of string representations of the data to print.
        delimeter (:obj:`str`, optional):
            If the first row in the table containts the column headers (see
            ``has_header``), this sets the delimeter between first table row and
            the column data. Use ``'print'`` for a simple line of hyphens,
            anything else results in an ``rst`` style table formatting.
        has_header (:obj:`bool`, optional):
            The first row in ``tbl`` contains the column headers.

    Returns:
        :obj:`str`: Single long string with the data table.
    """
    nrows, ncols = tbl.shape
    col_width = [numpy.amax([len(dij) for dij in dj]) for dj in tbl.T]

    _nrows = nrows
    start = 1
    if delimeter != 'print':
        _nrows += 2
        start += 1
    if has_header:
        _nrows += 1
        start += 1

    row_string = ['']*_nrows

    for i in range(start,nrows+start-1):
        row_string[i] = '  '.join([tbl[1+i-start,j].ljust(col_width[j]) for j in range(ncols)])
    if delimeter == 'print':
        # Heading row
        row_string[0] = '  '.join([tbl[0,j].ljust(col_width[j]) for j in range(ncols)])
        # Delimiter
        if has_header:
            row_string[1] = '-'*len(row_string[0])
        return '\n'.join(row_string)+'\n'

    # For an rst table
    row_string[0] = '  '.join([ '='*col_width[j] for j in range(ncols)])
    row_string[1] = '  '.join([tbl[0,j].ljust(col_width[j]) for j in range(ncols)])
    if has_header:
        row_string[2] = row_string[0]
    row_string[-1] = row_string[0]
    return '\n'.join(row_string)+'\n'

