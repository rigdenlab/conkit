"""I/O interface for file reading, writing and conversions"""

__author__ = 'Felix Simkovic'
__date__ = '13 Sep 2016'
__version__ = 0.1

from conkit.io.CaspIO import CaspParser
from conkit.io.CCMpredIO import CCMpredParser
from conkit.io.ComsatIO import ComsatParser
from conkit.io.BbcontactsIO import BbcontactsParser
from conkit.io.EVfoldIO import EVfoldParser
from conkit.io.FreeContactIO import FreeContactParser
from conkit.io.GremlinIO import GremlinParser
from conkit.io.MemBrainIO import MemBrainParser
from conkit.io.PconsIO import PconsParser
from conkit.io.PdbIO import PdbParser
from conkit.io.PlmDCAIO import PlmDCAParser
from conkit.io.PsicovIO import PsicovParser

from conkit.io.A3mIO import A3mIO
from conkit.io.FastaIO import FastaIO
from conkit.io.JonesIO import JonesIO
from conkit.io.StockholmIO import StockholmIO


CONTACT_FILE_PARSERS = {
    'casprr': CaspParser,
    'ccmpred': CCMpredParser,
    'comsat': ComsatParser,
    'bbcontacts': BbcontactsParser,
    'evfold': None,
    'freecontact': FreeContactParser,
    'gremlin': GremlinParser,
    'membrain': MemBrainParser,
    'metapsicov': PsicovParser,
    'pconsc': PconsParser,
    'pconsc2': PconsParser,
    'pconsc3': PconsParser,
    'pdb' : PdbParser,
    'plmdca': PlmDCAParser,
    'psicov': PsicovParser,
}

SEQUENCE_FILE_PARSERS = {
    'a3m': A3mIO,
    'a3m-inserts': A3mIO,
    'fasta': FastaIO,
    'jones': JonesIO,
    'stockholm': StockholmIO,
}


def convert(f_in, format_in, f_out, format_out):
    """Convert a file in format x to file in format y

    Parameters
    ----------
    f_in
       Open file handle for input file [read-permissions]
    format_in : str
       File format of f_in
    f_out
       Open file handle for output file [write-permissions]
    format_out : str
       File format of f_out

    Examples
    --------
    1) Convert a sequence file from A3M format to FASTA format:

    >>> from conkit import io
    >>> io.convert(open('example.a3m', 'r'), 'a3m', open('example.fas', 'w'), 'fasta')

    Notes
    -----
    A3M format comes by default WITHOUT insert states, these are removed. To obtain
    an alignment WITH insert states, use format ``a3m-inserts``.

    2) Convert a PconsC3 contact prediction file to the standard Casp RR format:

    >>> from conkit import io
    >>> io.convert(open('example.out', 'r'), 'pconsc3', open('example.rr', 'w'), 'casprr'))

    """

    if not (format_in in CONTACT_FILE_PARSERS or format_in in SEQUENCE_FILE_PARSERS):
        raise ValueError("Unrecognised input file format: '{selected}'".format(selected=format_in))

    elif not (format_out in CONTACT_FILE_PARSERS or format_out in SEQUENCE_FILE_PARSERS):
        raise ValueError("Unrecognised output file format: '{selected}'".format(selected=format_out))

    elif format_in in CONTACT_FILE_PARSERS and format_out in SEQUENCE_FILE_PARSERS:
        raise ValueError("Cannot convert contact file to sequence file")

    elif format_in in SEQUENCE_FILE_PARSERS and format_out in CONTACT_FILE_PARSERS:
        raise ValueError("Cannot convert sequence file to contact file")

    elif format_in in CONTACT_FILE_PARSERS:
        parser_in = CONTACT_FILE_PARSERS[format_in]()
        parser_out = CONTACT_FILE_PARSERS[format_out]()
        hierarchy = parser_in.read(f_in)
        parser_out.write(f_out, hierarchy)

    elif format_in in SEQUENCE_FILE_PARSERS:
        parser_in = SEQUENCE_FILE_PARSERS[format_in]()
        parser_out = SEQUENCE_FILE_PARSERS[format_out]()
        if format_in == 'a3m-inserts':
            hierarchy = parser_in.read(f_in, remove_inserts=False)
        else:
            hierarchy = parser_in.read(f_in)
        parser_out.write(f_out, hierarchy)

    else:
        raise Exception("Should never be here")

    return


def read(handle, format, f_id='conkit'):
    """Parse a file handle to read into structure

    Parameters
    ----------
    handle
       Open file handle for input file [read-permissions]
    format : str
       File format of handle
    f_id : str
       Identifier for the returned file

    Returns
    -------
    hierarchy
       The hierarchy instance of the requested file

    Examples
    --------
    1) Read a Multiple Sequence Alignment file into a ConKit hierarchy:

    >>> from conkit import io
    >>> hierarchy = io.read(open('example.a3m', 'r'), 'a3m')

    2) Read a contact prediction file into a conkit hierarchy:

    >>> from conkit import io
    >>> hierarchy = io.read(open('example.mat', 'r'), 'ccmpred')

    """

    if not (format in CONTACT_FILE_PARSERS or format in SEQUENCE_FILE_PARSERS):
        raise ValueError("Unrecognised format: '{selected}'".format(selected=format))

    elif format in CONTACT_FILE_PARSERS:
        parser = CONTACT_FILE_PARSERS[format]()

    elif format in SEQUENCE_FILE_PARSERS:
        parser = SEQUENCE_FILE_PARSERS[format]()

    else:
        raise Exception("Should never be here")

    return parser.read(handle, f_id=f_id)


def write(handle, format, hierarchy):
    """Parse a file handle to read into structure

    Parameters
    ----------
    handle
       Open file handle for output file [write-permissions]
    format : str
       File format of handle
    hierarchy
       ConKit hierarchy to write

    Examples
    --------
    1) Write a ConKit hierarchy into a Multiple Sequence Alignment file:

    >>> from conkit import io
    >>> fasta_hierarchy = io.read(open('example.fas', 'r'), 'fasta')
    >>> io.write(open('example.a3m', 'w'), 'a3m', fasta_hierarchy)

    2) Write a ConKit hierarchy into a contact prediction file:

    >>> from conkit import io
    >>> psicov_hierarchy = io.read(open('example.txt', 'r'), 'psicov')
    >>> io.write(open('example.rr', 'w'), 'casprr')

    """

    if not (format in CONTACT_FILE_PARSERS or format in SEQUENCE_FILE_PARSERS):
        raise ValueError("Unrecognised format: '{selected}'".format(selected=format))

    elif format in CONTACT_FILE_PARSERS:
        parser = CONTACT_FILE_PARSERS[format]()

    elif format in SEQUENCE_FILE_PARSERS:
        parser = SEQUENCE_FILE_PARSERS[format]()

    else:
        raise Exception("Should never be here")

    parser.write(handle, hierarchy)

    return
