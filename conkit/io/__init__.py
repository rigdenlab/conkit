"""I/O interface for file reading, writing and conversions"""

__author__ = 'Felix Simkovic'
__date__ = '13 Sep 2016'
__version__ = 0.1

from conkit.io.CaspIO import CaspParser
from conkit.io.CCMpredIO import CCMpredParser
from conkit.io.ComsatIO import ComsatParser
from conkit.io.BbcontactsIO import BbcontactsParser
from conkit.io.BCLContactIO import BCLContactParser
from conkit.io.EPCMapIO import EPCMapParser
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
from conkit.io._iotools import open_f_handle


CONTACT_FILE_PARSERS = {
    'casprr': CaspParser,
    'ccmpred': CCMpredParser,
    'comsat': ComsatParser,
    'bbcontacts': BbcontactsParser,
    'bclcontact': BCLContactParser,
    'epcmap': EPCMapParser,
    'evfold': EVfoldParser,
    'freecontact': FreeContactParser,
    'gremlin': GremlinParser,
    'membrain': MemBrainParser,
    'metapsicov': PsicovParser,
    'pconsc': PconsParser,
    'pconsc2': PconsParser,
    'pconsc3': PconsParser,
    'pdb': PdbParser,
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


def convert(fname_in, format_in, fname_out, format_out):
    """Convert a file in format x to file in format y

    Parameters
    ----------
    fname_in : filehandle, filename
    format_in : str
       File format of f_in
    fname_out : filehandle, filename
    format_out : str
       File format of f_out

    Examples
    --------
    1) Convert a sequence file from A3M format to FASTA format:

    >>> from conkit import io
    >>> with open('example.a3m', 'r') as f_in, open('example.fas', 'w') as f_out:
    ...     io.convert(f_in, 'a3m', f_out, 'fasta')

    Notes
    -----
    A3M format comes by default WITHOUT insert states, these are removed. To obtain
    an alignment WITH insert states, use format ``a3m-inserts``.

    2) Convert a PconsC3 contact prediction file to the standard Casp RR format:

    >>> from conkit import io
    >>> with open('example.out', 'r') as f_in, open('example.rr', 'w') as f_out:
    ...     io.convert(f_in, 'pconsc3', f_out, 'casprr'))

    """
    # Check for the correct format and values provided
    kwargs = {}
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
    elif format_in in SEQUENCE_FILE_PARSERS:
        parser_in = SEQUENCE_FILE_PARSERS[format_in]()
        parser_out = SEQUENCE_FILE_PARSERS[format_out]()
        if format_in == 'a3m-inserts':
            kwargs['remove_inserts'] = False
    else:
        raise Exception("Should never be here")

    with open_f_handle(fname_in, 'read') as f_in, open_f_handle(fname_out, 'write') as f_out: 
        hierarchy = parser_in.read(f_in, **kwargs)
        parser_out.write(f_out, hierarchy)

    return


def read(f_in, format, f_id='conkit'):
    """Parse a file handle to read into structure

    Parameters
    ----------
    f_in
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
    >>> with open('example.a3m', 'r') as f_in:
    ...     hierarchy = io.read(f_in, 'a3m')

    2) Read a contact prediction file into a conkit hierarchy:

    >>> from conkit import io
    >>> with open('example.mat', 'r') as f_in:
    ...     hierarchy = io.read(f_in, 'ccmpred')

    """
    # Check for the correct format and values provided
    if not (format in CONTACT_FILE_PARSERS or format in SEQUENCE_FILE_PARSERS):
        raise ValueError("Unrecognised format: '{selected}'".format(selected=format))
    elif format in CONTACT_FILE_PARSERS:
        parser = CONTACT_FILE_PARSERS[format]()
    elif format in SEQUENCE_FILE_PARSERS:
        parser = SEQUENCE_FILE_PARSERS[format]()
    else:
        raise Exception("Should never be here")
    
    with open_f_handle(f_in, 'read') as f_in:
        hierarchy = parser.read(f_in, f_id=f_id)

    return hierarchy


def write(fname, format, hierarchy):
    """Parse a file handle to read into structure

    Parameters
    ----------
    fname : filehandle, filename
    format : str
       File format of handle
    hierarchy
       ConKit hierarchy to write

    Examples
    --------
    1) Write a ConKit hierarchy into a Multiple Sequence Alignment file:

    >>> from conkit import io
    >>> with open('example.fas', 'r') as f_in, open('example.a3m', 'w') as f_out:
    ...     hierarchy = io.read(f_in, 'fasta')
    ...     io.write(f_out, 'a3m', hierarchy)

    2) Write a ConKit hierarchy into a contact prediction file:

    >>> from conkit import io
    >>> with open('example.txt', 'r') as f_in, open('example.rr', 'w') as f_out:
    ...     hierarchy = io.read(f_in, 'psicov')
    ...     io.write(f_out, 'casprr', hierarchy)

    """
    # Check for the correct format and values provided
    if not (format in CONTACT_FILE_PARSERS or format in SEQUENCE_FILE_PARSERS):
        raise ValueError("Unrecognised format: '{selected}'".format(selected=format))
    elif format in CONTACT_FILE_PARSERS:
        parser = CONTACT_FILE_PARSERS[format]()
    elif format in SEQUENCE_FILE_PARSERS:
        parser = SEQUENCE_FILE_PARSERS[format]()
    else:
        raise Exception("Should never be here")
    
    with open_f_handle(fname, 'write') as f_out:
        parser.write(f_out, hierarchy)

    return
