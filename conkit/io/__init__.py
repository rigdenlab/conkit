# BSD 3-Clause License
#
# Copyright (c) 2016-18, University of Liverpool
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""I/O interface for file reading, writing and conversions"""

__author__ = 'Felix Simkovic'
__date__ = '13 Sep 2016'
__version__ = "0.1"

import importlib
import os

from conkit.io._cache import PARSER_CACHE
from conkit.io._iotools import open_f_handle

# Accessed by some modules - might be deprecated in the future
CONTACT_FILE_PARSERS = PARSER_CACHE.contact_file_parsers
SEQUENCE_FILE_PARSERS = PARSER_CACHE.sequence_file_parsers


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

    Note
    ----
    A3M format comes by default WITHOUT insert states, these are removed. To obtain
    an alignment WITH insert states, use format ``a3m-inserts``.

    2) Convert a PconsC3 contact prediction file to the standard Casp RR format:

    >>> from conkit import io
    >>> with open('example.out', 'r') as f_in, open('example.rr', 'w') as f_out:
    ...     io.convert(f_in, 'pconsc3', f_out, 'casprr'))

    """
    if format_in in CONTACT_FILE_PARSERS and format_out in SEQUENCE_FILE_PARSERS:
        raise ValueError("Cannot convert contact file to sequence file")
    elif format_in in SEQUENCE_FILE_PARSERS and format_out in CONTACT_FILE_PARSERS:
        raise ValueError("Cannot convert sequence file to contact file")
    else:
        hierarchy = read(fname_in, format_in)
        write(fname_out, format_out, hierarchy)


def read(fname, format, f_id='conkit'):
    """Parse a file handle to read into structure

    Parameters
    ----------
    fname : filehandle, filename
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
    if format in PARSER_CACHE:
        parser_in = PARSER_CACHE.import_class(format)()
    else:
        raise ValueError("Unrecognised format: '{}'".format(format))

    kwargs = {"f_id": f_id}
    if format == "a3m-inserts":
        kwargs["remove_inserts"] = False

    with open_f_handle(fname, "read") as f_in:
        hierarchy = parser_in.read(f_in, **kwargs)

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
    if format in PARSER_CACHE:
        parser_out = PARSER_CACHE.import_class(format)()
    else:
        raise ValueError("Unrecognised format: '{}'".format(format))

    kwargs = {}
    if format in ["flib", "pconsc", "pconsc2"]:
        kwargs["write_header_footer"] = False

    with open_f_handle(fname, 'write') as f_out:
        parser_out.write(f_out, hierarchy, **kwargs)
