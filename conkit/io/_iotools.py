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
"""Internal utility functions"""

__author__ = "Felix Simkovic"
__date__ = "20 Nov 2016"
__version__ = "0.1"

import io
import sys
import tempfile


def create_tmp_f(content=None, mode='w'):
    """Create a temporary file

    Parameters
    ----------
    mode : str, optional
       Write-mode [default: 'w']
    content : str, optional
       The content to write to the file

       If `None`, leave the file empty

    Returns
    -------
    str
       The path to the filename

    """
    f_in = tempfile.NamedTemporaryFile(mode=mode, delete=False)
    if content:
        f_in.write(content)
    f_in.close()
    return f_in.name


def is_str_like(content):
    """Check if an instance is string-like

    Parameters
    ----------
    content : str
       The content to check

    Returns
    -------
    bool
       True if is str like else False

    Note
    ----
    Function identical to :func:`numpy.is_str_like`, credits to the author

    """
    try:
        content + ''
    except (TypeError, ValueError):
        return False
    return True


def open_f_handle(f_handle, mode):
    """Open a filehandle

    Parameters
    ----------
    f_handle : file_handle, file_name
       A file handle or a file name
    mode : str
       read, write or append

    Returns
    -------
    f_handle
       The opened file handle

    Raises
    ------
    :exc:`TypeError`
       f_handle must be str of filehandle
    :exc:`ValueError`
       Mode needs to be one of: append, read, write

    """
    if mode not in ['append', 'read', 'write']:
        raise ValueError('Mode needs to be one of: append, read, write')

    try:
        if is_str_like(f_handle) and sys.version_info.major >= 3:
            return io.open(f_handle, mode[0], encoding="utf-8")
        elif is_str_like(f_handle):
            return open(f_handle, mode[0])
        elif f_handle.mode == mode[0]:
            return f_handle
        else:
            raise TypeError("f_handle must be str or filehandle")
    except AttributeError:
        raise TypeError("f_handle must be str or filehandle")
