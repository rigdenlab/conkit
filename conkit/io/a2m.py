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
"""
Parser module specific to HH-suite A2M sequence files

"""

__author__ = "Felix Simkovic"
__date__ = "30 Jul 2018"
__version__ = "1.0"

import os

from conkit.io._parser import SequenceFileParser
from conkit.core.sequence import Sequence
from conkit.core.sequencefile import SequenceFile


class A2mParser(SequenceFileParser):
    """Parser class for A2M sequence files

    """

    def __init__(self):
        super(A2mParser, self).__init__()

    def read(self, f_handle, f_id='a2m'):
        """Read a sequence file

        Parameters
        ----------
        f_handle
           Open file handle [read permissions]
        f_id : str, optional
           Unique sequence file identifier

        Returns
        -------
        :obj:`~conkit.core.sequencefile.SequenceFile`

        """

        hierarchy = SequenceFile(f_id)

        for i, line in enumerate(f_handle):

            line = line.rstrip()
            if not line:
                continue

            for a, c in enumerate(line):
                if c.isalpha() or c == '-':
                    continue
                else:
                    indicator = ['-'] * len(line)
                    indicator[a] = '^'
                    msg = "Unknown character in line {0}:{1}{1}{2}{1}{3}"
                    msg = msg.format(i + 1, os.linesep, line, ''.join(indicator))
                    raise ValueError(msg)
            sequence = line
            sequence_entry = Sequence('seq_{i}'.format(i=i), sequence)

            hierarchy.add(sequence_entry)

        return hierarchy

    def write(self, f_handle, hierarchy):
        """Write a sequence file instance to to file

        Parameters
        ----------
        f_handle
           Open file handle [write permissions]
        hierarchy : :obj:`~conkit.core.sequencefile.SequenceFile` or :obj:`~conkit.core.sequence.Sequence`

        """
        # Double check the type of hierarchy and reconstruct if necessary
        sequence_file = self._reconstruct(hierarchy)

        content = ""

        for sequence_entry in sequence_file:
            content += sequence_entry.seq + os.linesep

        f_handle.write(content)
