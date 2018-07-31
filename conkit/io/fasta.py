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
Parser module specific to FASTA sequence files
"""

__author__ = "Felix Simkovic"
__date__ = "09 Sep 2016"
__version__ = "0.1"

import os

from conkit.io._parser import SequenceFileParser
from conkit.core.sequence import Sequence
from conkit.core.sequencefile import SequenceFile


class FastaParser(SequenceFileParser):
    """Parser class for FASTA sequence files
    """

    def __init__(self):
        super(FastaParser, self).__init__()

    def read(self, f_handle, f_id='fasta'):
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

        Raises
        ------
        :exc:`ValueError`
           FASTA record needs to start with >

        """
        hierarchy = SequenceFile(f_id)

        while True:
            line = f_handle.readline().rstrip()

            if not line:
                continue
            elif line.startswith('#'):
                hierarchy.remark = line[1:]
            elif line.startswith('>'):
                break

        while True:
            if not line.startswith('>'):
                raise ValueError("Fasta record needs to start with '>'")

            id = line[1:]  # Header without '>'

            chunks = []
            line = f_handle.readline().rstrip()
            while True:
                if not line:
                    break
                elif line.startswith('>'):
                    break
                chunks.append(line)
                line = f_handle.readline().rstrip()
            _seq_string = "".join(chunks)  # Sequence from chunks

            sequence_entry = Sequence(id, _seq_string)

            hierarchy.add(sequence_entry)

            if not line:
                break

        return hierarchy

    def write(self, f_handle, hierarchy):
        """Write a sequence file instance to to file

        Parameters
        ----------
        f_handle
           Open file handle [write permissions]
        hierarchy : :obj:`~conkit.core.sequencefile.SequenceFile`, :obj:`~conkit.core.sequence.Sequence`

        """
        # Double check the type of hierarchy and reconstruct if necessary
        hierarchy = self._reconstruct(hierarchy)

        content = ""

        # Write remarks
        for remark in hierarchy.remark:
            content += '#{remark}'.format(remark=remark) + os.linesep

        for sequence_entry in hierarchy:
            header = '>{id}'.format(id=sequence_entry.id)
            if len(sequence_entry.remark) > 0:
                header = '|'.join([header] + sequence_entry.remark)
            content += header + os.linesep

            # Cut the sequence into chunks [FASTA <= 60 chars per line]
            sequence_string = sequence_entry.seq.upper()  # UPPER CASE !!!
            for i in range(0, sequence_entry.seq_len, 60):
                content += sequence_string[i:i + 60] + os.linesep

        f_handle.write(content)
