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
Parser module specific to Stockholm sequence files
"""

__author__ = "Felix Simkovic"
__date__ = "09 Sep 2016"
__version__ = "0.1"

import os
import re

from conkit.io._parser import SequenceFileParser
from conkit.core.sequence import Sequence
from conkit.core.sequencefile import SequenceFile

V_RECORD = re.compile(r'^#(\s+STOCKHOLM.*)$')
GF_RECORD = re.compile(r'^#=GF\s+\S+\s+(.*)$')
GR_RECORD = re.compile(r'^#=GR\s+(\S+)\s+(\S+)\s+(.*)$')
GS_RECORD = re.compile(r'^#=GS\s+(\S+)\s+(\S+)\s+(.*)$')
SEQ_RECORD = re.compile(r'^(\S+)\s+([A-Z0-9~-]+)$')
END_RECORD = re.compile(r'^//$')


class StockholmParser(SequenceFileParser):
    """Parser class for Stockholm sequence files
    """

    def __init__(self):
        super(StockholmParser, self).__init__()

    def read(self, f_handle, f_id='stockholm'):
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

        # Create a new sequence file instance
        sequence_file = SequenceFile(f_id)

        # Read any possible comments and store in file remarks
        while True:
            line = f_handle.readline().rstrip()

            if not line:
                continue
            elif V_RECORD.match(line):
                _ = V_RECORD.match(line).group(1)
                # sequence_file.add_remark(version)
            elif GF_RECORD.match(line) or GS_RECORD.match(line):
                break

        # Read the sequence record(s)
        while True:
            if GF_RECORD.match(line):
                ident = GF_RECORD.match(line).group(1)[:-3]
                sequence_entry = Sequence(ident, "")
                sequence_file.add(sequence_entry)

            elif GS_RECORD.match(line):
                ident, _, desc = GS_RECORD.match(line).groups()
                sequence_entry = Sequence(ident, "")
                sequence_entry.remark = desc
                sequence_file.add(sequence_entry)

            elif GR_RECORD.match(line):
                pass

            elif len(line.split()) == 2 and line.split()[0] in sequence_file:
                ident, seq = line.replace('.', '-').split()
                sequence_file[ident].seq = sequence_file[ident].seq + seq

            line = f_handle.readline().rstrip()

            # // char in alignment reached
            if END_RECORD.match(line):
                break

        return sequence_file

    def write(self, f_handle, hierarchy):
        """Write a sequence file instance to to file

        Parameters
        ----------
        f_handle
           Open file handle [write permissions]
        hierarchy : :obj:`~conkit.core.sequencefile.SequenceFile` or :obj:`~conkit.core.sequence.Sequence`

        """
        sequence_file = self._reconstruct(hierarchy)

        content = "# STOCKHOLM 1.0" + os.linesep
        content += "#=GF ID {}".format(sequence_file.top_sequence.id) + os.linesep * 2

        chunks = []
        for i, sequence_entry in enumerate(sequence_file):

            if i != 0:
                content += '#=GS {:33} DE {}'.format(sequence_entry.id, " ".join(sequence_entry.remark)) + os.linesep

            # Cut the sequence into chunks [ <= 200 seq chars per line]
            chunk = []
            sequence_string = sequence_entry.seq
            sequence_string = sequence_string.upper()  # UPPER CASE !!!
            for j in range(0, sequence_entry.seq_len, 200):
                chunk.append(sequence_string[j:j + 200])
            chunks.append(tuple([sequence_entry.id, chunk]))

        # Write the sequence out in chunks
        for j in range(len(chunks[0][1])):
            content += os.linesep
            for i in range(len(chunks)):
                content += "{:41} {}".format(chunks[i][0], chunks[i][1][j]) + os.linesep

        content += "//" + os.linesep

        f_handle.write(content)
