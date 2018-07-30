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
Parser module specific to HH-suite A3M sequence files

Credits
-------
Stefan Seemayer and his A3MIO project [https://github.com/sseemayer/BioPython-A3MIO]

"""

__author__ = "Felix Simkovic"
__credits__ = "Stefan Seemayer"
__date__ = "11 Sep 2016"
__version__ = "0.1"

import numpy as np
import os
import re

from conkit.io._parser import SequenceFileParser
from conkit.core.sequence import Sequence
from conkit.core.sequencefile import SequenceFile


class A3mParser(SequenceFileParser):
    """Parser class for A3M sequence files

    """

    def __init__(self):
        super(A3mParser, self).__init__()

    def read(self, f_handle, f_id='a3m', remove_inserts=True):
        """Read a sequence file

        Parameters
        ----------
        f_handle
           Open file handle [read permissions]
        f_id : str, optional
           Unique sequence file identifier
        remove_inserts : bool, optional
           Remove insert states [default: True]

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
            elif line.startswith('#'):
                sequence_file.remark = line[1:]
            elif line.startswith('>'):
                break

        while True:
            if not line.startswith('>'):
                raise ValueError("Fasta record needs to start with '>'")
            id = line[1:]
            chunks = []
            line = f_handle.readline().rstrip()
            while True:
                if not line:
                    break
                elif line.startswith('>'):
                    break
                chunks.append(line)
                line = f_handle.readline().rstrip()
            seq_string = "".join(chunks)
            if remove_inserts:
                seq_string = self._remove_inserts(seq_string)
            sequence_entry = Sequence(id, seq_string)
            try:
                sequence_file.add(sequence_entry)
            except ValueError:
                while True:
                    new_id = sequence_entry.id + "_{0}".format(np.random.randint(0, 100000))
                    if new_id in sequence_file:
                        continue
                    else:
                        break
                sequence_entry.id = new_id
                sequence_file.add(sequence_entry)
            if not line:
                break

        if not remove_inserts:
            self._adjust_insert(sequence_file)

        return sequence_file

    def _adjust_insert(self, hierarchy):
        """Adjust insert states

        Credits
        -------
        This function was adapted from Stefan Seemayer's BioPython-A3MIO
        repository - https://github.com/sseemayer/BioPython-A3MIO

        """

        # Determine the insert states by splitting the sequences into chunks based
        # on the case of the letter.
        INSERT_STATE = re.compile(r'([A-Z0-9~-])')
        inserts = [INSERT_STATE.split(sequence_entry.seq) for sequence_entry in hierarchy]

        # Determine maximum insert length at each position
        insert_max_lengths = [max(len(inserts[i][j]) for i in range(len(inserts))) for j in range(len(inserts[0]))]

        # Add gaps where gaps are needed
        def pad(chunk, length, pad_char='-'):
            return chunk + pad_char * (length - len(chunk))

        # Manipulate each sequence to match insert states
        for sequence_entry, seq in zip(hierarchy, inserts):
            sequence_entry.seq = "".join(pad(insert, insert_len) for insert, insert_len in zip(seq, insert_max_lengths))
        return

    def _remove_inserts(self, seq):
        """Remove insert states"""
        return "".join([char for char in seq if not char.islower()])

    def write(self, f_handle, hierarchy):
        """Write a sequence file instance to to file

        Parameters
        ----------
        f_handle
           Open file handle [write permissions]
        hierarchy : :obj:`~conkit.core.sequencefile.SequenceFile`, :obj:`~conkit.core.sequence.Sequence`

        """
        # Double check the type of hierarchy and reconstruct if necessary
        sequence_file = self._reconstruct(hierarchy)

        content = ""

        # Write remarks
        for remark in sequence_file.remark:
            content += '#{remark}'.format(remark=remark) + os.linesep

        for sequence_entry in sequence_file:
            header = '>{id}'.format(id=sequence_entry.id)
            if len(sequence_entry.remark) > 0:
                header = '|'.join([header] + sequence_entry.remark)
            content += header + os.linesep
            content += sequence_entry.seq + os.linesep

        f_handle.write(content)
