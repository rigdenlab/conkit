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
Parser module specific to CLUSTAL sequence files
"""

__author__ = 'Felix Simkovic'
__date__ = '30 Jul 2018'
__version__ = '0.1'

import collections
import os

from conkit.io._parser import SequenceFileParser
from conkit.core.sequence import Sequence
from conkit.core.sequencefile import SequenceFile


class ClustalParser(SequenceFileParser):
    """Parser class for CLUSTAL sequence files
    """

    def __init__(self):
        super(ClustalParser, self).__init__()

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
        :exc:`TypeError` 
           Incorrect file format

        """

        header = f_handle.readline().rstrip()
        if header[:7].upper() != 'CLUSTAL':
            raise TypeError('Incorrect file format')

        cache = collections.OrderedDict()
        for line in f_handle:
            line = line.strip()
            line_set = set(line)
            if not line or any(char in line_set for char in ['*', ':', '.']):
                continue
            else:
                parts = line.strip().split()
                id_, chunk = parts[:2]
                if id_ in cache:
                    cache[id_] += chunk
                else:
                    cache[id_] = chunk

        hierarchy = SequenceFile(f_id)
        while len(cache) > 0:
            id_, seq = cache.popitem(last=False)  # FIFO
            sequence = Sequence(id_, seq)
            hierarchy.add(sequence)

        return hierarchy

    def write(self, f_handle, hierarchy):
        """Write a sequence file instance to to file

        Parameters
        ----------
        f_handle
           Open file handle [write permissions]
        hierarchy : :obj:`~conkit.core.sequencefile.SequenceFile`, :obj:`~conkit.core.sequence.Sequence`

        """
        hierarchy = self._reconstruct(hierarchy)
        chunker = []
        longest = 0
        for sequence in hierarchy:
            this = [sequence.id]
            if len(sequence.id) > longest:
                longest = len(sequence.id)
            for i in range(0, sequence.seq_len, 60):
                this += [sequence.seq[i:i + 60]]
            chunker += [this]

        content = 'CLUSTAL FORMAT written with ConKit\n'
        content += '\n'
        linetemplate = '%-{}s\t%s\n'.format(longest)
        while len(chunker) > 0:
            entry = chunker.pop(0)
            content += linetemplate % (entry[0], entry[1])
            entry.pop(1)
            if len(entry) > 1:
                chunker.append(entry)

        f_handle.write(content)
