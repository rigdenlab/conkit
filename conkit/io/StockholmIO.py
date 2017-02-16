"""
Parser module specific to Stockholm sequence files
"""

__author__ = "Felix Simkovic"
__date__ = "09 Sep 2016"
__version__ = 0.1

from conkit.core import Sequence
from conkit.core import SequenceFile
from conkit.io._ParserIO import _SequenceFileParser

import os
import re

V_RECORD = re.compile(r'^#(\s+STOCKHOLM.*)$')
GF_RECORD = re.compile(r'^#=GF\s+\S+\s+(.*)$')
GR_RECORD = re.compile(r'^#=GR\s+(\S+)\s+(\S+)\s+(.*)$')
GS_RECORD = re.compile(r'^#=GS\s+(\S+)\s+(\S+)\s+(.*)$')
SEQ_RECORD = re.compile(r'^(\S+)\s+([A-Z0-9~-]+)$')
END_RECORD = re.compile(r'^//$')


class StockholmIO(_SequenceFileParser):
    """Parser class for Stockholm sequence files
    """
    def __init__(self):
        super(StockholmIO, self).__init__()

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
        :obj:`SequenceFile <conkit.core.SequenceFile>`

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
        hierarchy : :obj:`SequenceFile <conkit.core.SequenceFile>` or :obj:`Sequence <conkit.core.Sequence>`

        """
        # Double check the type of sequence_file and reconstruct if necessary
        sequence_file = self._reconstruct(hierarchy)

        f_handle.write("# STOCKHOLM 1.0" + os.linesep)
        f_handle.write("#=GF ID {ident}{ls}{ls}".format(ident=sequence_file.top_sequence.id,
                                                        ls=os.linesep))
        chunks = []
        for i, sequence_entry in enumerate(sequence_file):

            if i != 0:
                f_handle.write(
                    '#=GS {ident:33} DE {remark}{ls}'.format(
                        ident=sequence_entry.id,
                        remark=" ".join(sequence_entry.remark),
                        ls=os.linesep,
                    )
                )

            # Cut the sequence into chunks [ <= 200 seq chars per line]
            chunk = []
            sequence_string = sequence_entry.seq
            sequence_string = sequence_string.upper()  # UPPER CASE !!!
            for j in range(0, sequence_entry.seq_len, 200):
                chunk.append(sequence_string[j:j + 200])
            chunks.append(tuple([sequence_entry.id, chunk]))

        # Write the sequence out in chunks
        for j in range(len(chunks[0][1])):
            f_handle.write(os.linesep)
            for i in range(len(chunks)):
                f_handle.write(
                    "{ident:41} {seq_chunk}{ls}".format(
                        ident=chunks[i][0],
                        seq_chunk=chunks[i][1][j],
                        ls=os.linesep,
                    )
                )

        f_handle.write("//" + os.linesep)

        return
