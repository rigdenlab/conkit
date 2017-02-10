"""
Parser module specific to FASTA sequence files
"""

__author__ = "Felix Simkovic"
__date__ = "09 Sep 2016"
__version__ = 0.1

from conkit.core import Sequence
from conkit.core import SequenceFile
from conkit.io._ParserIO import _SequenceFileParser

import os


class FastaIO(_SequenceFileParser):
    """Parser class for FASTA sequence files
    """
    def __init__(self):
        super(FastaIO, self).__init__()

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
        :obj:`SequenceFile <conkit.core.SequenceFile>`

        """

        # Create a new sequence file instance
        hierarchy = SequenceFile(f_id)

        # Read any possible comments and store in file remarks
        while True:
            line = f_handle.readline().rstrip()

            if not line:
                continue
            elif line.startswith('#'):
                hierarchy.remark = line[1:]
            elif line.startswith('>'):
                break

        # Read the sequence record(s) and store them
        while True:
            if not line.startswith('>'):
                raise ValueError("Fasta record needs to start with '>'")

            id = line[1:]   # Header without '>'

            chunks = []
            line = f_handle.readline().rstrip()
            while True:
                if not line:
                    break
                elif line.startswith('>'):
                    break
                chunks.append(line)
                line = f_handle.readline().rstrip()
            _seq_string = "".join(chunks)   # Sequence from chunks

            # Create the sequence record instance
            sequence_entry = Sequence(id, _seq_string)

            # Store the sequence in the file
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
        hierarchy : :obj:`SequenceFile <conkit.core.SequenceFile>`, :obj:`Sequence <conkit.core.Sequence>`

        """
        # Double check the type of hierarchy and reconstruct if necessary
        hierarchy = self._reconstruct(hierarchy)

        # Write remarks
        for remark in hierarchy.remark:
            f_handle.write('#{remark}'.format(remark=remark) + os.linesep)

        for sequence_entry in hierarchy:
            header = '>{id}'.format(id=sequence_entry.id)
            if len(sequence_entry.remark) > 0:
                header = '|'.join([header] + sequence_entry.remark)
            f_handle.write(header + os.linesep)

            # Cut the sequence into chunks [FASTA <= 60 chars per line]
            sequence_string = sequence_entry.seq.upper()       # UPPER CASE !!!
            for i in range(0, sequence_entry.seq_len, 60):
                f_handle.write(sequence_string[i:i+60] + os.linesep)

        return
