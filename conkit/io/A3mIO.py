"""
Parser module specific to HH-suite A3M sequence files

Credits
-------
Stefan Seemayer and his A3MIO project [https://github.com/sseemayer/BioPython-A3MIO]
"""

__author__ = "Felix Simkovic"
__credits__ = "Stefan Seemayer"
__date__ = "11 Sep 2016"
__version__ = 0.1

from conkit.core import Sequence
from conkit.core import SequenceFile
from conkit.io._ParserIO import _SequenceFileParser

import numpy
import os
import re


class A3mIO(_SequenceFileParser):
    """Parser class for A3M sequence files

    """
    def __init__(self):
        super(A3mIO, self).__init__()

    def read(self, f_handle, f_id='a3m', remove_insert=True):
        """Read a sequence file

        Parameters
        ----------
        f_handle
           Open file handle [read permissions]
        f_id : str, optional
           Unique sequence file identifier
        remove_insert : bool, optional
           Remove insert states [default: True]

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
            elif line.startswith('#'):
                sequence_file.remark = line[1:]
            elif line.startswith('>'):
                break

        # Read the sequence record(s) and store them
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
            seq_string = "".join(chunks)   # Sequence from chunks

            # Remove insert states
            if remove_insert:
                seq_string = self._remove_insert(seq_string)

            # Create the sequence record instance
            sequence_entry = Sequence(id, seq_string)

            # Store the sequence in the file
            try:
                sequence_file.add(sequence_entry)
            except ValueError:
                while True:
                    new_id = sequence_entry.id + "_{0}".format(numpy.random.randint(0, 100000))
                    if new_id in sequence_file:
                        continue
                    else:
                        break
                sequence_entry.id = new_id
                sequence_file.add(sequence_entry)

            if not line:
                break

        # Match the insert states of the sequence
        if not remove_insert:
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
        inserts = [
            INSERT_STATE.split(sequence_entry.seq)
            for sequence_entry in hierarchy
        ]

        # Determine maximum insert length at each position
        insert_max_lengths = [
            max(
                len(inserts[i][j])
                for i in range(len(inserts))
            )
            for j in range(len(inserts[0]))
        ]

        # Add gaps where gaps are needed
        def pad(chunk, length, pad_char='-'):
            return chunk + pad_char * (length - len(chunk))

        # Manipulate each sequence to match insert states
        for sequence_entry, seq in zip(hierarchy, inserts):
            sequence_entry.seq = "".join(pad(insert, insert_len)
                                         for insert, insert_len
                                         in zip(seq, insert_max_lengths))
        return

    def _remove_insert(self, seq):
        """Remove insert states"""
        return "".join([char for char in seq if not char.islower()])


    def write(self, f_handle, hierarchy):
        """Write a sequence file instance to to file

        Parameters
        ----------
        f_handle
           Open file handle [write permissions]
        hierarchy : :obj:`SequenceFile <conkit.core.SequenceFile>`, :obj:`Sequence <conkit.core.Sequence>`

        """
        # Double check the type of hierarchy and reconstruct if necessary
        sequence_file = self._reconstruct(hierarchy)

        # Write remarks
        for remark in sequence_file.remark:
            f_handle.write('#{remark}'.format(remark=remark) + os.linesep)

        for sequence_entry in sequence_file:
            header = '>{id}'.format(id=sequence_entry.id)
            if len(sequence_entry.remark) > 0:
                header = '|'.join([header] + sequence_entry.remark)
            f_handle.write(header + os.linesep)
            f_handle.write(sequence_entry.seq + os.linesep)

        return
