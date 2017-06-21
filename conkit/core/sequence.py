# coding=utf-8
#
# BSD 3-Clause License
#
# Copyright (c) 2016-17, University of Liverpool
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
"""Sequence container used throughout ConKit"""

from __future__ import division
from __future__ import print_function

__author__ = "Felix Simkovic"
__date__ = "03 Aug 2016"
__version__ = "1.0"

from Bio import pairwise2
from conkit.core._entity import _Entity

# One to three amino acid letter code conversions
ONE_TO_THREE = {'A': 'ALA', 'C': 'CYS', 'B': 'ASX', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY', 'F': 'PHE', 'I': 'ILE',
                'H': 'HIS', 'K': 'LYS', 'J': 'XLE', 'M': 'MET', 'L': 'LEU', 'O': 'PYL', 'N': 'ASN', 'Q': 'GLN',
                'P': 'PRO', 'S': 'SER', 'R': 'ARG', 'U': 'SEC', 'T': 'THR', 'W': 'TRP', 'V': 'VAL', 'Y': 'TYR',
                'X': 'XAA', 'Z': 'GLX'}

# Three to one amino acid letter code conversions
THREE_TO_ONE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CME': 'C', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E',
                'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'MSE': 'M', 'PHE': 'F',
                'PRO': 'P', 'PYL': 'O', 'SER': 'S', 'SEC': 'U', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V',
                'ASX': 'B', 'GLX': 'Z', 'XAA': 'X', 'UNK': 'X', 'XLE': 'J'}


class Sequence(_Entity):
    """A sequence template to store all associated information

    Attributes
    ----------
    id : str
       A unique identifier
    remark : list
       The :obj:`Sequence <conkit.core.Sequence>`-specific remarks
    seq : str
       The protein sequence as :obj:`str`
    seq_len : int
       The protein sequence length

    Examples
    --------
    >>> from conkit.core import Sequence
    >>> sequence_entry = Sequence("example", "ABCDEF")
    >>> print(sequence_entry)
    Sequence(id="example" seq="ABCDEF" seqlen=6)

    """
    __slots__ = ['_remark', '_seq']

    def __init__(self, id, seq):
        """Initialise a generic sequence

        Parameters
        ----------
        id : str
           A unique sequence identifier
        seq : str
           The protein sequence

        """
        self._remark = []
        self._seq = None

        # Assign values post creation to use setter/getter methods
        # Possibly very bad practice but no better alternative for now
        self.seq = seq

        super(Sequence, self).__init__(id)

    def __add__(self, other):
        """Concatenate two sequence instances to a new"""
        id = self.id + '_' + other.id
        seq = self.seq + other.seq
        return Sequence(id, seq)

    def __repr__(self):
        if self.seq_len > 12:
            seq_string = self.seq[:5] + '...' + self.seq[-5:]
        else:
            seq_string = self.seq
        return "{0}(id=\"{1}\" seq=\"{2}\" seq_len={3})".format(
            self.__class__.__name__, self.id, seq_string, len(self.seq)
        )

    @property
    def remark(self):
        """The :obj:`Sequence <conkit.core.Sequence>`-specific remarks"""
        return self._remark

    @remark.setter
    def remark(self, remark):
        """Set the :obj:`Sequence <conkit.core.Sequence>` remark

        Parameters
        ----------
        remark : str, list
           The remark will be added to the list of remarks

        """
        if isinstance(remark, list):
            self._remark += remark
        elif isinstance(remark, tuple):
            self._remark += list(remark)
        else:
            self._remark += [remark]

    @property
    def seq(self):
        """The protein sequence as :obj:`str`"""
        return self._seq

    @seq.setter
    def seq(self, seq):
        """Set the sequence

        Parameters
        ----------
        seq : str

        Raises
        ------
        ValueError
           One or more amino acids in the sequence are not recognised

        """
        if any(c not in list(ONE_TO_THREE.keys()) for c in seq.upper() if c != '-'):
            raise ValueError('Unrecognized amino acids in sequence')
        self._seq = seq

    @property
    def seq_len(self):
        """The protein sequence length"""
        return len(self.seq)

    def align_global(self, other, id_chars=2, nonid_chars=1, gap_open_pen=-0.5, gap_ext_pen=-0.1, inplace=False):
        """Generate a global alignment between two :obj:`Sequence <conkit.core.Sequence>` instances

        Parameters
        ----------
        other : :obj:`Sequence <conkit.core.Sequence>`
        id_chars : int, optional
        nonid_chars : int, optional
        gap_open_pen : float, optional
        gap_ext_pen : float, optional
        inplace : bool, optional
           Replace the saved order of residues [default: False]

        Returns
        -------
        obj
           The reference to the :obj:`Sequence`, regardless of inplace
        obj
           The reference to the :obj:`Sequence`, regardless of inplace

        """
        sequence1 = self._inplace(inplace)
        sequence2 = other._inplace(inplace)

        alignment = pairwise2.align.globalms(
            sequence1.seq, sequence2.seq, id_chars, nonid_chars, gap_open_pen, gap_ext_pen
        )

        sequence1.seq = alignment[-1][0]
        sequence2.seq = alignment[-1][1]

        return sequence1, sequence2

    def align_local(self, other, id_chars=2, nonid_chars=1, gap_open_pen=-0.5, gap_ext_pen=-0.1, inplace=False):
        """Generate a local alignment between two :obj:`Sequence <conkit.core.Sequence>` instances

        Parameters
        ----------
        other : :obj:`Sequence <conkit.core.Sequence>`
        id_chars : int, optional
        nonid_chars : int, optional
        gap_open_pen : float, optional
        gap_ext_pen : float, optional
        inplace : bool, optional
           Replace the saved order of residues [default: False]

        Returns
        -------
        obj
           The reference to the :obj:`Sequence <conkit.core.Sequence>`, regardless of inplace
        obj
           The reference to the :obj:`Sequence <conkit.core.Sequence>`, regardless of inplace

        """
        sequence1 = self._inplace(inplace)
        sequence2 = other._inplace(inplace)

        alignment = pairwise2.align.localms(
            sequence1.seq, sequence2.seq, id_chars, nonid_chars, gap_open_pen, gap_ext_pen
        )

        sequence1.seq = alignment[-1][0]
        sequence2.seq = alignment[-1][1]

        return sequence1, sequence2
