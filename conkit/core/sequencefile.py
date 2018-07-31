# coding=utf-8
#
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
"""SequenceFile container used throughout ConKit"""

from __future__ import division
from __future__ import print_function

__author__ = "Felix Simkovic"
__date__ = "03 Aug 2016"
__version__ = "1.0"

import numpy as np
import os
import sys

if sys.version_info.major < 3:
    from itertools import izip as zip

from conkit.core.entity import Entity
from conkit.core.mappings import AminoAcidMapping, SequenceAlignmentState
from conkit.misc import deprecate


class SequenceFile(Entity):
    """A sequence file object representing a single sequence file

    The :obj:`~conkit.core.sequencefile.SequenceFile` class represents a data structure to hold
    :obj:`~conkit.core.sequence.Sequence` instances in a single sequence file. It contains
    functions to store and analyze sequences.

    Attributes
    ----------
    id : str
       A unique identifier
    is_alignment : bool
       A boolean status for the alignment
    meff : int
       The number of effective sequences in the :obj:`~conkit.core.sequencefile.SequenceFile`
    nseq : int
       The number of sequences in the :obj:`~conkit.core.sequencefile.SequenceFile`
    remark : list
       The :obj:`~conkit.core.sequencefile.SequenceFile`-specific remarks
    status : int
       An indication of the sequence file, i.e alignment, no alignment, or unknown
    top_sequence : :obj:`~conkit.core.sequence.Sequence`, None
       The first :obj:`~conkit.core.sequence.Sequence` entry in the file


    Examples
    --------
    >>> from conkit.core import Sequence, SequenceFile
    >>> sequence_file = SequenceFile("example")
    >>> sequence_file.add(Sequence("foo", "ABCDEF"))
    >>> sequence_file.add(Sequence("bar", "ZYXWVU"))
    >>> print(sequence_file)
    SequenceFile(id="example" nseq=2)

    """
    __slots__ = ['_remark', '_status']

    def __init__(self, id):
        """Initialise a new :obj:`~conkit.core.sequencefile.SequenceFile`

        Parameters
        ----------
        id : str
           A unique identifier for the sequence file

        """
        self._remark = []
        self._status = SequenceAlignmentState.unknown
        super(SequenceFile, self).__init__(id)

    def __repr__(self):
        return '{}(id="{}" nseq={})'.format(self.__class__.__name__, self.id, self.nseq)

    @property
    def ascii_matrix(self):
        """The alignment encoded in a 2-D ASCII matrix"""
        return [list(seq.seq_ascii) for seq in self]

    @property
    def encoded_matrix(self):
        """The alignment encoded for contact prediction"""
        return [list(seq.seq_encoded) for seq in self]

    @property
    def is_alignment(self):
        """A boolean status for the alignment

        Returns
        -------
        bool
           A boolean status for the alignment

        """
        seq_length = self.top_sequence.seq_len
        self._status = SequenceAlignmentState.aligned
        for sequence in self:
            if sequence.seq_len != seq_length:
                self._status = SequenceAlignmentState.unaligned
                break
        return self._status == SequenceAlignmentState.aligned

    @property
    def diversity(self):
        """The diversity of an alignment defined by :math:`\\sqrt{N}/L`. 
        
        ``N`` equals the number of sequences in
        the alignment and ``L`` the sequence length
        
        """
        if self.empty:
            return 0.0
        elif self.is_alignment:
            return (np.sqrt(len(self)) / float(self.top.seq_len)).round(decimals=3).item()
        else:
            raise ValueError('This is not an alignment')

    @property
    def empty(self):
        """Status of emptiness of sequencefile"""
        return len(self) < 1

    @property
    @deprecate('0.11', msg='Use meff instead.')
    def neff(self):
        """The number of effective sequences"""
        return self.meff

    @property
    def meff(self):
        """The number of effective sequences"""
        return int(sum(self.get_weights()))

    @property
    def nseq(self):
        """The number of sequences"""
        return len(self)

    @property
    def remark(self):
        """The :obj:`~conkit.core.sequencefile.SequenceFile`-specific remarks"""
        return self._remark

    @remark.setter
    def remark(self, remark):
        """Set the :obj:`~conkit.core.sequencefile.SequenceFile` remark

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
    def status(self):
        """An indication of the residue status, i.e true positive, false positive, or unknown"""
        return self._status.value

    @status.setter
    def status(self, status):
        """Set the status

        Parameters
        ----------
        status : int
           [0] for `unknown`, 
           [-1] for `no alignment`, or 
           [1] for `alignment`

        Raises
        ------
        :exc:`ValueError`
           Cannot determine if your sequence file is an alignment or not

        """
        self._status = SequenceAlignmentState(status)

    @property
    def top_sequence(self):
        """The first :obj:`~conkit.core.sequence.Sequence` entry in :obj:`~conkit.core.sequencefile.SequenceFile`

        Returns
        -------
        :obj:`~conkit.core.sequence.Sequence` 
           The first :obj:`~conkit.core.sequence.Sequence` entry in :obj:`~conkit.core.sequencefile.SequenceFile`

        """
        return self.top

    @deprecate('0.11', msg='Use calculate_meff_with_identity instead.')
    def calculate_meff(self, identity=0.8):
        """Calculate the number of effective sequences"""
        return self.calculate_meff_with_identity(identity)

    @deprecate('0.11', msg='Use calculate_meff_with_identity instead.')
    def calculate_neff_with_identity(self, identity):
        """Calculate the number of effective sequences with specified sequence identity"""
        return self.calculate_meff_with_identity(identity)

    @deprecate('0.11', msg='Use get_weights instead.')
    def calculate_weights(self, identity=0.8):
        """Calculate the sequence weights"""
        return self.get_weights(identity=identity)

    @deprecate('0.11', msg='Use get_frequency instead.')
    def calculate_freq(self):
        """Calculate the gap frequency in each alignment column"""
        return self.get_frequency("X")

    def get_meff_with_id(self, identity):
        """Calculate the number of effective sequences with specified sequence identity
        
        See Also
        --------
        meff, get_weights

        """
        return int(sum(self.get_weights(identity=identity)))

    def get_weights(self, identity=0.8):
        """Calculate the sequence weights

        This function calculates the sequence weights in the
        the Multiple Sequence Alignment.

        The mathematical function used to calculate `Meff` is

        .. math::

           M_{eff}=\\sum_{i}\\frac{1}{\\sum_{j}S_{i,j}}

        Parameters
        ----------
        identity : float, optional
           The sequence identity to use for similarity decision [default: 0.8]

        Returns
        -------
        list
           A list of the sequence weights in the alignment 

        Raises
        ------
        :exc:`ValueError`
           :obj:`~conkit.core.sequencefile.SequenceFile` is not an alignment
        :exc:`ValueError`
           Sequence Identity needs to be between 0 and 1

        """
        if identity < 0 or identity > 1:
            raise ValueError("Sequence Identity needs to be between 0 and 1")

        if self.is_alignment:
            from conkit.core.ext.c_sequencefile import c_get_weights
            X = np.array(self.ascii_matrix, dtype=np.int64)
            hamming = np.zeros(X.shape[0], dtype=np.float64)
            c_get_weights(X, identity, hamming)
            return hamming.tolist()
        else:
            raise ValueError('This is not an alignment')

    def get_frequency(self, symbol):
        """Calculate the frequency of an amino acid (symbol) in each Multiple Sequence Alignment column

        Returns
        -------
        list
           A list containing the per alignment-column amino acid frequency count

        Raises
        ------
        :exc:`RuntimeError`
           :obj:`~conkit.core.sequencefile.SequenceFile` is not an alignment

        """
        if self.is_alignment:
            from conkit.core.ext.c_sequencefile import c_get_frequency
            X = np.array(self.encoded_matrix, dtype=np.int64)
            symbol = getattr(AminoAcidMapping, symbol, AminoAcidMapping["X"]).value
            frequencies = np.zeros(X.shape[1], dtype=np.int64)
            c_get_frequency(X, symbol, frequencies)
            return frequencies.tolist()
        else:
            raise ValueError('This is not an alignment')

    def filter(self, min_id=0.3, max_id=0.9, inplace=False):
        """Filter sequences from an alignment according to the minimum and maximum identity
        between the sequences

        Parameters
        ----------
        min_id : float, optional
           Minimum sequence identity
        max_id : float, optional
           Maximum sequence identity
        inplace : bool, optional
           Replace the saved order of sequences [default: False]

        Returns
        -------
        :obj:`~conkit.core.sequencefile.SequenceFile` 
           The reference to the :obj:`~conkit.core.sequencefile.SequenceFile`, regardless of inplace

        Raises
        ------
        :exc:`ValueError`
           :obj:`~conkit.core.sequencefile.SequenceFile` is not an alignment
        :exc:`ValueError`
           Minimum sequence identity needs to be between 0 and 1
        :exc:`ValueError`
           Maximum sequence identity needs to be between 0 and 1

        """
        if 0 > min_id > 1:
            raise ValueError("Minimum sequence identity needs to be between 0 and 1")
        elif 0 > max_id > 1:
            raise ValueError("Maximum sequence identity needs to be between 0 and 1")

        if self.is_alignment:
            from conkit.core.ext.c_sequencefile import c_filter
            X = np.array(self.ascii_matrix, dtype=np.int64)
            throwables = np.full(X.shape[0], False, dtype=np.bool)
            c_filter(X, min_id, max_id, throwables)
            filtered = self._inplace(inplace)
            for i, sequence in enumerate(self):
                if throwables[i]:
                    filtered.remove(sequence.id)
            return filtered
        else:
            raise ValueError('This is not an alignment')

    def filter_gapped(self, min_prop=0.0, max_prop=0.9, inplace=True):
        """Filter all sequences a gap proportion greater than the limit
        
        Parameters
        ----------
        min_prop : float, optional
           Minimum allowed gap proportion [default: 0.0]
        max_prop : float, optional
           Maximum allowed gap proportion [default: 0.9]
        inplace : bool, optional
           Replace the saved order of sequences [default: False]

        Returns
        -------
        :obj:`~conkit.core.sequencefile.SequenceFile`
           The reference to the :obj:`~conkit.core.sequencefile.SequenceFile`, regardless of inplace

        Raises
        ------
        :exc:`ValueError`
           :obj:`~conkit.core.sequencefile.SequenceFile` is not an alignment
        :exc:`ValueError`
           Minimum gap proportion needs to be between 0 and 1
        :exc:`ValueError`
           Maximum gap proportion needs to be between 0 and 1

        """
        if 0. > min_prop > 1.:
            raise ValueError("Minimum gap proportion needs to be between 0 and 1")
        elif 0. > max_prop > 1.:
            raise ValueError("Maximum gap proportion needs to be between 0 and 1")

        symbol = "X"

        if self.is_alignment:
            from conkit.core.ext.c_sequencefile import c_filter_symbol
            X = np.array(self.encoded_matrix, dtype=np.int64)
            symbol = getattr(AminoAcidMapping, symbol, AminoAcidMapping["X"]).value
            throwables = np.full(X.shape[0], False, dtype=np.bool)
            c_filter_symbol(X, min_prop, max_prop, symbol, throwables)
            filtered = self._inplace(inplace)
            for i, sequence in enumerate(self):
                if throwables[i]:
                    filtered.remove(sequence.id)
            return filtered
        else:
            raise ValueError('This is not an alignment')

    def sort(self, kword, reverse=False, inplace=False):
        """Sort the :obj:`~conkit.core.sequencefile.SequenceFile`

        Parameters
        ----------
        kword : str
           The dictionary key to sort sequences by
        reverse : bool, optional
           Sort the sequences in reverse order [default: False]
        inplace : bool, optional
           Replace the saved order of sequences [default: False]

        Returns
        -------
        :obj:`~conkit.core.sequencefile.SequenceFile` 
           The reference to the :obj:`~conkit.core.sequencefile.SequenceFile`, regardless of inplace

        Raises
        ------
        :exc:`ValueError`
           ``kword`` not in :obj:`~conkit.core.sequencefile.SequenceFile`

        """
        sequence_file = self._inplace(inplace)
        sequence_file._sort(kword, reverse)
        return sequence_file

    def to_string(self):
        """Return the :obj:`~conkit.core.sequencefile.SequenceFile` as :obj:`str`"""
        content = [s.seq for s in self]
        return os.linesep.join(content)

    def trim(self, start, end, inplace=False):
        """Trim the :obj:`~conkit.core.sequencefile.SequenceFile`

        Parameters
        ----------
        start : int
           First residue to include
        end : int
           Final residue to include
        inplace : bool, optional
           Replace the saved order of sequences [default: False]

        Returns
        -------
        :obj:`~conkit.core.sequencefile.SequenceFile` 
           The reference to the :obj:`~conkit.core.sequencefile.SequenceFile`, regardless of inplace

        """
        sequence_file = self._inplace(inplace)
        if self.is_alignment:
            i = start - 1
            j = end
            for sequence in sequence_file:
                sequence.seq = sequence.seq[i:j]
            return sequence_file
        else:
            raise ValueError("This is not an alignment")

    def summary(self):
        """Generate a summary for the :obj:`~conkit.core.sequencefile.SequenceFile`
        
        Returns
        -------
        str
        
        """
        sstream = 'Summary for %s\n' % self.id
        sstream += '-------------------------------\n'
        sstream += 'Alignment:\t\t%s\n' % self.is_alignment
        sstream += 'Number of sequences:\t%d\n' % self.nseq
        sstream += 'Alignment depth (0.8):\t%d\n' % self.meff
        return sstream
