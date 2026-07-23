# coding=utf-8
#
# BSD 3-Clause License
#
# Copyright (c) 2016-21, University of Liverpool
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
"""ConfidenceMap container used throughout ConKit"""

from __future__ import division

__author__ = "Aderik Voorspoels"
__date__ = "07 Apr 2025"
__version__ = "0.14.1"

import collections
import numpy as np
import operator
import sys

from conkit.core.entity import Entity
from conkit.core.struct import Gap, Residue
#from conkit.core.mappings import AminoAcidMapping, ContactMatchState
from conkit.core.sequence import Sequence
#from conkit.misc import normalize


class ConfidenceMap(Entity):
    """A confidence map object representing a single prediction

    The :obj:`~conkit.core.confidencemap.ConfidenceMap` class represents a data structure to hold a single
    confidence map prediction in one place. It contains functions to store, manipulate and organise
    :obj:`~conkit.core.confidence.Confidence` instances.

    Examples
    --------
    >>> from conkit.core import Confidence, ConfidenceMap
    >>> confidence_map = ConfidenceMap("example")
    >>> confidence_map.add(Confidence(1, 10, 0.333))
    >>> confidence_map.add(Confidence(5, 30, 0.667))
    >>> print(confidence_map)
    ConfidenceMap(id="example" nconfidences=2)

    Attributes
    ----------
    coverage : float
       The sequence coverage score
    id : str
       A unique identifier
    nconfidences : int
       The number of :obj:`~conkit.core.confidence.Confidence` instances in the :obj:`~conkit.core.confidencemap.ConfidenceMap`
    precision : float
       The precision (Positive Predictive Value) score
    repr_sequence : :obj:`~conkit.core.sequence.Sequence`
       The representative :obj:`~conkit.core.sequence.Sequence` associated with the :obj:`~conkit.core.confidencemap.ConfidenceMap`
    repr_sequence_altloc : :obj:`~conkit.core.sequence.Sequence`
       The representative altloc :obj:`~conkit.core.sequence.Sequence` associated with the :obj:`~conkit.core.confidencemap.ConfidenceMap`
    sequence : :obj:`~conkit.core.sequence.Sequence`
       The :obj:`~conkit.core.sequence.Sequence` associated with the :obj:`~conkit.core.confidencemap.ConfidenceMap`
    top_confidence : :obj:`~conkit.core.confidence.Confidence`
       The first :obj:`~conkit.core.confidence.Confidence` entry

    """

    __slots__ = ["_sequence","_type","_per"]

    def __init__(self, id):
        """Initialise a new confidence map"""
        self._sequence = None
        self._type = None
        self._per = None
        super(ConfidenceMap, self).__init__(id)  ## give the confidence map all atributes and methods of Entity
              
    def __repr__(self):
        return '{}(id="{}", contains={}, by={})'.format(self.__class__.__name__, self.id, self._type, self._per)


    @property
    def empty(self):
        """Empty confidence map"""
        return len(self) < 1

    @property
    def nconfidences(self):
        """The number of :obj:`~conkit.core.confidence.Confidence` instances

        Returns
        -------
        int
           The number of confidences in the :obj:`~conkit.core.confidencemap.ConfidenceMap`

        """
        return len(self)


    @property
    def repr_sequence(self):
        """The representative :obj:`~conkit.core.sequence.Sequence` associated
        with the :obj:`~conkit.core.confidencemap.ConfidenceMap`

        The peptide sequence constructed from the available
        confidences using the normal res_seq positions

        Returns
        -------
        :obj:`~conkit.core.sequence.Sequence`

        Raises
        ------
        :exc:`TypeError`
           Sequence undefined

        See Also
        --------
        repr_sequence_altloc, sequence

        """
        if isinstance(self.sequence, Sequence):
            res_seqs = np.unique(np.array(self.as_list()).flatten()).tolist()
            return self._construct_repr_sequence(res_seqs)
        else:
            raise TypeError("Define the sequence as Sequence() instance")

    @property
    def repr_sequence_altloc(self):
        """The representative altloc :obj:`~conkit.core.sequence.Sequence` associated
        with the :obj:`~conkit.core.confidencemap.ConfidenceMap`

        The peptide sequence constructed from the available
        confidences using the :attr:`~conkit.core.confidence.Confidence.res_altseq` positions

        Returns
        -------
        :obj:`~conkit.core.sequence.Sequence`

        Raises
        ------
        :exc:`TypeError`
           Sequence undefined

        See Also
        --------
        repr_sequence, sequence

        """
        if isinstance(self.sequence, Sequence):
            res_seqs = np.unique(np.array(self.as_list(altloc=True)).flatten()).tolist()
            return self._construct_repr_sequence(res_seqs)
        else:
            raise TypeError("Define the sequence as Sequence() instance")

    @property
    def highest_residue_number(self):
        """The highest residue sequence number among confidences in the :obj:`~conkit.core.confidencemap.ConfidenceMap`

        Returns
        -------
        int
           Highest residue sequence number in the confidence map

        """
        if len(self) == 0:
            return None
        else:
            return max([max(confidence.id) for confidence in self])

    @property
    def sequence(self):
        """The :obj:`~conkit.core.sequence.Sequence` associated with the :obj:`~conkit.core.confidencemap.ConfidenceMap`

        Returns
        -------
        :obj:`~conkit.core.sequence.Sequence`
           A :obj:`~conkit.core.sequence.Sequence` object

        See Also
        --------
        repr_sequence, repr_sequence_altloc

        """
        return self._sequence

    @sequence.setter
    def sequence(self, sequence):
        """Associate a :obj:`~conkit.core.sequence.Sequence` instance with the :obj:`~conkit.core.confidencemap.ConfidenceMap`

        Parameters
        ----------
        sequence : :obj:`~conkit.core.sequence.Sequence`

        Raises
        ------
        :exc:`TypeError`
           Incorrect hierarchy instance provided

        """
        if isinstance(sequence, Sequence):
            self._sequence = sequence
        else:
            raise TypeError("Instance of Sequence() required: {}".format(sequence))


    @property
    def top_confidence(self):
        """The first :obj:`~conkit.core.confidence.Confidence` entry

        Returns
        -------
        :obj:`~conkit.core.confidence.Confidence`
           The first :obj:`~conkit.core.confidence.Confidence` entry in :obj:`~conkit.core.confidencefile.ConfidenceFile`

        """
        return self.top

    def _construct_repr_sequence(self, res_seqs):
        """Construct the representative sequence"""
        representative_sequence = ""
        for i in np.arange(1, self.sequence.seq_len + 1):
            if i in res_seqs:
                representative_sequence += self.sequence.seq[i - 1]
            else:
                representative_sequence += "-"
        return Sequence(self.sequence.id + "_repr", representative_sequence)

    def as_dict(self, altloc=False):
        """The :obj:`~conkit.core.confidencemap.ConfidenceMap` as a dictionary where each key corresponds with the residue
        number and the values are sets of tuples with the :attr:`~conkit.core.confidence.Confidence.id`

        Parameters
        ----------
        altloc : bool
           Use the :attr:`~conkit.core.confidence.Confidence.res_altloc` positions [default: False]

        Returns
        -------
        dict
            A dictionary representation of the :obj:`~conkit.core.confidencemap.ConfidenceMap` instance
        """
        if self.sequence is None:
            seq_len = self.highest_residue_number
        else:
            seq_len = len(self.sequence)

        result = {}

        if altloc:
            for resn in range(1, seq_len + 1):
                result[resn] = {(c.res2_altseq, c.res2_altseq) for c in self if resn in (c.res1_altseq, c.res2_altseq)}
        else:
            for resn in range(1, seq_len + 1):
                result[resn] = {(c.res1_seq, c.res2_seq) for c in self if resn in (c.res1_seq, c.res2_seq)}

        return result

    def as_list(self, altloc=False):
        """The :obj:`~conkit.core.confidencemap.ConfidenceMap` as a 2D-list containing confidence-pair residue indexes

        Parameters
        ----------
        altloc : bool
           Use the :attr:`~conkit.core.confidence.Confidence.res_altloc` positions [default: False]

        """
        if altloc:
            return [[c.res1_altseq, c.res2_altseq] for c in self]
        else:
            return [[c.res1_seq, c.res2_seq] for c in self]

    def as_set(self, altloc=False):
        """The :obj:`~conkit.core.confidencemap.ConfidenceMap` as a 2D-set containing confidence-pair residue indexes

        Parameters
        ----------
        altloc : bool
           Use the :attr:`~conkit.core.confidence.Confidence.res_altloc` positions [default: False]

        """
        if altloc:
            return {(c.res1_altseq, c.res2_altseq) for c in self}
        else:
            return {c.id for c in self}

    def set_sequence_register(self, altloc=False):
        """Assign the amino acids from :obj:`~conkit.core.sequence.Sequence` to all :obj:`~conkit.core.confidence.Confidence` instances

        Parameters
        ----------
        altloc : bool
           Use the :attr:`~conkit.core.confidence.Confidence.res_altloc` positions [default: False]

        Raises
        ------
        :exc:`ValueError`
           Undefined sequence
        """
        if self.sequence is None:
            raise ValueError("No sequence defined")

        seq_len = len(self.sequence)

        for c in self:
            if altloc:
                res1_index = c.res1_altseq
                res2_index = c.res2_altseq
            else:
                res1_index = c.res1_seq
                res2_index = c.res2_seq
            if res1_index <= seq_len and res2_index <= seq_len:
                c.res1 = self.sequence.seq[res1_index - 1]
                c.res2 = self.sequence.seq[res2_index - 1]
            else:
                raise ValueError('Confidence {} is out of sequence bounds'.format(c.id))

    def reindex(self, index, altloc=False, inplace=False):
        """Re-index the :obj:`~conkit.core.confidencemap.ConfidenceMap`

        Parameters
        ----------
        index : int
           The new starting index [assigned to the lowest existing index in the confidence map]
        altloc : bool
           Use the res_altloc positions [default: False]
        inplace : bool
           Replace the saved order of confidences [default: False]

        Returns
        -------
        :obj:`~conkit.core.confidencemap.ConfidenceMap`
           The reference to the :obj:`~conkit.core.confidencemap.ConfidenceMap`, regardless of inplace

        Raises
        ------
        :exc:`ValueError`
           Index must be positive

        """
        if index < 0:
            raise ValueError("Index must be positive!")
        confidence_map = self._inplace(inplace)
        if confidence_map.empty:
            return confidence_map
        res1s, res2s = zip(*confidence_map.as_list(altloc=altloc))
        offset = min(res1s) - index
        for confidence in confidence_map:
            if altloc:
                confidence.res1_altseq -= offset
                confidence.res2_altseq -= offset
            else:
                confidence.res1_seq -= offset
                confidence.res2_seq -= offset
        for confidence in confidence_map:
            confidence.id = (confidence.res1_seq, confidence.res2_seq)
        return confidence_map

    def to_string(self):
        """Return the :obj:`ConfidenceMap <conkit.core.confidencemap.ConfidenceMap>` as :obj:`str`"""
        content = ["%d\t%d\t%.5f" % (c.res1_seq, c.res2_seq, c.raw_score) for c in self]
        return "\n".join(content)

    @staticmethod
    def _adjust(confidence_map, keymap):
        """Adjust res_altseq entries to insertions and deletions"""
        encoder = dict((x.res_seq, x.res_altseq) for x in keymap if isinstance(x, Residue))
        for confidence in confidence_map:
            if confidence.res1_seq in encoder:
                confidence.res1_altseq = encoder[confidence.res1_seq]
            if confidence.res2_seq in encoder:
                confidence.res2_altseq = encoder[confidence.res2_seq]
        return confidence_map

    @staticmethod
    def _create_keymap(confidence_map, altloc=False):
        """Create a simple keymap

        Parameters
        ----------
        altloc : bool
           Use the res_altloc positions [default: False]

        Returns
        -------
        list
           A list of residue mappings

        """
        confidence_map_keymap = collections.OrderedDict()
        for confidence in confidence_map:
            pos1 = Residue(confidence.res1_seq, confidence.res1_altseq, confidence.res1, confidence.res1_chain)
            pos2 = Residue(confidence.res2_seq, confidence.res2_altseq, confidence.res2, confidence.res2_chain)
            if altloc:
                res1_index, res2_index = confidence.res1_altseq, confidence.res2_altseq
            else:
                res1_index, res2_index = confidence.res1_seq, confidence.res2_seq
            confidence_map_keymap[res1_index] = pos1
            confidence_map_keymap[res2_index] = pos2
        confidence_map_keymap_sorted = sorted(list(confidence_map_keymap.items()), key=lambda x: int(x[0]))
        return list(zip(*confidence_map_keymap_sorted))[1]

    @staticmethod
    def _find_single(confidence_map, index):
        """Find all confidences associated with ``index`` based on id property"""
        for c in confidence_map:
            if c.id[0] == index or c.id[1] == index:
                yield c

    @staticmethod
    def _insert_states(sequence, keymap):
        """Create a sequence matching keymap including deletions and insertions"""
        it = iter(keymap)
        keymap_ = []
        for amino_acid in sequence:
            if amino_acid == ord("-"):
                keymap_.append(Gap())
            else:
                keymap_.append(next(it))
        return keymap_

    @staticmethod
    def _reindex_by_keymap(keymap):
        """Reindex a key map"""
        return [residue._replace(res_altseq=i + 1) for i, residue in enumerate(keymap)]

    @staticmethod
    def _renumber(confidence_map, self_keymap, other_keymap):
        """Renumber the confidence map based on the mapping of self and other keymaps"""
        for self_residue, other_residue in zip(self_keymap, other_keymap):
            if isinstance(self_residue, Gap):
                continue
            for confidence in ConfidenceMap._find_single(confidence_map, self_residue.res_seq):
                # Make sure we check with the ID, which doesn't change
                if confidence.id[0] == self_residue.res_altseq:
                    confidence.res1_seq = other_residue.res_seq
                    confidence.res1_chain = other_residue.res_chain
                elif confidence.id[1] == self_residue.res_altseq:
                    confidence.res2_seq = other_residue.res_seq
                    confidence.res2_chain = other_residue.res_chain
                else:
                    raise ValueError("Error renumbering confidence map - please report this bug")

        return confidence_map
