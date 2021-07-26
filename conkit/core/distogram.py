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
"""Distogram container used throughout ConKit"""

import numpy as np
from operator import attrgetter
from conkit.core.contactmap import ContactMap


class Distogram(ContactMap):
    """A distogram objecte to store all associated information. This class inherits methods and attributes
    from :obj:`~conkit.core.contactmap.ContactMap`

    Examples
    --------
    >>> import numpy as np
    >>> from conkit.core import Distance
    >>> from conkit.core import Distogram
    >>> distogram = Distogram("example")
    >>> distogram.add(Distance(1, 25, (0.25, 0.45, 0.25, 0.05), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
    >>> distogram.add(Distance(7, 19, (0.15, 0.15, 0.60, 0.1), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
    >>> print(distogram)
    ContactMap(id="example" ndistances=2)

    Attributes
    ----------
    id : str
       A unique identifier
    original_file_format : str
       The original file format used to create the :obj:`~conkit.core.distogram.Distogram` instance
    ndistances : int
       The number of :obj:`~conkit.core.distance.Distance` instances in the :obj:`~conkit.core.distogram.Distogram`

    """

    __slots__ = ["_original_file_format"]

    def __init__(self, id):
        self._original_file_format = None
        super(Distogram, self).__init__(id)

    def __repr__(self):
        return '{}(id="{}", ndistances={})'.format(self.__class__.__name__, self.id, self.ndistances)

    @property
    def ndistances(self):
        """The number of :obj:`~conkit.core.distance.Distance` instances

        Returns
        -------
        int
           The number of distance pairs in the :obj:`~conkit.core.distogram.Distogram`

        """
        return len(self)

    @property
    def original_file_format(self):
        """The original file format used to create the :obj:`~conkit.core.distogram.Distogram` instance"""
        return self._original_file_format

    @original_file_format.setter
    def original_file_format(self, value):
        self._original_file_format = value

    def get_unique_distances(self, inplace=False):
        """Filter the :obj:`~conkit.core.distance.Distance` instances so that each residue pairs is present only once

        Parameters
        ----------
        inplace : bool, optional
           Replace the saved order of contacts [default: False]

        Returns
        -------
        :obj:`~conkit.core.contactmap.ContactMap`
            :obj:`~conkit.core.contactmap.ContactMap` instance, regardless of inplace
        """
        distogram = self._inplace(inplace)
        unique_pairs = {tuple(sorted(el.id)): el for el in self}
        distogram.child_list = list(unique_pairs.values())
        distogram.child_dict = unique_pairs
        return distogram

    def get_absent_residues(self, seq_len=None):
        """Get residues not represented by any :obj:`~conkit.core.distance.Distance` instance

        Parameters
        ----------
        seq_len : int, optional
           Sequence length. If not provided, it will be pulled from :attr:`~conkit.core.contactmap.ContactMap.sequence`
           [default: None]

        Returns
        -------
        list
            A list of absent residues

        Raises
        ------
        :exc:`ValueError`
           No seq_len was provided and :attr:`~conkit.core.contactmap.ContactMap.sequence` is not defined
        """
        if seq_len is None:
            if self.sequence is None:
                raise ValueError('Need to define a sequence or provide seq_len')
            seq_len = self.sequence.seq_len

        absent_residues = []
        for residue in range(1, seq_len + 1):
            if not any([c.id for c in self if residue in c.id]):
                absent_residues.append(residue)
        return absent_residues

    def as_array(self, seq_len=None, get_weigths=False):
        """Transform the :obj:`~conkit.core.distogram.Distogram` instance into a :obj:numpy.array instance with shape
         (seq_len, seq_len) where each element represents the predicted distance between residues

        Parameters
        ----------
        seq_len : int, optional
           Sequence length. If not provided, it will be pulled from :attr:`~conkit.core.contactmap.ContactMap.sequence`
           [default: None]
        get_weigths : bool, optional
           If True the resulting array contains the confidence for each predicted distance rather than actual distances
           [default: False]

        Returns
        -------
        :obj:`numpy.array`
           :obj:`numpy.array` instance that represents the distogram. Note: change of residue indexing, now starts in 0

        Raises
        ------
        :exc:`ValueError`
           No seq_len was provided and :attr:`~conkit.core.contactmap.ContactMap.sequence` is not defined
        """
        if seq_len is None:
            if self.sequence is None:
                raise ValueError('Need to define a sequence or provide seq_len')
            seq_len = self.sequence.seq_len

        if get_weigths:
            getter = attrgetter('max_score')
        else:
            getter = attrgetter('predicted_distance')

        array = np.full((seq_len + 1, seq_len + 1), np.nan)
        for distance in self:
            array[distance.res1_seq, distance.res2_seq] = getter(distance)
            array[distance.res2_seq, distance.res1_seq] = getter(distance)
        array = np.delete(array, 0, axis=0)
        array = np.delete(array, 0, axis=1)

        return array
