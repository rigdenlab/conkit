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
from conkit.core.distance import Distance
from conkit.core.contact import Contact
from conkit.core.contactmap import ContactMap


class Distogram(ContactMap):
    """A distogram object to store all associated information. This class inherits methods and attributes
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
    Distogram(id="example" ndistances=2)

    Attributes
    ----------
    id : str
       A unique identifier
    original_file_format : str
       The original file format used to create the :obj:`~conkit.core.distogram.Distogram` instance
    ndistances : int
       The number of :obj:`~conkit.core.distance.Distance` instances in the :obj:`~conkit.core.distogram.Distogram`

    """

    __slots__ = ["_original_file_format", "_sequence"]

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

        if seq_len < self.highest_residue_number:
            raise ValueError('Sequence length does not match contact map')

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

    def reshape_bins(self, new_bins):
        """Reshape the predicted distance bins for all :obj:`~conkit.core.distance.Distance` instances. This will
        update :attr:`~conkit.core.distance.Distance.distance_scores` and
        :attr:`~conkit.core.distance.Distance.distance_bins` to fit the new bins.

        Parameters
        ----------
        new_bins : tuple
           A tuple of tuples, where each element corresponds with the upper and lower edges of the intervals for
           the new distance bins

        Raises
        ------
        :exc:`ValueError`
           The new distance bins are not valid
        """
        if self.original_file_format == 'PDB':
            raise ValueError('Cannot re-shape bins obtained from a PDB structure file')
        Distance._assert_valid_bins(new_bins)

        for distance in self:
            distance._reshape_bins(new_bins)

    def as_contactmap(self, distance_cutoff=8):
        """Create a :obj:`~conkit.core.contactmap.ContactMap` instance with the contacts present in this
        :obj:`~conkit.core.distogram.Distogram` instance.

        Parameters
        ----------
        distance_cutoff : int, float
           The distance cutoff used to consider a residue pair within contact of each other

        Returns
        -------
        :obj:`~conkit.core.contactmap.ContactMap`
            A contactmap with the contacts present in this distogram instance.
        """
        contactmap = ContactMap("map_1")
        for dist in self:
            if dist.predicted_distance <= distance_cutoff:
                contact = dist.as_contact(distance_cutoff)
                contactmap.add(contact)
        if self.sequence is not None:
            contactmap.sequence = self.sequence._inplace(False)

        return contactmap

    @staticmethod
    def calculate_rmsd(prediction, model, seq_len=None, calculate_wrmsd=False):
        """Calculate the RMSD between two :obj:`~conkit.core.distogram.Distogram` instances.

        Parameters
        ----------
        prediction: :obj:`~conkit.core.distogram.Distogram`
           A ConKit :obj:`~conkit.core.distogram.Distogram` used as the prediction for the RMSD
        model: :obj:`~conkit.core.distogram.Distogram`
           A ConKit :obj:`~conkit.core.distogram.Distogram` used as the model to calculate the RMSD
        seq_len: int, optional
           Sequence length. If not provided, it will be pulled from :attr:`~conkit.core.contactmap.ContactMap.sequence`
           [default: None]
        calculate_wrmsd: bool
           If set to True wRMSD is calculated using distance confidence scores [default: False]

        Returns
        -------
        list
            A list of floats with the RMSD values along the sequence

        Raises
        ------
        :exc:`ValueError`
           other is not a  :obj:`~conkit.core.distogram.Distogram` instance.
        """
        if not isinstance(model, Distogram) or not isinstance(prediction, Distogram):
            raise ValueError('Need to provide a conkit.core.distogram.Distogram instance')

        max_distance = prediction.top.distance_bins[-1][0]

        model_array = model.as_array(seq_len=seq_len)
        model_array[model_array > max_distance] = max_distance
        prediction_array = prediction.as_array(seq_len=seq_len)
        prediction_array[prediction_array > max_distance] = max_distance

        if prediction_array.shape != model_array.shape:
            raise ValueError('Distograms cannot be matched')

        difference = prediction_array - model_array
        squared_difference = difference ** 2

        if calculate_wrmsd:
            prediction_weights = prediction.as_array(seq_len=seq_len, get_weigths=True)
            squared_difference *= prediction_weights

        sum_squared_differences = np.nansum(squared_difference, axis=0)
        n_observations_array = np.sum(~np.isnan(squared_difference), axis=0)
        rmsd = np.sqrt(sum_squared_differences / n_observations_array)
        return rmsd

    def find_residues_within(self, resnum, distance_cutoff):
        """Find all residues within a given distance of a given residue

        Parameters
        ----------
        resnum: int
           The residue number of the residue of interest
        distance_cutoff: int, float
           The distance cutoff used to find residues

        Returns
        -------
        set
           A set with the residue numbers of residues within the given distance
        """
        result = []

        for distance in self:
            if distance.predicted_distance <= distance_cutoff and resnum in distance.id:
                result += list(distance.id)

        return set(result)

    @staticmethod
    def merge_arrays(distogram_1, distogram_2):
        """Take two :obj:`~conkit.core.distogram.Distogram` instances and merge them together into the same
        :obj:`numpy.array` instance. Each half square in this array will correspond with the predicted distances
        at each hierarchy

        Parameters
        ----------
        distogram_1: :obj:`~conkit.core.distogram.Distogram`
           First :obj:`~conkit.core.distogram.Distogram` instance, used to populate top half square of the array
        distogram_2: :obj:`~conkit.core.distogram.Distogram`
           Second :obj:`~conkit.core.distogram.Distogram` instance, used to populate lower half square of the array

        Returns
        -------
        :obj:`numpy.array`
           :obj:`numpy.array` instance that represents the combined distograms.

        Raises
        ------
        :exc:`ValueError`
           No sequence has been registered for one of the :obj:`~conkit.core.distogram.Distogram` instances
        :exc:`ValueError`
            The sequence length associated to the :obj:`~conkit.core.distogram.Distogram` instances is incompatible
        """

        if distogram_1.sequence is None or distogram_2.sequence is None:
            raise ValueError("All hierarchies must have a sequence registered")
        if distogram_1.sequence.seq_len != distogram_2.sequence.seq_len:
            raise ValueError("Sequence lengths are incompatible")

        array = np.full((distogram_1.sequence.seq_len + 1, distogram_1.sequence.seq_len + 1), np.nan)

        for distance in distogram_1:
            array[distance.res1_seq, distance.res2_seq] = distance.predicted_distance
        for distance in distogram_2:
            array[distance.res2_seq, distance.res1_seq] = distance.predicted_distance

        array = np.delete(array, 0, axis=0)
        array = np.delete(array, 0, axis=1)
        return array
