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
"""Residue distance prediction container used throughout ConKit"""

import math
import numpy as np
import statistics
from conkit.core.contact import Contact


class Distance(Contact):
    """A residue pair distance template to store all associated information. This class inherits methods and attributes
    from :obj:`~conkit.core.contact.Contact`

    Examples
    --------
    >>> import numpy as np
    >>> from conkit.core import Distance
    >>> distance = Distance(1, 25, (0.25, 0.45, 0.25, 0.05), ((0, 4), (4, 6), (6, 8), (8, np.inf)))
    >>> print(distance)
    Distance(id="(1, 25)" res1="A" res1_seq=1 res2="A" res2_seq=25 raw_score=0.95,
    distance_scores=(0.25, 0.45, 0.25, 0.05), distance_bins=((0, 4), (4, 6), (6, 8), (8, inf)))

    Attributes
    ----------
    distance_bins : tuple
       The distance boundaries of the bins associated to this residue pair in Ångstrom.
       Intervals are open on the left, i.e. a < d ≤ b
    distance_scores: tuple
        The prediction scores associated to each distance bin for this residue pair
    raw_score : float
       The prediction score for the residue pair to be within 8Å of each other

    """

    __slots__ = [
        "_distance_bound",
        "raw_score",
        "_res1",
        "_res2",
        "res1_chain",
        "res2_chain",
        "_res1_seq",
        "_res2_seq",
        "_res1_altseq",
        "_res2_altseq",
        "scalar_score",
        "_status",
        "weight",
        "distance_bins",
        "distance_scores"
    ]

    def __init__(self, res1_seq, res2_seq, distance_scores, distance_bins, raw_score=None, distance_bound=(0, 8)):
        """Initialize a generic distance residue pair

        Parameters
        ----------
        res1_seq : int
           The residue sequence number of residue 1
        res2_seq : int
           The residue sequence number of residue 2
        distance_scores: tuple
            The prediction score associated to the distance bins of this residue pair.
        distance_bins : tuple
           The lower and upper distance boundary values of the bins associated to this residue pair distance prediction.
        raw_score : float
           The covariance score for the contact pair
           Default is set to None, in which case the raw_score is calculated using distance_scores
        distance_bound : tuple, optional
           The lower and upper distance boundary values of a contact pair in Ångstrom.
           Default is set to between 0.0 and 8.0 Å.

        """
        self.distance_bins = distance_bins
        self.distance_scores = distance_scores
        self.parent = None
        if raw_score is None:
            raw_score = self.get_probability_within_distance(distance_bound[-1])

        super(Distance, self).__init__(res1_seq, res2_seq, raw_score, distance_bound)

    def __repr__(self):
        text = (
            "{name}(id={id} res1={_res1} res1_chain={res1_chain} res1_seq={_res1_seq} "
            "res2={_res2} res2_chain={res2_chain} res2_seq={_res2_seq} raw_score={raw_score} "
            "distance_bins={distance_bins} distance_scores={distance_scores})"
        )
        return text.format(name=self.__class__.__name__, **{k: getattr(self, k) for k in self.__dir__()})

    @property
    def max_score(self):
        """Maximum confidence score observed across the different distance bins"""
        return max(self.distance_scores)

    @property
    def predicted_distance_bin(self):
        """Distance bin with the highest confidence score"""
        bin_index = self.distance_scores.index(self.max_score)
        return self.distance_bins[bin_index]

    @property
    def predicted_distance(self):
        """Median distance associated with the distance bin with the highest confidence score"""
        return statistics.median(self.predicted_distance_bin)

    def get_probability_within_distance(self, distance):
        """Calculate the probability that the residue pair is within a given distance

        Parameters
        ----------
        distance : int, float

        Returns
        -------
        None, float
           The probability that the residue pair is within the specified distance

        Raises
        ------
        :exc:`ValueError`
           distance is not a positive number
        """
        if not self.distance_bins:
            raise ValueError('No distance bins have been defined')
        elif self.parent is not None and self.parent.original_file_format == 'PDB':
            if self.distance_bins[0][-1] < distance:
                return 1.0
            return 0.0
        elif distance == 0:
            return 0.0
        elif distance < 0:
            raise ValueError('Distance must be a positive value')

        probability = 0
        for distance_score, distance_bin in zip(self.distance_scores, self.distance_bins):
            # Last bin is special case because interval goes to Inf
            if np.isinf(distance_bin[1]):
                factor = math.e ** (-distance) / math.e ** (-distance_bin[0])
                probability += distance_score * (1 - factor)
                break
            # Assume other bins have continuous probability
            elif distance_bin[0] < distance <= distance_bin[1]:
                bin_diff = distance_bin[1] - distance_bin[0]
                distance_diff = distance - distance_bin[0]
                probability += distance_score / bin_diff * distance_diff
                break
            probability += distance_score

        return probability

    def reshape_bins(self, new_bins):
        """Reshape the predicted distance bins and update :attr:`~conkit.core.distance.Distance.distance_scores` and
        :attr:`~conkit.core.distance.Distance.distance_bins` accordingly

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
        if self.parent is not None and self.parent.original_file_format == 'PDB':
            raise ValueError('Cannot re-shape bins obtained from a PDB structure file')
        self._assert_valid_bins(new_bins)
        self._reshape_bins(new_bins)

    def _reshape_bins(self, new_bins):
        """Reshape the predicted distance bins and update :attr:`~conkit.core.distance.Distance.distance_scores` and
        :attr:`~conkit.core.distance.Distance.distance_bins` accordingly

        Parameters
        ----------
        new_bins : tuple
           A tuple of tuples, where each element corresponds with the upper and lower edges of the intervals for
           the new distance bins
        """
        new_distance_scores = []
        for current_new_bin in new_bins:
            probability_lower_bound = self.get_probability_within_distance(current_new_bin[0])
            probability_upper_bound = self.get_probability_within_distance(current_new_bin[1])
            new_probability = probability_upper_bound - probability_lower_bound
            new_distance_scores.append(new_probability)

        self.distance_bins = tuple(new_bins)
        self.distance_scores = tuple(new_distance_scores)

    def as_contact(self, distance_cutoff):
        """Create a :obj:`~conkit.core.contact.Contact` instance with the information in this
        :obj:`~conkit.core.distance.Distance` instance.

        Parameters
        ----------
        distance_cutoff : int, float
           The distance cutoff used to consider a residue pair within contact of each other

        Returns
        -------
        :obj:`~conkit.core.contact.Contact`
            A contact with the information present in this distance instance.
        """
        if self.predicted_distance > distance_cutoff:
            contact = Contact(self.res1_seq, self.res2_seq, 0, distance_bound=(0, distance_cutoff))
        else:
            contact = Contact(self.res1_seq, self.res2_seq, self.raw_score, distance_bound=(0, distance_cutoff))
        for attr in self.__slots__:
            if hasattr(contact, attr):
                setattr(contact, attr, getattr(self, attr))
        return contact

    @staticmethod
    def _assert_valid_bins(distance_bins):
        """Determine whether a set of distance bins is valid. Valid distance bins must follow these rules:
            - There is more than one bin
            - Lower limit of first bin is 0
            - Upper limit of last bin is Inf
            - Only Inf value is the upper limit of last bin
            - Upper and lower limits of bins are different
            - Upper limit is higher than lower limit across all bins

        Parameters
        ----------
        distance_bins : tuple
           The tuple with the distance bins to be tested

        Raises
        ------
        :exc:`ValueError`
           The distance bins are not valid
        """
        if len(distance_bins) <= 1:
            raise ValueError('New distance bins are invalid')
        elif not np.isinf(distance_bins[-1][1]):
            raise ValueError('New distance bins are invalid')
        elif distance_bins[0][0] != 0:
            raise ValueError('New distance bins are invalid')

        temp_list = []
        for dbin in distance_bins:
            if len(dbin) != 2 or dbin[0] >= dbin[1]:
                raise ValueError('New distance bins are invalid')
            temp_list += list(dbin)

        if len(np.unique(temp_list[1:-1])) != len(distance_bins) - 1 or temp_list.count(np.inf) != 1:
            raise ValueError('New distance bins are invalid')
