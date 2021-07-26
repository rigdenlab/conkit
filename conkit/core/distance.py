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

import statistics
from conkit.core.contact import Contact


class Distance(Contact):
    """A residue pair distance template to store all associated information. This class inherits methods and attributes
    from :obj:`~conkit.core.distance.Distance`

    Examples
    --------
import numpy as np    >>> from conkit.core import Distance
    >>> distance = Distance(1, 25, (0.25, 0.45, 0.25, 0.05), ((0, 4), (4, 6), (6, 8), (8, 10)))
    >>> print(distance)
    Distance(id="(1, 25)" res1="A" res1_seq=1 res2="A" res2_seq=25 raw_score=0.95,
    distance_scores=(0.25, 0.45, 0.25, 0.05), distance_bins=((0, 4), (4, 6), (6, 8), (8, 10)))

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
        if raw_score is None:
            raw_score = self.get_probability_within_distance(distance_bound[-1])

        super(Distance, self).__init__(res1_seq, res2_seq, raw_score, distance_bound)

    def __repr__(self):
        text = (
            "{name}(id={id} res1={_res1} res1_chain={res1_chain} res1_seq={_res1_seq} "
            "res2={_res2} res2_chain={res2_chain} res2_seq={_res2_seq} raw_score={raw_score}"
            "distance_bins={distance_bins} distance_scores={distance_scores})"
        )
        return text.format(name=self.__class__.__name__, id=self._id, **{k: getattr(self, k) for k in self.__slots__})

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
        None, int, float
           The probability that the residue pair is within the specified distance

        """

        if not self.distance_bins:
            return None
        elif len(self.distance_bins) == 1:
            if self.distance_bins[0][-1] < distance:
                return 0
            return 1
        elif self.distance_bins[-1][1] <= distance:
            return 1

        for bin_index, upper_limit in enumerate((d_bin[-1] for d_bin in self.distance_bins)):
            if upper_limit > distance:
                break

        return sum(self.distance_scores[:bin_index])
