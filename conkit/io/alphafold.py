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
"""
Parser module specific to AF2 distance predictions
"""

import numpy as np
from scipy.special import softmax
from conkit.io._parser import BinaryDistanceFileParser
from conkit.core.distance import Distance
from conkit.core.distogram import Distogram
from conkit.core.distancefile import DistanceFile


class AlphaFold2Parser(BinaryDistanceFileParser):
    """Parser class for AF2 distance prediction file"""

    def read(self, f_handle, f_id="alphafold2"):
        """Read a distance prediction file

        Parameters
        ----------
        f_handle
           Open file handle [read permissions]
        f_id : str, optional
           Unique contact file identifier

        Returns
        -------
        :obj:`~conkit.core.distancefile.DistanceFile`

        """

        hierarchy = DistanceFile(f_id)
        hierarchy.original_file_format = "ALPHAFOLD2"
        _map = Distogram("distogram_1")
        hierarchy.add(_map)

        prediction = np.load(f_handle, allow_pickle=True)
        predicted_distogram = prediction['distogram']
        probs = softmax(predicted_distogram['logits'], axis=-1)
        bin_edges = predicted_distogram['bin_edges']

        distance_bins = [(0, bin_edges[0])]
        distance_bins += [(bin_edges[idx], bin_edges[idx + 1]) for idx in range(len(bin_edges) - 1)]
        distance_bins.append((bin_edges[-1], np.inf))
        distance_bins = tuple(distance_bins)
        L = probs.shape[0]
        for i in range(L):
            for j in range(i, L):
                _distance = Distance(i + 1, j + 1, tuple(probs[i, j, :].tolist()), distance_bins)
                _map.add(_distance)

        return hierarchy

    def write(self, f_handle, hierarchy):
        """Write a distance file instance to a file

        Raises
        ------
        :exc:`NotImplementedError`
           Write function not available

        """
        raise NotImplementedError("Write function not available yet")
