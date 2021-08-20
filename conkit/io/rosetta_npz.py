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
Parser module specific to rosetta NPZ distance predictions
"""

import numpy as np
from conkit.io._parser import BinaryDistanceFileParser
from conkit.core.distance import Distance
from conkit.core.distogram import Distogram
from conkit.core.distancefile import DistanceFile

DISTANCE_BINS = ((0, 2), (2, 2.5), (2.5, 3), (3, 4), (4, 4.5), (4.5, 5), (5, 5.5), (5.5, 6), (6, 6.5), (6.5, 7),
                 (7, 7.5), (7.5, 8), (8, 8.5), (8.5, 9), (9, 9.5), (9.5, 10), (10, 10.5), (10.5, 11), (11, 11.5),
                 (11.5, 12), (12, 12.5), (12.5, 13), (13, 13.5), (13.5, 14), (14, 14.5), (14.5, 15), (15, 15.5),
                 (15.5, 16), (16, 16.5), (16.5, 17), (17, 17.5), (17.5, 18), (18, 18.5), (18.5, 19), (19, 19.5),
                 (19.5, 20), (20, np.inf))


class RosettaNpzParser(BinaryDistanceFileParser):
    """Parser class for rosetta NPZ distance prediction file"""

    def read(self, f_handle, f_id="rosettanpz"):
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
        hierarchy.original_file_format = "ROSETTA_NPZ"
        _map = Distogram("distogram_1")
        hierarchy.add(_map)

        prediction = np.load(f_handle, allow_pickle=True)
        probs = prediction['dist']
        # Bin #0 corresponds with d>20A & bins #1 ~ #36 correspond with 2A<d<20A in increments of 0.5A
        probs = probs[:, :, [x for x in range(1, 37)] + [0]]

        L = probs.shape[0]
        for i in range(L):
            for j in range(i, L):
                _distance = Distance(i + 1, j + 1, tuple(probs[i, j, :].tolist()), DISTANCE_BINS)
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
