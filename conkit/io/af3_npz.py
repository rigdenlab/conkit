# BSD 3-Clause License
#
# Copyright (c) 2016-25, University of Liverpool
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
Parser module specific to Alphafold3 NPZ distogram predictions
"""

import numpy as np
from conkit.io._parser import BinaryDistanceFileParser
from conkit.core.distance import Distance
from conkit.core.distogram import Distogram
from conkit.core.distancefile import DistanceFile

# magic that defines af3 distogram bin borders taken from af3 git
DISTANCE_BINS = ((0, 2.3125),(2.3125,2.625),(2.625,2.9375),(2.9375,3.25),(3.25,3.5625),(3.5625,3.875),
                 (3.875,4.1875),(4.1875,4.5),(4.5,4.8125),(4.8125,5.125),(5.125,5.4375),(5.4375,5.75),
                 (5.75,6.0625),(6.0625,6.375),(6.375,6.6875),(6.6875,7.0),(7.0,7.3125),(7.3125,7.625),
                 (7.625,7.9375),(7.9375,8.25),(8.25,8.5625),(8.5625,8.875),(8.875,9.1875),(9.1875,9.5),
                 (9.5,9.8125),(9.8125,10.125),(10.125,10.4375),(10.4375,10.75),(10.75,11.0625),(11.0625,11.375),
                 (11.375,11.6875),(11.6875,12.0),(12.0,12.3125),(12.3125,12.625),(12.625,12.9375),(12.9375,13.25),
                 (13.25,13.5625),(13.5625,13.875),(13.875,14.1875),(14.1875,14.5),(14.5,14.8125),(14.8125,15.125),
                 (15.125,15.4375),(15.4375,15.75),(15.75,16.0625),(16.0625,16.375),(16.375,16.6875),(16.6875,17.0),
                 (17.0,17.3125),(17.3125,17.625),(17.625,17.9375),(17.9375,18.25),(18.25,18.5625),(18.5625,18.875),
                 (18.875,19.1875),(19.1875,19.5),(19.5,19.8125),(19.8125,20.125),(20.125,20.4375),(20.4375,20.75),
                 (20.75,21.0625),(21.0625,21.375),(21.375,21.6875), (21.6875, np.inf)) 

class AF3NpzParser(BinaryDistanceFileParser):
    """Parser class for Alphafol 3 NPZ distogram file"""

    def read(self, f_handle, f_id="af3npz"):
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
        hierarchy.original_file_format = "af3npz"
        _map = Distogram("distogram_1")
        hierarchy.add(_map)

        prediction = np.load(f_handle, allow_pickle=True)
        distogram_array = prediction['distogram']

        L = distogram_array.shape[0]
        for i in range(L):
            for j in range(i, L):
                _distance = Distance(i + 1, j + 1, tuple(distogram_array[i, j, :].tolist()), DISTANCE_BINS)
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
