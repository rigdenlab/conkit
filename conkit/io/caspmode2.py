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
Parser module specific to CASP-RR MODE 2 distance predictions
"""

import numpy as np
from conkit.io._parser import DistanceFileParser
from conkit.core.distance import Distance
from conkit.core.distogram import Distogram
from conkit.core.distancefile import DistanceFile

DISNTACE_BINS = ((0, 4), (4, 6), (6, 8), (8, 10), (10, 12), (12, 14), (14, 16), (16, 18), (18, 20), (20, np.inf))


class CaspMode2Parser(DistanceFileParser):
    """Parser class for CASP RR MODE 2 distance prediction file"""

    def read(self, f_handle, f_id="casp2"):
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
        hierarchy.original_file_format = "CASPRR_MODE_2"
        _map = Distogram("distogram_1")
        hierarchy.add(_map)

        for idx, line in enumerate(f_handle.readlines()):
            line = line.lstrip().rstrip().split()
            if not line or len(line) != 13 or not line[0].isdigit() or not line[1].isdigit():
                continue

            res1_seq = int(line[0])
            res2_seq = int(line[1])
            raw_score = float(line[2])
            distance_scores = tuple([float(p) for p in line[3:]])
            _distance = Distance(res1_seq, res2_seq, distance_scores, DISNTACE_BINS, raw_score=raw_score)
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

