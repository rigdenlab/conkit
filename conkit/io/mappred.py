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
Parser module specific to MapPred distance predictions
"""

import numpy as np
from conkit.io._parser import DistanceFileParser
from conkit.core.distance import Distance
from conkit.core.distogram import Distogram
from conkit.core.distancefile import DistanceFile

DISTANCE_BINS = ((0, 4), (4, 4.5), (4.5, 5), (5, 5.5), (5.5, 6), (6, 6.5), (6.5, 7), (7, 7.5), (7.5, 8), (8, 8.5),
                 (8.5, 9), (9, 9.5), (9.5, 10), (10, 10.5), (10.5, 11), (11, 11.5), (11.5, 12), (12, 12.5), (12.5, 13),
                 (13, 13.5), (13.5, 14), (14, 14.5), (14.5, 15), (15, 15.5), (15.5, 16), (16, 16.5), (16.5, 17),
                 (17, 17.5), (17.5, 18), (18, 18.5), (18.5, 19), (19, 19.5), (19.5, 20), (20, np.inf))


class MapPredParser(DistanceFileParser):
    """Parser class for MapPred distance prediction file"""

    def read(self, f_handle, f_id="mappred"):
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
        hierarchy.original_file_format = "MAPPRED"
        _map = Distogram("distogram_1")
        hierarchy.add(_map)

        for line in f_handle.readlines():
            line = line.lstrip().rstrip().split()
            if not line or len(line) != 36 or not line[0].isdigit() or not line[1].isdigit():
                continue

            res1_seq = int(line[0])
            res2_seq = int(line[1])
            distance_scores = tuple([float(p) for p in line[2:]])
            _distance = Distance(res1_seq, res2_seq, distance_scores, DISTANCE_BINS)
            _map.add(_distance)

        return hierarchy

    def write(self, f_handle, hierarchy):
        """Write a distance prediction file instance to a file

        Parameters
        ----------
        f_handle
           Open file handle [write permissions]
        hierarchy : :obj:`~conkit.core.distancefile.DistanceFile`, :obj:`~conkit.core.distogram.Distogram`
                    or :obj:`~conkit.core.distance.Distance`

        Raises
        ------
        :exc:`RuntimeError`
           More than one contact map in the hierarchy

        """
        distancefile = self._reconstruct(hierarchy)
        if len(distancefile) > 1:
            raise RuntimeError("More than one distogram provided")
        distogram = distancefile.top_map

        content = "#REMARK MapPred 1.1\n#REMARK idx_i, idx_j, distance distribution of 34 bins\n#REMARK 34 bins " \
                  "consist of 32 normal bins (4-20A with a step of 0.5A) and two boundary bins ( [0,4) and [20, inf) " \
                  "), as follows: [0,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15," \
                  "15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,inf]\n"
        line_template = "{} {}" + " {:.6f}" * 34 + "\n"
        for distance in distogram:
            distance.reshape_bins(DISTANCE_BINS)
            content += line_template.format(distance.res1_seq, distance.res2_seq, *distance.distance_scores)
        f_handle.write(content)
