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
"""
Command line object for map_align contact map alignment application
"""

from Bio.Application import _Option
from Bio.Application import AbstractCommandline


class MapAlignCommandline(AbstractCommandline):
    """
    Command line object for map_align [#]_ [#]_

    https://github.com/sokrypton/map_align


    Examples
    --------
    >>> from conkit.applications import MapAlignCommandline
    >>> mapalign_cline = MapAlignCommandline(map_a='cmap_1.mapalign', map_b='cmap_2.mapalign')
    >>> print(mapalign_cline)
    map_align cmap_1.mapalign cmap_2.mapalign

    You would typically run the command line with :func:`mapalign_cline` or via
    the :mod:`~subprocess` module.

    """

    def __init__(self, cmd="map_align", **kwargs):
        self.parameters = [
            _Option(
                ["-a", "contact_map_a"],
                "contact map A",
                filename=True,
                equate=False,
                is_required=True,
            ),
            _Option(
                ["-b", "contact_map_b"],
                "contact map B",
                filename=True,
                equate=False,
                is_required=True,
            ),
            _Option(["-gap_o", "gap_opening_penalty"], "Gap opening penalty [default=-1]", equate=False),
            _Option(["-gap_e", "gap_extension_penalty"], "Gap extension penalty [default=-0.01]", equate=False),
            _Option(["-sep_cut", "seq_separation_cutoff"], "Sequence separation cutoff [default=3]", equate=False),
            _Option(["-iter", "n_iterations"], "Number of iterations [default=20]", equate=False),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
