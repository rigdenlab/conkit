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
"""A module to produce a distogram heatmap plot"""

import numpy as np

from conkit.core.struct import Gap
from conkit.plot.contactmapmatrix import ContactMapMatrixFigure
from conkit.plot.tools import _isinstance


class DistogramHeatmapFigure(ContactMapMatrixFigure):
    """A Figure object specifically for a :obj:`~conkit.core.distogram.Distogram`

    This figure will illustrate the predicted residues distances in a
    distogram. If two instances of :obj:`~conkit.core.distogram.Distogram`
    are used, then each half square of the heatmap corresponds with the distances
    at each instance

    Attributes
    ----------
    hierarchy : :obj:`~conkit.core.distogram.Distogram`
       The default distogram hierarchy
    other : :obj:`~conkit.core.distogram.Distogram`
       The second distogram hierarchy
    altloc : bool
       Use the :attr:`~conkit.core.distance.Distance.res_altloc` positions [default: False]

    Examples
    --------
    >>> import conkit
    >>> distogram = conkit.io.read('toxd/toxd.npz', 'rosettanpz').top_map
    >>> sequence  = conkit.io.read('toxd/toxd.fasta', 'fasta').top
    >>> distogram.sequence = sequence
    >>> conkit.plot.DistogramHeatmapFigure(distogram)

    """

    def __init__(self, *args, **kwargs):
        """A new distogram heatmap plot

        Parameters
        ----------
        hierarchy : :obj:`~conkit.core.distogram.Distogram`
           The default distogram hierarchy
        other : :obj:`~conkit.core.distogram.Distogram`, optional
           The second distogram hierarchy
        altloc : bool, optional
           Use the :attr:`~conkit.core.distance.Distance.res_altloc` positions [default: False]
        lim : tuple, list, optional
           The [min, max] residue numbers to show
        **kwargs
           General :obj:`~conkit.plot.figure.Figure` keyword arguments

        """
        super(DistogramHeatmapFigure, self).__init__(*args, **kwargs)

    @property
    def hierarchy(self):
        return self._hierarchy

    @hierarchy.setter
    def hierarchy(self, hierarchy):
        if hierarchy and _isinstance(hierarchy, "Distogram"):
            self._hierarchy = hierarchy
        else:
            raise TypeError("Invalid hierarchy type: %s" % hierarchy.__class__.__name__)

    @property
    def other(self):
        return self._other

    @other.setter
    def other(self, hierarchy):
        if hierarchy and _isinstance(hierarchy, "Distogram"):
            self._other = hierarchy
        else:
            raise TypeError("Invalid hierarchy type: %s" % hierarchy.__class__.__name__)

    def draw(self):
        if self._other:
            array = self._hierarchy.merge_arrays(self._hierarchy, self._other)
            self_data = np.array([c for c in self._hierarchy.as_list() if all(ci != Gap.IDENTIFIER for ci in c)])
            other_data = np.array([c for c in self._other.as_list() if all(ci != Gap.IDENTIFIER for ci in c)])
        else:
            array = self._hierarchy.as_array()
            self_data = np.array([c for c in self._hierarchy.as_list() if all(ci != Gap.IDENTIFIER for ci in c)])
            other_data = self_data

        im = self.ax.imshow(array, cmap='hot')
        self.fig.colorbar(im, orientation='vertical')
        self.define_axis_settings(self_data, other_data)

        # TODO: deprecate this in 0.10
        if self._file_name:
            self.savefig(self._file_name, dpi=self._dpi)
