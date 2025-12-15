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
"""A module to produce a contact map plot"""

from __future__ import division
from __future__ import print_function

__author__ = "Felix Simkovic"
__date__ = "10 Jan 2018"
__version__ = "0.13.3"

import matplotlib.collections as mcoll
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np

from conkit.core.struct import Gap
from conkit.misc import normalize
from conkit.plot.figure import Figure
from conkit.plot.tools import ColorDefinitions, _isinstance


class ContactMapMatrixFigure(Figure):
    """A Figure object specifically for a :obj:`~conkit.core.contactmap.ContactMap`

    This figure will illustrate the contacts in a contact
    map matrix. This plot is a very common representation of contacts.
    With this figure, you can illustrate either your contact
    map by itself, compared against a second contact map, and/or
    matched against contacts extracted from a contact map.

    Attributes
    ----------
    hierarchy : :obj:`~conkit.core.contactmap.ContactMap`
       The default contact map hierarchy
    other : :obj:`~conkit.core.contactmap.ContactMap`
       The second contact map hierarchy
    altloc : bool
       Use the :attr:`~conkit.core.contact.Contact.res_altloc` positions [default: False]

    Examples
    --------
    >>> import conkit
    >>> cmap = conkit.io.read('toxd/toxd.mat', 'ccmpred').top_map
    >>> conkit.plot.ContactMapMatrixFigure(cmap)

    """

    def __init__(self, hierarchy, overlay=None, other=None, altloc=False, lim=None, cmap="Greys", **kwargs):
        """A new contact map plot

        Parameters
        ----------
        hierarchy : :obj:`~conkit.core.contactmap.ContactMap`
           The default contact map hierarchy
        other : :obj:`~conkit.core.contactmap.ContactMap`, optional
           The second contact map hierarchy
        altloc : bool, optional
           Use the :attr:`~conkit.core.contact.Contact.res_altloc` positions [default: False]
        lim : tuple, list, optional
           The [min, max] residue numbers to show
        **kwargs
           General :obj:`~conkit.plot.figure.Figure` keyword arguments

        """
        super(ContactMapMatrixFigure, self).__init__(**kwargs)

        self._overlay = None
        self._hierarchy = None                  
        self._other = None
        self._lim = None

        self.altloc = altloc

        self.hierarchy = hierarchy
        if other:
            self.other = other
        if lim:
            self.lim = lim
        if overlay:
            self.overlay = overlay

        self.draw(cmap=cmap)

    def __repr__(self):
        return self.__class__.__name__

    @property
    def overlay(self):
        return self._overlay

    @overlay.setter
    def overlay(self, overlay):
        if overlay and _isinstance(overlay, "ContactMap"):
            self._overlay = overlay
        else:
            raise TypeError("Invalid overlay type: %s" % overlay.__class__.__name__)

    @property
    def hierarchy(self):
        return self._hierarchy

    @hierarchy.setter
    def hierarchy(self, hierarchy):
        if hierarchy and _isinstance(hierarchy, "ContactMap"):
            self._hierarchy = hierarchy
        else:
            raise TypeError("Invalid hierarchy type: %s" % hierarchy.__class__.__name__)

    @property
    def other(self):
        return self._other

    @other.setter
    def other(self, hierarchy):
        if hierarchy and _isinstance(hierarchy, "ContactMap"):
            self._other = hierarchy
        else:
            raise TypeError("Invalid hierarchy type: %s" % hierarchy.__class__.__name__)

    @property
    def lim(self):
        return self._lim

    @lim.setter
    def lim(self, lim):
        if isinstance(lim, (list, tuple)) and len(lim) == 2:
            self._lim = lim
        elif isinstance(lim, (list, tuple)):
            raise ValueError("A list with 2 entries is required!")
        else:
            raise TypeError("A list with [min, max] limits is required!")

    def draw(self,cmap="Greys"):
        _hierarchy = self._hierarchy.rescale()

        self_data = np.array([c for c in _hierarchy.as_list() if all(ci != Gap.IDENTIFIER for ci in c)])
        self_colors = ContactMapMatrixFigure._determine_color(_hierarchy,cmap=cmap)
        self_rawsc = np.array(
            [c.raw_score for c in _hierarchy if all(ci != Gap.IDENTIFIER for ci in [c.res1_seq, c.res2_seq])]
        )

        if self._other:
            _other = self._other.rescale()
            other_data = np.array([c for c in _other.as_list() if any(ci != Gap.IDENTIFIER for ci in c)])
            other_colors = ContactMapMatrixFigure._determine_color(_other,cmap=cmap)
            other_rawsc = np.array(
                [c.raw_score for c in _other if all(ci != Gap.IDENTIFIER for ci in [c.res1_seq, c.res2_seq])]
            )
        else:
            other_data = self_data
            other_colors = self_colors
            other_rawsc = self_rawsc

        self._patch_scatter(
            self_data[:, 1], self_data[:, 0], symbol="s", facecolor=self_colors, radius=1.0, linewidth=0
        )
        self._patch_scatter(
            other_data[:, 0], other_data[:, 1], symbol="s", facecolor=other_colors, radius=1.0, linewidth=0
        )

        if self._overlay:
            _overlay = self._overlay.rescale()
            overlay_data = np.array([c for c in _overlay.as_list() if any(ci != Gap.IDENTIFIER for ci in c)])
            self._patch_scatter( overlay_data[:, 0]+0.5, overlay_data[:, 1]+0.5, symbol="o", facecolor='black', edgecolor='white', radius=0.25, linewidth=0.2)
            self._patch_scatter( overlay_data[:, 1]+0.5, overlay_data[:, 0]+0.5, symbol="o", facecolor='black', edgecolor='white', radius=0.25, linewidth=0.2)

        self.fig.colorbar(cm.ScalarMappable(cmap=cmap), orientation='vertical',ax=self.ax, label='Contact occurrence rate')
        self.define_axis_settings(self_data, other_data)

        # TODO: deprecate this in 0.14
        if self._file_name:
            self.savefig(self._file_name, dpi=self._dpi)

    def define_axis_settings(self, data_1, data_2):
        if self.lim:
            min_max_data = np.arange(self.lim[0], self.lim[1] + 1)
            self.ax.set_xlim(self.lim[0] - 0.5, self.lim[1] + 0.5)
            self.ax.set_ylim(self.lim[0] - 0.5, self.lim[1] + 0.5)
        else:
            min_max_data = np.append(data_1[:, 0], data_1[:, 1])
            min_max_data = np.append(min_max_data, data_2[:, 0])
            min_max_data = np.append(min_max_data, data_2[:, 1])

        self.ax.set_xlim(min_max_data.min(), min_max_data.max() + 1.0)
        self.ax.set_ylim(min_max_data.min(), min_max_data.max() + 1.0)
        gap = 10 * (min_max_data.max() - min_max_data.min()) // 100
        if gap < 1:
            gap = 1
        tick_range = np.arange(min_max_data.min(), min_max_data.max(), gap, dtype=np.int64)

        self.ax.set_xticks(tick_range + 0.5)
        self.ax.set_xticklabels(tick_range)
        self.ax.set_yticks(tick_range + 0.5)
        self.ax.set_yticklabels(tick_range)

        self.ax.set_xlabel("Residue number")
        self.ax.set_ylabel("Residue number")

    @staticmethod
    def _determine_color(h,cmap="Greys"):
        """Determine the color of the contacts in order"""
        cmap = plt.get_cmap(cmap)
        return [cmap(contact.raw_score) for contact in h]
