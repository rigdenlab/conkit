# BSD 3-Clause License
#
# Copyright (c) 2016-18, University of Liverpool
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
__version__ = "1.0"

import matplotlib.collections as mcoll
import matplotlib.pyplot as plt
import numpy as np

from conkit.core.struct import Gap
from conkit.misc import deprecate, normalize
from conkit.plot.figure import Figure
from conkit.plot.tools import ColorDefinitions, _isinstance


class ContactMapFigure(Figure):
    """A Figure object specifically for a :obj:`~conkit.core.contactmap.ContactMap`

    This figure will illustrate the contacts in a contact
    map. This plot is a very common representation of contacts.
    With this figure, you can illustrate either your contact
    map by itself, compared against a second contact map, and/or
    matched against contacts extracted from a contact map.

    Attributes
    ----------
    hierarchy : :obj:`~conkit.core.contactmap.ContactMap`
       The default contact map hierarchy
    other : :obj:`~conkit.core.contactmap.ContactMap`
       The second contact map hierarchy
    reference : :obj:`~conkit.core.contactmap.ContactMap`
       The reference contact map hierarchy
    altloc : bool
       Use the :attr:`~conkit.core.contact.Contact.res_altloc` positions [default: False]
    use_conf : bool
       The marker size will correspond to the raw score [default: False]

    Examples
    --------
    >>> import conkit
    >>> cmap = conkit.io.read('toxd/toxd.mat', 'ccmpred').top_map
    >>> conkit.plot.ContactMapFigure(cmap)

    """

    def __init__(self, hierarchy, other=None, reference=None, altloc=False, use_conf=False, lim=None, **kwargs):
        """A new contact map plot

        Parameters
        ----------
        hierarchy : :obj:`~conkit.core.contactmap.ContactMap`
           The default contact map hierarchy
        other : :obj:`~conkit.core.contactmap.ContactMap`, optional
           The second contact map hierarchy
        reference : :obj:`~conkit.core.contactmap.ContactMap`, optional
           The reference contact map hierarchy
        altloc : bool, optional
           Use the :attr:`~conkit.core.contact.Contact.res_altloc` positions [default: False]
        use_conf : bool, optional
           The marker size will correspond to the raw score [default: False]
        lim : tuple, list, optional
           The [min, max] residue numbers to show
        **kwargs
           General :obj:`~conkit.plot.figure.Figure` keyword arguments

        """
        super(ContactMapFigure, self).__init__(**kwargs)

        self._hierarchy = None
        self._other = None
        self._reference = None
        self._lim = None

        self.altloc = altloc
        self.use_conf = use_conf

        self.hierarchy = hierarchy
        if other:
            self.other = other
        if reference:
            self.reference = reference
        if lim:
            self.lim = lim

        self.draw()

    def __repr__(self):
        return self.__class__.__name__

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
    def reference(self):
        return self._reference

    @reference.setter
    def reference(self, hierarchy):
        if hierarchy and _isinstance(hierarchy, "ContactMap"):
            self._reference = hierarchy
        else:
            raise TypeError("Invalid hierarchy type: %s" % hierarchy.__class__.__name__)

    @property
    def lim(self):
        return self._lim

    @lim.setter
    def lim(self, lim):
        if (isinstance(lim, list) or isinstance(lim, tuple)) and len(lim) == 2:
            self._lim = lim
        elif (isinstance(lim, list) or isinstance(lim, tuple)):
            raise ValueError("A list with 2 entries is required!")
        else:
            raise TypeError("A list with [min, max] limits is required!")

    @deprecate('0.11', msg='Use draw instead')
    def redraw(self):
        self.draw()

    def draw(self):

        if self._reference:
            if self.altloc:
                reference_data = np.array(self._reference.as_list(altloc=True))
            else:
                reference_data = np.array(self._reference.as_list())
            reference_colors = [ColorDefinitions.STRUCTURAL for _ in np.arange(reference_data.shape[0])]
        else:
            reference_data = np.empty((0, 2))
            reference_colors = []

        self_data = np.array([c for c in self._hierarchy.as_list() if all(ci != Gap.IDENTIFIER for ci in c)])
        self_colors = ContactMapFigure._determine_color(self._hierarchy)
        self_rawsc = np.array(
            [c.raw_score for c in self._hierarchy if all(ci != Gap.IDENTIFIER for ci in [c.res1_seq, c.res2_seq])])

        if self._other:
            other_data = np.array([c for c in self._other.as_list() if any(ci != Gap.IDENTIFIER for ci in c)])
            other_colors = ContactMapFigure._determine_color(self._other)
            other_rawsc = np.array(
                [c.raw_score for c in self._other if all(ci != Gap.IDENTIFIER for ci in [c.res1_seq, c.res2_seq])])
        else:
            other_data = self_data
            other_colors = self_colors
            other_rawsc = self_rawsc

        if self.use_conf:
            self_radius = normalize(self_rawsc, vmin=0.2, vmax=0.48)
            other_radius = normalize(other_rawsc, vmin=0.2, vmax=0.48)
        else:
            self_radius = other_radius = 0.48

        self._patch_scatter(
            reference_data[:, 0], reference_data[:, 1], symbol="o", facecolor=reference_colors, radius=0.5, linewidth=0)
        self._patch_scatter(
            reference_data[:, 1], reference_data[:, 0], symbol="o", facecolor=reference_colors, radius=0.5, linewidth=0)
        self._patch_scatter(
            self_data[:, 1], self_data[:, 0], symbol="o", facecolor=self_colors, radius=self_radius, linewidth=0)
        self._patch_scatter(
            other_data[:, 0], other_data[:, 1], symbol="o", facecolor=other_colors, radius=other_radius, linewidth=0)

        if self.lim:
            min_max_data = np.arange(self.lim[0], self.lim[1] + 1)
            self.ax.set_xlim(self.lim[0] - 0.5, self.lim[1] + 0.5)
            self.ax.set_ylim(self.lim[0] - 0.5, self.lim[1] + 0.5)
        else:
            min_max_data = np.append(self_data[:, 0], self_data[:, 1])
            if self._reference:
                min_max_data = np.append(min_max_data, reference_data[:, 0])
                min_max_data = np.append(min_max_data, reference_data[:, 1])
            if self._other:
                min_max_data = np.append(min_max_data, other_data[:, 0])
                min_max_data = np.append(min_max_data, other_data[:, 1])

        self.ax.set_xlim(min_max_data.min() - 0.5, min_max_data.max() + 0.5)
        self.ax.set_ylim(min_max_data.min() - 0.5, min_max_data.max() + 0.5)
        gap = int(10 * (min_max_data.max() - min_max_data.min()) / 100)
        tick_range = np.arange(min_max_data.min(), min_max_data.max(), gap, dtype=np.int64)

        self.ax.set_xticks(tick_range)
        self.ax.set_yticks(tick_range)

        self.ax.set_xlabel('Residue number')
        self.ax.set_ylabel('Residue number')

        line2d_kws = dict(marker="o", linestyle="")
        pseudo_data = [[0, 1], [0, 0]]
        if self._reference:
            tp_artist = plt.Line2D(*pseudo_data, color=ColorDefinitions.MATCH, label='Match', **line2d_kws)
            fp_artist = plt.Line2D(*pseudo_data, color=ColorDefinitions.MISMATCH, label='Mismatch', **line2d_kws)
            rf_artist = plt.Line2D(*pseudo_data, color=ColorDefinitions.STRUCTURAL, label='Structural', **line2d_kws)
            artists = [tp_artist, fp_artist, rf_artist]
        else:
            nt_artist = plt.Line2D(*pseudo_data, color=ColorDefinitions.GENERAL, label='Contact', **line2d_kws)
            artists = [nt_artist]

        if self.legend:
            self.ax.legend(
                handles=artists,
                numpoints=1,
                fontsize=10,
                bbox_to_anchor=(0., 1.02, 1., .102),
                loc=3,
                ncol=3,
                mode="expand",
                borderaxespad=0.)
        # TODO: deprecate this in 0.10
        if self._file_name:
            self.savefig(self._file_name, dpi=self._dpi)

    @staticmethod
    def _determine_color(h):
        """Determine the color of the contacts in order"""

        def color_selector(contact):
            if contact.true_positive:
                return ColorDefinitions.MATCH
            elif contact.false_positive:
                return ColorDefinitions.MISMATCH
            else:
                return ColorDefinitions.GENERAL

        return [color_selector(contact) for contact in h]
