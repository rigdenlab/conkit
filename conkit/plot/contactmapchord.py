# BSD 3-Clause License
#
# Copyright (c) 2016-17, University of Liverpool
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
A module to produce a contact map chord diagram
"""

from __future__ import division
from __future__ import print_function

__author__ = "Felix Simkovic"
__date__ = "13 Feb 2017"
__version__ = "0.1"

import matplotlib.pyplot as plt
import numpy as np

from conkit.core.contact import Contact
from conkit.plot._figure import Figure
from conkit.plot._plottools import ColorDefinitions, points_on_circle


class ContactMapChordFigure(Figure):
    """A Figure object specifically for a Contact Map chord diagram

    This figure will illustrate the contacts linking the residues
    in the target sequence. This plot is a very common representation
    of contacts. With this figure, you can illustrate intra-molecular.

    Color scheme: 
    
    ==========  ===========   ==========  ===========  ==========  ===========  ==========  ===========  ==========  ===========
    Amino acid  Hex code      Amino acid  Hex code     Amino acid  Hex code     Amino acid  Hex code     Amino acid  Hex code
    ==========  ===========   ==========  ===========  ==========  ===========  ==========  ===========  ==========  ===========
    Ala         ``#882D17``   Arg         ``#B3446C``  Asn         ``#F99379``  Asp         ``#875692``  Cys         ``#F3C300`` 
    Gln         ``#F6A600``   Glu         ``#F38400``  Gly         ``#BE0032``  His         ``#C2B280``  Ile         ``#848482``
    Leu         ``#E68FAC``   Lys         ``#008856``  Met         ``#0067A5``  Phe         ``#A1CAF1``  Pro         ``#604E97``
    Ser         ``#DCD300``   Thr         ``#8DB600``  Trp         ``#E25822``  Tyr         ``#2B3D26``  Val         ``#654522``
    Unk         ``#000000``
    ==========  ===========   ==========  ===========  ==========  ===========  ==========  ===========  ==========  ===========

    Attributes
    ----------
    hierarchy : :obj:`ContactMap <conkit.core.ContactMap>`
       The default contact map hierarchy

    Examples
    --------
    >>> import conkit
    >>> cmap = conkit.io.read('toxd/toxd.mat', 'ccmpred').top_map
    >>> conkit.plot.ContactMapChordFigure(cmap)

    """

    def __init__(self, hierarchy, use_conf=False, **kwargs):
        """A new contact map plot

        Parameters
        ----------
        hierarchy : :obj:`ContactMap <conkit.core.ContactMap>`
           The default contact map hierarchy
        use_conf : bool, optional
           The marker size will correspond to the raw score [default: False]
        **kwargs
           General :obj:`Figure <conkit.plot._Figure.Figure>` keyword arguments

        """
        super(ContactMapChordFigure, self).__init__(**kwargs)

        self._hierarchy = None
        self.hierarchy = hierarchy
        self.use_conf = use_conf

        self._draw()

    def __repr__(self):
        return "{0}(file_name=\"{1}\")".format(
                self.__class__.__name__, self.file_name
        )

    @property
    def hierarchy(self):
        """The default contact map hierarchy"""
        return self._hierarchy

    @hierarchy.setter
    def hierarchy(self, hierarchy):
        """Define the default contact map hierarchy"""
        if hierarchy:
            Figure._check_hierarchy(hierarchy, "ContactMap")
        if hierarchy.sequence:
            Figure._check_hierarchy(hierarchy.sequence, "Sequence")
        self._hierarchy = hierarchy

    def redraw(self):
        """Re-draw the plot with updated parameters"""
        self._draw()

    def _draw(self):
        """Draw the actual plot"""

        # Re-normalize the data for the lines
        hierarchy = self.hierarchy.rescale()

        # Obtain the data from the hierarchy
        self_data = np.asarray([(c.res1, c.res1_seq, c.res2, c.res2_seq, c.raw_score, c.status) for c in hierarchy])
        _drange = np.append(self_data[:, 1], self_data[:, 3]).astype(np.int64)
        self_data_range = np.arange(_drange.min(), _drange.max() + 1)

        # The number of points on the outer circle and their coordinates
        npoints = self_data_range.shape[0]
        coords = np.asarray(points_on_circle(npoints))

        # Instantiate the figure
        fig, ax = plt.subplots()
        
        # Calculate and plot the Bezier curves
        bezier_path = np.arange(0, 1.01, 0.01)
        for c in self_data:
            x1, y1 = coords[int(c[1]) - self_data_range.min()]
            x2, y2 = coords[int(c[3]) - self_data_range.min()]
            xb, yb = [0, 0]                  # Midpoint the curve is supposed to approach
            x = (1 - bezier_path) ** 2 * x1 + 2 * (1 - bezier_path) * bezier_path * xb + bezier_path ** 2 * x2
            y = (1 - bezier_path) ** 2 * y1 + 2 * (1 - bezier_path) * bezier_path * yb + bezier_path ** 2 * y2
            # 0.0 transparent through 1.0 opaque
            alpha = float(c[4]) if self.use_conf else 1.0
            color = {
                Contact._MISMATCH: ColorDefinitions.MISMATCH,
                Contact._MATCH: ColorDefinitions.MATCH,
            }.get(int(c[5]), ColorDefinitions.MATCH)
            ax.plot(x, y, color=color, alpha=alpha, linestyle="-", zorder=0)
            if int(c[5]) == Contact._MATCH:
                ax.plot(x, y, color=color, alpha=alpha, linestyle="-", zorder=1, linewidth=1)
            else:
                ax.plot(x, y, color=color, alpha=alpha, linestyle="-", zorder=0, linewidth=1)

        # Get the amino acids if available
        #     - get the residue data from the original data array
        residue_data = np.append(self_data[:, [1, 0]], self_data[:, [3, 2]])
        residue_data = residue_data.reshape(self_data[:, 0].shape[0] * 2, 2)

        #     - compute a default color list
        color_codes = dict([(k, ColorDefinitions.AA_ENCODING['X']) for k in self_data_range])

        #     - fill default dict with data we have
        for k, v in np.vstack({tuple(row) for row in residue_data}):
            color_codes[int(k)] = ColorDefinitions.AA_ENCODING[v]

        #     - create a color list
        colors = [color_codes[k] for k in sorted(color_codes.keys())]

        # Plot the residue points
        ax.scatter(coords[:, 0], coords[:, 1], marker='o', color=colors, edgecolors="none", zorder=1)

        # Annotate some residue
        # TODO: Use _plottools module to process this
        label_data = set([int(x) for x in zip(*residue_data)[0]])
        label_coords = np.zeros((npoints, 2))
        space = 2 * np.pi / npoints
        for i in np.arange(npoints):
            label_coords[i] = [
                (npoints + npoints / 10) * np.cos(space * i) - npoints / 20,
                (npoints + npoints / 10) * np.sin(space * i) - npoints / 40 
            ]
        for r in sorted(label_data)[::int(npoints / (npoints / 10))]:
            i = r - self_data_range.min()
            xy = x, y = coords[i]
            xytext = label_coords[i]
            ax.annotate(r, xy=xy, xytext=xytext)
            ax.scatter(x, y, marker='o', facecolors="none", edgecolors="#000000", zorder=2)

        # Arrow for the start
        arrow_x, arrow_y = (npoints + npoints / 5, 0)
        ax.arrow(arrow_x, arrow_y, 0, npoints / 10, head_width=1.5, color="#000000")

        # Prettify the plot
        ax.set_xlim(-arrow_x, arrow_x + 2)
        ax.set_ylim(-arrow_x, arrow_x)
        ax.axis("off")

        # Make both axes identical in length and remove whitespace around the plot
        aspectratio = Figure._correct_aspect(ax, 1.0)
        ax.set(aspect=aspectratio)
        fig.tight_layout()

        fig.savefig(self.file_name, bbox_inches='tight', dpi=self.dpi)
