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
A module to produce a sequence coverage plot
"""

from __future__ import division
from __future__ import print_function

__author__ = "Felix Simkovic"
__date__ = "07 Feb 2017"
__version__ = "0.1"

import matplotlib.pyplot as plt
import numpy as np

from conkit.plot._figure import Figure
from conkit.plot._plottools import ColorDefinitions


class SequenceCoverageFigure(Figure):
    """A Figure object specifically for a Sequence coverage illustration.

    This figure will illustrate the coverage of sequences at each position
    in the provided alignment. It counts the frequency of sequences at which
    a residue is present and plots it.

    This figure can be particularly useful in cases where domain boundaries
    could be redefined. It is also useful in cases where only parts of your
    alignment are well covered, and thus trimming the alignment might produce
    much more accurate Evoluationary Covariance predictions.

    Attributes
    ----------
    hierarchy : :obj:`SequenceFile <conkit.core.SequenceFile>`
       The Multiple Sequence Alignment hierarchy

    Examples
    --------
    >>> import conkit
    >>> msa = conkit.io.read('toxd/toxd.a3m', 'a3m')
    >>> conkit.plot.SequenceCoverageFigure(msa)
    
    """

    def __init__(self, hierarchy, **kwargs):
        """A new sequence coverage plot

        Parameters
        ----------
        hierarchy : :obj:`SequenceFile <conkit.core.SequenceFile>`
           The Multiple Sequence Alignment hierarchy
        **kwargs
           General :obj:`Figure <conkit.plot._Figure.Figure>` keyword arguments

        """
        super(SequenceCoverageFigure, self).__init__(**kwargs)
        self._hierarchy = None

        self.hierarchy = hierarchy

        self._draw()

    def __repr__(self):
        return "{0}(file_name=\"{1}\")".format(
                self.__class__.__name__, self.file_name
        )

    @property
    def hierarchy(self):
        """A ConKit :obj:`SequenceFile <conkit.core.SequenceFile>`"""
        return self._hierarchy

    @hierarchy.setter
    def hierarchy(self, hierarchy):
        """Define the ConKit :obj:`SequenceFile <conkit.core.SequenceFile>`

        Raises
        ------
        RuntimeError
           The hierarchy is not an alignment

        """
        if hierarchy:
            Figure._check_hierarchy(hierarchy, "SequenceFile")
            if not hierarchy.is_alignment:
                raise RuntimeError("Provided hierarchy does not show characteristics of an alignment")
        self._hierarchy = hierarchy

    def redraw(self):
        """Re-draw the plot with updated parameters"""
        self._draw()

    def _draw(self):
        """Draw the actual plot"""

        residues = np.arange(1, self._hierarchy.top_sequence.seq_len + 1)
        aa_counts = np.asarray(self._hierarchy.calculate_freq()) * self._hierarchy.nseqs

        fig, ax = plt.subplots()

        ax.plot(residues, aa_counts, color=ColorDefinitions.GENERAL, marker='o', markersize=5, linestyle='-',
                label='Amino acid count', zorder=1)

        # Add lines as quality indicators
        ax.axhline(self._hierarchy.top_sequence.seq_len * 5, color=ColorDefinitions.L5CUTOFF,
                   label='5 x Nresidues', zorder=0)
        if any(x >= self._hierarchy.top_sequence.seq_len * 20 for x in aa_counts):
            ax.axhline(self._hierarchy.top_sequence.seq_len * 20, color=ColorDefinitions.L20CUTOFF,
                       label='20 x Nresidues', zorder=0)


        # Prettify the plot
        ax.set_xlim(residues[0], residues[-1])
        xticks = ax.get_xticks().astype(np.int64) + residues[0]
        # Remove any excess xticks
        xticks = np.delete(xticks, [i for i, t in enumerate(xticks) if t > residues[-1]])
        ax.set_xticks(xticks)
        ax.set_xticklabels(xticks)

        ax.set_xlabel('Residue number')
        ax.set_ylabel('Sequence Count')
        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)

        # Make axes length proportional and remove whitespace around the plot
        aspectratio = Figure._correct_aspect(ax, 0.3)
        ax.set(aspect=aspectratio)
        fig.tight_layout()

        fig.savefig(self.file_name, bbox_inches='tight', dpi=self.dpi)
