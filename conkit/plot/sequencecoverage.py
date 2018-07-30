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
"""A module to produce a sequence coverage plot"""

from __future__ import division
from __future__ import print_function

__author__ = "Felix Simkovic"
__date__ = "07 Feb 2017"
__version__ = "0.1"

import matplotlib.pyplot as plt
import numpy as np

from conkit.misc import deprecate
from conkit.plot.figure import Figure
from conkit.plot.tools import ColorDefinitions, _isinstance


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
    hierarchy : :obj:`~conkit.core.sequencefile.SequenceFile`
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
        hierarchy : :obj:`~conkit.core.sequencefile.SequenceFile`
           The Multiple Sequence Alignment hierarchy
        **kwargs
           General :obj:`~conkit.plot.figure.Figure` keyword arguments

        """
        super(SequenceCoverageFigure, self).__init__(**kwargs)
        self._hierarchy = None

        self.hierarchy = hierarchy

        self.draw()

    def __repr__(self):
        return self.__class__.__name__

    @property
    def hierarchy(self):
        """A ConKit :obj:`~conkit.core.sequencefile.SequenceFile`"""
        return self._hierarchy

    @hierarchy.setter
    def hierarchy(self, hierarchy):
        """Define the ConKit :obj:`~conkit.core.sequencefile.SequenceFile` """
        if hierarchy and _isinstance(hierarchy, "SequenceFile") and hierarchy.is_alignment:
            self._hierarchy = hierarchy
        elif hierarchy and _isinstance(hierarchy, "SequenceFile"):
            raise TypeError("Provided hierarchy does not show characteristics of an alignment")
        else:
            raise TypeError("Invalid hierarchy type: %s" % hierarchy.__class__.__name__)

    @deprecate('0.11', msg='Use draw instead')
    def redraw(self):
        self.draw()

    def draw(self):
        residues = np.arange(1, self._hierarchy.top_sequence.seq_len + 1)
        aa_counts = (np.asarray(self._hierarchy.get_frequency("X")) - self._hierarchy.nseq) * (-1)

        self.ax.plot(
            residues,
            aa_counts,
            color=ColorDefinitions.GENERAL,
            marker=None,
            linestyle='-',
            label='Amino acid count',
            zorder=1)

        self.ax.axhline(
            self._hierarchy.top_sequence.seq_len * 5, color=ColorDefinitions.L5CUTOFF, label='5 x Nresidues', zorder=0)
        if np.any(aa_counts >= self._hierarchy.top_sequence.seq_len * 20):
            self.ax.axhline(
                self._hierarchy.top_sequence.seq_len * 20,
                color=ColorDefinitions.L20CUTOFF,
                label='20 x Nresidues',
                zorder=0)

        self.ax.set_xlim(residues[0], residues[-1])
        xticks = self.ax.get_xticks().astype(np.int64) + residues[0]
        xticks = np.delete(xticks, [i for i, t in enumerate(xticks) if t > residues[-1]])
        self.ax.set_xticks(xticks)
        self.ax.set_xticklabels(xticks)

        self.ax.set_xlabel('Residue number')
        self.ax.set_ylabel('Sequence Count')

        if self.legend:
            self.ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)
        # TODO: deprecate this in 0.10
        if self._file_name:
            self.savefig(self._file_name, dpi=self._dpi)
