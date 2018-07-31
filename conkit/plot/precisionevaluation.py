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
"""A module to produce a precision evaluation plot"""

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


class PrecisionEvaluationFigure(Figure):
    """A Figure object specifically for a Precision evaluation.

    This figure will illustrate the precision scores of a contact
    map at different precision scores. These can be determined at
    various start and end points with different stepwise increases
    in between. 

    Attributes
    ----------
    hierarchy : :obj:`~conkit.core.contactmap.ContactMap`
       The contact map hierarchy
    cutoff_step : float
       The cutoff step
    min_cutoff : float
       The minimum cutoff factor
    max_cutoff : float
       The maximum cutoff factor

    Examples
    --------
    >>> import conkit
    >>> cmap = conkit.io.read('toxd/toxd.mat', 'ccmpred').top_map
    >>> cmap.sequence = conkit.io.read('toxd/toxd.fasta', 'fasta').top_sequence
    >>> pdb = conkit.io.read('toxd/toxd.pdb', 'pdb').top_map
    >>> cmap.match(pdb, inplace=True)
    >>> conkit.plot.PrecisionEvaluationFigure(cmap)

    """

    def __init__(self, hierarchy, min_cutoff=0.0, max_cutoff=100.0, cutoff_step=0.2, **kwargs):
        """A precision evaluation figure

        Parameters
        ----------
        hierarchy : :obj:`~conkit.core.contactmap.ContactMap`
           The contact map hierarchy
        min_cutoff : float, optional
           The minimum factor
        max_cutoff : float, optional
           The maximum facotr
        cutoff_step : float, optional
           The cutoff step
        **kwargs
           General :obj:`~conkit.plot.figure.Figure` keyword arguments

        """
        super(PrecisionEvaluationFigure, self).__init__(**kwargs)
        self._hierarchy = None
        self._cutoff_boundaries = [0.0, 100.0]
        self._cutoff_step = 0.2

        self.hierarchy = hierarchy
        self.cutoff_step = cutoff_step
        self.min_cutoff = min_cutoff
        self.max_cutoff = max_cutoff

        self.draw()

    def __repr__(self):
        return self.__class__.__name__

    @property
    def cutoff_step(self):
        """The cutoff step"""
        return self._cutoff_step

    @cutoff_step.setter
    def cutoff_step(self, cutoff_step):
        """Define the cutoff step"""
        self._cutoff_step = cutoff_step

    @property
    def hierarchy(self):
        """A ConKit :obj:`~conkit.core.contactmap.ContactMap`"""
        return self._hierarchy

    @hierarchy.setter
    def hierarchy(self, hierarchy):
        """Define the ConKit :obj:`~conkit.core.contactmap.ContactMap`

        Raises
        ------
        :exc:`RuntimeError`
           The hierarchy is not an alignment

        """
        if hierarchy and _isinstance(hierarchy, "ContactMap"):
            self._hierarchy = hierarchy
        else:
            raise TypeError("Invalid hierarchy type: %s" % hierarchy.__class__.__name__)

    @property
    def min_cutoff(self):
        """The minimum cutoff factor

        Raises
        ------
        :obj:`ValueError`
           The minimum cutoff value is larger than or equal to the maximum

        """
        if self._cutoff_boundaries[0] >= self._cutoff_boundaries[1]:
            msg = "The minimum cutoff value is larger than or equal to the maximum"
            raise ValueError(msg)
        return self._cutoff_boundaries[0]

    @min_cutoff.setter
    def min_cutoff(self, min_cutoff):
        """Define the minimum cutoff factor"""
        if min_cutoff < 0.0:
            raise ValueError("Minimum factor cannot be negative")
        self._cutoff_boundaries[0] = min_cutoff

    @property
    def max_cutoff(self):
        """The maximum cutoff factor

        Raises
        ------
        :obj:`ValueError`
           The maximum cutoff value is smaller than the the minimum

        """
        if self._cutoff_boundaries[1] < self._cutoff_boundaries[0]:
            msg = "The maximum cutoff value is smaller than the the minimum"
            raise ValueError(msg)
        return self._cutoff_boundaries[1]

    @max_cutoff.setter
    def max_cutoff(self, max_cutoff):
        """Define the maximum cutoff factor"""
        if max_cutoff > 100:
            raise ValueError("Maximum factor cannot be greater than 100")
        self._cutoff_boundaries[1] = max_cutoff

    @deprecate('0.11', msg='Use draw instead')
    def redraw(self):
        self.draw()

    def draw(self):
        factors = np.arange(self.min_cutoff, self.max_cutoff + 0.1, self.cutoff_step)
        precisions = np.zeros(factors.shape[0])
        for i, factor in enumerate(factors):
            ncontacts = int(self._hierarchy.sequence.seq_len * factor)
            m = self._hierarchy[:ncontacts]
            precisions[i] = m.precision

        self.ax.plot(
            factors,
            precisions,
            color=ColorDefinitions.GENERAL,
            marker=None,
            linestyle='-',
            label='Precision score',
            zorder=1)

        self.ax.axhline(0.5, color=ColorDefinitions.PRECISION50, linestyle='-', label='50% Precision', zorder=0)
        if self.min_cutoff <= 1.0:
            self.ax.axvline(1.0, color=ColorDefinitions.FACTOR1, linestyle='-', label='Factor L', zorder=0)
        if self.min_cutoff <= 0.5:
            self.ax.axvline(0.5, color=ColorDefinitions.FACTOR1, linestyle='--', label='Factor L/2', zorder=0)
        if self.min_cutoff <= 0.2:
            self.ax.axvline(0.2, color=ColorDefinitions.FACTOR1, linestyle='-.', label='Factor L/5', zorder=0)
        if self.min_cutoff <= 0.1:
            self.ax.axvline(0.1, color=ColorDefinitions.FACTOR1, linestyle=':', label='Factor L/10', zorder=0)

        self.ax.set_xlim(self.min_cutoff, self.max_cutoff)
        xticks = (self.ax.get_xticks() * self._hierarchy.sequence.seq_len).astype(np.int64)
        self.ax.set_xticklabels(xticks)
        self.ax.set_ylim(0.0, 1.0)

        self.ax.set_xlabel('Number of Contacts')
        self.ax.set_ylabel('Precision')

        if self.legend:
            self.ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)
        # TODO: deprecate this in 0.10
        if self._file_name:
            self.savefig(self._file_name, dpi=self._dpi)
