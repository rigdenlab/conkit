"""
A module to produce a precision evaluation plot
"""

from __future__ import division

__author__ = "Felix Simkovic"
__date__ = "07 Feb 2017"
__version__ = 0.1

import matplotlib.pyplot
import numpy

from conkit.plot._Figure import Figure
from conkit.plot._plottools import ColorDefinitions


class PrecisionEvaluationFigure(Figure):
    """A Figure object specifically for a Precision evaluation.

    This figure will illustrate the precision scores of a contact
    map at different precision scores. These can be determined at
    various start and end points with different stepwise increases
    in between. 

    Attributes
    ----------
    hierarchy : :obj:`conkit.core.ContactMap`
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
        hierarchy : :obj:`conkit.core.ContactMap`
           The contact map hierarchy
        min_cutoff : float, optional
           The minimum factor
        max_cutoff : float, optional
           The maximum facotr
        cutoff_step : float, optional
           The cutoff step
        **kwargs
           General :obj:`Figure <conkit.plot._Figure.Figure>` keyword arguments

        """
        super(PrecisionEvaluationFigure, self).__init__(**kwargs)
        self._hierarchy = None
        self._cutoff_boundaries = [0.0, 100.0]
        self._cutoff_step = 0.2

        self.hierarchy = hierarchy
        self.cutoff_step = cutoff_step
        self.min_cutoff = min_cutoff
        self.max_cutoff = max_cutoff

        self._draw()

    def __repr__(self):
        return "PrecisionEvaluationFigure(file_name=\"{0}\")".format(self.file_name)

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
        """A ConKit :obj:`conkit.core.ContactMap`"""
        return self._hierarchy

    @hierarchy.setter
    def hierarchy(self, hierarchy):
        """Define the ConKit :obj:`conkit.core.ContactMap`

        Raises
        ------
        RuntimeError
           The hierarchy is not an alignment

        """
        if hierarchy:
            Figure._check_hierarchy(hierarchy, "ContactMap")
        self._hierarchy = hierarchy

    @property
    def min_cutoff(self):
        """The minimum cutoff factor

        Raises
        ------
        ValueError
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
        ValueError
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

    def redraw(self):
        """Re-draw the plot with updated parameters"""
        self._draw()

    def _draw(self):
        """Draw the actual plot"""

        factors = numpy.arange(self.min_cutoff, self.max_cutoff + 0.1, self.cutoff_step)
        precisions = numpy.zeros(factors.shape[0])
        for i, factor in enumerate(factors):
            ncontacts = int(self._hierarchy.sequence.seq_len * factor)
            m = self._hierarchy[:ncontacts]
            precisions[i] = m.precision

        fig, ax = matplotlib.pyplot.subplots()

        # Add indicator lines for clarity of data
        ax.axhline(0.5, color=ColorDefinitions.PRECISION50, linestyle='-',  label='50% Precision')
        if self.min_cutoff <= 1.0:
            ax.axvline(1.0, color=ColorDefinitions.FACTOR1, linestyle='--', label='Factor L')

        # Add data points itself
        ax.plot(factors, precisions, color=ColorDefinitions.GENERAL, marker='o', markersize=5, linestyle='-',
                label='Precision score')

        # Prettify the plot
        ax.set_xlim(self.min_cutoff, self.max_cutoff)
        xticks = (ax.get_xticks() * self._hierarchy.sequence.seq_len).astype(numpy.int64)
        ax.set_xticklabels(xticks)
        ax.set_ylim(0.0, 1.0)

        ax.set_xlabel('Number of Contacts')
        ax.set_ylabel('Precision')
        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)

        # Make axes length proportional and remove whitespace around the plot
        aspectratio = Figure._correct_aspect(ax, 0.3)
        ax.set(aspect=aspectratio)
        fig.tight_layout()

        fig.savefig(self.file_name, bbox_inches='tight', dpi=self.dpi)
