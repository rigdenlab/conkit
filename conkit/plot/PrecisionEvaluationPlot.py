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


class PrecisionEvaluationFigure(Figure):
    """A Figure object specifically for a Precision evaluation.

    Description
    -----------
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

    def __init__(self, hierarchy, cutoff_step=0.2, **kwargs):
        """A precision evaluation figure

        Parameters
        ----------
        hierarchy : :obj:`conkit.core.ContactMap`
           The contact map hierarchy
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

        ax.axhline(0.5, color='g', label='50% Precision')

        ax.plot(factors, precisions, color='#000000', marker='o', linestyle='-',
                markersize=2, label='Precision score')

        # Prettify the plot
        step = int(factors.shape[0] / 6)
        xticklabels = (factors * self._hierarchy.sequence.seq_len).astype(dtype=numpy.int64)
        ax.set_xticks(factors[::step])
        ax.set_xticklabels(xticklabels[::step])

        yticks = numpy.arange(0, 1.01, 0.2)
        ax.set_yticks(yticks)
        ax.set_yticklabels(yticks)

        ax.set_xlabel('Number of Contacts')
        ax.set_ylabel('Precision')
        ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)

        # Make axes length proportional and remove whitespace around the plot
        aspectratio = Figure._correct_aspect(ax, 0.3)
        ax.set(aspect=aspectratio)
        fig.tight_layout()

        fig.savefig(self.file_name, bbox_inches='tight', dpi=self.dpi)
