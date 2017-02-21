"""
A module to produce a sequence coverage plot
"""

from __future__ import division

__author__ = "Felix Simkovic"
__date__ = "07 Feb 2017"
__version__ = 0.1

import matplotlib.pyplot
import numpy

from conkit.plot._Figure import Figure
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
        return "SequenceCoverageFigure(file_name=\"{0}\")".format(self.file_name)

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

        residues = numpy.arange(1, self._hierarchy.top_sequence.seq_len + 1)
        aa_counts = numpy.asarray(self._hierarchy.calculate_freq()) * self._hierarchy.nseqs

        fig, ax = matplotlib.pyplot.subplots()

        # Add lines as quality indicators
        ax.axhline(self._hierarchy.top_sequence.seq_len * 5, color=ColorDefinitions.L5CUTOFF, label='5 x Nresidues')
        if any(x >= self._hierarchy.top_sequence.seq_len * 20 for x in aa_counts):
            ax.axhline(self._hierarchy.top_sequence.seq_len * 20, color=ColorDefinitions.L20CUTOFF, label='20 x Nresidues')

        ax.plot(residues, aa_counts, color=ColorDefinitions.GENERAL, marker='o', linestyle='-',
                markersize=5, label='Amino acid count')

        # Prettify the plot
        ax.set_xlim(residues[0], residues[-1])
        xticks = ax.get_xticks().astype(numpy.int64) + residues[0]
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
