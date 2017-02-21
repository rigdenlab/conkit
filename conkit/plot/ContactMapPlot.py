"""
A module to produce a contact map plot
"""


from __future__ import division

__author__ = "Felix Simkovic"
__date__ = "07 Feb 2017"
__version__ = 0.1

import matplotlib.pyplot
import numpy

from conkit.plot._Figure import Figure
from conkit.plot._plottools import ColorDefinitions


class ContactMapFigure(Figure):
    """A Figure object specifically for a Contact Map

    This figure will illustrate the contacts in a contact
    map. This plot is a very common representation of contacts.
    With this figure, you can illustrate either your contact 
    map by itself, compared against a second contact map, and/or
    matched against contacts extracted from a contact map.

    Attributes
    ----------
    hierarchy : :obj:`ContactMap <conkit.core.ContactMap>`
       The default contact map hierarchy
    other : :obj:`ContactMap <conkit.core.ContactMap>`
       The second contact map hierarchy
    reference : :obj:`ContactMap <conkit.core.ContactMap>`
       The reference contact map hierarchy
    altloc : bool
       Use the res_altloc positions [default: False]
    use_conf : bool
       The marker size will correspond to the raw score [default: False]
    
    Examples
    --------
    >>> import conkit
    >>> cmap = conkit.io.read('toxd/toxd.mat', 'ccmpred').top_map
    >>> conkit.plot.ContactMapFigure(cmap)

    """

    def __init__(self, hierarchy, other=None, reference=None, altloc=False, use_conf=False, **kwargs):
        """A new contact map plot

        Parameters
        ----------
        hierarchy : :obj:`ContactMap <conkit.core.ContactMap>`
           The default contact map hierarchy
        other : :obj:`ContactMap <conkit.core.ContactMap>`, optional
           The second contact map hierarchy
        reference : :obj:`ContactMap <conkit.core.ContactMap>`, optional
           The reference contact map hierarchy
        altloc : bool, optional
           Use the res_altloc positions [default: False]
        use_conf : bool, optional
           The marker size will correspond to the raw score [default: False]
        **kwargs
           General :obj:`Figure <conkit.plot._Figure.Figure>` keyword arguments

        """
        super(ContactMapFigure, self).__init__(**kwargs)

        self._hierarchy = None
        self._other = None
        self._reference = None

        self.altloc = altloc
        self.use_conf = use_conf

        self.hierarchy = hierarchy
        self.other = other
        self.reference = reference

        self._draw()

    def __repr__(self):
        return "ContactMapFigure(file_name=\"{0}\")".format(self.file_name)

    @property
    def hierarchy(self):
        """The default contact map hierarchy"""
        return self._hierarchy

    @hierarchy.setter
    def hierarchy(self, hierarchy):
        """Define the default contact map hierarchy"""
        if hierarchy:
            Figure._check_hierarchy(hierarchy, "ContactMap")
        self._hierarchy = hierarchy

    @property
    def other(self):
        """The second contact map hierarchy"""
        return self._other

    @other.setter
    def other(self, hierarchy):
        """Define the default contact map hierarchy"""
        if hierarchy:
            Figure._check_hierarchy(hierarchy, "ContactMap")
        self._other = hierarchy

    @property
    def reference(self):
        """The reference contact map hierarchy"""
        return self._reference

    @reference.setter
    def reference(self, hierarchy):
        """Define the reference contact map hierarchy"""
        if hierarchy:
            Figure._check_hierarchy(hierarchy, "ContactMap")
        self._reference = hierarchy

    def redraw(self):
        """Re-draw the plot with updated parameters"""
        self._draw()

    def _draw(self):
        """Draw the actual plot"""

        fig, ax = matplotlib.pyplot.subplots()

        # Plot the other_ref contacts
        if self._reference:
            if self.altloc:
                reference_data = numpy.asarray([(c.res1_altseq, c.res2_altseq)
                                                for c in self._reference if c.is_true_positive])
            else:
                reference_data = numpy.asarray([(c.res1_seq, c.res2_seq)
                                                for c in self._reference if c.is_true_positive])
            reference_colors = [ColorDefinitions.STRUCTURAL for _ in range(len(reference_data))]
            ax.scatter(reference_data.T[0], reference_data.T[1], color=reference_colors,
                       s=10, marker='o', edgecolor='none', linewidths=0.0)
            ax.scatter(reference_data.T[1], reference_data.T[0], color=reference_colors,
                       s=10, marker='o', edgecolor='none', linewidths=0.0)

        # Plot the self contacts
        self_data = numpy.asarray([(c.res1_seq, c.res2_seq, c.raw_score) for c in self._hierarchy])
        self_colors = ContactMapFigure._determine_color(self._hierarchy)
        if self.use_conf:
            self_sizes = (self_data.T[2] - self_data.T[2].min()) / (self_data.T[2].max() - self_data.T[2].min())
            self_sizes = self_sizes * 20 + 10
        else:
            self_sizes = [10] * len(self_data.T[2])

        # This is the bottom triangle
        ax.scatter(self_data.T[1], self_data.T[0], color=self_colors, marker='o',
                   s=self_sizes, edgecolor='none', linewidths=0.0)

        # Plot the other contacts
        if self._other:
            other_data = numpy.asarray([(c.res1_seq, c.res2_seq, c.raw_score) for c in self._other])
            other_colors = ContactMapFigure._determine_color(self._other)
            if self.use_conf:
                other_sizes = (other_data.T[2] - other_data.T[2].min()) / (other_data.T[2].max() - other_data.T[2].min())
                other_sizes = other_sizes * 20 + 10
            else:
                other_sizes = [10] * len(other_data.T[2])
            # This is the upper triangle
            ax.scatter(other_data.T[0], other_data.T[1], color=other_colors, marker='o',
                       s=other_sizes, edgecolor='none', linewidths=0.0)
        else:
            # This is the upper triangle
            ax.scatter(self_data.T[0], self_data.T[1], color=self_colors, marker='o',
                       s=self_sizes, edgecolor='none', linewidths=0.0)

        # Allow dynamic x and y limits
        min_max_data = numpy.append(self_data.T[0], self_data.T[1])
        if self._reference:
            min_max_data = numpy.append(min_max_data, reference_data.T[0])
            min_max_data = numpy.append(min_max_data, reference_data.T[1])
        if self._other:
            min_max_data = numpy.append(min_max_data, other_data.T[0])
            min_max_data = numpy.append(min_max_data, other_data.T[1])
        ax.set_xlim(min_max_data.min() - 0.5, min_max_data.max() + 0.5)
        ax.set_ylim(min_max_data.min() - 0.5, min_max_data.max() + 0.5)

        # Set the xticks and yticks dynamically
        gap = int(10 * (min_max_data.max() - min_max_data.min()) / 100)
        tick_range = numpy.arange(min_max_data.min(), min_max_data.max(), gap, dtype=numpy.int64)
        ax.set_xticks(tick_range)
        ax.set_yticks(tick_range)

        # Prettify the plot
        ax.set_xlabel('Residue number')
        ax.set_ylabel('Residue number')

        # Create a custom legend
        if self._reference:
            tp_artist = matplotlib.pyplot.Line2D((0, 1), (0, 0), color=ColorDefinitions.MATCH,
                                                 marker='o', linestyle='', label='Match')
            fp_artist = matplotlib.pyplot.Line2D((0, 1), (0, 0), color=ColorDefinitions.MISMATCH,
                                                 marker='o', linestyle='', label='Mismatch')
            rf_artist = matplotlib.pyplot.Line2D((0, 1), (0, 0), color=ColorDefinitions.STRUCTURAL,
                                                 marker='o', linestyle='', label='Structural')
            artists = [tp_artist, fp_artist, rf_artist]
        else:
            nt_artist = matplotlib.pyplot.Line2D((0, 1), (0, 0), color=ColorDefinitions.GENERAL,
                                                 marker='o', linestyle='', label='Contact')
            artists = [nt_artist]
        ax.legend(handles=artists, numpoints=1, fontsize=10, bbox_to_anchor=(0., 1.02, 1., .102),
                  loc=3, ncol=3, mode="expand", borderaxespad=0.)

        # Make both axes identical in length and remove whitespace around the plot
        aspectratio = Figure._correct_aspect(ax, 1.0)
        ax.set(aspect=aspectratio)
        fig.tight_layout()

        fig.savefig(self.file_name, bbox_inches='tight', dpi=self.dpi)

    @staticmethod
    def _determine_color(h):
        """Determine the color of the contacts in order"""
        return [
            ColorDefinitions.MATCH if contact.is_true_positive
            else ColorDefinitions.MISMATCH if contact.is_false_positive
            else ColorDefinitions.GENERAL for contact in h
        ]
