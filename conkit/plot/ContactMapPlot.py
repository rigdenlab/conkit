"""
A module to produce a contact map plot
"""


from __future__ import division

__author__ = "Felix Simkovic"
__date__ = "07 Feb 2017"
__version__ = 0.1

import matplotlib.pyplot
import numpy

from conkit import constants
from conkit.core import ContactMap
from conkit.plot._Figure import Figure


class ContactMapFigure(Figure):

    def __init__(self, hierarchy, other=None, reference=None, altloc=False, use_conf=False, **kwargs):
        """A new contact map plot

        Parameters
        ----------
        hierarchy : :obj:`conkit.core.ContactMap`
           A ConKit :obj:`conkit.core.ContactMap`
        other : :obj:`conkit.core.ContactMap`, optional
           A ConKit :obj:`conkit.core.ContactMap`
        reference : :obj:`conkit.core.ContactMap`, optional
           A ConKit :obj:`conkit.core.ContactMap` [this map refers to the reference contacts]
        altloc : bool, optional
           Use the res_altloc positions [default: False]
        use_conf : bool, optional
           The marker size will correspond to the raw score [default: False]
        **kwargs
           General :obj:`conkit.plot._Figure.Figure` keyword arguments

        Raises
        ------
        RuntimeError
           The hierarchy is not a :obj:`conkit.core.SequenceFile` object

        """
        super(ContactMapFigure, self).__init__(**kwargs)

        if not isinstance(hierarchy, ContactMap):
            raise RuntimeError("Provided hierarchy is not a ContactMap")
        elif other and not isinstance(other, ContactMap):
            raise RuntimeError("Provided other hierarchy is not a ContactMap")
        elif reference and not isinstance(reference, ContactMap):
            raise RuntimeError("Provided reference hierarchy is not a ContactMap")

        self._hierarchy = hierarchy
        self._other = other
        self._reference = reference
        self._altloc = altloc
        self._use_conf = use_conf

        self._draw()

    def draw(self):
        """Draw the actual plot"""

        fig, ax = matplotlib.pyplot.subplots(dpi=self.dpi)

        # Plot the other_ref contacts
        if self._reference:
            if self._altloc:
                reference_data = numpy.asarray([(c.res1_altseq, c.res2_altseq)
                                                for c in self._reference if c.is_true_positive])
            else:
                reference_data = numpy.asarray([(c.res1_seq, c.res2_seq)
                                                for c in self._reference if c.is_true_positive])
            reference_colors = [constants.RFCOLOR for _ in range(len(reference_data))]
            ax.scatter(reference_data.T[0], reference_data.T[1], color=reference_colors,
                       marker='.', edgecolor='none', linewidths=0.0)
            ax.scatter(reference_data.T[1], reference_data.T[0], color=reference_colors,
                       marker='.', edgecolor='none', linewidths=0.0)

        # Plot the self contacts
        self_data = numpy.asarray([(c.res1_seq, c.res2_seq, c.raw_score) for c in self._hierarchy])
        self_colors = [
            constants.TPCOLOR if contact.is_true_positive
            else constants.FPCOLOR if contact.is_false_positive
            else constants.NTCOLOR for contact in self._hierarchy
            ]
        if self._use_conf:
            # TODO: Find a better scaling algorithm
            self_sizes = (self_data.T[2] - self_data.T[2].min()) / (self_data.T[2].max() - self_data.T[2].min()) * 500
        else:
            self_sizes = None

        # This is the bottom triangle
        ax.scatter(self_data.T[1], self_data.T[0], color=self_colors, marker='.',
                   s=self_sizes, edgecolor='none', linewidths=0.0)

        # Plot the other contacts
        if self._other:
            other_data = numpy.asarray([(c.res1_seq, c.res2_seq, c.raw_score) for c in self._other])
            other_colors = [
                constants.TPCOLOR if contact.is_true_positive
                else constants.FPCOLOR if contact.is_false_positive
                else constants.NTCOLOR for contact in self._other
                ]
            if self._use_conf:
                # TODO: Find a better scaling algorithm
                other_sizes = (other_data.T[2] - other_data.T[2].min()) / (other_data.T[2].max() - other_data.T[2].min()) * 500
            else:
                other_sizes = None
            # This is the upper triangle
            ax.scatter(other_data.T[0], other_data.T[1], color=other_colors, marker='.',
                       s=other_sizes, edgecolor='none', linewidths=0.0)
        else:
            # This is the upper triangle
            ax.scatter(self_data.T[0], self_data.T[1], color=self_colors, marker='.',
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
        tick_range = numpy.arange(min_max_data.min(), min_max_data.max(), 10, dtype=numpy.int64)
        ax.set_xticks(tick_range)
        ax.set_yticks(tick_range)

        # Prettify the plot
        ax.set_xlabel('Residue number')
        ax.set_ylabel('Residue number')

        # Create a custom legend
        if self._reference:
            tp_artist = matplotlib.pyplot.Line2D((0, 1), (0, 0), color=constants.TPCOLOR,
                                                 marker='o', linestyle='', label='True positive')
            fp_artist = matplotlib.pyplot.Line2D((0, 1), (0, 0), color=constants.FPCOLOR,
                                                 marker='o', linestyle='', label='False positive')
            rf_artist = matplotlib.pyplot.Line2D((0, 1), (0, 0), color=constants.RFCOLOR,
                                                 marker='o', linestyle='', label='Reference')
            artists = [tp_artist, fp_artist, rf_artist]
        else:
            nt_artist = matplotlib.pyplot.Line2D((0, 1), (0, 0), color=constants.NTCOLOR,
                                                 marker='o', linestyle='', label='Contact')
            artists = [nt_artist]
        ax.legend(handles=artists, numpoints=1, fontsize=10, bbox_to_anchor=(0., 1.02, 1., .102),
                  loc=3, ncol=3, mode="expand", borderaxespad=0.)

        # Make both axes identical in length and remove whitespace around the plot
        ax.set(aspect=1.0)
        fig.tight_layout()

        fig.savefig(self.file_name, bbox_inches='tight')
