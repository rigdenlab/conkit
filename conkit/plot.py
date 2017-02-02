"""
Plot functions collected in single space
"""

from __future__ import division
from __future__ import print_function

__author__ = "Felix Simkovic"
__date__ = "02 Feb 2017"
__version__ = 0.1

from conkit import constants
from conkit.core import ContactMap
from conkit.core import SequenceFile

import numpy

try:
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot
    MATPLOTLIB = True
except ImportError:
    MATPLOTLIB = False


def contact_map(hierarchy, other=None, reference=None, altloc=False, dpi=300,
                file_name='contactmap.png', file_format='png'):
    """Produce a 2D contact map plot

    Parameters
    ----------
    hierarchy : :obj:`conkit.core.ContactMap`
       A ConKit :obj:`conkit.core.ContactMap`
    other : :obj:`conkit.core.ContactMap`, optional
       A ConKit :obj:`conkit.core.ContactMap`
    reference : :obj:`conkit.core.ContactMap`, optional
       A ConKit :obj:`conkit.core.ContactMap` [this map refers to the reference contacts]
    altloc : bool
       Use the res_altloc positions [default: False]
    dpi : int, optional
       The resolution of the plot [default: 300]
    file_format : str, optional
       Plot figure format. See :func:`matplotlib.pyplot.savefig` for options  [default: png]
    file_name : str, optional
       File name to which the contact map will be printed  [default: contactmap.png]

    Warnings
    --------
    If the ``file_name`` variable is not changed, the current file will be
    continuously overwritten.

    Raises
    ------
    RuntimeError
       Matplotlib not installed
    RuntimeError
       Hierarchy is not :obj:`conkit.core.ContactMap`

    """
    if not MATPLOTLIB:
        raise RuntimeError("Matplotlib not installed")

    if not isinstance(hierarchy, ContactMap):
        raise RuntimeError("Provided hierarchy is not a ContactMap")
    elif other and not isinstance(other, ContactMap):
        raise RuntimeError("Provided other hierarchy is not a ContactMap")
    elif reference and not isinstance(reference, ContactMap):
        raise RuntimeError("Provided reference hierarchy is not a ContactMap")

    fig, ax = matplotlib.pyplot.subplots(dpi=dpi)

    # Plot the other_ref contacts
    if reference:
        if altloc:
            reference_data = numpy.asarray([(c.res1_altseq, c.res2_altseq)
                                            for c in reference if c.is_true_positive])
        else:
            reference_data = numpy.asarray([(c.res1_seq, c.res2_seq)
                                            for c in reference if c.is_true_positive])
        reference_colors = [constants.RFCOLOR for _ in range(len(reference_data))]
        ax.scatter(reference_data.T[0], reference_data.T[1], color=reference_colors,
                   marker='.', edgecolor='none', linewidths=0.0)
        ax.scatter(reference_data.T[1], reference_data.T[0], color=reference_colors,
                   marker='.', edgecolor='none', linewidths=0.0)

    # Plot the self contacts
    self_data = numpy.asarray([(c.res1_seq, c.res2_seq) for c in hierarchy])
    self_colors = [
        constants.TPCOLOR if contact.is_true_positive
        else constants.FPCOLOR if contact.is_false_positive
        else constants.NTCOLOR for contact in hierarchy
    ]
    # This is the bottom triangle
    ax.scatter(self_data.T[1], self_data.T[0], color=self_colors,
               marker='.', edgecolor='none', linewidths=0.0)

    # Plot the other contacts
    if other:
        other_data = numpy.asarray([(c.res1_seq, c.res2_seq) for c in other])
        other_colors = [
            constants.TPCOLOR if contact.is_true_positive
            else constants.FPCOLOR if contact.is_false_positive
            else constants.NTCOLOR for contact in other
        ]
        # This is the upper triangle
        ax.scatter(other_data.T[0], other_data.T[1], color=other_colors,
                   marker='.', edgecolor='none', linewidths=0.0)
    else:
        # This is the upper triangle
        ax.scatter(self_data.T[0], self_data.T[1], color=self_colors,
                   marker='.', edgecolor='none', linewidths=0.0)

    # Allow dynamic x and y limits
    min_res_seq = numpy.min(self_data.ravel())
    max_res_seq = numpy.max(self_data.ravel())
    if other:
        min_res_seq = numpy.min(numpy.append(self_data.ravel(), other_data.ravel()))
        max_res_seq = numpy.max(numpy.append(self_data.ravel(), other_data.ravel()))
    ax.set_xlim(min_res_seq - 0.5, max_res_seq + 0.5)
    ax.set_ylim(min_res_seq - 0.5, max_res_seq + 0.5)

    # Set the xticks and yticks dynamically
    tick_range = numpy.arange(min_res_seq, max_res_seq, 10)
    ax.set_xticks(tick_range)
    ax.set_yticks(tick_range)

    # Prettify the plot
    ax.set_xlabel('Residue number')
    ax.set_ylabel('Residue number')

    # Create a custom legend
    if reference:
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

    _, file_extension = file_name.rsplit('.', 1)
    if file_extension != file_format:
        raise ValueError('File extension and file format have to be identical: '
                         '{0} - {1} are not'.format(file_extension, file_format))
    fig.savefig(file_name, format=file_format.lower(), bbox_inches='tight')


def sequence_coverage(hierarchy, dpi=300, file_name='seqcov.png', file_format='png'):
    """Plot the sequence coverage of each alignment position

    This function calculates and plots the coverage of sequences at
    each position in the Multiple Sequence Alignment.

    Parameters
    ----------
    hierarchy : :obj:`conkit.core.SequenceFile`
       The Multiple Sequence Alignment hierarchy
    dpi : int, optional
       The resolution of the plot [default: 300]
    file_format : str, optional
       Plot figure format. See :func:`matplotlib.pyplot.savefig` for options  [default: png]
    file_name : str, optional
       File name to which the plot will be printed  [default: seqcov.png]

    Warnings
    --------
    If the ``file_name`` variable is not changed, the current file will be
    continuously overwritten.

    Raises
    ------
    RuntimeError
       :obj:`conkit.core.SequenceFile` is not an alignment
    RuntimeError
       Matplotlib not installed
    RuntimeError
       Hierarchy is not :obj:`conkit.core.SequenceFile`

    """
    if not MATPLOTLIB:
        raise RuntimeError("Matplotlib not installed")

    if not isinstance(hierarchy, SequenceFile):
        raise RuntimeError("Provided hierarchy is not a SequenceFile")

    if not hierarchy.is_alignment:
        raise RuntimeError("Provided hierarchy does not show characteristics of an alignment")

    residues = numpy.arange(1, hierarchy.top_sequence.seq_len + 1)
    aa_frequencies = numpy.asarray(hierarchy.calculate_freq()) * hierarchy.top_sequence.seq_len

    fig, ax = matplotlib.pyplot.subplots(dpi=dpi)
    ax.plot(residues, aa_frequencies, color='#000000', marker='.', linestyle='-',
            label='Amino acid count')

    ax.axhline(hierarchy.top_sequence.seq_len * 0.3, color='r', label='30% Coverage')
    ax.axhline(hierarchy.top_sequence.seq_len * 0.6, color='g', label='60% Coverage')

    # Prettify the plot
    ax.set_xlabel('Residue number')
    ax.set_ylabel('Sequence Count')
    ax.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=3, mode="expand", borderaxespad=0.)

    # Make axes length proportional and remove whitespace around the plot
    ax.set(aspect=0.3)
    fig.tight_layout()

    _, file_extension = file_name.rsplit('.', 1)
    if file_extension != file_format:
        raise ValueError('File extension and file format have to be identical: '
                         '{0} - {1} are not'.format(file_extension, file_format))
    fig.savefig(file_name, format=file_format.lower(), bbox_inches='tight')
