"""
Storage space for a contact map
"""

from __future__ import division

__author__ = "Felix Simkovic"
__credits__ = "Jens Thomas"
__date__ = "03 Aug 2016"
__version__ = 0.1

from conkit import constants
from conkit.core.Entity import Entity
from conkit.core.Sequence import Sequence

import collections
import itertools
import matplotlib.pyplot
import numpy
import warnings


class ContactMap(Entity):
    """A contact map object representing a single prediction

    Description
    -----------
    The :obj:`ContactMap` class represents a data structure to hold a single
    contact map prediction in one place. It contains functions to store,
    manipulate and organise :obj:`Contact` instances.

    Attributes
    ----------
    coverage : float
       The sequence coverage score
    id : str
       A unique identifier
    ncontacts : int
       The number of :obj:`Contact` instances in the :obj:`ContactMap`
    precision : float
       The precision (Positive Predictive Value) score
    repr_sequence : :obj:`Sequence`
       The representative :obj:`Sequence` associated with the :obj:`ContactMap`
    repr_sequence_altloc : :obj:`Sequence`
       The representative altloc :obj:`Sequence` associated with the :obj:`ContactMap`
    sequence : :obj:`Sequence`
       The :obj:`Sequence` associated with the :obj:`ContactMap`
    top_contact : :obj:`Contact`
       The first :obj:`Contact` entry in :obj:`ContactMap`

    Examples
    --------
    >>> from conkit.core import Contact, ContactMap
    >>> contact_map = ContactMap("example")
    >>> contact_map.add(Contact(1, 10, 0.333))
    >>> contact_map.add(Contact(5, 30, 0.667))
    >>> print(contact_map)
    ContactMap(id="example" ncontacts=2)

    """
    __slots__ = ['_sequence']

    def __init__(self, id):
        """Initialise a new contact map"""
        self._sequence = None
        super(ContactMap, self).__init__(id)

    def __repr__(self):
        return "ContactMap(id=\"{0}\", ncontacts={1})".format(self.id, self.ncontacts)

    @property
    def coverage(self):
        """The sequence coverage score

        The coverage score is calculated by analysing the number of residues
        covered by the predicted contact pairs.

        .. math::

           Coverage=\\frac{x_{cov}}{L}

        The coverage score is calculated by dividing the number of contacts
        :math:`x_{cov}` by the number of residues in the sequence :math:`L`.

        Returns
        -------
        cov : float
           The calculated coverage score

        See Also
        --------
        precision

        """
        seq_array = numpy.fromstring(self.repr_sequence.seq, dtype='uint8')
        gaps = numpy.where(seq_array == ord('-'), 1, 0)
        cov = (seq_array.size - numpy.sum(gaps, axis=0)) / seq_array.size
        return cov

    @property
    def ncontacts(self):
        """The number of :obj:`Contact` instances in the :obj:`ContactMap`

        Returns
        -------
        ncontacts : int
           The number of sequences in the :obj:`ContactMap`

        """
        return len(self)

    @property
    def precision(self):
        """The precision (Positive Predictive Value) score

        The precision value is calculated by analysing the true and false
        postive contacts.

        .. math::

           Precision=\\frac{TruePositives}{TruePositives - FalsePositives}

        The status of each contact, i.e true or false positive status, can be
        determined by running the :obj:`match()` function providing a reference
        structure.

        Returns
        -------
        ppv : float
           The calculated precision score

        See Also
        --------
        coverage

        """
        status_array = numpy.asarray([c.status for c in self])

        fp_count = numpy.sum(numpy.where(status_array == -1, 1, 0))
        uk_count = numpy.sum(numpy.where(status_array == 0, 1, 0))
        tp_count = numpy.sum(numpy.where(status_array == 1, 1, 0))

        if fp_count == 0.0 and tp_count == 0.0:
            warnings.warn('Contact status unknown. Match two contact maps first!')
            return 0.0
        elif uk_count > 0:
            warnings.warn('Some contacts in your map are unmatched - Precision value might be inaccurate')

        return tp_count / (tp_count + fp_count)

    @property
    def repr_sequence(self):
        """The representative :obj:`Sequence` associated with the :obj:`ContactMap`

        The peptide sequence constructed from the available
        contacts using the normal res_seq positions

        Returns
        -------
        sequence : :obj:`Sequence`

        Raises
        ------
        TypeError
           Sequence undefined

        See Also
        --------
        repr_sequence_altloc, sequence

        """
        if not isinstance(self.sequence, Sequence):
            raise TypeError('Define the sequence as Sequence() instance')
        # Get all resseqs that are the contact map
        res1_seqs, res2_seqs = zip(*[contact.id for contact in self])
        res_seqs = set(
            sorted(res1_seqs + res2_seqs)
        )
        return self._construct_repr_sequence(list(res_seqs))

    @property
    def repr_sequence_altloc(self):
        """The representative altloc :obj:`Sequence` associated with the :obj:`ContactMap`

        The peptide sequence constructed from the available
        contacts using the altloc res_seq positions

        Returns
        -------
        sequence : :obj:`Sequence`

        Raises
        ------
        ValueError
           Sequence undefined

        See Also
        --------
        repr_sequence, sequence

        """
        if not isinstance(self.sequence, Sequence):
            raise TypeError('Define the sequence as Sequence() instance')
        # Get all resseqs that are the contact map
        res1_seqs, res2_seqs = zip(*[(contact.res1_altseq, contact.res2_altseq) for contact in self])
        res_seqs = set(
            sorted(res1_seqs + res2_seqs)
        )
        return self._construct_repr_sequence(list(res_seqs))

    @property
    def sequence(self):
        """The :obj:`Sequence` associated with the :obj:`ContactMap`

        Returns
        -------
        :obj:`Sequence`

        See Also
        --------
        repr_sequence, repr_sequence_altloc

        """
        return self._sequence

    @sequence.setter
    def sequence(self, sequence):
        """Associate a :obj:`Sequence` instance with the :obj:`ContactMap`

        Parameters
        ----------
        sequence : :obj:`Sequence`

        Raises
        ------
        ValueError
           Incorrect hierarchy instance provided

        """
        if not isinstance(sequence, Sequence):
            raise TypeError('Instance of Sequence() required: {0}'.format(sequence))
        self._sequence = sequence

    @property
    def top_contact(self):
        """The first :obj:`Contact` entry in :obj:`ContactMap`

        Returns
        -------
        top_contact : :obj:`Contact`, None
           The first :obj:`Contact` entry in :obj:`ContactFile`

        """
        if len(self) > 0:
            return self[0]
        else:
            return None

    def _construct_repr_sequence(self, res_seqs):
        """Construct the representative sequence"""
        # Determine which are present and which are not
        representative_sequence = ''
        for i in xrange(1, self.sequence.seq_len + 1):
            if i in res_seqs:
                representative_sequence += self.sequence.seq[i - 1]
            else:
                representative_sequence += '-'
        return Sequence(self.sequence.id + '_repr', representative_sequence)

    def assign_sequence_register(self, altloc=False):
        """Assign the amino acids from :obj:`Sequence` to all :obj:`Contact` instances

        Parameters
        ----------
        altloc : bool
           Use the res_altloc positions [default: False]

        """
        for c in self:
            if altloc:
                res1_index, res2_index = c.res1_altseq, c.res2_altseq
            else:
                res1_index, res2_index = c.res1_seq, c.res2_seq
            c.res1 = self.sequence.seq[res1_index - 1]
            c.res2 = self.sequence.seq[res2_index - 1]

    def calculate_scalar_score(self):
        """Calculate a scaled score for the :obj:`ContactMap`

        This score is a scaled score for all raw scores in a contact
        map. It is defined by the formula

        .. math::

           {x}'=\\frac{x}{\\overline{d}}

        where :math:`x` corresponds to the raw score of each predicted
        contact and :math:`\overline{d}` to the mean of all raw scores.

        The score is saved in a separate :obj:`Contact` attribute called
        ``scalar_score``

        This score is described in more detail in Ovchinnikov et al. (2015).

        """
        raw_scores = numpy.asarray([c.raw_score for c in self])
        sca_scores = raw_scores / numpy.mean(raw_scores)
        for contact, sca_score in zip(self, sca_scores):
            contact.scalar_score = sca_score
        return

    def create_keymap(self, altloc=False):
        """Create a simple keymap

        Paramters
        ---------
        altloc : bool
           Use the res_altloc positions [default: False]

        """
        # Template namedtuple for temporary storage
        Template = collections.namedtuple('Template', ['res_seq', 'res_altseq', 'res_matchseq', 'res_name', 'status'])

        # Faster and remembers order?!
        cmap_mapping = collections.OrderedDict()
        for c in self:
            pos1 = Template(res_seq=c.res1_seq, res_altseq=c.res1_altseq, res_matchseq=constants.UNKNOWN,
                            res_name=c.res1, status=constants.UNREGISTERED)
            pos2 = Template(res_seq=c.res2_seq, res_altseq=c.res2_altseq, res_matchseq=constants.UNKNOWN,
                            res_name=c.res2, status=constants.UNREGISTERED)

            if altloc:
                res1_index, res2_index = c.res1_altseq, c.res2_altseq
            else:
                res1_index, res2_index = c.res1_seq, c.res2_seq

            if res1_index in cmap_mapping and cmap_mapping[res1_index]:
                assert cmap_mapping[res1_index] == pos1
            if res2_index in cmap_mapping and cmap_mapping[res2_index]:
                assert cmap_mapping[res2_index] == pos2

            cmap_mapping[res1_index] = pos1
            cmap_mapping[res2_index] = pos2

        cmap_sorted = sorted(cmap_mapping.items(), key=lambda x: int(x[0]))
        return [nt._asdict() for nt in zip(*cmap_sorted)[1]]

    def match(self, other, remove_unmatched=False, renumber=False, inplace=False):
        """Modify both hierarchies so residue numbers match one another.

        This function is key when plotting contact maps or visualising
        contact maps in 3-dimensional space. In particular, when residue
        numbers in the structure do not start at count 0 or when peptide
        chain breaks are present.

        Parameters
        ----------
        other : :obj:`ContactMap`
           A ConKit :obj:`ContactMap`
        remove_unmatched : bool, optional
           Remove all unmatched contacts [default: False]
        renumber : bool, optional
           Renumber the res_seq entries [default: False]

           If ``True``, ``res1_seq`` and ``res2_seq`` changes
           but ``id`` remains the same
        inplace : bool, optional
           Replace the saved order of contacts [default: False]

        Returns
        -------
        hierarchy_mod
           A ConKit :obj:`ContactMap`, regardless of inplace
        hierarchy_ref
           A ConKit :obj:`ContactMap`, regardless of inplace

        Warnings
        --------
        This function is highly unstable for inter-molecular contacts. Please remove all
        intra-molecular contacts from both maps manually before matching!

        """
        warnings.warn('Unstable for inter-molecular contacts!')

        contact_map1 = self._inplace(inplace)
        contact_map2 = other._inplace(inplace)

        # Align both full sequences against each other
        contact_map1.sequence.align_local(contact_map2.sequence,
            id_chars=2, nonid_chars=1, gap_open_pen=-0.5, gap_ext_pen=-0.1,
            inplace=True
        )
        # Align the repr sequence of map 1 against full sequence of map 1
        contact_map1.sequence.align_local(contact_map1.repr_sequence,
            id_chars=2, nonid_chars=1, gap_open_pen=-0.5, gap_ext_pen=-0.2,
            inplace=True
        )
        # Align the __altloc__ repr sequence of map 2 against full sequence of map 2
        contact_map2.sequence.align_local(contact_map2.repr_sequence_altloc,
            id_chars=2, nonid_chars=1, gap_open_pen=-0.5, gap_ext_pen=-0.2,
            inplace=True
        )
        # Align the repr sequence of map 1 against the __altloc__ repr sequence of map 2
        contact_map1.repr_sequence.align_local(
            contact_map2.repr_sequence_altloc,
            id_chars=2, nonid_chars=1, gap_open_pen=-1.0, gap_ext_pen=-0.5,
            inplace=True
        )

        # Make sure all of our sequences are identical in length
        if not (contact_map1.sequence.seq_len
                == contact_map1.repr_sequence.seq_len
                == contact_map2.sequence.seq_len
                == contact_map2.repr_sequence_altloc.seq_len):
            raise ValueError('Sequences do not match in length')

        if all(c == '-' for c in contact_map1.repr_sequence.seq):
            raise ValueError('No representative sequence available for - {0}'.format(contact_map1))
        elif all(c == '-' for c in contact_map2.repr_sequence_altloc.seq):
            raise ValueError('No representative sequence available for - {0}'.format(contact_map2))

        # Create default residue key mappings
        cmap1_keymap = contact_map1.create_keymap()
        cmap2_keymap = contact_map2.create_keymap(altloc=True)

        # ------------------------------------------------------------------------------------------------------------------
        # Match the mappings using their res_seq entries
        cursor1 = cursor2 = 0
        for pos1, pos2 in zip(contact_map1.repr_sequence.seq, contact_map2.repr_sequence_altloc.seq):
            if pos1 == pos2 and pos1 == '-':  # Gap
                continue
            elif pos1 != pos2 and pos1 == '-':  # Gap in cmap1 sequence
                cmap2_keymap[cursor2]['status'] = constants.UNMATCHED
                cursor2 += 1
            elif pos1 != pos2 and pos2 == '-':  # Gap in cmap2 sequence
                cmap1_keymap[cursor1]['status'] = constants.UNMATCHED
                cursor1 += 1
            else:  # Match
                if cmap1_keymap[cursor1]['res_name'] != cmap2_keymap[cursor2]['res_name']:
                    raise RuntimeError('Matched residues with different residue names')
                cmap1_keymap[cursor1]['res_matchseq'] = cmap2_keymap[cursor2]['res_seq']
                cmap1_keymap[cursor1]['status'] = constants.MATCHED
                cmap2_keymap[cursor2]['res_matchseq'] = cmap1_keymap[cursor1]['res_seq']
                cmap2_keymap[cursor2]['status'] = constants.MATCHED
                cursor1, cursor2 = cursor1 + 1, cursor2 + 1

        # ------------------------------------------------------------------------------------------------------------------
        # Modify contact_map1 to match contact_map2
        for combo in itertools.combinations(cmap1_keymap, 2):

            # Note: contact id is tuple of strings not ints
            contact_id_orig = tuple([combo[0]['res_seq'], combo[1]['res_seq']])
            contact_id_matched = tuple([combo[0]['res_matchseq'], combo[1]['res_matchseq']])
            contact_status = tuple([combo[0]['status'], combo[1]['status']])

            # Skip non-existent combinations
            if contact_id_orig not in contact_map1:
                continue

            # THIS FIXES IT ALL - WORK OUT THE LIMITS FOR INTER-CHAINS
            # elif contact_id_orig[0] < 128 and contact_id_orig[1] < 128:
            #     cmap1.remove(contact_id_orig)
            #     continue
            # elif contact_id_orig[0] > 127 and contact_id_orig[1] > 127:
            #     cmap1.remove(contact_id_orig)
            #     continue

            # Define the TRUE and FALSE positive status of the contact
            if sum(contact_status) == constants.MATCHED and contact_id_matched in contact_map2:
                contact_map1[contact_id_orig].define_true_positive()
            elif sum(contact_status) == constants.MATCHED:
                contact_map1[contact_id_orig].define_false_positive()

            # Note: It is essential to renumber the contact first!
            #       To avoid indexing issues or formatting problems, we do __NOT__ change the
            #       key in our contact map of this particular contact -
            #       only the res1_seq, res2_seq, res1_altseq and res2_altseq attributes.
            if renumber:
                contact_map1[contact_id_orig].res1_seq, contact_map1[contact_id_orig].res2_seq = contact_id_matched
                contact_map1[contact_id_orig].res1_altseq, contact_map1[contact_id_orig].res2_altseq = contact_id_orig

            if remove_unmatched and sum(contact_status) != 0:
                contact_map1.remove(contact_id_orig)

        return contact_map1, contact_map2

    def remove_neighbors(self, min_distance=5, inplace=False):
        """Remove contacts between neighboring residues

        Parameters
        ----------
        min_distance : int, optional
           The minimum number of residues between contacts  [default: 5]
        inplace : bool, optional
           Replace the saved order of contacts [default: False]

        Returns
        -------
        contact_map : :obj:`ContactMap`
           The reference to the :obj:`ContactMap`, regardless of inplace

        """
        contact_map = self._inplace(inplace)

        for contact in reversed(contact_map):
            if abs(contact.res1_seq - contact.res2_seq) < min_distance:
                contact_map.remove(contact.id)

        return contact_map

    def plot_map(self, other=None, altloc=False, file_format='png', file_name='contactmap.png'):
        """Produce a 2D contact map plot

        Parameters
        ----------
        other : :obj:`ContactMap`, optional
           A ConKit :obj:`ContactMap`
        altloc : bool
           Use the res_altloc positions [default: False]
        file_format : str, optional
           Plot figure format. See matplotlib :obj:`savefig()` for options  [default: png]
        file_name : str, optional
           File name to which the contact map will be printed  [default: contactmap.png]

        Warnings
        --------
        If the ``file_name`` variable is not changed, the current file will be
        continuously overwritten.

        """

        fig, ax = matplotlib.pyplot.subplots(figsize=(5, 5), dpi=600)

        def draw(x, y, c):
            ax.scatter(x, y, color=c, marker='.', s=10, edgecolor='none', linewidths=0.0)

        # Plot the other_ref contacts
        if other:
            if altloc:
                other_data = numpy.asarray([(c.res1_altseq, c.res2_altseq) for c in other])
            else:
                other_data = numpy.asarray([(c.res1_seq, c.res2_seq) for c in other])
            other_colors = [constants.RFCOLOR for _ in xrange(len(other))]
            draw(other_data.T[0], other_data.T[1], other_colors)  # upper
            draw(other_data.T[1], other_data.T[0], other_colors)  # lower

        # Plot the self contacts
        if altloc:
            self_data = numpy.asarray([(c.res1_altseq, c.res2_altseq) for c in self])
        else:
            self_data = numpy.asarray([(c.res1_seq, c.res2_seq) for c in self])
        self_colors = [
            constants.TPCOLOR if contact.is_true_positive else
            constants.FPCOLOR if contact.is_false_positive else
            constants.NTCOLOR
            for contact in self
        ]
        draw(self_data.T[0], self_data.T[1], self_colors)  # upper
        draw(self_data.T[1], self_data.T[0], self_colors)  # lower

        # Prettify the plot
        label = 'Residue number'
        ax.set_xlabel(label)
        ax.set_ylabel(label)

        if self.sequence is not None:
            lim = self.sequence.seq_len
        else:
            lim = self_data.flatten().max()
        ax.set_xlim(1, lim)
        ax.set_ylim(1, lim)

        fig.tight_layout()

        _, file_extension = file_name.rsplit('.', 1)
        if file_extension != file_format:
            raise ValueError('File extension and file format have to be identical: '
                             '{0} - {1} are not'.format(file_extension, file_format))
        fig.savefig(file_name, format=file_format.lower())

    def rescale(self, inplace=False):
        """Rescale the raw scores in :obj:`ContactMap`

        Rescaling of the data is done to normalize the raw scores
        to be in the range [0, 1]. The formula to rescale the data is:

        .. math::

           {x}'=\\frac{x-min(d)}{max(d)-min(d)}

        :math:`x` is the original value and :math:`d` are all values to be
        rescaled.

        Parameters
        ----------
        inplace : bool, optional
           Replace the saved order of contacts [default: False]

        Returns
        -------
        contact_map : :obj:`ContactMap`
           The reference to the :obj:`ContactMap`, regardless of inplace

        """
        contact_map = self._inplace(inplace)

        raw_scores = numpy.asarray([c.raw_score for c in contact_map])
        norm_raw_scores = (raw_scores - raw_scores.min()) / (raw_scores.max() - raw_scores.min())

        for contact, norm_raw_score in zip(contact_map, norm_raw_scores):
            contact.raw_score = norm_raw_score

        return contact_map

    def sort(self, kword, reverse=False, inplace=False):
        """Sort the :obj:`ContactMap`

        Parameters
        ----------
        kword : str
           The dictionary key to sort contacts by
        reverse : bool, optional
           Sort the contact pairs in descending order [default: False]
        inplace : bool, optional
           Replace the saved order of contacts [default: False]

        Returns
        -------
        contact_map : :obj:`ContactMap`
           The reference to the :obj:`ContactMap`, regardless of inplace

        Raises
        ------
        ValueError
           ``kword`` not in :obj:`ContactMap`

        """
        contact_map = self._inplace(inplace)
        contact_map._sort(kword, reverse)
        return contact_map
