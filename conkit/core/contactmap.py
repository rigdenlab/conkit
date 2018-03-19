# coding=utf-8
#
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
"""ContactMap container used throughout ConKit"""

from __future__ import division
from __future__ import print_function

__author__ = "Felix Simkovic"
__date__ = "03 Aug 2016"
__version__ = "1.0"

import collections
import numpy as np
import os
import sys

if sys.version_info.major < 3:
    from itertools import izip as zip

from conkit.core._entity import _Entity
from conkit.core._struct import _Gap, _Residue
from conkit.core.mappings import ContactMatchState
from conkit.core.sequence import Sequence
from conkit.misc import normalize


class ContactMap(_Entity):
    """A contact map object representing a single prediction

    The :obj:`ContactMap <conkit.core.contactmap.ContactMap>` class represents a data structure to hold a single
    contact map prediction in one place. It contains functions to store,
    manipulate and organise :obj:`Contact <conkit.core.contact.Contact>` instances.

    Attributes
    ----------
    coverage : float
       The sequence coverage score
    id : str
       A unique identifier
    ncontacts : int
       The number of :obj:`Contact <conkit.core.contact.Contact>` instances in the :obj:`ContactMap 
       <conkit.core.contactmap.ContactMap>`
    precision : float
       The precision (Positive Predictive Value) score
    repr_sequence : :obj:`Sequence <conkit.core.sequence.Sequence>`
       The representative :obj:`Sequence <conkit.core.sequence.Sequence>` associated with the :obj:`ContactMap
       <conkit.core.contactmap.ContactMap>`
    repr_sequence_altloc : :obj:`Sequence <conkit.core.sequence.Sequence>`
       The representative altloc :obj:`Sequence <conkit.core.sequence.Sequence>` associated with the :obj:`ContactMap 
       <conkit.core.contactmap.ContactMap>`
    sequence : :obj:`Sequence <conkit.core.sequence.Sequence>`
       The :obj:`Sequence <conkit.core.sequence.Sequence>` associated with the :obj:`ContactMap 
       <conkit.core.contactmap.ContactMap>`
    top_contact : :obj:`Contact <conkit.core.contact.Contact>`
       The first :obj:`Contact <conkit.core.contact.Contact>` entry in :obj:`ContactMap 
       <conkit.core.contactmap.ContactMap>`

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
        return "{0}(id=\"{1}\", ncontacts={2})".format(self.__class__.__name__, self.id, self.ncontacts)

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
        float
           The calculated coverage score

        See Also
        --------
        precision

        """
        seq_array = np.array(list(self.repr_sequence.seq_ascii))
        gaps = np.where(seq_array == ord('-'), 1, 0)
        return (seq_array.size - np.sum(gaps, axis=0)) / seq_array.size

    @property
    def empty(self):
        """Empty contact map"""
        return len(self) < 1

    @property
    def ncontacts(self):
        """The number of :obj:`Contact <conkit.core.contact.Contact>` instances in the :obj:`ContactMap <conkit.core.contactmap.ContactMap>`

        Returns
        -------
        int
           The number of contacts in the :obj:`ContactMap <conkit.core.contactmap.ContactMap>`

        """
        return len(self)

    @property
    def short_range_contacts(self):
        """The short range contacts found :obj:`ContactMap <conkit.core.contactmap.ContactMap>`

        Short range contacts are defined as 6 <= x <= 11 residues apart

        Returns
        -------
        obj
           A copy of the :obj:`ContactMap <conkit.core.contactmap.ContactMap>` with short-range contacts only

        See Also
        --------
        medium_range_contacts, long_range_contacts

        """
        return self.remove_neighbors(min_distance=6, max_distance=11)

    @property
    def medium_range_contacts(self):
        """The medium range contacts found :obj:`ContactMap <conkit.core.contactmap.ContactMap>`

        Medium range contacts are defined as 12 <= x <= 23 residues apart

        Returns
        -------
        obj
           A copy of the :obj:`ContactMap <conkit.core.contactmap.ContactMap>` with medium-range contacts only

        See Also
        --------
        short_range_contacts, long_range_contacts

        """
        return self.remove_neighbors(min_distance=12, max_distance=23)

    @property
    def long_range_contacts(self):
        """The long range contacts found :obj:`ContactMap <conkit.core.contactmap.ContactMap>`

        long range contacts are defined as 24 <= x residues apart

        Returns
        -------
        obj
           A copy of the :obj:`ContactMap <conkit.core.contactmap.ContactMap>` with long-range contacts only

        See Also
        --------
        short_range_contacts, medium_range_contacts

        """
        return self.remove_neighbors(min_distance=24)

    @property
    def precision(self):
        """The precision (Positive Predictive Value) score

        The precision value is calculated by analysing the true and false
        postive contacts.

        .. math::

           Precision=\\frac{TruePositives}{TruePositives - FalsePositives}

        The status of each contact, i.e true or false positive status, can be
        determined by running the :func:`match` function providing a reference
        structure.

        Returns
        -------
        float
           The calculated precision score

        See Also
        --------
        coverage

        """
        if self.empty:
            return 0.0

        import warnings

        s = np.array([c.status for c in self])
        cdict = dict(zip(np.unique(s), np.array([s[s == i].shape[0] for i in np.unique(s)])))
        fp_count = cdict[ContactMatchState.mismatched.value] if ContactMatchState.mismatched.value in cdict else 0.0
        uk_count = cdict[ContactMatchState.unknown.value] if ContactMatchState.unknown.value in cdict else 0.0
        tp_count = cdict[ContactMatchState.matched.value] if ContactMatchState.matched.value in cdict else 0.0

        if fp_count == 0.0 and tp_count == 0.0:
            warnings.warn("No matches or mismatches found in your contact map. " "Match two ContactMaps first.")
            return 0.0
        elif uk_count > 0:
            warnings.warn("Some contacts between the ContactMaps are unmatched due to non-identical "
                          "sequences. The precision value might be inaccurate.")

        return tp_count / (tp_count + fp_count)

    @property
    def repr_sequence(self):
        """The representative :obj:`Sequence <conkit.core.sequence.Sequence>` associated with the :obj:`ContactMap
        <conkit.core.contactmap.ContactMap>`

        The peptide sequence constructed from the available
        contacts using the normal res_seq positions

        Returns
        -------
        obj
           A :obj:`conkit.coreSequence` object

        Raises
        ------
        TypeError
           Sequence undefined

        See Also
        --------
        repr_sequence_altloc, sequence

        """
        if isinstance(self.sequence, Sequence):
            res_seqs = np.unique(np.array(self.as_list()).flatten()).tolist()
            return self._construct_repr_sequence(res_seqs)
        else:
            raise TypeError('Define the sequence as Sequence() instance')

    @property
    def repr_sequence_altloc(self):
        """The representative altloc :obj:`Sequence <conkit.core.sequence.Sequence>` associated with the :obj:`ContactMap
        <conkit.core.contactmap.ContactMap>`

        The peptide sequence constructed from the available
        contacts using the altloc res_seq positions

        Returns
        -------
        obj
           A :obj:`Sequence <conkit.core.sequence.Sequence>` object

        Raises
        ------
        ValueError
           Sequence undefined

        See Also
        --------
        repr_sequence, sequence

        """
        if isinstance(self.sequence, Sequence):
            res_seqs = np.unique(np.array(self.as_list(altloc=True)).flatten()).tolist()
            return self._construct_repr_sequence(res_seqs)
        else:
            raise TypeError('Define the sequence as Sequence() instance')

    @property
    def singletons(self):
        """Singleton contact pairs in the current :obj:`ContactMap <conkit.core.contactmap.ContactMap>`
        
        Contacts are identified by a distance-based grouping analysis. A :obj:`Contact <conkit.core.contact.Contact>` is
        classified as singleton if not other contacts are found within 2 residues.

        Returns
        -------
        :obj:`ContactMap <conkit.core.contactmap.ContactMap>`

        """
        contacts = np.array(self.as_list(), dtype=np.int64)
        removables = set()
        for i in np.arange(contacts.shape[0]):
            for j in np.arange(i + 1, contacts.shape[0]):
                dist = np.abs(contacts[i] - contacts[j])
                if np.all(dist <= 2):
                    removables.update({tuple(contacts[i]), tuple(contacts[j])})
        singletons = self.deepcopy()
        [singletons.remove(i) for i in removables]
        return singletons

    @property
    def sequence(self):
        """The :obj:`Sequence <conkit.core.sequence.Sequence>` associated with the :obj:`ContactMap <conkit.core.contactmap.ContactMap>`

        Returns
        -------
        obj
           A :obj:`Sequence <conkit.core.sequence.Sequence>` object

        See Also
        --------
        repr_sequence, repr_sequence_altloc

        """
        return self._sequence

    @sequence.setter
    def sequence(self, sequence):
        """Associate a :obj:`Sequence <conkit.core.sequence.Sequence>` instance with the :obj:`ContactMap
        <conkit.core.contactmap.ContactMap>`

        Parameters
        ----------
        sequence : :obj:`Sequence <conkit.core.sequence.Sequence>`

        Raises
        ------
        ValueError
           Incorrect hierarchy instance provided

        """
        if isinstance(sequence, Sequence):
            self._sequence = sequence
        else:
            msg = "Instance of Sequence() required: {}".format(sequence)
            raise TypeError(msg)

    @property
    def top_contact(self):
        """The first :obj:`Contact <conkit.core.contact.Contact>` entry in :obj:`ContactMap <conkit.core.contactmap.ContactMap>`

        Returns
        -------
        obj, None
           The first :obj:`Contact <conkit.core.contact.Contact>` entry in :obj:`ContactFile <conkit.core.contactfile.ContactFile>`

        """
        return self.top

    def _construct_repr_sequence(self, res_seqs):
        """Construct the representative sequence"""
        representative_sequence = ''
        for i in np.arange(1, self.sequence.seq_len + 1):
            if i in res_seqs:
                representative_sequence += self.sequence.seq[i - 1]
            else:
                representative_sequence += '-'
        return Sequence(self.sequence.id + '_repr', representative_sequence)

    def as_list(self, altloc=False):
        """The :obj:`ContactMap <conkit.core.contactmap.ContactMap>` as a 2D-list containing contact-pair residue indexes

        Parameters
        ----------
        altloc : bool
           Use the res_altloc positions [default: False]

        """
        if altloc:
            return [[c.res1_altseq, c.res2_altseq] for c in self]
        else:
            return [[c.res1_seq, c.res2_seq] for c in self]

    def assign_sequence_register(self, altloc=False):
        """Assign the amino acids from :obj:`Sequence <conkit.core.sequence.Sequence>` to all :obj:`Contact <conkit.core.contact.Contact>` instances

        Parameters
        ----------
        altloc : bool
           Use the res_altloc positions [default: False]

        """
        for c in self:
            if altloc:
                res1_index = c.res1_altseq
                res2_index = c.res2_altseq
            else:
                res1_index = c.res1_seq
                res2_index = c.res2_seq
            c.res1 = self.sequence.seq[res1_index - 1]
            c.res2 = self.sequence.seq[res2_index - 1]

    def calculate_jaccard_index(self, other):
        """Calculate the Jaccard index between two :obj:`ContactMap <conkit.core.contactmap.ContactMap>` instances

        This score analyzes the difference of the predicted contacts from two maps,

        .. math::

           J_{x,y}=\\frac{\\left|x \\cap y\\right|}{\\left|x \\cup y\\right|}

        where :math:`x` and :math:`y` are the sets of predicted contacts from two
        different predictors, :math:`\\left|x \\cap y\\right|` is the number of
        elements in the intersection of :math:`x` and :math:`y`, and the
        :math:`\\left|x \\cup y\\right|` represents the number of elements in the
        union of :math:`x` and :math:`y`.

        The J-score has values in the range of :math:`[0, 1]`, with a value of :math:`1`
        corresponding to identical contact maps and :math:`0` to dissimilar ones.

        Parameters
        ----------
        other : :obj:`ContactMap <conkit.core.contactmap.ContactMap>`
           A ConKit :obj:`ContactMap <conkit.core.contactmap.ContactMap>`

        Returns
        -------
        float
           The Jaccard index

        See Also
        --------
        match, precision

        Warnings
        --------
        The Jaccard distance ranges from :math:`[0, 1]`, where :math:`1` means
        the maps contain identical contacts pairs.

        Notes
        -----
        The Jaccard index is different from the Jaccard distance mentioned in [#]_. The
        Jaccard distance corresponds to :math:`1-Jaccard_{index}`.

        .. [#] Q. Wuyun, W. Zheng, Z. Peng, J. Yang (2016). A large-scale comparative assessment
           of methods for residue-residue contact prediction. *Briefings in Bioinformatics*,
           [doi: 10.1093/bib/bbw106].

        """
        union = len(self) + np.sum([1 for contact in other if contact.id not in self])
        if union == 0:
            return 1.0
        intersection = np.sum([1 for contact in self if contact.id in other])
        return float(intersection) / union

    def calculate_kernel_density(self, *args, **kwargs):
        """Calculate the contact density in the contact map using Gaussian kernels"""
        import warnings
        warnings.warn("This function will be deprecated in a future release! Use calculate_contact_density instead!")
        return self.calculate_contact_density(*args, **kwargs)

    def calculate_contact_density(self, bw_method="amise"):
        """Calculate the contact density in the contact map using Gaussian kernels

        Various algorithms can be used to estimate the bandwidth. To calculate the
        bandwidth for an 1D data array ``X`` with ``n`` data points and ``d`` dimensions,
        the listed algorithms have been implemented. Please note, in rules 2 and 3, the
        value of :math:`\\sigma` is the smaller of the standard deviation of ``X`` or
        the normalized interquartile range.

        Parameters
        ----------
        bw_method : str, optional
           The bandwidth estimator to use [default: amise]

        Returns
        -------
        list
           The list of per-residue density estimates

        Raises
        ------
        RuntimeError
           Cannot find SciKit package
        ValueError
           Undefined bandwidth method

        """
        try:
            import sklearn.neighbors
        except ImportError:
            raise RuntimeError("Cannot find SciKit package")

        if self.empty:
            raise ValueError("ContactMap is empty")

        # TODO: Chunan suggested to fix this bug - results are usually marginally better
        # REM: Bug in Sadowski's algorithm, res2 is excluded from list to train KDE
        # REM: Remember to change test cases when corrected implementation benchmarked
        #  x = np.array([i for c in self for i in np.arange(c.res1_seq, c.res2_seq + 1)])[:, np.newaxis]
        x = np.array([i for c in self for i in np.arange(c.res1_seq, c.res2_seq)])[:, np.newaxis]
        x_fit = np.arange(x.min(), x.max() + 1)[:, np.newaxis]
        from conkit.misc.bandwidth import bandwidth_factory
        bandwidth = bandwidth_factory(bw_method)(x).bw
        kde = sklearn.neighbors.KernelDensity(bandwidth=bandwidth).fit(x)
        return np.exp(kde.score_samples(x_fit)).tolist()

    def calculate_scalar_score(self):
        """Calculate a scaled score for the :obj:`ContactMap <conkit.core.contactmap.ContactMap>`

        This score is a scaled score for all raw scores in a contact
        map. It is defined by the formula

        .. math::

           {x}'=\\frac{x}{\\overline{d}}

        where :math:`x` corresponds to the raw score of each predicted
        contact and :math:`\overline{d}` to the mean of all raw scores.

        The score is saved in a separate :obj:`Contact <conkit.core.contact.Contact>` attribute called ``scalar_score``

        This score is described in more detail in [#]_.

        .. [#] S. Ovchinnikov, L. Kinch, H. Park, Y. Liao, J. Pei, D.E. Kim,
           H. Kamisetty, N.V. Grishin, D. Baker (2015). Large-scale determination
           of previously unsolved protein structures using evolutionary information.
           *Elife* **4**, e09248.

        """
        raw_scores = np.array([c.raw_score for c in self])
        sca_scores = raw_scores / np.mean(raw_scores)
        for contact, sca_score in zip(self, sca_scores):
            contact.scalar_score = sca_score

    def find(self, register, altloc=False, strict=False):
        """Find all contacts with one or both residues in ``register``

        Parameters
        ----------
        register : int, list, tuple
           A list of residue register to find
        altloc : bool
           Use the res_altloc positions [default: False]
        strict : bool
           Both residues of :obj:`Contact <conkit.core.contact.Contact>` in register [default: False]

        Returns
        -------
        obj
           A modified version of the :obj:`ContactMap <conkit.core.contactmap.ContactMap>` containing
           the found contacts

        """
        if isinstance(register, int):
            register = [register]
        register = set(register)
        comparison_operator = _AND if strict else _OR
        contact_map = self.deepcopy()
        for contactid in self.as_list(altloc=altloc):
            if not comparison_operator(contactid[0] in register, contactid[1] in register):
                contact_map.remove(tuple(contactid))
        return contact_map

    def match(self, other, match_other=False, remove_unmatched=False, renumber=False, inplace=False):
        """Modify both hierarchies so residue numbers match one another.

        This function is key when plotting contact maps or visualising
        contact maps in 3-dimensional space. In particular, when residue
        numbers in the structure do not start at count 0 or when peptide
        chain breaks are present.

        Parameters
        ----------
        other : :obj:`ContactMap <conkit.core.contactmap.ContactMap>`
           A ConKit :obj:`ContactMap <conkit.core.contactmap.ContactMap>`
        match_other: bool, optional
           Match `other` to `self` [default: False]
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
        obj
            :obj:`ContactMap <conkit.core.contactmap.ContactMap>` instance, regardless of inplace

        Raises
        ------
        ValueError
           Error creating reliable keymap matching the sequence in :obj:`ContactMap <conkit.core.contactmap.ContactMap>`

        """
        contact_map1 = self._inplace(inplace)
        if match_other:
            contact_map2 = other._inplace(inplace)
        else:
            contact_map2 = other._inplace(False)

        # ================================================================
        # 1. Align all sequences
        # ================================================================

        # Align both full sequences against each other
        aligned_sequences_full = contact_map1.sequence.align_local(
            contact_map2.sequence, id_chars=2, nonid_chars=1, gap_open_pen=-0.5, gap_ext_pen=-0.1)
        contact_map1_full_sequence, contact_map2_full_sequence = aligned_sequences_full

        # Align contact map 1 full sequences with representative sequence
        aligned_sequences_map1 = contact_map1_full_sequence.align_local(
            contact_map1.repr_sequence, id_chars=2, nonid_chars=1, gap_open_pen=-0.5, gap_ext_pen=-0.2, inplace=True)
        contact_map1_repr_sequence = aligned_sequences_map1[-1]

        # Align contact map 2 full sequences with __ALTLOC__ representative sequence
        aligned_sequences_map2 = contact_map2_full_sequence.align_local(
            contact_map2.repr_sequence_altloc,
            id_chars=2,
            nonid_chars=1,
            gap_open_pen=-0.5,
            gap_ext_pen=-0.2,
            inplace=True)
        contact_map2_repr_sequence = aligned_sequences_map2[-1]

        # Align both aligned representative sequences
        aligned_sequences_repr = contact_map1_repr_sequence.align_local(
            contact_map2_repr_sequence, id_chars=2, nonid_chars=1, gap_open_pen=-1.0, gap_ext_pen=-0.5, inplace=True)
        contact_map1_repr_sequence, contact_map2_repr_sequence = aligned_sequences_repr

        # ================================================================
        # 2. Identify TPs in other, map them, and match them to self
        # ================================================================

        # Encode the sequences to uint8 character arrays for easier and faster handling
        encoded_repr = np.array(
            [list(contact_map1_repr_sequence.seq_ascii),
             list(contact_map2_repr_sequence.seq_ascii)])

        # Create mappings for both contact maps
        contact_map1_keymap = ContactMap._create_keymap(contact_map1)
        contact_map2_keymap = ContactMap._create_keymap(contact_map2, altloc=True)

        # Some checks
        msg = "Error creating reliable keymap matching the sequence in ContactMap: "
        if len(contact_map1_keymap) != np.where(encoded_repr[0] != ord('-'))[0].shape[0]:
            raise ValueError(msg + contact_map1.id)
        elif len(contact_map2_keymap) != np.where(encoded_repr[1] != ord('-'))[0].shape[0]:
            raise ValueError(msg + contact_map2.id)

        # Create a sequence matching keymap including deletions and insertions
        contact_map1_keymap = ContactMap._insert_states(encoded_repr[0], contact_map1_keymap)
        contact_map2_keymap = ContactMap._insert_states(encoded_repr[1], contact_map2_keymap)

        # Reindex the altseq positions to account for insertions/deletions
        contact_map1_keymap = ContactMap._reindex_by_keymap(contact_map1_keymap)
        contact_map2_keymap = ContactMap._reindex_by_keymap(contact_map2_keymap)

        # Adjust the res_altseq based on the insertions and deletions
        contact_map2 = ContactMap._adjust(contact_map2, contact_map2_keymap)

        # Get the residue list for matching UNKNOWNs
        residues_map2 = tuple(i + 1 for i, a in enumerate(aligned_sequences_full[1].seq) if a != '-')

        # Adjust true and false positive statuses
        for contact in contact_map1:
            _id = (contact.res1_seq, contact.res2_seq)
            _id_alt = tuple(r.res_seq for r in contact_map2_keymap for i in _id if i == r.res_altseq)

            if any(i == _Gap.IDENTIFIER for i in _id_alt) and any(j not in residues_map2 for j in _id):
                contact_map1[_id].define_unknown()
            elif all(i in residues_map2 for i in _id):
                if _id_alt in contact_map2:
                    contact_map1[_id].define_match()
                else:
                    contact_map1[_id].define_mismatch()
            else:
                msg = "Error matching two contact maps - this should never happen"
                raise RuntimeError(msg)

        # ================================================================
        # 3. Remove unmatched contacts
        # ================================================================
        if remove_unmatched:
            for contactid in contact_map1.as_list():
                if contact_map1[tuple(contactid)].is_unknown:
                    contact_map1.remove(tuple(contactid))

        # ================================================================
        # 4. Renumber the contact map 1 based on contact map 2
        # ================================================================
        if renumber:
            contact_map1 = ContactMap._renumber(contact_map1, contact_map1_keymap, contact_map2_keymap)

        return contact_map1

    def reindex(self, index, altloc=False, inplace=False):
        """Re-index the :obj:`ContactMap <conkit.core.contactmap.ContactMap>`

        Parameters
        ----------
        index : int
           The new starting index [assigned to the lowest existing index in the contact map]
        altloc : bool
           Use the res_altloc positions [default: False]
        inplace : bool
           Replace the saved order of contacts [default: False]

        Returns
        -------
        obj
           The reference to the :obj:`ContactMap <conkit.core.contactmap.ContactMap>`, regardless of inplace

        Raises
        ------
        ValueError
           Index must be positive

        """
        if index < 0:
            raise ValueError("Index must be positive!")
        contact_map = self._inplace(inplace)
        if contact_map.empty:
            return contact_map
        res1s, res2s = zip(*contact_map.as_list(altloc=altloc))
        offset = min(res1s) - index
        for contact in contact_map:
            if altloc:
                contact.res1_altseq -= offset
                contact.res2_altseq -= offset
            else:
                contact.res1_seq -= offset
                contact.res2_seq -= offset
        for contact in contact_map:
            contact.id = (contact.res1_seq, contact.res2_seq)
        return contact_map

    def remove_neighbors(self, min_distance=5, max_distance=sys.maxsize, inplace=False):
        """Remove contacts between neighboring residues

        The algorithm works by keeping contact pairs that satisfy

            ``min_distance`` <= ``x`` <= ``max_distance``

        Parameters
        ----------
        min_distance : int, optional
           The minimum number of residues between contacts [default: 5]
        max_distance : int, optional
           The maximum number of residues between contacts [default: :func:``sys.maxsize``]
        inplace : bool, optional
           Replace the saved order of contacts [default: False]

        Returns
        -------
        obj
           The reference to the :obj:`ContactMap <conkit.core.contactmap.ContactMap>`, regardless of inplace

        """
        contact_map = self._inplace(inplace)
        for contactid in contact_map.as_list():
            if min_distance <= abs(contactid[1] - contactid[0]) <= max_distance:
                continue
            else:
                contact_map.remove(tuple(contactid))
        return contact_map

    def rescale(self, inplace=False):
        """Rescale the raw scores in :obj:`ContactMap <conkit.core.contactmap.ContactMap>`

        Parameters
        ----------
        inplace : bool, optional
           Replace the saved order of contacts [default: False]

        Returns
        -------
        obj
           The reference to the :obj:`ContactMap <conkit.core.contactmap.ContactMap>`, regardless of inplace

        """
        contact_map = self._inplace(inplace)

        raw_scores = np.array([c.raw_score for c in contact_map])
        norm_raw_scores = normalize(raw_scores)

        if np.isnan(norm_raw_scores).all():
            norm_raw_scores = np.where(norm_raw_scores == np.isnan, 0, 1)

        for contact, norm_raw_score in zip(contact_map, norm_raw_scores):
            contact.raw_score = norm_raw_score

        return contact_map

    def sort(self, kword, reverse=False, inplace=False):
        """Sort the :obj:`ContactMap <conkit.core.contactmap.ContactMap>`

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
        obj
           The reference to the :obj:`ContactMap <conkit.core.contactmap.ContactMap>`, regardless of inplace

        Raises
        ------
        ValueError
           ``kword`` not in :obj:`ContactMap <conkit.core.contactmap.ContactMap>`

        """
        contact_map = self._inplace(inplace)
        contact_map._sort(kword, reverse)
        return contact_map

    def to_string(self):
        """Return the :obj:`ContactMap <conkit.core.contactmap.ContactMap>` as :obj:`str`"""
        content = ["%d\t%d\t%.5f" % (c.res1_seq, c.res2_seq, c.raw_score) for c in self]
        return os.linesep.join(content)

    @staticmethod
    def _adjust(contact_map, keymap):
        """Adjust res_altseq entries to insertions and deletions"""
        encoder = dict((x.res_seq, x.res_altseq) for x in keymap if isinstance(x, _Residue))
        for contact in contact_map:
            if contact.res1_seq in encoder:
                contact.res1_altseq = encoder[contact.res1_seq]
            if contact.res2_seq in encoder:
                contact.res2_altseq = encoder[contact.res2_seq]
        return contact_map

    @staticmethod
    def _create_keymap(contact_map, altloc=False):
        """Create a simple keymap

        Parameters
        ----------
        altloc : bool
           Use the res_altloc positions [default: False]

        Returns
        -------
        list
           A list of residue mappings

        """
        contact_map_keymap = collections.OrderedDict()
        for contact in contact_map:
            pos1 = _Residue(contact.res1_seq, contact.res1_altseq, contact.res1, contact.res1_chain)
            pos2 = _Residue(contact.res2_seq, contact.res2_altseq, contact.res2, contact.res2_chain)
            if altloc:
                res1_index, res2_index = contact.res1_altseq, contact.res2_altseq
            else:
                res1_index, res2_index = contact.res1_seq, contact.res2_seq
            contact_map_keymap[res1_index] = pos1
            contact_map_keymap[res2_index] = pos2
        contact_map_keymap_sorted = sorted(list(contact_map_keymap.items()), key=lambda x: int(x[0]))
        return list(zip(*contact_map_keymap_sorted))[1]

    @staticmethod
    def _find_single(contact_map, index):
        """Find all contacts associated with ``index`` based on id property"""
        for c in contact_map:
            if c.id[0] == index or c.id[1] == index:
                yield c

    @staticmethod
    def _insert_states(sequence, keymap):
        """Create a sequence matching keymap including deletions and insertions"""
        it = iter(keymap)
        keymap_ = []
        for amino_acid in sequence:
            if amino_acid == ord('-'):
                keymap_.append(_Gap())
            else:
                keymap_.append(next(it))
        return keymap_

    @staticmethod
    def _reindex_by_keymap(keymap):
        """Reindex a key map"""
        for i, residue in enumerate(keymap):
            residue.res_altseq = i + 1
        return keymap

    @staticmethod
    def _renumber(contact_map, self_keymap, other_keymap):
        """Renumber the contact map based on the mapping of self and other keymaps"""
        for self_residue, other_residue in zip(self_keymap, other_keymap):
            if isinstance(self_residue, _Gap):
                continue
            for contact in ContactMap._find_single(contact_map, self_residue.res_seq):
                # Make sure we check with the ID, which doesn't change
                if contact.id[0] == self_residue.res_altseq:
                    contact.res1_seq = other_residue.res_seq
                    contact.res1_chain = other_residue.res_chain
                elif contact.id[1] == self_residue.res_altseq:
                    contact.res2_seq = other_residue.res_seq
                    contact.res2_chain = other_residue.res_chain
                else:
                    raise ValueError('Should never get here')

        return contact_map


def _AND(a, b):
    return a and b


def _OR(a, b):
    return a or b
