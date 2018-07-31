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
"""Contact container used throughout ConKit"""

from __future__ import division
from __future__ import print_function

__author__ = "Felix Simkovic"
__date__ = "03 Aug 2016"
__version__ = "1.0"

from enum import Enum, unique
from conkit.core.entity import Entity
from conkit.core.mappings import AminoAcidOneToThree, AminoAcidThreeToOne, ContactMatchState
from conkit.misc import deprecate


class Contact(Entity):
    """A contact pair template to store all associated information

    Examples
    --------
    >>> from conkit.core import Contact
    >>> contact = Contact(1, 25, 1.0)
    >>> print(contact)
    Contact(id="(1, 25)" res1="A" res1_seq=1 res2="A" res2_seq=25 raw_score=1.0)

    Attributes
    ----------
    distance_bound : tuple
       The lower and upper distance boundary values of a contact pair in Ångstrom [Default: 0-8Å].
    id : str
       A unique identifier
    true_positive : bool
       A boolean status for the contact
    true_negative : bool
       A boolean status for the contact
    false_positive : bool
       A boolean status for the contact
    false_negative : bool
       A boolean status for the contact
    status_unknown : bool
       A boolean status for the contact
    lower_bound : int
       The lower distance boundary value
    raw_score : float
       The prediction score for the contact pair
    res1 : str
       The amino acid of residue 1 [default: X]
    res2 : str
       The amino acid of residue 2 [default: X]
    res1_chain : str
       The chain for residue 1
    res2_chain : str
       The chain for residue 2
    res1_seq : int
       The residue sequence number of residue 1
    res2_seq : int
       The residue sequence number of residue 2
    res1_altseq : int
       The alternative residue sequence number of residue 1
    res2_altseq : int
       The alternative residue sequence number of residue 2
    scalar_score : float
       The :attr:`~conkit.core.contact.Contact.raw_score` scaled according to its average
    status : int
       An indication of the residue status
    upper_bound : int
       The upper distance boundary value
    weight : float
       A separate internal weight factor for the contact pair

    """
    __slots__ = [
        '_distance_bound', 'raw_score', '_res1', '_res2', 'res1_chain', 'res2_chain', '_res1_seq', '_res2_seq',
        '_res1_altseq', '_res2_altseq', 'scalar_score', '_status', 'weight'
    ]

    def __init__(self, res1_seq, res2_seq, raw_score, distance_bound=(0, 8)):
        """Initialize a generic contact pair

        Parameters
        ----------
        distance_bound : tuple, optional
           The lower and upper distance boundary values of a contact pair in Ångstrom.
           Default is set to between 0.0 and 8.0 Å.
        raw_score : float
           The covariance score for the contact pair
        res1_seq : int
           The residue sequence number of residue 1
        res2_seq : int
           The residue sequence number of residue 2

        """
        self.raw_score = raw_score
        self.res1_chain = ''
        self.res2_chain = ''
        self.scalar_score = 0.0
        self.weight = 1.0

        self._distance_bound = [0.0, 8.0]
        self._res1 = 'X'
        self._res2 = 'X'
        self._res1_seq = 0
        self._res2_seq = 0
        self._res1_altseq = 0
        self._res2_altseq = 0
        self._status = ContactMatchState.unknown

        self.distance_bound = distance_bound
        self.res1_seq = res1_seq
        self.res2_seq = res2_seq

        super(Contact, self).__init__((res1_seq, res2_seq))

    def __repr__(self):
        text = "{name}(id={id} res1={_res1} res1_chain={res1_chain} res1_seq={_res1_seq} " \
               "res2={_res2} res2_chain={res2_chain} res2_seq={_res2_seq} raw_score={raw_score})"
        return text.format(
            name=self.__class__.__name__, id=self._id, **{k: getattr(self, k)
                                                          for k in self.__class__.__slots__})

    @property
    def distance_bound(self):
        """The lower and upper distance boundary values of a contact pair in Ångstrom [Default: 0-8Å]."""
        return tuple(self._distance_bound)

    @distance_bound.setter
    def distance_bound(self, distance_bound):
        """Define the lower and upper distance boundary value

        Parameters
        ----------
        distance_bound : list, tuple
           A 2-element list/tuple with a lower and upper distance boundary value

        """
        if isinstance(distance_bound, tuple) or isinstance(distance_bound, list):
            self._distance_bound = list(map(float, distance_bound))
        else:
            raise TypeError("Data of type list or tuple required")

    @property
    @deprecate('0.11', msg='Use true_positive instead')
    def is_match(self):
        return self._status == ContactMatchState.true_positive

    @property
    @deprecate('0.11', msg='Use false_positive instead')
    def is_mismatch(self):
        return self._status == ContactMatchState.false_positive

    @property
    @deprecate('0.11', msg='Use status_unknown instead')
    def is_unknown(self):
        return self._status == ContactMatchState.unknown

    @property
    def lower_bound(self):
        """The lower distance boundary value"""
        return self.distance_bound[0]

    @lower_bound.setter
    def lower_bound(self, value):
        """Set the lower distance boundary value

        Parameters
        ----------
        value : int, float 

        Raises
        ------
        :exc:`ValueError`
           :attr:`~conkit.core.contact.Contact.lower_bound` must be positive
        :exc:`ValueError`
           :attr:`~conkit.core.contact.Contact.lower_bound` must be smaller than
           :attr:`~conkit.core.contact.Contact.upper_bound`

        """
        if 0 < value < self.upper_bound:
            self._distance_bound[0] = float(value)
        else:
            raise ValueError('Lower bound must be positive and smaller than upper bound')

    @property
    def upper_bound(self):
        """The upper distance boundary value"""
        return self.distance_bound[1]

    @upper_bound.setter
    def upper_bound(self, value):
        """Set the upper distance boundary value

        Parameters
        ----------
        value : int, float

        Raises
        ------
        :exc:`ValueError`
           :attr:`~conkit.core.contact.Contact.upper_bound` must be positive
        :exc:`ValueError`
           :attr:`~conkit.core.contact.Contact.upper_bound` must be larger than
           :attr:`~conkit.core.contact.Contact.lower_bound`

        """
        if 0 < value > self.lower_bound:
            self._distance_bound[1] = float(value)
        else:
            raise ValueError('Upper bound must be positive and larger than lower bound')

    @property
    def res1(self):
        """The amino acid of residue 1 [default: X]"""
        return self._res1

    @res1.setter
    def res1(self, amino_acid):
        """Define the amino acid of residue 1

        Parameters
        ----------
        amino_acid : str
           The one- or three-letter code of an amino acid

        """
        self._res1 = Contact._set_residue(amino_acid)

    @property
    def res2(self):
        """The amino acid of residue 2 [default: X]"""
        return self._res2

    @res2.setter
    def res2(self, amino_acid):
        """Define the amino acid of residue 2

        Parameters
        ----------
        amino_acid : str
           The one- or three-letter code of an amino acid

        """
        self._res2 = Contact._set_residue(amino_acid)

    @property
    def res1_altseq(self):
        """The alternative residue sequence number of residue 1"""
        return self._res1_altseq

    @res1_altseq.setter
    def res1_altseq(self, index):
        """Define the alternative residue 1 sequence index

        Parameters
        ----------
        index : int

        """
        if isinstance(index, int):
            self._res1_altseq = index
        else:
            raise TypeError('Data type int required for res_seq')

    @property
    def res2_altseq(self):
        """The alternative residue sequence number of residue 2"""
        return self._res2_altseq

    @res2_altseq.setter
    def res2_altseq(self, index):
        """Define the alternative residue 2 sequence index

        Parameters
        ----------
        index : int

        """
        if isinstance(index, int):
            self._res2_altseq = index
        else:
            raise TypeError('Data type int required for res_seq')

    @property
    def res1_seq(self):
        """The residue sequence number of residue 1"""
        return self._res1_seq

    @res1_seq.setter
    def res1_seq(self, index):
        """Define residue 1 sequence index

        Parameters
        ----------
        index : int

        Raises
        ------
        :exc:`TypeError`
           Data type :obj:`int` required for :attr:`~conkit.core.contact.Contact.res1_seq`

        """
        if isinstance(index, int):
            self._res1_seq = index
        else:
            raise TypeError('Data type int required for res_seq')

    @property
    def res2_seq(self):
        """The residue sequence number of residue 2"""
        return self._res2_seq

    @res2_seq.setter
    def res2_seq(self, index):
        """Define residue 2 sequence index

        Parameters
        ----------
        index : int

        Raises
        ------
        :exc:`TypeError`
           Data type :obj:`int` required for :attr:`~conkit.core.contact.Contact.res2_seq`

        """
        if isinstance(index, int):
            self._res2_seq = index
        else:
            raise TypeError('Data type int required for res_seq')

    @property
    def status(self):
        """An indication of the residue status"""
        return self._status.value

    @status.setter
    def status(self, status):
        """Set the status

        Parameters
        ----------
        status : int, :obj:`~conkit.core.mappings.ContactMatchState`
           [0] for :attr:`~conkit.core.mappings.ContactMatchState.unknown`,
           [1] for :attr:`~conkit.core.mappings.ContactMatchState.true_positive`,
           [2] for :attr:`~conkit.core.mappings.ContactMatchState.true_negative`,
           [3] for :attr:`~conkit.core.mappings.ContactMatchState.false_positive`,
           [4] for :attr:`~conkit.core.mappings.ContactMatchState.false_negative`,

        Raises
        ------
        :exc:`ValueError`
           Not a valid :obj:`~conkit.core.mappings.ContactMatchState`

        """
        self._status = ContactMatchState(status)

    @property
    def true_positive(self):
        return self._status == ContactMatchState.true_positive

    @true_positive.setter
    def true_positive(self, is_tp):
        if is_tp:
            self._status = ContactMatchState.true_positive
        else:
            self.status_unknown = True

    @property
    def true_negative(self):
        return self._status == ContactMatchState.true_negative

    @true_negative.setter
    def true_negative(self, is_tn):
        if is_tn:
            self._status = ContactMatchState.true_negative
        else:
            self.status_unknown = True

    @property
    def false_positive(self):
        return self._status == ContactMatchState.false_positive

    @false_positive.setter
    def false_positive(self, is_fp):
        if is_fp:
            self._status = ContactMatchState.false_positive
        else:
            self.status_unknown = True

    @property
    def false_negative(self):
        return self._status == ContactMatchState.false_negative

    @false_negative.setter
    def false_negative(self, is_fn):
        if is_fn:
            self._status = ContactMatchState.false_negative
        else:
            self.status_unknown = True

    @property
    def status_unknown(self):
        return self._status == ContactMatchState.unknown

    @status_unknown.setter
    def status_unknown(self, is_unknown):
        if is_unknown:
            self._status = ContactMatchState.unknown
        else:
            raise ValueError("Choose one of true_positive, false_positive, true_negative, false_negative instead!")

    @deprecate('0.11', msg='Use true_positive instead')
    def define_match(self):
        """Define a contact as matching contact"""
        self._status = ContactMatchState.true_positive

    @deprecate('0.11', msg='Use false_positive instead')
    def define_mismatch(self):
        """Define a contact as mismatching contact"""
        self._status = ContactMatchState.false_positive

    @deprecate('0.11', msg='Use status_unknown instead')
    def define_unknown(self):
        """Define a contact with unknown status"""
        self._status = ContactMatchState.unknown

    def _to_dict(self):
        """Convert the object into a dictionary"""
        keys = ['id', 'true_positive', 'false_positive', 'status_unknown', 'lower_bound', 'upper_bound'] \
                + [k for k in self.__slots__]

        dict_ = {}
        for k in keys:
            if k[0] == '_':
                k = k[1:]
            dict_[k] = getattr(self, k)
        return dict_

    @staticmethod
    def _set_residue(amino_acid):
        """Assign the residue to the corresponding amino_acid"""
        a_a = amino_acid.upper()
        if a_a in AminoAcidOneToThree.__members__:
            return a_a
        elif a_a in AminoAcidThreeToOne.__members__:
            return AminoAcidThreeToOne[a_a].value
        else:
            raise ValueError("Unknown amino acid: {} (assert all is uppercase!)".format(amino_acid))
