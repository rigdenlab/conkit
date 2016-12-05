# coding=utf-8
"""
Storage space for a contact pair
"""

__author__ = "Felix Simkovic"
__date__ = "03 Aug 2016"
__version__ = 0.1

from conkit import constants
from conkit.core.Entity import Entity


class Contact(Entity):
    """A contact pair template to store all associated information

    Attributes
    ----------
    distance_bound : tuple
       The lower and upper distance boundary values of a contact pair in Ångstrom [Default: 0-8Å].
    id : str
       A unique identifier
    is_false_positive : bool
       A boolean status for the contact
    is_true_positive : bool
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
       The raw_score scaled according to the average ``raw_score``
    status : int
       An indication of the residue status, i.e true positive, false positive, or unknown
    upper_bound : int
       The upper distance boundary value
    weight : float
       A separate internal weight factor for the contact pair

    Examples
    --------
    >>> from conkit.core import Contact
    >>> contact = Contact(1, 25, 1.0)
    >>> print(contact)
    Contact(id="(1, 25)" res1="A" res1_seq=1 res2="A" res2_seq=25 raw_score=1.0)

    """
    __slots__ = ['_distance_bound', '_raw_score', '_res1', '_res2', '_res1_chain', '_res2_chain',
                 '_res1_seq', '_res2_seq', '_res1_altseq', '_res2_altseq', '_scalar_score',
                 '_status', '_weight']

    _UNKNOWN = 0
    _FALSE_POSITIVE = -1
    _TRUE_POSITIVE = 1

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
        self._distance_bound = [0, 8]
        self._raw_score = 1.0
        self._res1 = 'X'
        self._res2 = 'X'
        self._res1_chain = ''
        self._res2_chain = ''
        self._res1_seq = 0
        self._res2_seq = 0
        self._res1_altseq = 0
        self._res2_altseq = 0
        self._scalar_score = 0.0
        self._status = Contact._UNKNOWN
        self._weight = 1.0

        # Assign values post creation to use setter/getter methods
        # Possibly very bad practice but no better alternative for now
        self.distance_bound = distance_bound
        self.raw_score = raw_score
        self.res1_seq = res1_seq
        self.res2_seq = res2_seq

        super(Contact, self).__init__((res1_seq, res2_seq))

    def __repr__(self):
        return "Contact(id=\"{0}\" res1=\"{1}\" res1_chain=\"{2}\" res1_seq={3} " \
               "res2=\"{4}\" res2_chain=\"{5}\" res2_seq={6} raw_score={7})".format(
            self.id, self.res1, self.res1_chain, self.res1_seq,
            self.res2, self.res2_chain, self.res2_seq, self.raw_score
        )

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
        if isinstance(distance_bound, tuple):
            self._distance_bound = list(distance_bound)
        elif isinstance(distance_bound, list):
            self._distance_bound = distance_bound
        else:
            raise TypeError("Data of type list or tuple required")

    @property
    def is_false_positive(self):
        """A boolean status for the contact

        Parameters
        ----------
        is_false_positive : bool
           True / False

        See Also
        --------
        define_false_positive, define_true_positive, is_true_positive

        """
        return True if self.status == Contact._FALSE_POSITIVE else False

    @property
    def is_true_positive(self):
        """A boolean status for the contact

        Parameters
        ----------
        is_true_positive : bool
           True / False

        See Also
        --------
        define_true_positive, define_false_positive, is_false_positive

        """
        return True if self.status == Contact._TRUE_POSITIVE else False

    @property
    def lower_bound(self):
        """The lower distance boundary value"""
        return self.distance_bound[0]

    @lower_bound.setter
    def lower_bound(self, value):
        """Set the lower distance boundary value

        Parameters
        ----------
        value : int

        Raises
        ------
        ValueError
           Lower bound must be positive
        ValueError
           Lower bound must be smaller than upper bound

        """
        if value < 0:
            raise ValueError('Lower bound must be positive')
        elif value >= self.upper_bound:
            raise ValueError('Lower bound must be smaller than upper bound')
        self._distance_bound[0] = value

    @property
    def upper_bound(self):
        """The upper distance boundary value"""
        return self.distance_bound[1]

    @upper_bound.setter
    def upper_bound(self, value):
        """Set the upper distance boundary value

        Parameters
        ----------
        value : int

        Raises
        ------
        ValueError
           Upper bound must be positive
        ValueError
           Upper bound must be larger than lower bound

        """
        if value < 0:
            raise ValueError('Upper bound must be positive')
        elif value <= self.lower_bound:
            raise ValueError('Upper bound must be larger than lower bound')
        self._distance_bound[1] = value

    @property
    def raw_score(self):
        """The prediction score for the contact pair"""
        return self._raw_score

    @raw_score.setter
    def raw_score(self, score):
        """Define the raw score

        Parameters
        ----------
        score : float

        """
        self._raw_score = float(score)

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
        # Keep this statement in case we get a float
        if not isinstance(index, int):
            raise TypeError('Data type int required for res_seq')
        self._res1_altseq = index

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
        # Keep this statement in case we get a float
        if not isinstance(index, int):
            raise TypeError('Data type int required for res_seq')
        self._res2_altseq = index

    @property
    def res1_chain(self):
        """The chain for residue 1"""
        return self._res1_chain

    @res1_chain.setter
    def res1_chain(self, chain):
        """Define the chain for residue 1

        Parameters
        ----------
        chain : str

        """
        self._res1_chain = chain

    @property
    def res2_chain(self):
        """The chain for residue 2"""
        return self._res2_chain

    @res2_chain.setter
    def res2_chain(self, chain):
        """Define the chain for residue 2

        Parameters
        ----------
        chain : str

        """
        self._res2_chain = chain

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

        """
        # Keep this statement in case we get a float
        if not isinstance(index, int):
            raise TypeError('Data type int required for res_seq')
        self._res1_seq = index

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

        """
        # Keep this statement in case we get a float
        if not isinstance(index, int):
            raise TypeError('Data type int required for res_seq')
        self._res2_seq = index

    @property
    def scalar_score(self):
        """The raw_score scaled according to the average :obj:`raw_score`"""
        return self._scalar_score

    @scalar_score.setter
    def scalar_score(self, score):
        """Set the scalar score

        Parameters
        ----------
        score : float

        """
        self._scalar_score = float(score)

    @property
    def status(self):
        """An indication of the residue status, i.e true positive, false positive, or unknown"""
        return self._status

    @status.setter
    def status(self, status):
        """Set the status

        Parameters
        ----------
        status : int
           [0] for `unknown`, [-1] for `false positive`, or [1] for `true positive`

        Raises
        ------
        ValueError
           Unknown status

        """
        if any(i == status for i in [Contact._UNKNOWN, Contact._FALSE_POSITIVE, Contact._TRUE_POSITIVE]):
            self._status = status
        else:
            raise ValueError("Unknown status")

    @property
    def weight(self):
        """A separate internal weight factor for the contact pair"""
        return self._weight

    @weight.setter
    def weight(self, weight):
        """Define a separate internal weight factor for the contact pair

        Parameters
        ----------
        weight : float, int

        """
        self._weight = float(weight)

    def define_false_positive(self):
        """Define a contact as false positive

        See Also
        --------
        define_true_positive, is_true_positive, is_false_positive

        """
        self._status = Contact._FALSE_POSITIVE

    def define_true_positive(self):
        """Define a contact as true positive

        See Also
        --------
        define_false_positive, is_false_positive, is_true_positive

        """
        self._status = Contact._TRUE_POSITIVE

    def _to_dict(self):
        """Convert the object into a dictionary"""
        return {
            'id': self.id,
            'is_false_positive': self.is_false_positive,
            'is_true_positive': self.is_true_positive,
            'distance_bound': self.distance_bound,
            'lower_bound': self.lower_bound,
            'upper_bound': self.upper_bound,
            'raw_score': self.raw_score,
            'res1': self.res1,
            'res2': self.res2,
            'res1_chain': self.res1_chain,
            'res2_chain': self.res2_chain,
            'res1_seq': self.res1_seq,
            'res2_seq': self.res2_seq,
            'res1_altseq': self.res1_altseq,
            'res2_altseq': self.res2_altseq,
            'scalar_score': self.scalar_score,
            'status': self.status,
            'weight': self.weight,
        }

    @staticmethod
    def _set_residue(amino_acid):
        """Assign the residue to the corresponding amino_acid"""

        # Check that the amino acid exists
        msg = "Unknown amino acid: {0}".format(amino_acid)

        # Keep if statements separate to avoid type error for int and str.upper()
        if not isinstance(amino_acid, str):
            raise ValueError(msg)

        _amino_acid = amino_acid.upper()
        if not (len(_amino_acid) == 1 or len(_amino_acid) == 3):
            raise ValueError(msg)
        elif len(_amino_acid) == 1 and _amino_acid not in list(constants.THREE_TO_ONE.values()):
            raise ValueError(msg)
        elif len(_amino_acid) == 3 and _amino_acid not in list(constants.THREE_TO_ONE.keys()):
            raise ValueError(msg)

        # Save the one-letter-code
        if len(_amino_acid) == 3:
            _amino_acid = constants.THREE_TO_ONE[_amino_acid]

        return _amino_acid
