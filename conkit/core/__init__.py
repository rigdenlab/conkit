# coding=utf-8
"""Core entities for hierarchy construction"""

from __future__ import division
from __future__ import print_function

__author__ = "Felix Simkovic"
__contributing_authors__ = "Jens Thomas"
__date__ = "03 Aug 2016"
__version__ = "0.2"

from Bio import pairwise2

import collections
import copy
import numpy as np
import operator

try:
    import scipy.spatial
    SCIPY = True
except ImportError:
    SCIPY = False

try:
    import sklearn.neighbors
    SKLEARN = True
except ImportError:
    SKLEARN = False


# ================================================
# Amino acid conversions
# ================================================
ONE_TO_THREE = {'A': 'ALA', 'C': 'CYS', 'B': 'ASX', 'E': 'GLU', 'D': 'ASP', 'G': 'GLY', 'F': 'PHE', 'I': 'ILE',
                'H': 'HIS', 'K': 'LYS', 'J': 'XLE', 'M': 'MET', 'L': 'LEU', 'O': 'PYL', 'N': 'ASN', 'Q': 'GLN',
                'P': 'PRO', 'S': 'SER', 'R': 'ARG', 'U': 'SEC', 'T': 'THR', 'W': 'TRP', 'V': 'VAL', 'Y': 'TYR',
                'X': 'XAA', 'Z': 'GLX'}

THREE_TO_ONE = {'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CME': 'C', 'CYS': 'C', 'GLN': 'Q', 'GLU': 'E', 
                'GLY': 'G', 'HIS': 'H', 'ILE': 'I', 'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'MSE': 'M', 'PHE': 'F', 
                'PRO': 'P', 'PYL': 'O', 'SER': 'S', 'SEC': 'U', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 'VAL': 'V', 
                'ASX': 'B', 'GLX': 'Z', 'XAA': 'X', 'UNK': 'X', 'XLE': 'J'}



class _BandwidthEstimators(object):
    """A collection of bandwidth estimators for Kernel Density Estimation"""
    @staticmethod
    def amise(p, niterations=25, eps=1e-3):

        def curvature(p, x, w):
            z = (x - p) / w
            y = (1 * (z ** 2 - 1.0) * np.exp(-0.5 * z * z) / (w * np.sqrt(2. * np.pi)) / w ** 2).sum()
            return y / p.shape[0]

        def extended_range(mn, mx, bw, ext=3):
            return mn - ext * bw, mx + ext * bw

        def optimal_bandwidth_equation(p, default_bw):
            alpha = 1. / (2. * np.sqrt(np.pi))
            sigma = 1.0
            n = p.shape[0]
            q = stiffness_integral(p, default_bw)
            return default_bw - ((n * q * sigma ** 4) / alpha) ** (-1.0 / (p.shape[1] + 4))

        def stiffness_integral(p, default_bw, eps=1e-4):
            mn, mx = extended_range(p.min(), p.max(), default_bw, ext=3)
            n = 1
            dx = (mx - mn) / n
            yy = 0.5 * dx * (curvature(p, mn, default_bw) ** 2 +
                             curvature(p, mx, default_bw) ** 2)
            # The trapezoidal rule guarantees a relative error of better than eps
            # for some number of steps less than maxn.
            maxn = (mx - mn) / np.sqrt(eps)
            # Cap the total computation spent
            maxn = 2048 if maxn > 2048 else maxn
            n = 2
            while n <= maxn:
                dx /= 2.
                y = 0
                for i in np.arange(1, n, 2):
                    y += curvature(p, mn + i * dx, default_bw) ** 2
                yy = 0.5 * yy + y * dx
                if n > 8 and abs(y * dx - 0.5 * yy) < eps * yy:
                    break
                n *= 2
            return yy

        x0 = _BandwidthEstimators.bowman(p)
        y0 = optimal_bandwidth_equation(p, x0)

        x = 0.8 * x0
        y = optimal_bandwidth_equation(p, x)

        for _ in np.arange(niterations):
            x -= y * (x0 - x) / (y0 - y)
            y = optimal_bandwidth_equation(p, x)
            if abs(y) < (eps * y0):
                break
        return x

    @staticmethod
    def bowman(p):
        sigma = np.sqrt((p ** 2).sum() / p.shape[0] - (p.sum() / p.shape[0]) ** 2)
        return sigma * ((((p.shape[1] + 2) * p.shape[0]) / 4.) ** (-1. / (p.shape[1] + 4)))

    @staticmethod
    def scott(p):
        sigma = np.minimum(
            np.std(p, axis=0, ddof=1),
            (np.percentile(p, 75) - np.percentile(p, 25)) / 1.349
        ).astype(np.float64)
        return 1.059 * sigma * p.shape[0] ** (-1. / (p.shape[1] + 4))

    @staticmethod
    def silverman(p):
        sigma = np.minimum(
            np.std(p, axis=0, ddof=1),
            (np.percentile(p, 75) - np.percentile(p, 25)) / 1.349
        ).astype(np.float64)
        return 0.9 * sigma * (p.shape[0] * (p.shape[1] + 2) / 4.) ** (-1. / (p.shape[1] + 4))


class _Gap(object):
    """A basic class representing a gap residue"""
    __slots__ = ('res_seq', 'res_altseq', 'res_name', 'res_chain')

    _IDENTIFIER = -999999

    def __init__(self):
        self.res_seq = _Gap._IDENTIFIER
        self.res_altseq = _Gap._IDENTIFIER
        self.res_name = 'X'
        self.res_chain = ''

    def __repr__(self):
        string = "{0}(res_seq='{1}' res_altseq='{2}' res_name='{3}' res_chain='{4}')"
        return string.format(
                self.__class__.__name__, self.res_seq, 
                self.res_altseq, self.res_name, self.res_chain
        )


class _Residue(object):
    """A basic class representing a residue"""
    __slots__ = ('res_seq', 'res_altseq', 'res_name', 'res_chain')

    def __init__(self, res_seq, res_altseq, res_name, res_chain):
        self.res_seq = res_seq
        self.res_altseq = res_altseq
        self.res_name = res_name
        self.res_chain = res_chain

    def __repr__(self):
        string = "{0}(res_seq='{1}' res_altseq='{2}' res_name='{3}' res_chain='{4}')"
        return string.format(
                self.__class__.__name__, self.res_seq, 
                self.res_altseq, self.res_name, self.res_chain
        )


class _Entity(object):
    """Base class for all entities used in this interface.

    It handles the storage of data. It also provides a high-efficiency
    methods to allow fast lookup and iterations of each entity. It also
    provides a hierarchical structure to remember parent and child
    entities.

    Attributes
    ----------
    id : str, list, tuple
       The ID of the selected entity
    full_id : tuple
       A traceback id including all parent classes
    parent : :obj:`Entity <conkit.core.Entity>`
       An attribute to store the reference to the parent :obj:`Entity <conkit.core.Entity>`
    child_list : list
       A list storing the child entities
    child_dict : dict
       A dictionary storing the child entities

    Notes
    -----
    It is strongly advised against the use of the :obj:`Entity <conkit.core.Entity>` class directly.
    Instead, use one or more of the the remaining data models.

    """
    __slots__ = ['_id', '_parent', '_child_list', '_child_dict']

    def __init__(self, id):
        """Initialise a generic :obj:`Entity <conkit.core.Entity>`

        Parameters
        ----------
        id : str, list, tuple
           The ID of the selected entity

        """
        self._id = None
        self._parent = None
        self._child_list = []
        self._child_dict = {}

        # Assign values post creation to use setter/getter methods
        # Possibly very bad practice but no better alternative for now
        self.id = id

    def __contains__(self, id):
        """True if there is a child element with the given id"""
        return id in self.child_dict

    def __delitem__(self, id):
        """Remove a child with given id"""
        child = self[id]
        child.parent = None
        self.child_dict.pop(id)
        self.child_list.remove(child)

    def __getitem__(self, id):
        """Return the child with the given id"""
        if isinstance(id, slice):
            indexes_to_keep = list(range(*id.indices(len(self))))
            copy_to_return = self.copy()
            for i, child in enumerate(self):
                if i not in indexes_to_keep:
                    copy_to_return.remove(child.id)
            return copy_to_return
        elif isinstance(id, int):
            return self.child_list[id]
        else:
            return self.child_dict[id]

    def __iter__(self):
        """Iterate over children"""
        for child in self.child_list:
            yield child

    def __len__(self):
        """Return the number of children"""
        return len(self.child_list)

    def __reversed__(self):
        """Reversed list of the children"""
        for child in reversed(self.child_list):
            yield child

    @property
    def child_dict(self):
        """A dictionary storing the child entities"""
        return self._child_dict

    @child_dict.setter
    def child_dict(self, child_dict):
        """Define a dictionary storing the child entities

        Parameters
        ----------
        child_dict : dict

        """
        self._child_dict = child_dict

    @property
    def child_list(self):
        """A list storing the child entities"""
        return self._child_list

    @child_list.setter
    def child_list(self, child_list):
        """Define a list storing the child entities

        Parameters
        ----------
        child_list : dict

        """
        self._child_list = child_list

    @property
    def full_id(self):
        """A traceback id including all parent classes

        The full id is a tuple containing all id's starting from
        the top object (:obj:`ContactFile <conkit.core.ContactFile>`) down to the current object.
        A full id for a :obj:`Contact <conkit.core.Contact>` e.g. is something like:
        ('1aa', 1, (1, 10))

        This corresponds to:

        :obj:`ContactFile <conkit.core.ContactFile>` identifier => 1aaa
        :obj:`ContactMap <conkit.core.ContactMap>` identifier => 1
        :obj:`Contact <conkit.core.Contact>` identifier => (1, 10)

        """
        traceback = [self.id]
        mother = self.parent
        while mother is not None:
            traceback.append(mother.id)
            mother = mother.parent
        return tuple(reversed(traceback))

    @property
    def id(self):
        """The ID of the selected entity"""
        return self._id

    @id.setter
    def id(self, id):
        """Set the ID of the selected entity

        Parameters
        ----------
        id : str, list, tuple
           The unique ID for an :obj:`Entity <conkit.core.Entity>`

        Warnings
        --------
        You cannot provide an :obj:`int` or :obj:`float` as ID.

        """
        if isinstance(id, int) or isinstance(id, float):
            raise TypeError('Please provide data type of str, list, or tuple')
        elif isinstance(id, list):
            id = tuple(id)
        self._id = id

    @property
    def parent(self):
        """An attribute to store the reference to the parent :obj:`Entity <conkit.core.Entity>`"""
        return self._parent

    @parent.setter
    def parent(self, parent):
        """Define the reference to the parent :obj:`Entity <conkit.core.Entity>`

        Parameters
        ----------
        parent : :obj:`Entity <conkit.core.Entity>`

        """
        self._parent = parent

    def _inplace(self, inplace):
        """Modify the current version using a copy

        Parameters
        ----------
        inplace : bool

        """
        if inplace:
            return self
        else:
            return self.copy()

    def _sort(self, kword, reverse):
        """Sort the :obj:`Entity <conkit.core.Entity>`"""
        if any(not hasattr(e, kword) for e in self._child_list):
            raise ValueError('Attribute not defined')
        self.child_list.sort(key=operator.attrgetter(kword), reverse=reverse)

    def add(self, entity):
        """Add a child to the :obj:`Entity <conkit.core.Entity>`

        Parameters
        ----------
        entity : :obj:`Entity <conkit.core.Entity>`

        """
        if entity.id in self:
            raise ValueError("%s defined twice" % str(entity.id))
        entity.parent = self
        self.child_list.append(entity)
        self.child_dict[entity.id] = entity

    def copy(self):
        """Create a shallow copy of :obj:`Entity <conkit.core.Entity>`"""
        shallow = copy.copy(self)

        shallow.child_list = []
        shallow.child_dict = {}
        shallow.parent = None

        for child in self:
            shallow.add(child.copy())
        return shallow

    def deepcopy(self):
        """Create a deep copy of :obj:`Entity <conkit.core.Entity>`"""
        deep = copy.deepcopy(self)

        deep.child_list = []
        deep.child_dict = {}
        deep.parent = None

        for child in self:
            deep.add(child.copy())
        return deep

    def remove(self, id):
        """Remove a child

        Parameters
        ----------
        id : str, int, list, tuple

        Warnings
        --------
        If ``id`` is of type :obj:`int`, then the :obj:`Entity <conkit.core.Entity>`
        in the ``child_list`` at index ``id`` will be deleted

        """
        del self[id]


class Contact(_Entity):
    """A contact pair template to store all associated information

    Attributes
    ----------
    distance_bound : tuple
       The lower and upper distance boundary values of a contact pair in Ångstrom [Default: 0-8Å].
    id : str
       A unique identifier
    is_match : bool
       A boolean status for the contact
    is_mismatch : bool
       A boolean status for the contact
    is_unknown : bool
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
    _MISMATCH = -1
    _MATCH = 1

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
        return "{0}(id=\"{1}\" res1=\"{2}\" res1_chain=\"{3}\" res1_seq={4} " \
               "res2=\"{5}\" res2_chain=\"{6}\" res2_seq={7} raw_score={8})".format(
            self.__class__.__name__, self.id, self.res1, self.res1_chain, 
            self.res1_seq, self.res2, self.res2_chain, self.res2_seq, self.raw_score
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
    def is_match(self):
        """A boolean status for the contact"""
        return True if self.status == Contact._MATCH else False

    @property
    def is_mismatch(self):
        """A boolean status for the contact"""
        return True if self.status == Contact._MISMATCH else False

    @property
    def is_unknown(self):
        """A boolean status for the contact"""
        return True if self.status == Contact._UNKNOWN else False

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
        if any(i == status for i in [Contact._UNKNOWN, Contact._MISMATCH, Contact._MATCH]):
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

    def define_match(self):
        """Define a contact as matching contact"""
        self._status = Contact._MATCH

    def define_mismatch(self):
        """Define a contact as mismatching contact"""
        self._status = Contact._MISMATCH

    def define_unknown(self):
        """Define a contact with unknown status"""
        self._status = Contact._UNKNOWN

    def _to_dict(self):
        """Convert the object into a dictionary"""
        keys =  ['id', 'is_match', 'is_mismatch', 'is_unknown', 'lower_bound', 'upper_bound'] \
                + [k[1:] for k in self.__slots__]
        return {k: getattr(self, k) for k in keys}

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
        elif len(_amino_acid) == 1 and _amino_acid not in list(THREE_TO_ONE.values()):
            raise ValueError(msg)
        elif len(_amino_acid) == 3 and _amino_acid not in list(THREE_TO_ONE.keys()):
            raise ValueError(msg)

        # Save the one-letter-code
        if len(_amino_acid) == 3:
            _amino_acid = THREE_TO_ONE[_amino_acid]

        return _amino_acid


class ContactMap(_Entity):
    """A contact map object representing a single prediction

    The :obj:`ContactMap <conkit.core.ContactMap>` class represents a data structure to hold a single
    contact map prediction in one place. It contains functions to store,
    manipulate and organise :obj:`Contact <conkit.core.Contact>` instances.

    Attributes
    ----------
    coverage : float
       The sequence coverage score
    id : str
       A unique identifier
    ncontacts : int
       The number of :obj:`Contact <conkit.core.Contact>` instances in the :obj:`ContactMap <conkit.core.ContactMap>`
    precision : float
       The precision (Positive Predictive Value) score
    repr_sequence : :obj:`Sequence <conkit.core.Sequence>`
       The representative :obj:`Sequence <conkit.core.Sequence>` associated with the :obj:`ContactMap <conkit.core.ContactMap>`
    repr_sequence_altloc : :obj:`Sequence <conkit.core.Sequence>`
       The representative altloc :obj:`Sequence <conkit.core.Sequence>` associated with the :obj:`ContactMap <conkit.core.ContactMap>`
    sequence : :obj:`Sequence <conkit.core.Sequence>`
       The :obj:`Sequence <conkit.core.Sequence>` associated with the :obj:`ContactMap <conkit.core.ContactMap>`
    top_contact : :obj:`Contact <conkit.core.Contact>`
       The first :obj:`Contact <conkit.core.Contact>` entry in :obj:`ContactMap <conkit.core.ContactMap>`

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
        return "{0}(id=\"{1}\", ncontacts={2})".format(
                self.__class__.__name__, self.id, self.ncontacts
        )

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
        seq_array = np.fromstring(self.repr_sequence.seq, dtype='uint8')
        gaps = np.where(seq_array == ord('-'), 1, 0)
        cov = (seq_array.size - np.sum(gaps, axis=0)) / seq_array.size
        return cov

    @property
    def ncontacts(self):
        """The number of :obj:`Contact <conkit.core.Contact>` instances in the :obj:`ContactMap <conkit.core.ContactMap>`

        Returns
        -------
        ncontacts : int
           The number of sequences in the :obj:`ContactMap <conkit.core.ContactMap>`

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
        determined by running the :func:`match` function providing a reference
        structure.

        Returns
        -------
        ppv : float
           The calculated precision score

        See Also
        --------
        coverage

        """
        # ContactMap is empty
        if len(self) == 0:
            print('ContactMap is empty')
            return 0.0
        
        unique, counts = np.unique(np.asarray([c.status for c in self]), return_counts=True)
        cdict = dict(zip(unique, counts))
        fp_count = cdict[Contact._MISMATCH] if Contact._MISMATCH in cdict else 0.0
        uk_count = cdict[Contact._UNKNOWN] if Contact._UNKNOWN in cdict else 0.0
        tp_count = cdict[Contact._MATCH] if Contact._MATCH in cdict else 0.0
        
        if fp_count == 0.0 and tp_count == 0.0:
            print("No matches or mismatches found in your contact map. Match two ContactMaps first.")
            return 0.0
        elif uk_count > 0:
            print("Some contacts between the ContactMaps are unmatched due to non-identical "
                  "sequences. The precision value might be inaccurate.")

        return tp_count / (tp_count + fp_count)

    @property
    def repr_sequence(self):
        """The representative :obj:`Sequence <conkit.core.Sequence>` associated with the :obj:`ContactMap <conkit.core.ContactMap>`

        The peptide sequence constructed from the available
        contacts using the normal res_seq positions

        Returns
        -------
        sequence : :obj:`conkit.coreSequence`

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
        res1_seqs, res2_seqs = list(zip(*[contact.id for contact in self]))
        res_seqs = set(
            sorted(res1_seqs + res2_seqs)
        )
        return self._construct_repr_sequence(list(res_seqs))

    @property
    def repr_sequence_altloc(self):
        """The representative altloc :obj:`Sequence <conkit.core.Sequence>` associated with the :obj:`ContactMap <conkit.core.ContactMap>`

        The peptide sequence constructed from the available
        contacts using the altloc res_seq positions

        Returns
        -------
        sequence : :obj:`Sequence <conkit.core.Sequence>`

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
        res1_seqs, res2_seqs = list(zip(*[(contact.res1_altseq, contact.res2_altseq) for contact in self]))
        res_seqs = set(
            sorted(res1_seqs + res2_seqs)
        )
        return self._construct_repr_sequence(list(res_seqs))

    @property
    def sequence(self):
        """The :obj:`Sequence <conkit.core.Sequence>` associated with the :obj:`ContactMap <conkit.core.ContactMap>`

        Returns
        -------
        :obj:`Sequence <conkit.core.Sequence>`

        See Also
        --------
        repr_sequence, repr_sequence_altloc

        """
        return self._sequence

    @sequence.setter
    def sequence(self, sequence):
        """Associate a :obj:`Sequence <conkit.core.Sequence>` instance with the :obj:`ContactMap <conkit.core.ContactMap>`

        Parameters
        ----------
        sequence : :obj:`Sequence <conkit.core.Sequence>`

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
        """The first :obj:`Contact <conkit.core.Contact>` entry in :obj:`ContactMap <conkit.core.ContactMap>`

        Returns
        -------
        top_contact : :obj:`Contact <conkit.core.Contact>`, None
           The first :obj:`Contact <conkit.core.Contact>` entry in :obj:`ContactFile <conkit.core.ContactFile>`

        """
        if len(self) > 0:
            return self[0]
        else:
            return None

    def _construct_repr_sequence(self, res_seqs):
        """Construct the representative sequence"""
        # Determine which are present and which are not
        representative_sequence = ''
        for i in range(1, self.sequence.seq_len + 1):
            if i in res_seqs:
                representative_sequence += self.sequence.seq[i - 1]
            else:
                representative_sequence += '-'
        return Sequence(self.sequence.id + '_repr', representative_sequence)

    def assign_sequence_register(self, altloc=False):
        """Assign the amino acids from :obj:`Sequence <conkit.core.Sequence>` to all :obj:`Contact <conkit.core.Contact>` instances

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

    def calculate_jaccard_index(self, other):
        """Calculate the Jaccard index between two :obj:`ContactMap <conkit.core.ContactMap>` instances

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
        other : :obj:`ContactMap <conkit.core.ContactMap>`
           A ConKit :obj:`ContactMap <conkit.core.ContactMap>`
        
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
        intersection = np.sum([1 for contact in self if contact.id in other])
        union = len(self) + np.sum([1 for contact in other if contact.id not in self])
        # If self and other are both empty, we define J(x,y) = 1
        if union == 0:
            return 1.0
        return float(intersection) / union

    def calculate_kernel_density(self, bw_method="amise"):
        """Calculate the contact density in the contact map using Gaussian kernels

        Various algorithms can be used to estimate the bandwidth. To calculate the
        bandwidth for an 1D data array ``X`` with ``n`` data points and ``d`` dimensions,
        the listed algorithms have been implemented. Please note, in rules 2 and 3, the
        value of :math:`\\sigma` is the smaller of the standard deviation of ``X`` or
        the normalized interquartile range.

        1. Asymptotic Mean Integrated Squared Error (AMISE)

           This particular choice of bandwidth recovers all the important features whilst maintaining smoothness. 
           It is a direct implementation of the method used by [#]_.


        2. Bowman & Azzalini [#]_ implementation

        .. math::

           \\sqrt{\\frac{\\sum{X}^2}{n}-(\\frac{\\sum{X}}{n})^2}*(\\frac{(d+2)*n}{4})^\\frac{-1}{d+4}

        3. Scott's [#]_ implementation

        .. math::

           1.059*\\sigma*n^\\frac{-1}{d+4}
        4. Silverman's [#]_ implementation

        .. math::

           0.9*\\sigma*(n*\\frac{d+2}{4})^\\frac{-1}{d+4}

        .. [#] Sadowski, M.I. (2013). Prediction of protein domain boundaries from inverse covariances.
        .. [#] Bowman, A.W. & Azzalini, A. (1997). Applied Smoothing Techniques for Data Analysis.
        .. [#] Scott, D.W. (1992). Multivariate Density Estimation: Theory, Practice, and Visualization.
        .. [#] Silverman, B.W. (1986). Density Estimation for Statistics and Data Analysis.

        Parameters
        ----------
        bw_method : str, optional
           The bandwidth estimator to use [default: amise_sekant]

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
        if not SKLEARN:
            raise RuntimeError('Cannot find SciKit package')

        # Compute the relevant data we need
        x = np.asarray([i for c in self for i in np.arange(c.res1_seq, c.res2_seq)])[:, np.newaxis]
        x_fit = np.linspace(x.min(), x.max() + 1, x.max() - x.min() + 1)[:, np.newaxis]

        # Obtain the bandwidth as defined by user method
        if bw_method == "amise":
            bandwidth = _BandwidthEstimators.amise(x)
        elif bw_method == "bowman":
            bandwidth = _BandwidthEstimators.bowman(x)
        elif bw_method == "scott":
            bandwidth = _BandwidthEstimators.scott(x)
        elif bw_method == "silverman":
            bandwidth = _BandwidthEstimators.silverman(x)
        else:
            msg = "Undefined bandwidth method: {0}".format(bw_method)
            raise ValueError(msg)

        # Estimate the Kernel Density using original data and fit random sample
        kde = sklearn.neighbors.KernelDensity(bandwidth=bandwidth).fit(x)
        return np.exp(kde.score_samples(x_fit)).tolist()

    def calculate_scalar_score(self):
        """Calculate a scaled score for the :obj:`ContactMap <conkit.core.ContactMap>`

        This score is a scaled score for all raw scores in a contact
        map. It is defined by the formula

        .. math::

           {x}'=\\frac{x}{\\overline{d}}

        where :math:`x` corresponds to the raw score of each predicted
        contact and :math:`\overline{d}` to the mean of all raw scores.

        The score is saved in a separate :obj:`Contact <conkit.core.Contact>` attribute called
        ``scalar_score``

        This score is described in more detail in [#]_.

        .. [#] S. Ovchinnikov, L. Kinch, H. Park, Y. Liao, J. Pei, D.E. Kim,
           H. Kamisetty, N.V. Grishin, D. Baker (2015). Large-scale determination
           of previously unsolved protein structures using evolutionary information.
           *Elife* **4**, e09248.

        """
        raw_scores = np.asarray([c.raw_score for c in self])
        sca_scores = raw_scores / np.mean(raw_scores)
        for contact, sca_score in zip(self, sca_scores):
            contact.scalar_score = sca_score
        return

    def find(self, indexes, altloc=False):
        """Find all contacts associated with ``index``

        Parameters
        ----------
        index : list, tuple
           A list of residue indexes to find
        altloc : bool
           Use the res_altloc positions [default: False]

        Returns
        -------
        :obj:`ContactMap <conkit.core.ContactMap>`
           A modified version of the contact map containing
           the found contacts

        """
        contact_map = self.copy()
        for contact in self:
            if altloc and (contact.res1_altseq in indexes or contact.res2_altseq in indexes):
                continue
            elif contact.res1_seq in indexes or contact.res2_seq in indexes:
                continue
            else:
                contact_map.remove(contact.id)
        return contact_map

    def match(self, other, remove_unmatched=False, renumber=False, inplace=False):
        """Modify both hierarchies so residue numbers match one another.

        This function is key when plotting contact maps or visualising
        contact maps in 3-dimensional space. In particular, when residue
        numbers in the structure do not start at count 0 or when peptide
        chain breaks are present.

        Parameters
        ----------
        other : :obj:`ContactMap <conkit.core.ContactMap>`
           A ConKit :obj:`ContactMap <conkit.core.ContactMap>`
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
            :obj:`ContactMap <conkit.core.ContactMap>` instance, regardless of inplace

        Raises
        ------
        ValueError
           Error creating reliable keymap matching the sequence in :obj:`ContactMap <conkit.core.ContactMap>`

        """
        contact_map1 = self._inplace(inplace)
        contact_map2 = other._inplace(inplace)

        # ================================================================
        # 1. Align all sequences
        # ================================================================

        # Align both full sequences against each other
        aligned_sequences_full = contact_map1.sequence.align_local(contact_map2.sequence, id_chars=2,
                                                                   nonid_chars=1, gap_open_pen=-0.5,
                                                                   gap_ext_pen=-0.1)
        contact_map1_full_sequence, contact_map2_full_sequence = aligned_sequences_full

        # Align contact map 1 full sequences with representative sequence
        aligned_sequences_map1 = contact_map1_full_sequence.align_local(contact_map1.repr_sequence,
                                                                        id_chars=2, nonid_chars=1, gap_open_pen=-0.5,
                                                                        gap_ext_pen=-0.2, inplace=True)
        contact_map1_repr_sequence = aligned_sequences_map1[-1]

        # Align contact map 2 full sequences with __ALTLOC__ representative sequence
        aligned_sequences_map2 = contact_map2_full_sequence.align_local(contact_map2.repr_sequence_altloc,
                                                                        id_chars=2, nonid_chars=1, gap_open_pen=-0.5,
                                                                        gap_ext_pen=-0.2, inplace=True)
        contact_map2_repr_sequence = aligned_sequences_map2[-1]

        # Align both aligned representative sequences
        aligned_sequences_repr = contact_map1_repr_sequence.align_local(contact_map2_repr_sequence,
                                                                        id_chars=2, nonid_chars=1, gap_open_pen=-1.0,
                                                                        gap_ext_pen=-0.5, inplace=True)
        contact_map1_repr_sequence, contact_map2_repr_sequence = aligned_sequences_repr

        # ================================================================
        # 2. Identify TPs in other, map them, and match them to self
        # ================================================================

        # Encode the sequences to uint8 character arrays for easier and faster handling
        encoded_repr = np.asarray([
            np.fromstring(contact_map1_repr_sequence.seq, dtype='uint8'),
            np.fromstring(contact_map2_repr_sequence.seq, dtype='uint8')
        ])

        # Create mappings for both contact maps
        contact_map1_keymap = ContactMap._create_keymap(contact_map1)
        contact_map2_keymap = ContactMap._create_keymap(contact_map2, altloc=True)

        # Some checks
        msg = "Error creating reliable keymap matching the sequence in ContactMap: "
        if len(contact_map1_keymap) != np.where(encoded_repr[0] != ord('-'))[0].shape[0]:
            msg += contact_map1.id
            raise ValueError(msg)
        elif len(contact_map2_keymap) != np.where(encoded_repr[1] != ord('-'))[0].shape[0]:
            msg += contact_map2.id
            raise ValueError(msg)

        # Create a sequence matching keymap including deletions and insertions
        contact_map1_keymap = ContactMap._insert_states(encoded_repr[0], contact_map1_keymap)
        contact_map2_keymap = ContactMap._insert_states(encoded_repr[1], contact_map2_keymap)

        # Reindex the altseq positions to account for insertions/deletions
        contact_map1_keymap = ContactMap._reindex(contact_map1_keymap)
        contact_map2_keymap = ContactMap._reindex(contact_map2_keymap)

        # Adjust the res_altseq based on the insertions and deletions
        contact_map2 = ContactMap._adjust(contact_map2, contact_map2_keymap)

        # Get the residue list for matching UNKNOWNs
        residues_map2 = tuple(i+1 for i, a in enumerate(aligned_sequences_full[1].seq) if a != '-')

        # Adjust true and false positive statuses
        for contact in contact_map1:
            id = (contact.res1_seq, contact.res2_seq)
            id_alt = tuple(r.res_seq for r in contact_map2_keymap for i in id if i == r.res_altseq)

            if any(i == _Gap._IDENTIFIER for i in id_alt) and any(j not in residues_map2 for j in id):
                contact_map1[id].define_unknown()
            elif all(i in residues_map2 for i in id):
                if id_alt in contact_map2:
                    contact_map1[id].define_match()
                else:
                    contact_map1[id].define_mismatch()
            else:
                msg = "Error matching two contact maps - this should never happen"
                raise RuntimeError(msg)

        # ================================================================
        # 3. Remove unmatched contacts
        # ================================================================
        if remove_unmatched:
            for contact in contact_map1.copy():
                if contact.is_unknown:
                    contact_map1.remove(contact.id)

        # ================================================================
        # 4. Renumber the contact map 1 based on contact map 2
        # ================================================================
        if renumber:
            contact_map1 = ContactMap._renumber(contact_map1, contact_map1_keymap, contact_map2_keymap)

        return contact_map1

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
        contact_map : :obj:`ContactMap <conkit.core.ContactMap>`
           The reference to the :obj:`ContactMap <conkit.core.ContactMap>`, regardless of inplace

        """
        contact_map = self._inplace(inplace)
        for contact in contact_map.copy():
            if abs(contact.res1_seq - contact.res2_seq) < min_distance:
                contact_map.remove(contact.id)
        return contact_map

    def rescale(self, inplace=False):
        """Rescale the raw scores in :obj:`ContactMap <conkit.core.ContactMap>`

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
        contact_map : :obj:`ContactMap <conkit.core.ContactMap>`
           The reference to the :obj:`ContactMap <conkit.core.ContactMap>`, regardless of inplace

        """
        contact_map = self._inplace(inplace)

        raw_scores = np.asarray([c.raw_score for c in contact_map])
        norm_raw_scores = (raw_scores - raw_scores.min()) / (raw_scores.max() - raw_scores.min())

        # Important to not end up with raw scores == np.nan
        if np.isnan(norm_raw_scores).all():
            norm_raw_scores = np.where(norm_raw_scores == np.isnan, 0, 1)

        for contact, norm_raw_score in zip(contact_map, norm_raw_scores):
            contact.raw_score = norm_raw_score

        return contact_map

    def sort(self, kword, reverse=False, inplace=False):
        """Sort the :obj:`ContactMap <conkit.core.ContactMap>`

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
        contact_map : :obj:`ContactMap <conkit.core.ContactMap>`
           The reference to the :obj:`ContactMap <conkit.core.ContactMap>`, regardless of inplace

        Raises
        ------
        ValueError
           ``kword`` not in :obj:`ContactMap <conkit.core.ContactMap>`

        """
        contact_map = self._inplace(inplace)
        contact_map._sort(kword, reverse)
        return contact_map

    @staticmethod
    def _adjust(contact_map, keymap):
        """Adjust res_altseq entries to insertions and deletions"""
        encoder = dict((x.res_seq, x.res_altseq) for x in keymap if isinstance(x, _Residue))
        for contact in contact_map:
            if contact.res1_seq in list(encoder.keys()):
                contact.res1_altseq = encoder[contact.res1_seq]
            if contact.res2_seq in list(encoder.keys()):
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
    def _reindex(keymap):
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


class ContactFile(_Entity):
    """A contact file object representing a single prediction file

    The contact file class represents a data structure to hold all predictions
    with a single contact map file. It contains functions to store,
    manipulate and organise contact maps.

    Attributes
    ----------
    author : str
       The author of the :obj:`ContactFile <conkit.core.ContactFile>`
    method : list, str
       The :obj:`ContactFile <conkit.core.ContactFile>`-specific method
    remark : list, str
       The :obj:`ContactFile <conkit.core.ContactFile>`-specific remarks
    target : str
       The target name
    top_map : :obj:`ContactMap <conkit.core.ContactMap>`
       The first :obj:`ContactMap <conkit.core.ContactMap>` entry in :obj:`ContactFile <conkit.core.ContactFile>`

    Examples
    --------
    >>> from conkit.core import ContactMap, ContactFile
    >>> contact_file = ContactFile("example")
    >>> contact_file.add(ContactMap("foo"))
    >>> contact_file.add(ContactMap("bar"))
    >>> print(contact_file)
    ContactFile(id="example" nseqs=2)

    """
    __slots__ = ['_author', '_method', '_remark', '_target']

    def __init__(self, id):
        """Initialise a new contact map

        Parameters
        ----------
        id : str
           A unique identifier for the contact file

        """
        self._author = None
        self._method = []
        self._remark = []
        self._target = None
        super(ContactFile, self).__init__(id)

    def __repr__(self):
        return "{0}(id=\"{1}\" nmaps={2})".format(
                self.__class__.__name__, self.id, len(self)
        )

    @property
    def author(self):
        """The author of the :obj:`ContactFile <conkit.core.ContactFile>`"""
        return self._author

    @author.setter
    def author(self, author):
        """Define the author of the :obj:`ContactFile <conkit.core.ContactFile>`

        Parameters
        ----------
        author : str

        """
        self._author = author

    @property
    def method(self):
        """The :obj:`ContactFile <conkit.core.ContactFile>`-specific method"""
        return self._method

    @method.setter
    def method(self, method):
        """Set the :obj:`ContactFile <conkit.core.ContactFile>` method

        Parameters
        ----------
        method : str, list
           The method will be added to the list of methods

        """
        if isinstance(method, list):
            self._method += method
        elif isinstance(method, tuple):
            self._method += list(method)
        else:
            self._method += [method]

    @property
    def remark(self):
        """The :obj:`ContactFile <conkit.core.ContactFile>`-specific remarks"""
        return self._remark

    @remark.setter
    def remark(self, remark):
        """Set the :obj:`ContactFile <conkit.core.ContactFile>` remark

        Parameters
        ----------
        remark : str, list
           The remark will be added to the list of remarks

        """
        if isinstance(remark, list):
            self._remark += remark
        elif isinstance(remark, tuple):
            self._remark += list(remark)
        else:
            self._remark += [remark]

    @property
    def target(self):
        """The target name"""
        return self._target

    @target.setter
    def target(self, target):
        """Define the target name

        Parameters
        ----------
        target : str

        """
        self._target = target

    @property
    def top_map(self):
        """The first :obj:`ContactMap <conkit.core.ContactMap>` entry in :obj:`ContactFile <conkit.core.ContactFile>`

        Returns
        -------
        top_map : :obj:`ContactMap <conkit.core.ContactMap>`, None
           The first :obj:`ContactMap <conkit.core.ContactMap>` entry in :obj:`ContactFile <conkit.core.ContactFile>`

        """
        if len(self) > 0:
            return self[0]
        else:
            return None

    def sort(self, kword, reverse=False, inplace=False):
        """Sort the :obj:`ContactFile <conkit.core.ContactFile>`

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
        contact_map : :obj:`ContactMap <conkit.core.ContactMap>`
           The reference to the :obj:`ContactMap <conkit.core.ContactMap>`, regardless of inplace

        Raises
        ------
        ValueError
           ``kword`` not in :obj:`ContactFile <conkit.core.ContactFile>`

        """
        contact_file = self._inplace(inplace)
        contact_file._sort(kword, reverse)
        return contact_file


class Sequence(_Entity):
    """A sequence template to store all associated information

    Attributes
    ----------
    id : str
       A unique identifier
    remark : list
       The :obj:`Sequence <conkit.core.Sequence>`-specific remarks
    seq : str
       The protein sequence as :obj:`str`
    seq_len : int
       The protein sequence length

    Examples
    --------
    >>> from conkit.core import Sequence
    >>> sequence_entry = Sequence("example", "ABCDEF")
    >>> print(sequence_entry)
    Sequence(id="example" seq="ABCDEF" seqlen=6)

    """
    __slots__ = ['_remark', '_seq']

    def __init__(self, id, seq):
        """Initialise a generic sequence 

        Parameters
        ----------
        id : str
           A unique sequence identifier
        seq : str
           The protein sequence

        """
        self._remark = []
        self._seq = None

        # Assign values post creation to use setter/getter methods
        # Possibly very bad practice but no better alternative for now
        self.seq = seq

        super(Sequence, self).__init__(id)

    def __add__(self, other):
        """Concatenate two sequence instances to a new"""
        id = self.id + '_' + other.id
        seq = self.seq + other.seq
        return Sequence(id, seq)

    def __repr__(self):
        if self.seq_len > 12:
            seq_string = self.seq[:5] + '...' + self.seq[-5:]
        else:
            seq_string = self.seq
        return "{0}(id=\"{1}\" seq=\"{2}\" seq_len={3})".format(
            self.__class__.__name__, self.id, seq_string, len(self.seq)
        )

    @property
    def remark(self):
        """The :obj:`Sequence <conkit.core.Sequence>`-specific remarks"""
        return self._remark

    @remark.setter
    def remark(self, remark):
        """Set the :obj:`Sequence <conkit.core.Sequence>` remark

        Parameters
        ----------
        remark : str, list
           The remark will be added to the list of remarks

        """
        if isinstance(remark, list):
            self._remark += remark
        elif isinstance(remark, tuple):
            self._remark += list(remark)
        else:
            self._remark += [remark]

    @property
    def seq(self):
        """The protein sequence as :obj:`str`"""
        return self._seq

    @seq.setter
    def seq(self, seq):
        """Set the sequence

        Parameters
        ----------
        seq : str

        Raises
        ------
        ValueError
           One or more amino acids in the sequence are not recognised

        """
        if any(c not in list(ONE_TO_THREE.keys()) for c in seq.upper() if c != '-'):
            raise ValueError('Unrecognized amino acids in sequence')
        self._seq = seq

    @property
    def seq_len(self):
        """The protein sequence length"""
        return len(self.seq)

    def align_global(self, other, id_chars=2, nonid_chars=1, gap_open_pen=-0.5, gap_ext_pen=-0.1, inplace=False):
        """Generate a global alignment between two :obj:`Sequence <conkit.core.Sequence>` instances

        Parameters
        ----------
        other : :obj:`Sequence <conkit.core.Sequence>`
        id_chars : int, optional
        nonid_chars : int, optional
        gap_open_pen : float, optional
        gap_ext_pen : float, optional
        inplace : bool, optional
           Replace the saved order of residues [default: False]

        Returns
        -------
        :obj:`Sequence <conkit.core.Sequence>`
           The reference to the :obj:`Sequence`, regardless of inplace
        :obj:`Sequence <conkit.core.Sequence>`
           The reference to the :obj:`Sequence`, regardless of inplace

        """
        sequence1 = self._inplace(inplace)
        sequence2 = other._inplace(inplace)

        alignment = pairwise2.align.globalms(
            sequence1.seq, sequence2.seq, id_chars, nonid_chars, gap_open_pen, gap_ext_pen
        )

        sequence1.seq = alignment[-1][0]
        sequence2.seq = alignment[-1][1]

        return sequence1, sequence2

    def align_local(self, other, id_chars=2, nonid_chars=1, gap_open_pen=-0.5, gap_ext_pen=-0.1, inplace=False):
        """Generate a local alignment between two :obj:`Sequence <conkit.core.Sequence>` instances

        Parameters
        ----------
        other : :obj:`Sequence <conkit.core.Sequence>`
        id_chars : int, optional
        nonid_chars : int, optional
        gap_open_pen : float, optional
        gap_ext_pen : float, optional
        inplace : bool, optional
           Replace the saved order of residues [default: False]

        Returns
        -------
        :obj:`Sequence <conkit.core.Sequence>`
           The reference to the :obj:`Sequence <conkit.core.Sequence>`, regardless of inplace
        :obj:`Sequence <conkit.core.Sequence>`
           The reference to the :obj:`Sequence <conkit.core.Sequence>`, regardless of inplace

        """
        sequence1 = self._inplace(inplace)
        sequence2 = other._inplace(inplace)

        alignment = pairwise2.align.localms(
            sequence1.seq, sequence2.seq, id_chars, nonid_chars, gap_open_pen, gap_ext_pen
        )

        sequence1.seq = alignment[-1][0]
        sequence2.seq = alignment[-1][1]

        return sequence1, sequence2


class SequenceFile(_Entity):
    """A sequence file object representing a single sequence file

    The :obj:`SequenceFile <conkit.core.SequenceFile>` class represents a data structure to hold
    :obj:`Sequence <conkit.core.Sequence>` instances in a single sequence file. It contains
    functions to store and analyze sequences.

    Attributes
    ----------
    id : str
       A unique identifier
    is_alignment : bool
       A boolean status for the alignment
    nseqs : int
       The number of sequences in the :obj:`SequenceFile <conkit.core.SequenceFile>`
    remark : list
       The :obj:`SequenceFile <conkit.core.SequenceFile>`-specific remarks
    status : int
       An indication of the sequence file, i.e alignment, no alignment, or unknown
    top_sequence : :obj:`Sequence <conkit.core.Sequence>`, None
       The first :obj:`Sequence <conkit.core.Sequence>` entry in the file


    Examples
    --------
    >>> from conkit.core import Sequence, SequenceFile
    >>> sequence_file = SequenceFile("example")
    >>> sequence_file.add(Sequence("foo", "ABCDEF"))
    >>> sequence_file.add(Sequence("bar", "ZYXWVU"))
    >>> print(sequence_file)
    SequenceFile(id="example" nseqs=2)

    """
    __slots__ = ['_remark', '_status']

    _UNKNOWN = 0
    _NO_ALIGNMENT = -1
    _YES_ALIGNMENT = 1

    def __init__(self, id):
        """Initialise a new :obj:`SequenceFile <conkit.core.SequenceFile>`

        Parameters
        ----------
        id : str
           A unique identifier for the sequence file

        """
        self._remark = []
        self._status = SequenceFile._UNKNOWN
        super(SequenceFile, self).__init__(id)

    def __repr__(self):
        return "{0}(id=\"{1}\" nseqs={2})".format(
                self.__class__.__name__, self.id, self.nseqs
        )

    @property
    def is_alignment(self):
        """A boolean status for the alignment

        Returns
        -------
        bool
           A boolean status for the alignment

        """
        seq_length = self.top_sequence.seq_len
        self.status = SequenceFile._YES_ALIGNMENT
        for sequence in self:
            if sequence.seq_len != seq_length:
                self.status = SequenceFile._NO_ALIGNMENT
                break
        return True if self.status == SequenceFile._YES_ALIGNMENT else False

    @property
    def nseqs(self):
        """The number of :obj:`Sequence <conkit.core.Sequence>` instances
        in the :obj:`SequenceFile <conkit.core.SequenceFile>`

        Returns
        -------
        int
           The number of sequences in the :obj:`SequenceFile <conkit.core.SequenceFile>`

        """
        return len(self)

    @property
    def remark(self):
        """The :obj:`SequenceFile <conkit.core.SequenceFile>`-specific remarks"""
        return self._remark

    @remark.setter
    def remark(self, remark):
        """Set the :obj:`SequenceFile <conkit.core.SequenceFile>` remark

        Parameters
        ----------
        remark : str, list
           The remark will be added to the list of remarks

        """
        if isinstance(remark, list):
            self._remark += remark
        elif isinstance(remark, tuple):
            self._remark += list(remark)
        else:
            self._remark += [remark]

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
           [0] for `unknown`, [-1] for `no alignment`, or [1] for `alignment`

        Raises
        ------
        ValueError
           Cannot determine if your sequence file is an alignment or not

        """
        if any(i == status for i in [SequenceFile._UNKNOWN, SequenceFile._NO_ALIGNMENT, SequenceFile._YES_ALIGNMENT]):
            self._status = status
        else:
            raise ValueError("Cannot determine if your sequence file is an alignment or not")

    @property
    def top_sequence(self):
        """The first :obj:`Sequence <conkit.core.Sequence>` entry in :obj:`SequenceFile <conkit.core.SequenceFile>`

        Returns
        -------
        :obj:`Sequence <conkit.core.Sequence>`, None
           The first :obj:`Sequence <conkit.core.Sequence>` entry in :obj:`SequenceFile <conkit.core.SequenceFile>`

        """
        if len(self) > 0:
            return self[0]
        else:
            return None

    def calculate_meff(self, identity=0.7):
        """Calculate the number of effective sequences

        This function calculates the number of effective
        sequences (`Meff`) in the Multiple Sequence Alignment.

        The mathematical function used to calculate `Meff` is

        .. math::

           M_{eff}=\\sum_{i}\\frac{1}{\\sum_{j}S_{i,j}}

        Parameters
        ----------
        identity : float, optional
           The sequence identity to use for similarity decision [default: 0.7]

        Returns
        -------
        int
           The number of effective sequences

        Raises
        ------
        MemoryError
           Too many sequences in the alignment for Hamming distance calculation
        RuntimeError
           SciPy package not installed
        ValueError
           :obj:`SequenceFile <conkit.core.SequenceFile>` is not an alignment
        ValueError
           Sequence Identity needs to be between 0 and 1

        """
        if not SCIPY:
            raise RuntimeError('Cannot find SciPy package')

        if not self.is_alignment:
            raise ValueError('This is not an alignment')

        if identity < 0 or identity > 1:
            raise ValueError("Sequence Identity needs to be between 0 and 1")

        # Alignment to unsigned integer matrix
        msa_mat = np.asarray(
            [np.fromstring(sequence_entry.seq, dtype=np.uint8) for sequence_entry in self], dtype=np.uint8
        )

        # Pre-define some variables
        n = msa_mat.shape[0]                        # size of the data
        batch_size = min(n, 250)                    # size of the batches
        hamming = np.zeros(n, dtype=np.int)   # storage for data

        # Separate the distance calculations into batches to avoid MemoryError exceptions.
        # This answer was provided by a StackOverflow user. The corresponding suggestion by
        # user @WarrenWeckesser: http://stackoverflow.com/a/41090953/3046533
        num_full_batches, last_batch = divmod(n, batch_size)
        batches = [batch_size] * num_full_batches
        if last_batch != 0:
            batches.append(last_batch)
        for k, batch in enumerate(batches):
            i = batch_size * k
            dists = scipy.spatial.distance.cdist(msa_mat[i:i+batch], msa_mat, metric='hamming')
            hamming[i:i+batch] = (dists < (1 - identity)).sum(axis=1)

        return (1. / hamming).sum().astype(int).item()

    def calculate_freq(self):
        """Calculate the gap frequency in each alignment column

        This function calculates the frequency of gaps at each
        position in the Multiple Sequence Alignment.

        Returns
        -------
        list
           A list containing the per alignment-column amino acid
           frequency count

        Raises
        ------
        MemoryError
           Too many sequences in the alignment
        RuntimeError
           :obj:`SequenceFile <conkit.core.SequenceFile>` is not an alignment

        """
        if not self.is_alignment:
            raise ValueError('This is not an alignment')

        msa_mat = np.asarray(
            [np.fromstring(sequence_entry.seq, dtype='uint8') for sequence_entry in self], dtype='uint8'
        )
        # matrix of 0s and 1s; 1 if char is '-'
        aa_frequencies = np.where(msa_mat != 45, 1, 0)
        # sum all values per row
        aa_counts = np.sum(aa_frequencies, axis=0)
        # divide all by sequence length
        return (aa_counts / len(msa_mat[:, 0])).tolist()

    def sort(self, kword, reverse=False, inplace=False):
        """Sort the :obj:`SequenceFile <conkit.core.SequenceFile>`

        Parameters
        ----------
        kword : str
           The dictionary key to sort sequences by
        reverse : bool, optional
           Sort the sequences in reverse order [default: False]
        inplace : bool, optional
           Replace the saved order of sequences [default: False]

        Returns
        -------
        :obj:`SequenceFile <conkit.core.SequenceFile>`
           The reference to the :obj:`SequenceFile <conkit.core.SequenceFile>`, regardless of inplace

        Raises
        ------
        ValueError
           ``kword`` not in :obj:`SequenceFile <conkit.core.SequenceFile>`

        """
        sequence_file = self._inplace(inplace)
        sequence_file._sort(kword, reverse)
        return sequence_file

    def trim(self, start, end, inplace=False):
        """Trim the :obj:`SequenceFile <conkit.core.SequenceFile>`

        Parameters
        ----------
        start : int
           First residue to include
        end : int
           Final residue to include
        inplace : bool, optional
           Replace the saved order of sequences [default: False]

        Returns
        -------
        :obj:`SequenceFile <conkit.core.SequenceFile>`
           The reference to the :obj:`SequenceFile <conkit.core.SequenceFile>`, regardless of inplace

        """
        sequence_file = self._inplace(inplace)

        if not self.is_alignment:
            raise ValueError('This is not an alignment')

        i, j = start-1, end
        for sequence in sequence_file:
            sequence.seq = sequence.seq[i:j]

        return sequence_file

