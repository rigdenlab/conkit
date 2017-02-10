"""
Storage space for a protein sequence
"""

__author__ = "Felix Simkovic"
__date__ = "07 Sep 2016"
__version__ = 0.1

from Bio import pairwise2
from conkit.constants import ONE_TO_THREE
from conkit.core.Entity import Entity


class Sequence(Entity):
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
        return "Sequence(id=\"{0}\" seq=\"{1}\" seq_len={2})".format(
            self.id, seq_string, len(self.seq))

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
