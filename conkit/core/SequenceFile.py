"""
Storage space for a sequence file
"""

from __future__ import division
from __future__ import print_function

__author__ = "Felix Simkovic"
__date__ = "07 Sep 2016"
__version__ = 0.1

from conkit.core.Entity import Entity

import numpy
try:
    import scipy.spatial
    SCIPY = True
except ImportError:
    SCIPY = False


class SequenceFile(Entity):
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
        return "SequenceFile(id=\"{0}\" nseqs={1})".format(self.id, self.nseqs)

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
        msa_mat = numpy.asarray(
            [numpy.fromstring(sequence_entry.seq, dtype=numpy.uint8) for sequence_entry in self], dtype=numpy.uint8
        )

        # Pre-define some variables
        n = msa_mat.shape[0]                        # size of the data
        batch_size = min(n, 250)                    # size of the batches
        hamming = numpy.zeros(n, dtype=numpy.int)   # storage for data

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

        msa_mat = numpy.asarray(
            [numpy.fromstring(sequence_entry.seq, dtype='uint8') for sequence_entry in self], dtype='uint8'
        )
        # matrix of 0s and 1s; 1 if char is '-'
        aa_frequencies = numpy.where(msa_mat != 45, 1, 0)
        # sum all values per row
        aa_counts = numpy.sum(aa_frequencies, axis=0)
        # divide all by sequence length
        return (aa_counts / len(msa_mat.T[0])).tolist()

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

