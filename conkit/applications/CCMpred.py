# coding=utf-8
"""
Command line object for CCMpred contact prediction application
"""

__author__ = "Felix Simkovic"
__date__ = "04 Aug 2016"
__version__ = 0.1

from Bio.Application import _Argument
from Bio.Application import _Option
from Bio.Application import _Switch
from Bio.Application import AbstractCommandline


class CCMpredCommandLine(AbstractCommandline):
    """
    Command line object for CCMpred [#]_ contact prediction application

    https://github.com/soedinglab/CCMpred

    The CCMpred program is a very fast pseudo-likelihood maximisation
    implementation of covariance detection in a Multiple Sequence
    Alignment. This wrapper allows for easy-to-use Python implementation.
    
    .. [#] Seemayer S, Gruber M, Söding J (2014). CCMpred--fast and precise
       prediction of protein residue-residue contacts from correlated mutations.
       Bioinformatics 30(21), 3128-3130.

    Examples
    --------
    To predict a contact map using a Multiple Sequence Alignment in
    JONES format, use:

    >>> from conkit.applications import CCMpredCommandLine
    >>> ccmpred_cline = CCMpredCommandLine(
    ...    alnfile="test.aln", matfile="output.mat"
    ... )
    >>> print(ccmpred_cline)
    ccmpred test.aln output.mat

    You would typically run the command line with :func:`ccmpred_cline` or via
    the Python subprocess module.

    """
    
    def __init__(self, cmd="ccmpred", **kwargs):
        self.parameters = [
            _Option(['-n', "numiter"],
                    'Compute a maximum of NUMITER operations [default: 50]',
                    equate=False),
            _Option(['-e', 'epsilon'],
                    "Set convergence criterion for minimum decrease in the "
                    "last K iterations to EPSILON [default: 0.01]",
                    equate=False),
            _Option(['-k', 'lastk'],
                    "Set K parameter for convergence criterion to LASTK [default: 5]",
                    equate=False),

            _Option(['-i', 'inifile'],
                    "Read initial weights from INIFILE",
                    filename=True,
                    equate=False),
            _Option(['-r', 'rawfile'],
                    "Store raw prediction matrix in RAWFILE",
                    filename=True,
                    equate=False),
            _Option(['-t', 'threads'],
                    'Calculate using THREADS threads on the CPU (automatically disables CUDA if available) [default: 1]',
                    equate=False),

            _Option(['-w', 'idthres'],
                     "Set sequence reweighting identity threshold to IDTHRES [default: 0.8]",
                     equate=False),
            _Option(['-l', 'lfactor'],
                     "Set pairwise regularization coefficients to LFACTOR * (L-1) [default: 0.2]",
                     equate=False),

            _Switch(['-A', 'apc'],
                    "Disable average product correction (APC)"),
            _Switch(['-R', 'renormalize'],
                     "Re-normalize output matrix to [0,1]"),

            _Argument(['alnfile'],
                      "Input alignment file [JONES format]",
                      filename=True,
                      is_required=True),
            _Argument(['matfile'],
                       "Output matrix file",
                       filename=True,
                       is_required=True),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
