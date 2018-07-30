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
"""
Command line object for CCMpred contact prediction application
"""

__author__ = "Felix Simkovic"
__date__ = "04 Aug 2016"
__version__ = "0.1"

from Bio.Application import _Argument
from Bio.Application import _Option
from Bio.Application import _Switch
from Bio.Application import AbstractCommandline


class CCMpredCommandline(AbstractCommandline):
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

    >>> from conkit.applications import CCMpredCommandline
    >>> ccmpred_cline = CCMpredCommandline(
    ...    alnfile="test.aln", matfile="output.mat"
    ... )
    >>> print(ccmpred_cline)
    ccmpred test.aln output.mat

    You would typically run the command line with :func:`ccmpred_cline` or via
    the :mod:`~subprocess` module.

    """

    def __init__(self, cmd="ccmpred", **kwargs):
        self.parameters = [
            _Option(['-n', "numiter"], 'Compute a maximum of NUMITER operations [default: 50]', equate=False),
            _Option(
                ['-e', 'epsilon'],
                "Set convergence criterion for minimum decrease in the "
                "last K iterations to EPSILON [default: 0.01]",
                equate=False),
            _Option(['-k', 'lastk'], "Set K parameter for convergence criterion to LASTK [default: 5]", equate=False),
            _Option(['-i', 'inifile'], "Read initial weights from INIFILE", filename=True, equate=False),
            _Option(['-r', 'rawfile'], "Store raw prediction matrix in RAWFILE", filename=True, equate=False),
            _Option(
                ['-t', 'threads'],
                'Calculate using THREADS threads on the CPU (automatically disables CUDA if available) [default: 1]',
                equate=False),
            _Option(
                ['-w', 'idthres'],
                "Set sequence reweighting identity threshold to IDTHRES [default: 0.8]",
                equate=False),
            _Option(
                ['-l', 'lfactor'],
                "Set pairwise regularization coefficients to LFACTOR * (L-1) [default: 0.2]",
                equate=False),
            _Switch(['-A', 'apc'], "Disable average product correction (APC)"),
            _Switch(['-R', 'renormalize'], "Re-normalize output matrix to [0,1]"),
            _Argument(['alnfile'], "Input alignment file [JONES format]", filename=True, is_required=True),
            _Argument(['matfile'], "Output matrix file", filename=True, is_required=True),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
