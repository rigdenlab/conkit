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
Command line object for PSICOV contact prediction application
"""

__author__ = "Felix Simkovic"
__date__ = "04 Aug 2016"
__version__ = "0.1"

from Bio.Application import _Argument
from Bio.Application import _Option
from Bio.Application import _Switch
from Bio.Application import AbstractCommandline


class PsicovCommandline(AbstractCommandline):
    """
    Command line object for PSICOV [#]_ contact prediction application

    http://bioinfadmin.cs.ucl.ac.uk/downloads/PSICOV/

    The PSICOV program is a Accurate Contact Prediction from large
    protein alignments.

    .. [#] Jones, D.T., Buchan, D.W., Cozzetto, D. & Pontil, M. (2012). PSICOV:
       Precise structural contact prediction using sparse inverse covariance
       estimation on large multiple sequence alignments. Bioinformatics. 28, 184-190.

    Examples
    --------
    To predict a contact map using a Multiple Sequence Alignment in
    JONES format, use:

    >>> from conkit.applications import PsicovCommandline
    >>> psicov_cline = PsicovCommandline(alnfile="test.aln")
    >>> print(ccmpred_cline)
    psicov test.aln

    You would typically run the command line with :func:`psicov_cline` or via
    the :mod:`~subprocess` module.

    """

    def __init__(self, cmd='psicov', **kwargs):
        self.parameters = [
            _Switch(['-a', 'lasso'], "use approximate Lasso algorithm"),
            _Switch(['-n', 'noshrink'], "don't pre-shrink the sample covariance matrix"),
            _Switch(['-f', 'filter'], "filer low-scoring contacts"),
            _Switch(['-p', 'ppv_output'], "output PPV estimates rather than raw scores"),
            _Switch(['-l', 'noapc'], "don't apply APC to Lasso output"),
            _Option(['-r', 'rho'], "set initial rho paramter", equate=False),
            _Option(
                ['-d', 'sparsity'], "set target precision matrix sparsity [default: 0; not specified]", equate=False),
            _Option(['-t', 'convergence_threshold'], "set Lasso convergence threshold [default: 1e-4]", equate=False),
            _Option(
                ['-i', 'blosum_weighting'],
                "select BLOSUM-like weighting with given identity threshold "
                "[default selects threshold automatically]",
                equate=False),
            _Option(['-c', 'pseudocount'], "set pseudocount value [default: 1]", equate=False),
            _Option(['-j', 'sequence_separation'], "set minimum sequence sparation [default: 5]", equate=False),
            _Option(['-g', 'gap_fraction'], "set maximum fraction of gaps [default: 0.9]", equate=False),
            _Option(['-z', 'nr_threads'], "set maximum number of threads", equate=False),
            _Option(['-b', 'rho_parameter_file'], "read rho parameter file", filename=True, equate=False),
            _Argument(['alnfile'], "Input alignment file [JONES format]", filename=True, is_required=True),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
