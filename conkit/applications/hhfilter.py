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
Command line object for HHfilter Multiple Sequence Alignment filtering application
"""

__author__ = "Felix Simkovic"
__date__ = "05 Aug 2016"
__version__ = "0.1"

from Bio.Application import _Option
from Bio.Application import AbstractCommandline


class HHfilterCommandline(AbstractCommandline):
    """
    Command line object for HHfilter [#]_ [#]_ alignment filter application

    https://toolkit.tuebingen.mpg.de/hhfilter

    Filter an alignment by maximum sequence identity of match states and minimum coverage.

    .. [#] Alva V., Nam SZ., Söding J., Lupas AN. (2016). The MPI bioinformatics Toolkit as an
       integrative platform for advanced protein sequence and structure analysis. Nucleic Acids Res. pii: gkw348.

    .. [#] Remmert M., Biegert A., Hauser A., Söding J. (2011). HHblits: Lightning-fast iterative
       protein sequence searching by HMM-HMM alignment. Nat Methods. 9(2):173-5.

    Examples
    --------
    To generate a Multiple Sequence Alignment, use:

    >>> from conkit.applications import HHfilterCommandline
    >>> hhfilter_cline = HHfilterCommandline(
    ...     input='test.a3m', output='test.filtered.a3m'
    ... )
    >>> print(hhfilter_cline)
    hhfilter -i test.a3m -o test.filtered.a3m

    You would typically run the command line with :func:`hhfilter_cline` or via
    the :mod:`~subprocess` module.

    """

    def __init__(self, cmd='hhfilter', **kwargs):
        self.parameters = [
            _Option(
                ['-i', 'input'],
                'read input file in A3M/A2M or FASTA format',
                filename=True,
                is_required=True,
                equate=False),
            _Option(
                ['-o', 'output'], 'write to output file in A3M format', filename=True, is_required=True, equate=False),
            _Option(['-a', 'append_output'], 'append to output file in A3M format', filename=True, equate=False),

            # Options
            _Option(
                ['-v', 'verbose'],
                'verbose mode: 0:no screen output  1:only warings  2: verbose [default: 2]',
                equate=False),
            _Option(['-id', 'pairwise_identity'], 'maximum pairwise sequence identity [default: 90]', equate=False),
            _Option(
                ['-diff', 'diversity'],
                'filter MSAs by selecting most diverse set of sequences, keeping '
                'at least this many seqs in each MSA block of length 50 [default: 1000]',
                equate=False),
            _Option(['-cov', 'coverage'], 'minimum coverage with master sequence (%) [default: 0]', equate=False),
            _Option(
                ['-qid', 'query_identity'],
                'minimum sequence identity with master sequence (%) [default: 0]',
                equate=False),
            _Option(
                ['-qsc', 'per_column_score'],
                'minimum score per column with master sequence [default: -20.0]',
                equate=False),
            _Option(
                ['-neff', 'nr_effective_sequences'],
                'target diversity of multiple sequence alignment [default: off]',
                equate=False),

            # # Input alignment options
            # _Option(['-M', 'a2m'],
            #         'use A2M/A3M input alignment format',
            #         equate=False),
            # _Option(['-M', 'fasta'],
            #          'use FASTA input alignment format',
            #         equate=False),
            # _Option(['-M', 'match_states'],
            #         'use FASTA: columns with fewer than X% gaprs are match states',
            #         equate=False),
        ]

        AbstractCommandline.__init__(self, cmd, **kwargs)
