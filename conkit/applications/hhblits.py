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
Command line object for HHblits Multiple Sequence Alignment application
"""

__author__ = "Felix Simkovic"
__date__ = "05 Aug 2016"
__version__ = "0.1"

from Bio.Application import _Option
from Bio.Application import _Switch
from Bio.Application import AbstractCommandline


class HHblitsCommandline(AbstractCommandline):
    """
    Command line object for HHblits [#]_ [#]_ alignment generation

    https://toolkit.tuebingen.mpg.de/hhblits

    The HHblits program is a homology detection tool by iterative HMM-HMM comparison.

    .. [#] Alva V., Nam SZ., Söding J., Lupas AN. (2016). The MPI bioinformatics Toolkit as an
       integrative platform for advanced protein sequence and structure analysis. Nucleic Acids Res. pii: gkw348.

    .. [#] Remmert M., Biegert A., Hauser A., Söding J. (2011). HHblits: Lightning-fast iterative
       protein sequence searching by HMM-HMM alignment. Nat Methods. 9(2):173-5.

    Examples
    --------
    To generate a Multiple Sequence Alignment, use:

    >>> from conkit.applications import HHblitsCommandline
    >>> hhblits_cline = HHblitsCommandline(
    ...     input="test.fasta", database="uniprot20_29Feb2012"
    ... )
    >>> print(hhblits_cline)
    hhblits -i test.fasta -d uniprot20_29Feb2012

    You would typically run the command line with :func:`hhblits_cline` or via
    the :mod:`~subprocess` module.

    """

    def __init__(self, cmd="hhblits", **kwargs):

        # TODO: Figure out how to do mutual groups
        if 'local' in list(kwargs.keys()) and 'global' in list(kwargs.keys()):
            raise ValueError("Use only one of \"global_aln/local_aln\" alignment modes")

        self.parameters = [
            _Option(
                ['-i', 'input'],
                'single sequence or multiple sequence alignment in '
                'a3m, a2m, or FASTA format, or HMM in hmm format',
                filename=True,
                is_required=True,
                equate=False),

            # Options
            _Option(['-d', 'database'], 'database name (e.g. uniprot20_29Feb2012)', is_required=True, equate=False),
            _Option(['-n', 'niterations'], 'number of iterations [default: 2]', equate=False),
            _Option(
                ['-e', 'evalue'], 'E-value cutoff for inclusion in result alignment [default: 0.001]', equate=False),

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

            # Output options
            _Option(
                ['-o', 'output'],
                'write results in standard format to file [default: <infile.hhr>]',
                filename=True,
                equate=False),
            _Option(
                ['-oa3m', 'oa3m'],
                'write result MSA with significant matches in a3m format',
                filename=True,
                equate=False),
            _Option(
                ['-ohhm', 'ohhm'],
                'write result MSA with significant matches in hmm format',
                filename=True,
                equate=False),
            _Option(
                ['-opsi', 'opsi'],
                'write result MSA with significant matches in psi format',
                filename=True,
                equate=False),
            _Option(['-oalis', 'oalis'], 'write MSAs in A3M format after each iteration', filename=True, equate=False),

            # Filter options applied to query MSA, database MSAs, and result MSA
            _Switch(['-all', 'show_all'], 'show all sequences in result MSA; do not filter result MSA'),
            _Option(['-id', 'id'], 'maximum pairwise sequence identity [default: 90]', equate=False),
            _Option(
                ['-diff', 'diff'],
                'filter MSAs by selecting most diverse set of sequences, keeping '
                'at least this many seqs in each MSA block of length 50 [default: 1000]',
                equate=False),
            _Option(['-cov', 'cov'], 'minimum coverage with master sequence (%) [default: 0]', equate=False),
            _Option(['-qid', 'qid'], 'minimum sequence identity with master sequence (%) [default: 0]', equate=False),
            _Option(['-qsc', 'qsc'], 'minimum score per column with master sequence [default: -20.0]', equate=False),
            _Option(['-neff', 'neff'], 'target diversity of multiple sequence alignment [default: off]', equate=False),

            # HMM-HMM alignment options
            _Switch(['-norealign', 'norealign'], 'do NOT realign displayed hits with MAC algorithm [default: realign]'),
            _Option(
                ['-mact', 'mac_realignment_threshold'],
                'posterior probability threshold for MAC re-alignment [default: 0.350], '
                'Parameter controls alignment greediness: 0:global >0.1:local',
                equate=False),
            _Switch(['-glob', 'global_aln'], 'use global alignment mode for searching/ranking [default: local]'),
            _Switch(['-loc', 'loca_alnl'], 'use local alignment mode for searching/ranking [default: local]'),

            # Other options
            _Option(
                ['-v', 'verbose'],
                'verbose mode: 0:no screen output  1:only warings  2: verbose [default: 2]',
                equate=False),
            _Option(
                ['-neffmax', 'neffmax'],
                'skip further search iterations when diversity Neff of query '
                'MSA becomes larger than neffmax [default: 10.0]',
                equate=False),
            _Option(['-cpu', 'cpu'], 'number of CPUs to use (for shared memory SMPs) [default: 2]'),

            # Extra options from `-h all`
            _Option(
                ['-maxfilt', 'maxfilt'],
                'max number of hits allowed to pass 2nd prefilter (default=20000)',
                equate=False),
        ]

        AbstractCommandline.__init__(self, cmd, **kwargs)
