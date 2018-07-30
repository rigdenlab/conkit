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

from Bio.Application import _Option
from Bio.Application import AbstractCommandline


class CdhitCommandline(AbstractCommandline):
    """
    Command line object for Cd-hit [#]_ [#]_

    http://cd-hit.org

    CD-HIT is a very widely used program for clustering and comparing
    protein or nucleotide sequences. CD-HIT was originally developed
    by Dr. Weizhong Li at Dr. Adam Godzik's Lab at the Burnham Institute
    (now Sanford-Burnham Medical Research Institute).

    CD-HIT is very fast and can handle extremely large databases. CD-HIT
    helps to significantly reduce the computational and manual efforts in
    many sequence analysis tasks and aids in understanding the data
    structure and correct the bias within a dataset.
    
    .. [#] Li W, Jaroszewski L, Godzik A(2001). Clustering of highly homologous sequences
       to reduce thesize of large protein database. Bioinformatics 17, 282-283.

    .. [#] Li W, Jaroszewski L, Godzik A (2002). Tolerating some redundancy significantly
       speeds up clustering of large protein databases. Bioinformatics 18, 77-82.


    Examples
    --------
    >>> from conkit.applications import CdhitCommandline
    >>> cdhit_cline = CdhitCommandline()
    >>> print(cdhit_cline)

    You would typically run the command line with :func:`cdhit_cline` or via
    the :mod:`~subprocess` module.

    """

    def __init__(self, cmd="cd-hit", **kwargs):
        self.parameters = [
            _Option(
                ['-i', 'input'],
                'input filename in fasta format, required',
                filename=True,
                equate=False,
                is_required=True),
            _Option(['-o', 'output'], 'output filename, required', filename=True, equate=False, is_required=True),
            _Option(
                ['-c', 'seq_id_thres'],
                "sequence identity threshold, default 0.9 "
                "this is the default cd-hit's 'global sequence identity' calculated as: "
                "number of identical amino acids in alignment divided by "
                "the full length of the shorter sequence",
                equate=False),
            _Option(
                ['-G', 'global_seq_id'],
                "use global sequence identity, default 1 "
                "if set to 0, then use local sequence identity, calculated as : "
                "number of identical amino acids in alignment "
                "divided by the length of the alignment "
                "NOTE!!! don't use -G 0 unless you use alignment coverage controls "
                "see options -aL (kwarg: `cov_alignment_long`), -AL (kwarg: `cov_alignment_long_control`),"
                "            -aS (kwarg: `cov_alignment_short`), -AS (kwarg: `cov_alignment_short_control`)",
                equate=False),
            _Option(['-b', 'band_width'], 'band_width of alignment, default 20', equate=False),
            _Option(
                ['-M', 'memory_limit'],
                'memory limit (in MB) for the program, default 800; 0 for unlimited',
                equate=False),
            _Option(['-T', 'num_threads'], 'number of threads, default 1; with 0, all CPUs will be used', equate=False),
            _Option(['-n', 'word_length'], "word_length, default 5, see user's guide for choosing it", equate=False),
            _Option(['-l', 'len_throw_away_seqs'], "length of throw_away_sequences, default 10", equate=False),
            _Option(['-t', 'tol_4_redundance'], "tolerance for redundance, default 2", equate=False),
            _Option(
                ['-d', 'len_desc'],
                "length of description in .clstr file, default 20 "
                "if set to 0, it takes the fasta defline and stops at first space "
                "-s	length difference cutoff, default 0.0",
                equate=False),
            _Option(
                ['-s', 'len_diff_cutoff'],
                "length difference cutoff, default 0.0 "
                "if set to 0.9, the shorter sequences need to be "
                "at least 90% length of the representative of the cluster",
                equate=False),
            _Option(
                ['-S', 'len_diff_cutoff_aa'],
                "length difference cutoff in amino acid, default 999999 "
                "if set to 60, the length difference between the shorter sequences "
                "and the representative of the cluster can not be bigger than 60",
                equate=False),
            _Option(
                ['-aL', 'cov_alignment_long'],
                "alignment coverage for the longer sequence, default 0.0 "
                "if set to 0.9, the alignment must covers 90% of the sequence",
                equate=False),
            _Option(
                ['-AL', 'cov_alignment_long_control'],
                "alignment coverage control for the longer sequence, default 99999999 "
                "if set to 60, and the length of the sequence is 400, "
                "then the alignment must be >= 340 (400-60) residues",
                equate=False),
            _Option(
                ['-aS', 'cov_alignment_short'],
                "alignment coverage for the shorter sequence, default 0.0 "
                "if set to 0.9, the alignment must covers 90% of the sequence",
                equate=False),
            _Option(
                ['-AS', 'cov_alignment_short_control'],
                "alignment coverage control for the shorter sequence, default 99999999 "
                "if set to 60, and the length of the sequence is 400, "
                "then the alignment must be >= 340 (400-60) residues",
                equate=False),
            _Option(
                ['-A', 'cov_alignment'],
                "minimal alignment coverage control for the both sequences, default 0 "
                "alignment must cover >= this value for both sequences",
                equate=False),
            _Option(
                ['-uL', 'max_unmatched_percentage_long'],
                "maximum unmatched percentage for the longer sequence, default 1.0 "
                "if set to 0.1, the unmatched region (excluding leading and tailing gaps) "
                "must not be more than 10% of the sequence",
                equate=False),
            _Option(
                ['-uS', 'max_unmatched_percentage_short'],
                "maximum unmatched percentage for the shorter sequence, default 1.0 "
                "if set to 0.1, the unmatched region (excluding leading and tailing gaps) "
                "must not be more than 10% of the sequence",
                equate=False),
            _Option(
                ['-U', 'len_max_unmatched'],
                "maximum unmatched length, default 99999999 "
                "if set to 10, the unmatched region (excluding leading and tailing gaps) "
                "must not be more than 10 bases",
                equate=False),
            _Option(
                ['-B', 'hdd_storage'],
                "1 or 0, default 0, by default, sequences are stored in RAM "
                "if set to 1, sequence are stored on hard drive "
                "it is recommended to use -B 1 for huge databases",
                equate=False),
            _Option(
                ['-p', 'aln_overlap_2_file'],
                "1 or 0, default 0 "
                "if set to 1, print alignment overlap in .clstr file",
                equate=False),
            _Option(
                ['-g', 'accurate_mode'],
                "1 or 0, default 0 "
                "by cd-hit's default algorithm, a sequence is clustered to the first "
                "cluster that meet the threshold (fast cluster). If set to 1, the program "
                "will cluster it into the most similar cluster that meet the threshold "
                "(accurate but slow mode) "
                "but either 1 or 0 won't change the representatives of final clusters",
                equate=False),
            _Option(['-bak', 'backup'], "write backup cluster file (1 or 0, default 0)", equate=False),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
