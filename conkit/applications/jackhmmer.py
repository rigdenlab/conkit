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
Command line object for Jackhmmer Multiple Sequence Alignment generation
"""

__author__ = "Felix Simkovic"
__date__ = "01 June 2016"
__version__ = "0.1"

from Bio.Application import _Argument
from Bio.Application import _Option
from Bio.Application import _Switch
from Bio.Application import AbstractCommandline


class JackhmmerCommandline(AbstractCommandline):
    """
    Command line object for Jackhmmer [#]_ alignment generation

    http://hmmer.org/

    Jackhmmer is an algorithm that uses iterative searches a protein sequence
    against a protein sequence database to find sequence homologs.

    .. [#] Johnson L. S., Eddy S. R., Portugaly E. (2010). Hidden Markov
       Model Speed Heuristic and Iterative HMM Search Procedure. BMC Bioinformatics 11, 431.

    Examples
    --------
    To generate a Multiple Sequence Alignment, use:

    >>> from conkit.applications import JackhmmerCommandline
    >>> jackhmmer_cline = JackhmmerCommandline(
    ...     input="test.fasta", database="uniref100.fasta"
    ... )
    >>> print(jackhmmer_cline)
    jackhmmer test.fasta uniref100.fasta

    You would typically run the command line with :func:`jackhmmer_cline` or via
    the :mod:`~subprocess` module.

    """

    def __init__(self, cmd='jackhmmer', **kwargs):
        self.parameters = [
            _Option(['-N', 'niterations'], 'set maximum number of iterations [default: 5]', equate=False),

            # Options directing output
            _Option(['-o', 'output'], 'direct output to file <f>, not stdout', filename=True, equate=False),
            _Option(['-A', 'msa_hits'], 'save multiple alignment of hits to file <f>', filename=True, equate=False),
            _Option(
                ['--tblout', 'per_sequence_hits'],
                'save parseable table of per-sequence hits to file <f>',
                filename=True,
                equate=False),
            _Option(
                ['--domtblout', 'per_domain_hits'],
                'save parseable table of per-domain hits to file <f>',
                filename=True,
                equate=False),
            _Option(
                ['--chkhmm', 'hmm_checkpoints'],
                'save HMM checkpoints to files <f>-<iteration>.hmm',
                filename=True,
                equate=False),
            _Option(
                ['--chkali', 'alignment_checkpoints'],
                'save alignment checkpoints to files <f>-<iteration>.sto',
                filename=True,
                equate=False),
            _Switch(['--acc', 'accession'], 'prefer accessions over names in output'),
            _Switch(['--noali', 'no_alignment'], 'don\'t output alignments, so output is smaller'),
            _Switch(['--notextw', 'notextw'], 'unlimit ASCII text output line width'),
            _Switch(['--textw', 'textw'], 'set max width of ASCII text output lines [default: 120] (n>=120)'),

            # Options controlling scoring system in first iteration
            _Option(['--popen', 'gap_open_probability'], 'gap open probability', equate=False),
            _Option(['--pextend', 'gap_extend_probability'], 'gap extend probability', equate=False),
            _Option(
                ['--mx', 'matrix_choice'], 'substitution score matrix choice (of some built-in matrices)',
                equate=False),
            _Option(
                ['--mxfile', 'matrix_option'],
                'read substitution score matrix from file <f>',
                filename=True,
                equate=False),

            # Options controlling reporting thresholds
            _Option(
                ['-E', 'evalue'],
                'report sequences <= this E-value threshold in output [default: 10.0] (x>0)',
                equate=False),
            _Option(['-T', 'score_threshold'], 'report sequences >= this score threshold in output', equate=False),
            _Option(
                ['--domE', 'domain_evalue'],
                'report domains <= this E-value threshold in output [default: 10.0] (x>0)',
                equate=False),
            _Option(
                ['--domT', 'domain_score_threshold'], 'report domains >= this score cutoff in output', equate=False),

            # Options controlling significance thresholds for inclusion in next round
            _Option(
                ['--incE', 'inclusion_evalue'],
                'consider sequences <= this E-value threshold as significant',
                equate=False),
            _Option(
                ['--incT', 'inclusion_score_threshold'],
                'consider sequences >= this score threshold as significant',
                equate=False),
            _Option(
                ['--incdomE', 'inclusion_domain_evalue'],
                'consider domains <= this E-value threshold as significant',
                equate=False),
            _Option(
                ['--incdomT', 'inclusion_domain_score_threshold'],
                'consider domains >= this score threshold as significant',
                equate=False),

            # Options controlling acceleration heuristics
            _Switch(['--max', 'no_heuristics'], 'Turn all heuristic filters off (less speed, more power)'),
            _Option(
                ['--F1', 'stage1_threshold'],
                'Stage 1 (MSV) threshold: promote hits w/ P <= F1 [default: 0.02]',
                equate=False),
            _Option(
                ['--F2', 'stage2_threshold'],
                'Stage 2 (Vit) threshold: promote hits w/ P <= F2 [default: 1e-3]',
                equate=False),
            _Option(
                ['--F3', 'stage3_threshold'],
                'Stage 3 (Fwd) threshold: promote hits w/ P <= F3 [default: 1e-5]',
                equate=False),
            _Switch(['--nobias', 'nobias'], 'turn off composition bias filter'),

            # Options controlling model construction after first iteration
            _Switch(['--fast', 'fast'], 'assign cols w/ >= symfrac residues as consensus'),
            _Switch(['--hand', 'hand'], 'manual construction (requires reference annotation)'),
            _Option(['--symfrac', 'symfrac'], 'sets sym fraction controlling --fast construction', equate=False),
            _Option(['--fragthres', 'fragthres'], 'if L <= x*alen, tag sequence as a fragment', equate=False),

            # Options controlling relative weights in models after first iteration
            _Switch(['--wpb', 'henikoff_pb_weights'], 'Henikoff position-based weights  [default]'),
            _Switch(['--wgsc', 'GSC_weights'], 'Gerstein/Sonnhammer/Chothia tree weights'),
            _Switch(['--wblosum', 'henikoff_sf_weights'], 'Henikoff simple filter weights'),
            _Switch(['--wnone', 'no_weight'], 'don\'t do any relative weighting; set all to 1'),
            _Option(
                ['--wid', 'wblosum_cutoff'],
                'for --wblosum: set identity cutoff [default: 0.62] (0<=x<=1)',
                equate=False),

            # Options controlling effective seq number in models after first iteration
            _Switch(['--eent', 'eent'], 'adjust eff seq # to achieve relative entropy target [default]'),
            _Switch(['--eclust', 'ecluse'], 'eff seq # is # of single linkage clusters'),
            _Switch(['--enone', 'enone'], 'no effective seq # weighting: just use nseq'),
            _Option(['--eset', 'eset'], 'set eff seq # for all models to <x>', equate=False),
            _Option(['--ere', 'ere'], 'for --eent: set minimum rel entropy/position to <x>', equate=False),
            _Option(['--esigma', 'esigma'], 'for --eent: set sigma param to <x> [default: 45.0]', equate=False),
            _Option(
                ['--eid', 'eid'], 'for --eclust: set fractional identity cutoff to <x> [default: 0.62]', equate=False),

            # Options controlling prior strategy in models after first iteration
            _Switch(['--pnone', 'pnone'], 'don\'t use any prior; parameters are frequencies'),
            _Switch(['--plaplace', 'plaplace'], 'use a Laplace +1 prior'),

            # Options controlling E value calibration
            _Option(['--EmL', 'eml'], 'length of sequences for MSV Gumbel mu fit [default: 200] (n>0)', equate=False),
            _Option(['--EmN', 'emn'], 'number of sequences for MSV Gumbel mu fit [default: 200] (n>0)', equate=False),
            _Option(
                ['--EvL', 'evl'], 'length of sequences for Viterbi Gumbel mu fit [default: 200] (n>0)', equate=False),
            _Option(
                ['--EvN', 'evn'], 'number of sequences for Viterbi Gumbel mu fit [default: 200] (n>0)', equate=False),
            _Option(
                ['--EfL', 'efl'], 'length of sequences for Forward exp tail tau fit [default: 100] (n>0)',
                equate=False),
            _Option(
                ['--EfN', 'efn'], 'number of sequences for Forward exp tail tau fit [default: 200] (n>0)',
                equate=False),
            _Option(
                ['--Eft', 'eft'],
                'tail mass for Forward exponential tail tau fit [default: 0.04] (0<x<1)',
                equate=False),

            # Other expert options
            _Switch(['--nonull2', 'nonull2'], 'turn off biased composition score corrections'),
            _Option(['-Z', 'ncomparison'], 'set # of comparisons done, for E-value calculation', equate=False),
            _Option(['--domZ', 'domz'], 'set # of significant seqs, for domain E-value calculation', equate=False),
            _Option(
                ['--seed', 'seed'], 'set RNG seed to <n> (if 0: one-time arbitrary seed) [default: 42]', equate=False),
            _Option(
                ['--qformat', 'qformat'], 'assert query <seqfile> is in format <s>: no autodetection', equate=False),
            _Option(
                ['--tformat', 'tformat'], 'assert target < seqdb > is in format < s >>: no autodetection',
                equate=False),
            _Option(['--cpu', 'cpu'], 'number of parallel CPU workers to use for multithreads', equate=False),

            # Required arguments
            _Argument(['input'], 'sequence containing file', filename=True, is_required=True),
            _Argument(['database'], 'sequence database', filename=True, is_required=True),
        ]

        AbstractCommandline.__init__(self, cmd, **kwargs)
