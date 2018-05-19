#!/usr/bin/env python
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
"""This script provides a simplified contact prediction pipeline. It will
take either a sequence or alignment as starting point and predict, analyse
and evaluate the final prediction and any intermediate results.

It uses two external programs to perform this task.

   - HHblits for Sequence Alignment generation, and
   - CCMpred for Direct Coupling Analysis.

*** The two programs above need to be installed separately ***

"""

__author__ = "Felix Simkovic"
__date__ = "01 June 2016"
__version__ = "0.1"

import argparse
import os
import time

import conkit.applications
import conkit.command_line
import conkit.io
import conkit.plot
import conkit.plot.tools

logger = None


def add_default_args(parser):
    """Define default arguments"""
    parser.add_argument('-prefix', default='conkit', help='Job ID')
    parser.add_argument('-wdir', default=os.getcwd(), help='Working directory')
    parser.add_argument('--demo', default=False, action="store_true", help=argparse.SUPPRESS)


def add_alignment_args(subparsers):
    """Parser with alignment as starting point"""
    from_alignment_subparser = subparsers.add_parser('aln')
    add_default_args(from_alignment_subparser)
    from_alignment_subparser.add_argument('ccmpred', help='Path to the CCMpred executable')
    from_alignment_subparser.add_argument('alignment_file', help='Path to alignment file')
    from_alignment_subparser.add_argument('alignment_format', help='Alignment format')
    from_alignment_subparser.set_defaults(which='alignment')


def add_sequence_args(subparsers):
    """Parser with sequence as starting point"""
    from_sequence_subparser = subparsers.add_parser('seq')
    add_default_args(from_sequence_subparser)
    from_sequence_subparser.add_argument('--nodca', default=False, action='store_true', help=argparse.SUPPRESS)
    from_sequence_subparser.add_argument('ccmpred', help='Path to the CCMpred executable')
    from_sequence_subparser.add_argument('hhblits', help='Path to the HHblits executable')
    from_sequence_subparser.add_argument('hhblitsdb', help='Path to HHblits database')
    from_sequence_subparser.add_argument('sequence_file', help='Path to sequence file')
    from_sequence_subparser.add_argument('sequence_format', help='Sequence format')
    from_sequence_subparser.set_defaults(which='sequence')


def main(argl=None):
    """The main routine for conkit-predict functionality

    Parameters
    ----------
    argl : list, tuple, optional
       A list containing the command line flags

    """
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers()
    add_alignment_args(subparsers)
    add_sequence_args(subparsers)
    args = parser.parse_args(argl)

    global logger
    logger = conkit.command_line.setup_logging(level='info')

    logger.info('Prefix: %s', args.prefix)
    logger.info('Working dir: %s', args.wdir)

    if args.which == 'alignment':
        aln_fname = os.path.abspath(args.alignment_file)
        aln_format = args.alignment_format.lower()

        logger.info('Input alignment file: %s', aln_fname)
        logger.info('Input alignment format: %s', aln_format)

        if aln_format not in conkit.io.SEQUENCE_FILE_PARSERS:
            msg = 'Sequence format not yet implemented: %s' % aln_fname
            logger.critical(msg)
            raise ValueError(msg)

        if aln_format != 'jones':
            jon_fname = os.path.join(args.wdir, args.prefix + '.jones')
            conkit.io.convert(aln_fname, aln_format, jon_fname, 'jones')
        else:
            jon_fname = aln_fname

    elif args.which == 'sequence':
        hhblits = os.path.abspath(args.hhblits)
        hhblitsdb = os.path.abspath(args.hhblitsdb)
        seq_fname = os.path.abspath(args.sequence_file)
        seq_format = args.sequence_format.lower()
        a3m_fname = os.path.join(args.wdir, args.prefix + '.a3m')
        hhr_fname = os.path.join(args.wdir, args.prefix + '.hhr')

        logger.info('HHblits DB: %s', hhblitsdb)
        logger.info('Input sequence file: %s', seq_fname)
        logger.info('Input sequence format: %s', seq_format)

        # Check that we can handle the sequence file
        if seq_format not in conkit.io.SEQUENCE_FILE_PARSERS:
            msg = 'Sequence format not yet implemented: %s' % seq_format
            logger.critical(msg)
            raise ValueError(msg)

        # Convert our sequence file to FASTA format
        if seq_format != 'fasta':
            seq_fname_tmp = seq_fname.rsplit('.', 1)[0] + '.fasta'
            conkit.io.convert(seq_fname, seq_format, seq_fname_tmp, 'fasta')
            seq_fname = seq_fname_tmp

        # Generate a multiple sequence alignment
        hhblits_cline = conkit.applications.HHblitsCommandline(
            cmd=hhblits,
            input=seq_fname,
            output=hhr_fname,
            database=hhblitsdb,
            oa3m=a3m_fname,
            niterations=3,
            id=99,
            show_all=True,
            cov=60,
            diff='inf',
            maxfilt=500000)
        logger.info('Executing: %s', hhblits_cline)
        if args.demo:
            assert os.path.isfile(a3m_fname)
            time.sleep(5)
        else:
            hhblits_cline()

        jon_fname = os.path.join(args.wdir, args.prefix + '.jones')
        conkit.io.convert(a3m_fname, 'a3m', jon_fname, 'jones')

    else:
        raise RuntimeError('Should never get to here')

    # CCMpred requires alignments to be in the *jones* format - i.e. the format created
    # and used by David Jones in PSICOV
    msa_h = conkit.io.read(jon_fname, 'jones')
    freq_plot_fname = os.path.join(args.wdir, args.prefix + 'freq.png')
    figure = conkit.plot.SequenceCoverageFigure(msa_h, legend=True)
    figure.ax.set_aspect(conkit.plot.tools.get_adjusted_aspect(figure.ax, 0.3))
    figure.savefig(freq_plot_fname)

    logger.info('Final alignment file: %s', jon_fname)
    logger.info('|- Total Number of sequences: %d', msa_h.nseq)
    logger.info('|- Number of effective sequences: %d', msa_h.meff)
    logger.info('|- Plotted sequence coverage: %s', freq_plot_fname)

    if args.which == 'sequence' and args.nodca:
        return

    ccmpred = args.ccmpred
    matrix_fname = os.path.join(args.wdir, args.prefix + '.mat')
    ccmpred_cline = conkit.applications.CCMpredCommandline(
        cmd=ccmpred, alnfile=jon_fname, matfile=matrix_fname, threads=2, renormalize=True)
    logger.info('Executing: %s', ccmpred_cline)
    if args.demo:
        assert os.path.isfile(matrix_fname)
        time.sleep(5)
    else:
        ccmpred_cline()

    dtn = 5
    dfactor = 1.
    cmap = conkit.io.read(matrix_fname, 'ccmpred').top_map
    cmap.sequence = conkit.io.read(jon_fname, 'jones').top_sequence
    cmap.remove_neighbors(min_distance=dtn, inplace=True)
    cmap.sort('raw_score', reverse=True, inplace=True)
    cmap = cmap[:cmap.sequence.seq_len]

    contact_map_fname = os.path.join(args.wdir, args.prefix + 'cmap.png')
    figure = conkit.plot.ContactMapFigure(cmap, legend=True)
    figure.ax.set_aspect(conkit.plot.tools.get_adjusted_aspect(figure.ax, 1.0))
    figure.savefig(contact_map_fname)

    logger.info('Plotted contact map: %s', contact_map_fname)
    logger.info('|- Min sequence separation for contacting residues: %d', dtn)
    logger.info('|- Contact list cutoff factor: %f * L', dfactor)

    casprr_fname = os.path.join(args.wdir, args.prefix + '.rr')
    conkit.io.convert(matrix_fname, 'ccmpred', casprr_fname, 'casprr')
    logger.info('Final prediction file: %s', casprr_fname)


if __name__ == "__main__":
    import sys
    import traceback
    try:
        main()
        sys.exit(0)
    except Exception as e:
        if not isinstance(e, SystemExit):
            msg = "".join(traceback.format_exception(*sys.exc_info()))
            logger.critical(msg)
        sys.exit(1)
