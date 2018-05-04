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
"""This script provides a command-line interface to ConKit's plotting functionality.

You are provided with a single access point to many different kinds of plots.

For more specific descriptions, call each subcommand's help menu directly.

"""

__author__ = "Felix Simkovic"
__date__ = "08 Feb 2017"
__version__ = "0.1"

import argparse
import inspect

import conkit.command_line
import conkit.io
import conkit.plot
import conkit.plot.tools

logger = None


def _add_default_args(parser):
    parser.add_argument('-dpi', dest='dpi', default=300, type=int, help="the resolution in DPI")
    parser.add_argument('-o', dest='output', default=None, type=str, help="output file")
    parser.add_argument(
        "--overwrite", dest="overwrite", default=False, action="store_true", help="overwrite output file if exists")


def _add_contact_default_args(parser):
    parser.add_argument('confile', help="Path to contact file")
    parser.add_argument('conformat', help="Format of contact file")


def _add_sequence_default_args(parser):
    parser.add_argument('seqfile', help="Path to sequence file")
    parser.add_argument('seqformat', help="Format of sequence file")


def _add_msa_default_args(parser):
    parser.add_argument("msafile", help="Path to MSA file")
    parser.add_argument("msaformat", help="Format of MSA file")


def _add_structure_default_args(parser):
    parser.add_argument('pdbfile', help="Path to structure file")
    parser.add_argument('pdbformat', help="Format of structure file")


def add_contact_map_args(subparsers):
    description = u"""
This command will plot a contact map using the provided contacts
alongside any additional information.

If you provide a reference contact map, it is assumed that the sequence
is identical for both contact predictions.

If you provide a reference structure, the true positive contacts
are identified by a distance of <8\u212B between C\u03B2-C\u03B2 atoms.

!!! IMPORTANT
=============
If contacts cannot be matched between your prediction and the
reference structure, they will not be plotted.

"""
    subparser = subparsers.add_parser(
        'cmap',
        help="Plot a contact map",
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    subparser.add_argument(
        '-c',
        dest='refid',
        default=None,
        help='Reference identifier to use [default: first in file]. '
        'Inter-molecular predictions use two letter convention, '
        'i.e AD for contacts between A and D.')
    subparser.add_argument('-d', dest='dtn', default=5, type=int, help='Minimum sequence separation')
    subparser.add_argument('-e', dest='otherfile', default=None, help='a second contact map to plot for comparison')
    subparser.add_argument('-ef', dest='otherformat', default=None, help='the format of the second contact map')
    subparser.add_argument(
        '-f', dest='dfactor', default=1.0, type=float, help='number of contacts to include relative to sequence length')
    subparser.add_argument('-p', dest='reffile', default=None, type=str, help="A reference file")
    subparser.add_argument('-pf', dest='refformat', default=None, type=str, help="A reference file")
    subparser.add_argument('--confidence', action="store_true", default=False, help='Plot the confidence scores')
    subparser.add_argument('--interchain', action="store_true", default=False, help='Plot inter-chain contacts')
    _add_default_args(subparser)
    _add_sequence_default_args(subparser)
    _add_contact_default_args(subparser)
    subparser.set_defaults(which='contact_map')


def add_contact_map_chord_args(subparsers):
    description = u"""
This command will plot a contact map using the provided contacts
in a Chord diagram. This will illustrate your sequence in circular
style with residues being connected by their contacts
"""
    subparser = subparsers.add_parser(
        'chord',
        help="Plot a contact map chord diagram",
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    subparser.add_argument('-d', dest='dtn', default=5, type=int, help='Minimum sequence separation')
    subparser.add_argument(
        '-f', dest='dfactor', default=1.0, type=float, help='number of contacts to include relative to sequence length')
    subparser.add_argument('--confidence', action="store_true", default=False, help='Plot the confidence scores')
    _add_default_args(subparser)
    _add_sequence_default_args(subparser)
    _add_contact_default_args(subparser)
    subparser.set_defaults(which='contact_map_chord')


def add_contact_map_matrix_args(subparsers):
    description = u"""
This command will plot a contact map matrix using the provided contacts
alongside any additional information.

"""
    subparser = subparsers.add_parser(
        'cmat',
        help="Plot a contact map matrix",
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    subparser.add_argument('-e', dest='otherfile', default=None, help='a second contact map to plot for comparison')
    subparser.add_argument('-ef', dest='otherformat', default=None, help='the format of the second contact map')
    _add_default_args(subparser)
    _add_sequence_default_args(subparser)
    _add_contact_default_args(subparser)
    subparser.set_defaults(which='contact_map_matrix')


def add_contact_density_args(subparsers):
    description = u"""
This command will plot a contact density plot using the provided
contacts. It will illustrate the density to help define domain
boundaries better.
"""
    subparser = subparsers.add_parser(
        'cdens',
        help="Plot a contact density plot",
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    subparser.add_argument('-b', dest='bw_method', default='amise', help='The bandwidth estimation method')
    subparser.add_argument('-d', dest='dtn', default=5, type=int, help='Minimum sequence separation')
    subparser.add_argument(
        '-f', dest='dfactor', default=5.0, type=float, help='number of contacts to include relative to sequence length')
    _add_default_args(subparser)
    _add_sequence_default_args(subparser)
    _add_contact_default_args(subparser)
    subparser.set_defaults(which='contact_density')


def add_precision_evaluation_args(subparsers):
    description = u"""
This command will plot an evaluation plot illustrating the precision score of
the provided contact prediction compared against a protein structure at different
cutoff thresholds.

"""
    subparser = subparsers.add_parser(
        'peval',
        help="Plot the precision evaluation plot",
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    subparser.add_argument(
        '-c',
        dest='pdbchain',
        default=None,
        help=
        'PDB chain to use [default: first in file]. Inter-molecular predictions use two letter convention, i.e AD for contacts between A and D.'
    )
    subparser.add_argument('-d', dest='dtn', default=5, type=int, help='Minimum sequence separation')
    subparser.add_argument(
        '-j', dest='cutoff_step', default=0.2, type=float, help='The cutoff step for contact selection')
    subparser.add_argument(
        '-min', dest='min_cutoff', default=0.0, type=float, help='The minimum factor for contact selection')
    subparser.add_argument(
        '-max', dest='max_cutoff', default=100.0, type=float, help='The maximum factor for contact selection')
    subparser.add_argument('--interchain', action="store_true", default=False, help='Plot inter-chain contacts')
    _add_default_args(subparser)
    _add_structure_default_args(subparser)
    _add_sequence_default_args(subparser)
    _add_contact_default_args(subparser)
    subparser.set_defaults(which='precision_evaluation')


def add_sequence_coverage_args(subparsers):
    description = u"""
This command will plot a coverage plot for every position in your alignment.

"""
    subparser = subparsers.add_parser(
        'scov',
        help="Plot the sequence coverage",
        description=description,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    subparser.add_argument('-id', dest='identity', default=0.7, type=float, help='sequence identity')
    _add_default_args(subparser)
    _add_msa_default_args(subparser)
    subparser.set_defaults(which='sequence_coverage')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    subparsers = parser.add_subparsers()

    # Add the subparsers
    #
    # Note:
    #     New subparsers can be added automatically as long as the convention is followed.
    #     I.e., To add a new function automatically, call it add_*_args and it will be
    #     picked up here.
    functions = [
        k for k in globals()
        if k.startswith("add") and k.endswith("args") and (
            inspect.ismethod(globals()[k]) or inspect.isfunction(globals()[k]))
    ]
    for f_name in functions:
        globals()[f_name](subparsers)

    args = parser.parse_args()

    global logger
    logger = conkit.command_line.setup_logging(level='info')

    if args.which == 'contact_map':
        if args.interchain:
            logger.info('This script is experimental for inter-chain contact plotting')

        logger.info('Min sequence separation for contacting residues: %d', args.dtn)
        logger.info('Contact list cutoff factor: %f * L', args.dfactor)

        seq = conkit.io.read(args.seqfile, args.seqformat)[0]
        con = conkit.io.read(args.confile, args.conformat)[0]

        con.sequence = seq
        con.set_sequence_register()
        con.remove_neighbors(min_distance=args.dtn, inplace=True)
        ncontacts = int(seq.seq_len * args.dfactor)
        con.sort('raw_score', reverse=True, inplace=True)
        con_sliced = con[:ncontacts]

        if args.otherfile:
            other = conkit.io.read(args.otherfile, args.otherformat)[0]
            other.sequence = seq
            other.set_sequence_register()
            other.remove_neighbors(min_distance=args.dtn, inplace=True)
            other.sort('raw_score', reverse=True, inplace=True)
            other_sliced = other[:ncontacts]
        else:
            other_sliced = None

        if args.reffile:
            if args.refid:
                reference = conkit.io.read(args.reffile, args.refformat)[args.refid]
            else:
                reference = conkit.io.read(args.reffile, args.refformat)[0]

            if args.refformat not in ['pdb', 'mmcif']:
                msg = "The provided format {0} is not yet implemented for the reference flag".format(args.refformat)
                raise RuntimeError(msg)

            con_matched = con_sliced.match(reference, match_other=True, renumber=True, remove_unmatched=True)
            if other_sliced:
                other_matched = other_sliced.match(reference, renumber=True, remove_unmatched=True)
            else:
                other_matched = other_sliced
        else:
            reference = None
            con_matched = con_sliced
            other_matched = other_sliced

        def altloc_remove(cmap):
            """Remove alternative locations"""
            altloc = False
            for contact in cmap.copy():
                # For now we need this args.interchain check to account for gapped residues
                # where res chain was not assigned
                if contact.res1_chain != contact.res2_chain and args.interchain:
                    altloc = True
                    break
                if contact.res1_chain == contact.res2_chain and args.interchain:
                    cmap.remove(contact.id)
                elif contact.res1_chain != contact.res2_chain and not args.interchain:
                    cmap.remove(contact.id)
            return altloc

        altloc = altloc_remove(con_matched)
        if other_matched:
            altloc = altloc_remove(other_matched)

        figure = conkit.plot.ContactMapFigure(
            con_matched,
            other=other_matched,
            reference=reference,
            altloc=altloc,
            use_conf=args.confidence,
            legend=True,
        )
        figure_aspect_ratio = 1.0

    elif args.which == 'contact_map_chord':
        logger.info('Min sequence separation for contacting residues: %d', args.dtn)
        logger.info('Contact list cutoff factor: %f * L', args.dfactor)

        seq = conkit.io.read(args.seqfile, args.seqformat)[0]
        con = conkit.io.read(args.confile, args.conformat)[0]

        con.sequence = seq
        con.set_sequence_register()
        con.remove_neighbors(min_distance=args.dtn, inplace=True)
        con.sort('raw_score', reverse=True, inplace=True)
        ncontacts = int(seq.seq_len * args.dfactor)
        con_sliced = con[:ncontacts]

        figure = conkit.plot.ContactMapChordFigure(con_sliced, use_conf=args.confidence, legend=True)
        figure_aspect_ratio = 1.0

    elif args.which == 'contact_map_matrix':

        seq = conkit.io.read(args.seqfile, args.seqformat)[0]
        con = conkit.io.read(args.confile, args.conformat)[0]

        con.sequence = seq
        con.set_sequence_register()

        if args.otherfile:
            other = conkit.io.read(args.otherfile, args.otherformat)[0]
            other.sequence = seq
            other.set_sequence_register()
        else:
            other = con

        figure = conkit.plot.ContactMapMatrixFigure(con, other=other, legend=True)
        figure_aspect_ratio = 1.0

    elif args.which == 'contact_density':
        logger.info('Min sequence separation for contacting residues: %d', args.dtn)
        logger.info('Contact list cutoff factor: %f * L', args.dfactor)
        logger.info('Bandwidth estimator: %s', args.bw_method)

        seq = conkit.io.read(args.seqfile, args.seqformat)[0]
        con = conkit.io.read(args.confile, args.conformat)[0]

        con.sequence = seq
        con.set_sequence_register()
        con.remove_neighbors(min_distance=args.dtn, inplace=True)
        con.sort('raw_score', reverse=True, inplace=True)
        ncontacts = int(seq.seq_len * args.dfactor)
        con_sliced = con[:ncontacts]

        figure = conkit.plot.ContactDensityFigure(con_sliced, bw_method=args.bw_method, legend=True)
        figure_aspect_ratio = 0.3

    elif args.which == 'precision_evaluation':
        if args.interchain:
            logger.info('This script is experimental for inter-chain contact plotting')

        logger.info('Min sequence separation for contacting residues: %d', args.dtn)
        logger.info('Min contact list cutoff factor: %f * L', args.min_cutoff)
        logger.info('Max contact list cutoff factor: %f * L', args.max_cutoff)
        logger.info('Contact list cutoff factor step: %f', args.cutoff_step)

        seq = conkit.io.read(args.seqfile, args.seqformat)[0]
        con = conkit.io.read(args.confile, args.conformat)[0]

        con.sequence = seq
        con.set_sequence_register()
        con.remove_neighbors(min_distance=args.dtn, inplace=True)
        con.sort('raw_score', reverse=True, inplace=True)

        if args.pdbchain:
            pdb = conkit.io.read(args.pdbfile, 'pdb')[args.pdbchain]
        else:
            pdb = conkit.io.read(args.pdbfile, 'pdb')[0]
        con_matched = con.match(pdb, renumber=True, remove_unmatched=True)

        figure = conkit.plot.PrecisionEvaluationFigure(
            con_matched,
            cutoff_step=args.cutoff_step,
            min_cutoff=args.min_cutoff,
            max_cutoff=args.max_cutoff,
            legend=True)
        figure_aspect_ratio = 0.3

    elif args.which == 'sequence_coverage':
        hierarchy = conkit.io.read(args.msafile, args.msaformat)
        figure = conkit.plot.SequenceCoverageFigure(hierarchy, legend=True)
        figure_aspect_ratio = 0.3

    if not args.output:
        args.output = "conkit.png"
    figure.ax.set_aspect(conkit.plot.tools.get_adjusted_aspect(figure.ax, figure_aspect_ratio))
    figure.savefig(args.output, dpi=args.dpi, overwrite=args.overwrite)
    logger.info('Final plot written to %s', args.output)


if __name__ == "__main__":
    main()
