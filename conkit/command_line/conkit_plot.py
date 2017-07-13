#!/usr/bin/env python
#
# BSD 3-Clause License
#
# Copyright (c) 2016-17, University of Liverpool
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
 
logger = None

# Note:
#     New subparsers can be added automatically as long as the convention is followed.
#     I.e., To add a new function automatically, call it add_*_args and it will be
#     picked up here.

def _add_default_args(parser):
    """Define default arguments"""
    parser.add_argument('-dpi', dest='dpi', default=300, type=int,
                        help="the resolution in DPI [default: 300]")
    parser.add_argument('-o', dest='output', default=None, type=str,
                        help="the figure file. Note, the format is determined by the suffix")


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
    contact_map_subparser = subparsers.add_parser('cmap', help="Plot a contact map", description=description,
                                                  formatter_class=argparse.RawDescriptionHelpFormatter)
    _add_default_args(contact_map_subparser)
    contact_map_subparser.add_argument('-c', dest='refid', default=None,
                                       help='Reference identifier to use [default: first in file]. '
                                            'Inter-molecular predictions use two letter convention, '
                                            'i.e AD for contacts between A and D.')
    contact_map_subparser.add_argument('-d', dest='dtn', default=5, type=int,
                                       help='Minimum sequence separation [default: 5]')
    contact_map_subparser.add_argument('-e', dest='otherfile', default=None,
                                       help='a second contact map to plot for comparison')
    contact_map_subparser.add_argument('-ef', dest='otherformat', default=None,
                                       help='the format of the second contact map')
    contact_map_subparser.add_argument('-f', dest='dfactor', default=1.0, type=float,
                                       help='number of contacts to include relative to sequence length [default: 1.0]')
    contact_map_subparser.add_argument('-p', dest='reffile', default=None, type=str,
                                       help="A reference file")
    contact_map_subparser.add_argument('-pf', dest='refformat', default=None, type=str,
                                       help="A reference file")
    contact_map_subparser.add_argument('--confidence', action="store_true", default=False,
                                       help='Plot the confidence scores')
    contact_map_subparser.add_argument('--interchain', action="store_true", default=False,
                                       help='Plot inter-chain contacts')
    contact_map_subparser.add_argument('seqfile',
                                       help="Path to the sequence file")
    contact_map_subparser.add_argument('seqformat',
                                       help="Format of the sequence file")
    contact_map_subparser.add_argument('confile',
                                       help="Path to the contact file")
    contact_map_subparser.add_argument('conformat',
                                       help="Format of the contact file")
    contact_map_subparser.set_defaults(which='contact_map')


def add_contact_map_chord_args(subparsers):
    description = u"""
This command will plot a contact map using the provided contacts
in a Chord diagram. This will illustrate your sequence in circular
style with residues being connected by their contacts
"""
    contact_map_chord_subparser = subparsers.add_parser('chord', help="Plot a contact map chord diagram",
                                                        description=description,
                                                        formatter_class=argparse.RawDescriptionHelpFormatter)
    _add_default_args(contact_map_chord_subparser)
    contact_map_chord_subparser.add_argument('-d', dest='dtn', default=5, type=int,
                                             help='Minimum sequence separation [default: 5]')
    contact_map_chord_subparser.add_argument('-f', dest='dfactor', default=1.0, type=float,
                                             help='number of contacts to include relative to '
                                                  'sequence length [default: 1.0]')
    contact_map_chord_subparser.add_argument('--confidence', action="store_true", default=False,
                                             help='Plot the confidence scores')
    contact_map_chord_subparser.add_argument('seqfile',
                                             help="Path to the sequence file")
    contact_map_chord_subparser.add_argument('seqformat',
                                             help="Format of the sequence file")
    contact_map_chord_subparser.add_argument('confile',
                                             help="Path to the contact file")
    contact_map_chord_subparser.add_argument('conformat',
                                             help="Format of the contact file")
    contact_map_chord_subparser.set_defaults(which='contact_map_chord')


def add_contact_density_args(subparsers):
    description = u"""
This command will plot a contact density plot using the provided
contacts. It will illustrate the density to help define domain
boundaries better.
"""
    contact_density_subparser = subparsers.add_parser('cdens', help="Plot a contact density plot",
                                                        description=description,
                                                        formatter_class=argparse.RawDescriptionHelpFormatter)
    _add_default_args(contact_density_subparser)
    contact_density_subparser.add_argument('-b', dest='bw_method', default='amise',
                                           help='The bandwidth estimation method [default: amise]')
    contact_density_subparser.add_argument('-d', dest='dtn', default=5, type=int,
                                           help='Minimum sequence separation [default: 5]')
    contact_density_subparser.add_argument('-f', dest='dfactor', default=5.0, type=float,
                                           help='number of contacts to include relative to '
                                                'sequence length [default: 10.0]')
    contact_density_subparser.add_argument('seqfile',
                                           help="Path to the sequence file")
    contact_density_subparser.add_argument('seqformat',
                                           help="Format of the sequence file")
    contact_density_subparser.add_argument('confile',
                                           help="Path to the contact file")
    contact_density_subparser.add_argument('conformat',
                                           help="Format of the contact file")
    contact_density_subparser.set_defaults(which='contact_density')


def add_precision_evaluation_args(subparsers):
    description = u"""
This command will plot an evaluation plot illustrating the precision score of
the provided contact prediction compared against a protein structure at different
cutoff thresholds.

"""
    precision_evaluation_subparser = subparsers.add_parser('peval', help="Plot the precision evaluation plot",
                                                           description=description,
                                                           formatter_class=argparse.RawDescriptionHelpFormatter)
    _add_default_args(precision_evaluation_subparser)
    precision_evaluation_subparser.add_argument('-c', dest='pdbchain', default=None,
                                                help='PDB chain to use [default: first in file]. Inter-molecular '
                                                     'predictions use two letter convention, i.e AD for contacts '
                                                     'between A and D.')
    precision_evaluation_subparser.add_argument('-d', dest='dtn', default=5, type=int,
                                                help='Minimum sequence separation [default: 5]')
    precision_evaluation_subparser.add_argument('-j', dest='cutoff_step', default=0.2, type=float,
                                                help='The cutoff step for contact selection [default: 0.2]')
    precision_evaluation_subparser.add_argument('-min', dest='min_cutoff', default=0.0, type=float,
                                                help='The minimum factor for contact selection [default: 0.0]')
    precision_evaluation_subparser.add_argument('-max', dest='max_cutoff', default=100.0, type=float,
                                                help='The maximum factor for contact selection [default: 100.0]')
    precision_evaluation_subparser.add_argument('--interchain', action="store_true", default=False,
                                                help='Plot inter-chain contacts')
    precision_evaluation_subparser.add_argument('pdbfile',
                                                help="A reference PDB file")
    precision_evaluation_subparser.add_argument('pdbformat',
                                                help="A reference PDB file")
    precision_evaluation_subparser.add_argument('seqfile',
                                                help="Path to the sequence file")
    precision_evaluation_subparser.add_argument('seqformat',
                                                help="Format of the sequence file")
    precision_evaluation_subparser.add_argument('confile',
                                                help="Path to the contact file")
    precision_evaluation_subparser.add_argument('conformat',
                                                help="Format of the contact file")
    precision_evaluation_subparser.set_defaults(which='precision_evaluation')


def add_sequence_coverage_args(subparsers):
    description = u"""
This command will plot a coverage plot for every position in your alignment.

"""
    sequence_coverage_subparser = subparsers.add_parser('scov', help="Plot the sequence coverage",
                                                        description=description,
                                                        formatter_class=argparse.RawDescriptionHelpFormatter)
    _add_default_args(sequence_coverage_subparser)
    sequence_coverage_subparser.add_argument('-id', dest='identity', default=0.7, type=float,
                                             help='sequence identity [default: 0.7]')
    sequence_coverage_subparser.add_argument('msafile',
                                             help='Multiple Sequence Alignment file')
    sequence_coverage_subparser.add_argument('msaformat',
                                             help='Multiple Sequence Alignment format')
    sequence_coverage_subparser.set_defaults(which='sequence_coverage')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    subparsers = parser.add_subparsers()

    # Add the subparsers
    #
    # Note:
    #     New subparsers can be added automatically as long as the convention is followed.
    #     I.e., To add a new function automatically, call it add_*_args and it will be
    #     picked up here.
    functions = [k for k in globals() if k.startswith("add") and k.endswith("args")
                 and (inspect.ismethod(globals()[k]) or inspect.isfunction(globals()[k]))]
    for f_name in functions:
        globals()[f_name](subparsers)

    # Parse all arguments
    args = parser.parse_args()

    # Setup the logger
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
        con.assign_sequence_register()
        con.remove_neighbors(min_distance=args.dtn, inplace=True)
        ncontacts = int(seq.seq_len * args.dfactor)
        con.sort('raw_score', reverse=True, inplace=True)
        con_sliced = con[:ncontacts]

        if args.otherfile:
            other = conkit.io.read(args.otherfile, args.otherformat)[0]
            other.sequence = seq
            other.assign_sequence_register()
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

            con_matched = con_sliced.match(reference, renumber=True, remove_unmatched=True)
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

        outformat = 'png'
        outfile = args.output if args.output else args.confile.rsplit('.', 1)[0] + '.' + outformat
        plot = conkit.plot.ContactMapFigure(con_matched, other=other_matched, reference=reference,
                                            file_name=outfile, altloc=altloc, use_conf=args.confidence,
                                            dpi=args.dpi)

    elif args.which == 'contact_map_chord':
        logger.info('Min sequence separation for contacting residues: %d', args.dtn)
        logger.info('Contact list cutoff factor: %f * L', args.dfactor)

        seq = conkit.io.read(args.seqfile, args.seqformat)[0]
        con = conkit.io.read(args.confile, args.conformat)[0]

        con.sequence = seq
        con.assign_sequence_register()
        con.remove_neighbors(min_distance=args.dtn, inplace=True)
        con.sort('raw_score', reverse=True, inplace=True)
        ncontacts = int(seq.seq_len * args.dfactor)
        con_sliced = con[:ncontacts]
        outformat = 'png'
        outfile = args.output if args.output else args.confile.rsplit('.', 1)[0] + '.' + outformat
        plot = conkit.plot.ContactMapChordFigure(con_sliced, use_conf=args.confidence, file_name=outfile, 
                                                 dpi=args.dpi)

    elif args.which == 'contact_density':
        logger.info('Min sequence separation for contacting residues: %d', args.dtn)
        logger.info('Contact list cutoff factor: %f * L', args.dfactor)
        logger.info('Bandwidth estimator: %s', args.bw_method)

        seq = conkit.io.read(args.seqfile, args.seqformat)[0]
        con = conkit.io.read(args.confile, args.conformat)[0]

        con.sequence = seq
        con.assign_sequence_register()
        con.remove_neighbors(min_distance=args.dtn, inplace=True)
        con.sort('raw_score', reverse=True, inplace=True)
        ncontacts = int(seq.seq_len * args.dfactor)
        con_sliced = con[:ncontacts]
        outformat = 'png'
        outfile = args.output if args.output else args.confile.rsplit('.', 1)[0] + '.' + outformat
        plot = conkit.plot.ContactDensityFigure(con_sliced, bw_method=args.bw_method, file_name=outfile, dpi=args.dpi)

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
        con.assign_sequence_register()
        con.remove_neighbors(min_distance=args.dtn, inplace=True)
        con.sort('raw_score', reverse=True, inplace=True)

        if args.pdbchain:
            pdb = conkit.io.read(args.pdbfile, 'pdb')[args.pdbchain]
        else:
            pdb = conkit.io.read(args.pdbfile, 'pdb')[0]
        con_matched = con.match(pdb, renumber=True, remove_unmatched=True)

        outformat = 'png'
        outfile = args.output if args.output else args.confile.rsplit('.', 1)[0] + '.' + outformat
        plot = conkit.plot.PrecisionEvaluationFigure(con_matched, cutoff_step=args.cutoff_step, file_name=outfile,
                                                     min_cutoff=args.min_cutoff, max_cutoff=args.max_cutoff,
                                                     dpi=args.dpi)

    elif args.which == 'sequence_coverage':
        hierarchy = conkit.io.read(args.msafile, args.msaformat)
        outformat = 'png'
        outfile = args.output if args.output else args.msafile.rsplit('.', 1)[0] + '.' + outformat
        plot = conkit.plot.SequenceCoverageFigure(hierarchy, file_name=outfile, dpi=args.dpi)

    logger.info('Final plot written in %s format to %s', plot.format.upper(), plot.file_name)

    return


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

