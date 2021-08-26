#!/usr/bin/env python
#
# BSD 3-Clause License
#
# Copyright (c) 2016-21, University of Liverpool
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
"""This script provides a model validation pipeline. It will
take a sequence, a predicted distogram and a PDB model, then
compare the distances observed in the model and those predicted.
It will then report regions in the model where an outlier was
detected, and use map_align to provide a solution for potential
register errors.

It uses one external programs to perform this task:

   - map_align for contact map alignment

*** This program needs to be installed separately ***

"""

import argparse
import os
from prettytable import PrettyTable

import conkit.applications
import conkit.command_line
import conkit.command_line.cli_tools
import conkit.io
import conkit.plot

logger = None


def create_argument_parser():
    """Create a parser for the command line arguments used in conkit-validate"""

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("seqfile", type=conkit.command_line.cli_tools.check_file_exists, help="Path to sequence file")
    parser.add_argument("seqformat", type=str, help="Sequence format")
    parser.add_argument("distfile", type=conkit.command_line.cli_tools.check_file_exists,
                        help="Path to distance prediction file")
    parser.add_argument("distformat", type=str, help="Format of distance prediction file")
    parser.add_argument("pdbfile", type=conkit.command_line.cli_tools.check_file_exists,
                        help="Path to structure file")
    parser.add_argument("pdbformat", type=str, help="Format of structure file")
    parser.add_argument("--outdir", dest="outdir", type=str, help="Output directory [default: ./conkit-validate]",
                        default=os.path.join(os.getcwd(), 'conkit-validate'))
    parser.add_argument("--overwrite", dest="overwrite", default=False, action="store_true",
                        help="overwrite output directory if exists [default: False]")
    parser.add_argument("--map_align_exe", dest="map_align_exe", default="map_align",
                        type=conkit.command_line.cli_tools.check_file_exists,
                        help="Path to the map_align executable [default: map_align]")
    parser.add_argument("--gap_opening_penalty", dest="gap_opening_penalty", default=-1, type=float,
                        help="Gap opening penalty [default: -1]")
    parser.add_argument("--gap_extension_penalty", dest="gap_extension_penalty", default=-0.01, type=float,
                        help="Gap extension penalty [default: -0.01]")
    parser.add_argument("--seq_separation_cutoff", dest="seq_separation_cutoff", default=3, type=int,
                        help="Sequence separation cutoff [default: 3]"),
    parser.add_argument("--n_iterations", dest="n_iterations", default=20, type=int,
                        help="Number of iterations [default=20]")

    return parser


def main():
    """The main routine for conkit-validate functionality

    Parameters
    ----------
    argv : dict, optional
       A list containing the command line flags

    """
    parser = create_argument_parser()
    args = parser.parse_args()

    global logger
    logger = conkit.command_line.setup_logging(level="info")

    logger.info("Output directory:                           %s", args.outdir)

    conkit.command_line.cli_tools.prepare_output_directory(args.outdir, args.overwrite)

    sequence = conkit.io.read(args.seqfile, args.seqformat).top
    prediction = conkit.io.read(args.distfile, args.distformat).top
    model = conkit.io.read(args.pdbfile, args.pdbformat).top
    logger.info("Input PDB model:                            %s", args.pdbfile)
    logger.info("Input distance prediction:                  %s", args.distfile)
    logger.info("Length of the sequence:                     %d", len(sequence))

    fig_fname = os.path.join(args.outdir, 'conkit.png')
    figure = conkit.plot.ModelValidationFigure(model, prediction, sequence, use_weights=True)
    figure.savefig(fig_fname, overwrite=args.overwrite)
    logger.info("Validation plot written to %s", fig_fname)

    if any(figure.outliers):
        table = PrettyTable()
        table.field_names = ["Outlier no.", "Residue Range", "wRMSD", "FN Count"]
        for idx, outlier in enumerate(figure.outliers, 1):
            start_outlier = outlier - 10 if outlier > 10 else 0
            stop_outlier = outlier + 10 if outlier + 10 < len(sequence) else len(sequence)
            table.add_row([idx, '{} - {}'.format(start_outlier, stop_outlier),
                           '{0:.2f}'.format(figure.rmsd_profile[outlier]),
                           '{0:.2f}'.format(figure.fn_profile[outlier])])
        logger.info("\nList of detected outliers:")
        logger.info(table)

        contact_map_a = os.path.join(args.outdir, 'contact_map_a.mapalign')
        contact_map_b = os.path.join(args.outdir, 'contact_map_b.mapalign')
        conkit.io.write(contact_map_a, 'mapalign', prediction)
        conkit.io.write(contact_map_b, 'mapalign', model)

        map_align_cline = conkit.applications.MapAlignCommandline(
            cmd=args.map_align_exe,
            contact_map_a=contact_map_a,
            contact_map_b=contact_map_b,
            gap_opening_penalty=args.gap_opening_penalty,
            gap_extension_penalty=args.gap_extension_penalty,
            seq_separation_cutoff=args.seq_separation_cutoff,
            n_iterations=args.n_iterations

        )

        logger.info("\nExecuting: %s", map_align_cline)
        stdout, stderr = map_align_cline()
        map_align_log = os.path.join(args.outdir, 'map_align.log')
        with open(map_align_log, 'w') as fhandle:
            fhandle.write(stdout)
        alignment = conkit.command_line.cli_tools.parse_map_align_stdout(stdout)
        for idx, outlier in enumerate(figure.outliers, 1):
            table = PrettyTable()
            table.field_names = ["Current Residue", "New Residue"]
            start_outlier = outlier - 10 if outlier > 10 else 0
            stop_outlier = outlier + 10 if outlier + 10 < len(sequence) else len(sequence)
            found_fix = False
            for resnum in range(start_outlier, stop_outlier + 1):
                if resnum in alignment.keys() and alignment[resnum] != resnum:
                    found_fix = True
                    table.add_row(['{} ({})'.format(sequence.seq[resnum-1], resnum),
                                   '{} ({})'.format(sequence.seq[alignment[resnum]-1], alignment[resnum])])
            if found_fix:
                logger.info("\nList of proposed changes to fix outlier no. {}:".format(idx))
                logger.info(table)


    else:
        logger.info("No outliers were detected, finishing now.")


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