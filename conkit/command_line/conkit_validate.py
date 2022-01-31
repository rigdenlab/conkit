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

It uses one external program to perform this task:

   map_align for contact map alignment

*** This program needs to be installed separately from https://github.com/sokrypton/map_align***

"""

import argparse
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import DSSP
import os
from prettytable import PrettyTable

import conkit.applications
import conkit.command_line
import conkit.io
import conkit.plot
from conkit.plot.tools import is_executable

logger = None


def create_argument_parser():
    """Create a parser for the command line arguments used in conkit-validate"""

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("seqfile", type=check_file_exists, help="Path to sequence file")
    parser.add_argument("seqformat", type=str, help="Sequence format")
    parser.add_argument("distfile", type=check_file_exists, help="Path to distance prediction file")
    parser.add_argument("distformat", type=str, help="Format of distance prediction file",
                        choices=list(conkit.io.DISTANCE_FILE_PARSERS.keys()))
    parser.add_argument("pdbfile", type=check_file_exists, help="Path to structure file")
    parser.add_argument("pdbformat", type=str, help="Format of structure file", choices=['pdb', 'mmcif'])
    parser.add_argument("-dssp_exe", dest="dssp", default='mkdssp', help="path to dssp executable", type=is_executable)
    parser.add_argument("-output", dest="output", default="conkit.png", help="path to output figure png file", type=str)
    parser.add_argument("--overwrite", dest="overwrite", default=False, action="store_true",
                        help="overwrite output figure png file if it already exists")
    parser.add_argument("--map_align_exe", dest="map_align_exe", default=None,
                        type=is_executable, help="Path to the map_align executable")
    parser.add_argument("--gap_opening_penalty", dest="gap_opening_penalty", default=-1, type=float,
                        help="Gap opening penalty")
    parser.add_argument("--gap_extension_penalty", dest="gap_extension_penalty", default=-0.01, type=float,
                        help="Gap extension penalty")
    parser.add_argument("--seq_separation_cutoff", dest="seq_separation_cutoff", default=3, type=int,
                        help="Sequence separation cutoff"),
    parser.add_argument("--n_iterations", dest="n_iterations", default=20, type=int,
                        help="Number of iterations")

    return parser


def check_file_exists(input_path):
    """Check if a given path exists

    Parameters
    ----------
    input_path : str, None
       Location of the file to be tested

    Returns
    -------
    abspath : str, None
       The absolute path of the file if it exists, None if the input is None

    Raises
    ------
    :exc:`FileNotFoundError`
        The file doesn't exist
    """

    if input_path is None:
        return None
    if os.path.isfile(os.path.abspath(input_path)):
        return os.path.abspath(input_path)
    else:
        raise FileNotFoundError("{} cannot be found".format(input_path))


def main():
    """The main routine for conkit-validate functionality"""
    parser = create_argument_parser()
    args = parser.parse_args()

    global logger
    logger = conkit.command_line.setup_logging(level="info")

    if os.path.isfile(args.output) and not args.overwrite:
        raise FileExistsError('The output file {} already exists!'.format(args.output))

    logger.info(os.linesep + "Working directory:                           %s", os.getcwd())
    logger.info("Reading input sequence:                      %s", args.seqfile)
    sequence = conkit.io.read(args.seqfile, args.seqformat).top

    if len(sequence) < 5:
        raise ValueError('Cannot validate model with less than 5 residues')

    logger.info("Length of the sequence:                      %d", len(sequence))
    logger.info("Reading input distance prediction:           %s", args.distfile)
    prediction = conkit.io.read(args.distfile, args.distformat).top
    logger.info("Reading input PDB model:                     %s", args.pdbfile)
    model = conkit.io.read(args.pdbfile, args.pdbformat).top
    p = PDBParser()
    structure = p.get_structure('structure', args.pdbfile)[0]
    dssp = DSSP(structure, args.pdbfile, dssp=args.dssp, acc_array='Wilke')

    logger.info(os.linesep + "Validating model.")

    if len(sequence) > 500:
        logger.info("Input model has more than 500 residues, this might take a while...")

    figure = conkit.plot.ModelValidationFigure(model, prediction, sequence, dssp, map_align_exe=args.map_align_exe)
    figure.savefig(args.output, overwrite=args.overwrite)
    logger.info(os.linesep + "Validation plot written to %s", args.output)

    residue_info = figure.data.loc[:, ['RESNUM', 'SCORE', 'MISALIGNED']]
    table = PrettyTable()
    table.field_names = ["Residue", "Predicted score", "Suggested register"]

    _resnum_template = '{} ({})'
    _error_score_template = '*** {0:.2f} ***'
    _correct_score_template = '    {0:.2f}    '
    _register_template = '*** {} ({}) ***'
    _empty_register = '               '

    for residue in residue_info.values:
        resnum, score, misalignment = residue
        current_residue = _resnum_template.format(sequence.seq[resnum - 1], resnum)
        score = _error_score_template.format(score) if score > 0.5 else _correct_score_template.format(score)

        if misalignment and resnum in figure.alignment.keys():
            register = _register_template.format(sequence.seq[figure.alignment[resnum] - 1], figure.alignment[resnum])
        else:
            register = _empty_register

        table.add_row([current_residue, score, register])

    logger.info(os.linesep)
    logger.info(table)


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
