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
import subprocess
import json
import numpy as np
from prettytable import PrettyTable

import conkit.applications
import conkit.command_line
import conkit.io
import conkit.io.tools import set_contact_definition
import conkit.plot
from conkit.plot.tools import is_executable, areaimol_ACC
from conkit.misc.renumbering_tools import write_renumbered_version_of_chain_in_struct

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
    parser.add_argument("-output_json", dest="output_json", default=None, help="path to output json file", type=str)
    parser.add_argument("--overwrite", dest="overwrite", default=False, action="store_true",
                        help="overwrite output figure png file if it already exists")
    parser.add_argument("--map_align_exe", dest="map_align_exe", default=None,
                        type=is_executable, help="Path to the map_align executable")
    parser.add_argument("--gesamt_exe", dest="gesamt_exe", default=None,
                        type=is_executable, help="Path to the gesamt executable to check structural alignment")
    parser.add_argument("--areaimol_exe", dest="areaimol_exe", default="areaimol",
                        type=is_executable, help="Path to areaimol executable to calculate solvent accesibility for RNA")
    parser.add_argument("--gemmi_exe", dest="gemmi_exe", default="gemmi",
                        type=is_executable, help="Path to the gemmi executable for converting mmcif to legacy pdb required for areaimol")
    parser.add_argument("--gap_opening_penalty", dest="gap_opening_penalty", default=-1, type=float,
                        help="Gap opening penalty")
    parser.add_argument("--gap_extension_penalty", dest="gap_extension_penalty", default=-0.01, type=float,
                        help="Gap extension penalty")
    parser.add_argument("--seq_separation_cutoff", dest="seq_separation_cutoff", default=3, type=int,
                        help="Sequence separation cutoff"),
    parser.add_argument("--n_iterations", dest="n_iterations", default=20, type=int,
                        help="Number of iterations")
    parser.add_argument("--moltype", dest="moltype", default="Protein", type=str,
                        help="Type of molecule")
    parser.add_argument("--chain", dest="selected_chain", default="", type=str,
                        help="Type of molecule")
    parser.add_argument("--renumber_model", dest="RENUMBER", default='yes', type=str,
                        help="Whether to crops-style create a version of the input model renumbered to match the sequence")
    parser.add_argument("--run_svm", dest="RUN_SVM", default='yes if prediction not pdb or mmcif', type=str,
                        help="Whether to run the support vector machine validation")
    parser.add_argument("--run_map_align", dest="RUN_MAP_ALIGN", default='yes', type=str,
                        help="Whether to run the contactmap alignment validation")
    parser.add_argument("--run_filters", dest="RUN_FILTERS", default='yes', type=str,
                        help="Whether to run the filters against false positives(if possible given the provided info)")
    parser.add_argument("--contact_dist", dest="contact_distance_cutoff", default=None, type=float,
                        help="distance cutoff for contacts when using Custom moltype")
    parser.add_argument("--rep_atom", dest="rep_atom", default=None, type=str,
                        help="representative atom for contacts when using Custom moltype")
    parser.add_argument("--min_error_length", dest="min_err_size", default=6, type=int,
                        help="minimum number of consecutive residues in a error before the svm labels it")
    parser.add_argument("--svm_threshold", dest="score_threshold", default=0.9, type=float,
                        help="the svm probability of error threshold for calling errors")
    parser.add_argument("--confidence_file", dest="conf_file", default=None, type=check_file_exists,
                        help="File containing confidences of prediction")
    parser.add_argument("--confidence_file_type", dest="conf_file_type", default=None, type=str,
                        help="type of file containing confidences of prediction")
    parser.add_argument("--take_plddt_from_distance_prediction", dest="PLDDT_IN_DISTFILE", default='yes', type=str,
                        help="Whether the predicted distfile supplies plddts (for example a pdb or mmciff file might have plddt scores in the bfactor collumn), if not but filters are asked for the confidence file could supply plddts")

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
    

def touch(fname, content='', mode='wb'):
    with open(fname, mode) as fhandle:
        fhandle.write(content)
    fhandle.close()


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
    sequencefile = conkit.io.read(args.seqfile, args.seqformat)
    sequence = sequencefile.top

    if len(sequence) < 5:
        raise ValueError('Cannot validate model with less than 5 residues')

    rep_atom, cutoff = set_contact_definition(args.moltype,rep_atom=args.rep_atom,cutoff=args.contact_distance_cutoff)

    logger.info("Length of the sequence:                      %d", len(sequence))
    logger.info("Reading input distance prediction:           %s", args.distfile)

    if args.distformat in ['pdb', 'mmcif']:
        prediction_file = conkit.io.read(args.distfile, args.distformat, distance_cutoff=cutoff, atom_type=rep_atom)
        prediction = prediction_file.top
    elif args.distformat in ['rosettanpz']:
        prediction_file = conkit.io.read(args.distfile, args.distformat, atom_type=rep_atom)
        prediction = prediction_file.top
    else: 
        prediction_file = conkit.io.read(args.distfile, args.distformat)
        prediction = prediction_file.top

    logger.info("Reading input PDB model:                     %s", args.pdbfile)

    if args.RENUMBER == 'yes':
        try:
            usable_model, alignment_dict, reverse_alignment_dict = write_renumbered_version_of_chain_in_struct(args.pdbfile,args.pdbformat,sequence,selected_chain=args.selected_chain,moltype=args.moltype)
        except:
            logger.critical("No sufficient sequence alignment was found between chains in: %s and %s check whether these are the right files and consider specifying the chain by setting --chain", args.pdbfile, args.seqfile)
    else: 
        usable_model = args.pdbfile
    
    model = conkit.io.read(usable_model, args.pdbformat, distance_cutoff=cutoff, atom_type=rep_atom).top

    if len(sequence) > 500:
        logger.info("Input model has more than 500 residues, this might take a while...")

    logger.info(os.linesep + "Validating model.")

    validation = conkit.plot.ModelValidationFigure(model, prediction, sequence)

    if args.RUN_SVM=='yes if prediction not pdb or mmcif': #don't run the svm if prediction is a structure by default
        if args.distformat in ['pdb', 'mmcif']:
            args.RUN_SVM='no'
        else:
            args.RUN_SVM='yes'

    if args.RUN_SVM=='yes':
        logger.info(os.linesep + "Running Support Vector Machine.")

        if args.moltype=='Protein':
            p = PDBParser()
            structure = p.get_structure('structure', usable_model)[0]
            ext_info = DSSP(structure, usable_model, dssp=args.dssp, acc_array='Wilke')
            validation.calculate_features()
        elif args.moltype=='RNA': 
            ext_info = areaimol_ACC(usable_model, args.pdbformat, args.areaimol_exe, tempfile_instructions_name='areaimol_acc_instructions.txt', tempfile_out_name='areaimol_log.log', gemmi_exe=args.gemmi_exe)
            validation.calculate_features(z_radius=20)
        else:
            ext_info = None

        if args.distformat in ['pdb', 'mmcif']:
            validation.svm(ext_info,moltype=args.moltype,prediction_type='STRUCT')
        else:
            validation.svm(ext_info,moltype=args.moltype,prediction_type='DIST')
        
        validation.svm_error_calling(min_err_size=args.min_err_size,score_threshold=args.score_threshold)
        

    if args.RUN_MAP_ALIGN=='yes':
        logger.info(os.linesep + "Running Map Align.")
        validation.map_align(map_align_exe=args.map_align_exe)


    if args.RUN_FILTERS=='yes':
        logger.info(os.linesep + "Running Filters.")

        validation.count_contacts(cutoff=cutoff)

        if (prediction.plddt != None) and (args.PLDDT_IN_DISTFILE == 'yes'): ##turn into check for plddt

            ## add a check to see if any external plddts were suplied
            validation.add_plddt()

        elif args.conf_file: #replace with a check to see if plddts can be taken from conf_file in future
            plddts = conkit.io.read(args.conf_file, args.conf_file_type)
            # validation.add_plddt(externally_supplied_plddts = A_list_from_conf_file)
            logger.info(os.linesep + "now plddts would be added.")
            
        if args.gesamt_exe and (args.distformat in ['pdb', 'mmcif']):

            validation.Run_gesamt_filter(usable_model, args.distfile, args.gesamt_exe, moltype=args.moltype, experimentfiletype=args.pdbformat)
            # identify potential errors
            logger.info(os.linesep + "added Q-scores")            
   
    logger.info(os.linesep + "Creating Figure.")
    validation.draw(RUN_SVM=(args.RUN_SVM=='yes'), RUN_MAP_ALIGN=(args.RUN_MAP_ALIGN=='yes'), RUN_FILTERS=(args.RUN_FILTERS=='yes'), svm_threshold=args.score_threshold)

    validation.savefig(args.output, overwrite=args.overwrite)
    logger.info(os.linesep + "Validation plot written to %s", args.output)

    residue_info = validation.data.loc[:, ['RESNUM', 'SCORE', 'MISALIGNED', 'PLDDT', 'CONTACTS', 'Q_IN_ERROR']]
    residue_info['NEW_REGISTER'] = ''

    table = PrettyTable()
    table.field_names = ["Residue", "Predicted score", "Suggested register", "plddt", "predicted contacts", "Q in error"]

    _resnum_template = '{} ({})'
    _error_score_template = '*** {0:.2f} ***'
    _correct_score_template = '    {0:.2f}    '
    _register_template = '*** {} ({}) ***'
    _empty_register = '               '

    for residue in residue_info.values:
        resnum, score, misalignment, plddt, contacts, Qs, register = residue
        current_residue = _resnum_template.format(sequence.seq[resnum - 1], resnum)
        score = _error_score_template.format(score) if score > 0.5 else _correct_score_template.format(score)

        if misalignment and resnum in validation.alignment.keys():
            register = _register_template.format(sequence.seq[validation.alignment[resnum] - 1], validation.alignment[resnum])
            residue_info.loc[residue_info['RESNUM'] == resnum, 'NEW_REGISTER'] = register
        else:
            register = _empty_register
            residue_info.loc[residue_info['RESNUM'] == resnum, 'NEW_REGISTER'] = register

        table.add_row([current_residue, score, register, plddt, contacts, Qs])

    ### add json format report ###

    if args.output_json:
        residue_info_json = residue_info.to_dict(orient='list')
        with open(args.output_json+".json", "w") as outfile:
            json.dump(residue_info_json, outfile)


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
