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
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.DSSP import DSSP
import os
import re
import subprocess
import tempfile
import json
import numpy as np
import pandas as pd
from prettytable import PrettyTable

import logging

import conkit.applications
import conkit.command_line
import conkit.io
from conkit.io.tools import set_contact_definition
import conkit.plot
from conkit.plot.tools import is_executable, areaimol_ACC
from conkit.misc import DNATCO_CATEGORIES
from conkit.misc.renumbering_tools import write_renumbered_version_of_chain_in_struct, detect_numbering_anomalies

logger = logging.getLogger(__name__)

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
    parser.add_argument("-dssp_exe", dest="dssp", default=None, help="path to dssp executable")#, type=is_executable)
    parser.add_argument("-output", dest="output", default="conkit.png", help="path to output figure png file", type=str)
    parser.add_argument("-output_json", dest="output_json", default=None, help="path to output json file", type=str)
    parser.add_argument("-outdir", dest="outdir", default=None, help="path to write created contactmaps to for debugging, if not specified maps get deleted", type=str)
    
    parser.add_argument("--overwrite", dest="overwrite", default=False, action="store_true",
                        help="overwrite output figure png file if it already exists")

    parser.add_argument("--map_align_exe", dest="map_align_exe", default=None,
                        type=is_executable, help="Path to the map_align executable")
    parser.add_argument("--gesamt_exe", dest="gesamt_exe", default=None,
                        type=is_executable, help="Path to the gesamt executable to check structural alignment")
    parser.add_argument("--areaimol_exe", dest="areaimol_exe", default=None,
                        type=is_executable, help="Path to areaimol executable to calculate solvent accesibility for RNA")
    parser.add_argument("--dnatco_exe", dest="dnatco_exe", default=None,
                        type=is_executable, help="Path to dnatco executable to calculate CANA categories for RNA")
    parser.add_argument("--gemmi_exe", dest="gemmi_exe", default=None,
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
    parser.add_argument("--cmo_filter", dest="cmo_filter_threshold", default=0.54, type=float,
                        help="the threshold value for the cmo filter (residues with a lower score are considered false positives)")
    parser.add_argument("--rf_filter", dest="rf_filter_threshold", default=0.76, type=float,
                        help="the threshold value for the RF filter (residues with a lower score are considered false positives)")

    parser.add_argument("--contact_dist", dest="contact_distance_cutoff", default=None, type=float,
                        help="distance cutoff for contacts when using Custom moltype (in angstrom)")
    parser.add_argument("--rep_atom", dest="rep_atom", default=None, type=str,
                        help="representative atom for contacts when using Custom moltype")

    parser.add_argument("--min_error_length", dest="min_err_size", default=6, type=int,
                        help="minimum number of consecutive residues in a error before the svm labels it")
    parser.add_argument("--svm_threshold", dest="score_threshold", default=0.5, type=float,
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
    

def calculate_dnatco(structfile, filetype, dnatco_exe):
    from collections import defaultdict

    if filetype != 'mmcif':
        if filetype == 'pdb':
            subprocess.run(['pdb2cif',structfile,structfile.replace('.cif','.pdb')])  # step to try to rescue pdb file entered should be changed because it silently introduces a dependency perhaps this can/should be replaced by a gemmin based call? also the swapping of the extensions is likely not a great way to do this
            structfile = structfile.replace('.cif','.pdb')
        else:
            logger.warning("DNATCO: unrecognised structure file type %r (expected 'pdb' or 'mmcif'); skipping DNATCO calculation.", filetype)
            return

    structfile = os.path.abspath(structfile)

    with tempfile.TemporaryDirectory() as tmpdir:
        subprocess.run(['node', dnatco_exe, '--coords', structfile, '--reportText'], cwd=tmpdir)

        with open(os.path.join(tmpdir, 'custom_report.txt'), 'r') as f:
            lines = f.readlines()

    for l in range(len(lines)):
        if  '|                               All dinucleotides                              |' in lines[l]:
            start = l+4
        if  "|                             Dinucleotide outliers                            |" in lines[l]:
            end = l-6


    records = []
    for l in range(start,end):
        parts = lines[l].split()
        nt1, nt2 = parts[3], parts[4]
        cana = parts[6]
        rmsd = float(parts[7])

        r1 = int(re.search(r"\d+", nt1).group())
        r2 = int(re.search(r"\d+", nt2).group())

        records.append((r1, r2, cana, rmsd))

    steps_df = pd.DataFrame(records, columns=["res1", "res2", "CANA", "RMSD"])

    # -----------------------------
    # Accumulate RMSDs per (residue, CANA)
    # -----------------------------
    rmsd_map = defaultdict(list)

    for _, row in steps_df.iterrows():
        rmsd_map[(row["res1"], row["CANA"])].append(row["RMSD"])
        rmsd_map[(row["res2"], row["CANA"])].append(row["RMSD"])

    # all residues present
    residues = sorted(
        set(steps_df["res1"]).union(steps_df["res2"])
    )

    # -----------------------------
    # Build final dataframe
    # -----------------------------
    res_df = pd.DataFrame(0.0, index=residues, columns=DNATCO_CATEGORIES)
    res_df.index.name = "RESNUM"

    for (res, cana), rmsds in rmsd_map.items():
        avg_rmsd = sum(rmsds) / len(rmsds)
        res_df.loc[res, cana] = 1.0 / avg_rmsd
        res_df.loc[res,'DNATCO_TOT_RMSD'] += sum(rmsds)

    return res_df


def touch(fname, content='', mode='wb'):
    with open(fname, mode) as fhandle:
        fhandle.write(content)
    fhandle.close()


def main():
    """The main routine for conkit-validate functionality"""
    parser = create_argument_parser()
    args = parser.parse_args()

    include_hetatms=(args.moltype == 'RNA')  #modified bases are very common in RNA strcutures, removing all hetatm records is a bad idea in this case

    conkit.command_line.setup_logging(level="info")

    if os.path.isfile(args.output) and not args.overwrite:
        raise FileExistsError('The output file {} already exists!'.format(args.output))

    logger.info(os.linesep + "Working directory:                           %s", os.getcwd())
    logger.info("Reading input sequence:                      %s", args.seqfile)
    sequencefile = conkit.io.read(args.seqfile, args.seqformat)
    sequence = sequencefile.top

    if len(sequence) < 5:
        raise ValueError('Cannot validate model with less than 5 residues')

    rep_atom, cutoff = set_contact_definition(args.moltype, rep_atom=args.rep_atom, cutoff=args.contact_distance_cutoff)

    logger.info("Length of the sequence:                      %d", len(sequence))
    logger.info("Reading input distance prediction:           %s", args.distfile)

    if args.distformat in ['pdb', 'mmcif']:
        prediction_file = conkit.io.read(args.distfile, args.distformat, distance_cutoff=cutoff, atom_type=rep_atom, include_hetatms=include_hetatms)
        prediction = prediction_file.top
    elif args.distformat in ['rosettanpz']:
        prediction_file = conkit.io.read(args.distfile, args.distformat, atom_type=rep_atom)
        prediction = prediction_file.top
        prediction.distance_cutoff = cutoff
    else: 
        prediction_file = conkit.io.read(args.distfile, args.distformat)
        prediction = prediction_file.top
        prediction.distance_cutoff = cutoff

    logger.info("Reading input PDB model:                     %s", args.pdbfile)

    if args.RENUMBER == 'yes':
        try:
            usable_model, alignment_dict, reverse_alignment_dict, original_map = write_renumbered_version_of_chain_in_struct(args.pdbfile, args.pdbformat, sequence, selected_chain=args.selected_chain, moltype=args.moltype)
        except Exception as e:
            logger.critical("Renumbering failed: %s", e)
            logger.critical("No sufficient sequence alignment was found between chains in %s and %s. "
                            "Check these are the right files; consider specifying --chain or disabling "
                            "renumbering with --renumber_model no.", args.pdbfile, args.seqfile)
            raise SystemExit(1)
    else:
        usable_model = args.pdbfile
        original_map = {}

    numbering_anomalies = detect_numbering_anomalies(original_map)

    model_file = conkit.io.read(usable_model, args.pdbformat, distance_cutoff=cutoff, atom_type=rep_atom, include_hetatms=include_hetatms)
    model = model_file.top
    model.distance_cutoff = cutoff

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

        if args.distformat in ['pdb','mmcif']:    # the maximal distance to account for in calculating the wRMSD when inputting a distogram is set by the lower bound of the highest bin, to mirror this with predicted structure we put it at 25 (slightly above where it would be with af2 distograms)
            max_distance = 25
        else: 
            max_distance = None

        if args.moltype=='Protein':
            if args.pdbformat == 'pdb':
                p = PDBParser()
            elif args.pdbformat == 'mmcif':
                p = MMCIFParser()
            else:
                logger.warning("Unrecognized structure file type %r being passed to DSSP.", args.pdbformat)
            try:
                structure = p.get_structure('structure', usable_model)
                ext_info = DSSP(structure[0], usable_model, dssp=args.dssp, acc_array='Wilke')
                secondary_structure_determination = 'DSSP'
            except Exception as e:
                logger.warning("DSSP failed; proceeding without secondary structure features. Error: %s", e)
                ext_info = None
                secondary_structure_determination = None
            validation.calculate_features(max_distance = max_distance)

        elif args.moltype=='RNA':
            try:
                ext_info = areaimol_ACC(usable_model, args.pdbformat, args.areaimol_exe, tempfile_instructions_name='areaimol_acc_instructions.txt', tempfile_out_name='areaimol_log.log', gemmi_exe=args.gemmi_exe)
            except Exception as e:
                logger.warning("areaimol ACC calculation failed; proceeding with ACC=0. Error: %s", e)
                ext_info = None
            if args.dnatco_exe:
                try:
                    dnatco_cats = calculate_dnatco(usable_model, args.pdbformat, args.dnatco_exe)
                    secondary_structure_determination = 'DNATCO'
                    ext_info = pd.merge(ext_info, dnatco_cats, how='outer', on='RESNUM')
                except Exception as e:
                    logger.warning("DNATCO calculation failed; proceeding without DNATCO features. Error: %s", e)
                    secondary_structure_determination = None
            else:
                secondary_structure_determination = None
            validation.calculate_features(z_radius = 20, max_distance = max_distance)
        else:
            ext_info = None

        if args.distformat in ['pdb', 'mmcif']:
            validation.svm(ext_info,moltype=args.moltype,prediction_type='STRUCT',sec_struc_info=secondary_structure_determination)
        else:
            validation.svm(ext_info,moltype=args.moltype,prediction_type='DIST',sec_struc_info=secondary_structure_determination)
        
        validation.svm_error_calling(min_err_size=args.min_err_size,score_threshold=args.score_threshold)
        

    if args.RUN_MAP_ALIGN=='yes':
        logger.info(os.linesep + "Running Map Align.")
        validation.map_align(map_align_exe=args.map_align_exe,temp_dir_name=args.outdir)


    if args.RUN_FILTERS=='yes':
        logger.info(os.linesep + "Running Filters.")

        validation.count_contacts()

        if (prediction.plddt != None) and (args.PLDDT_IN_DISTFILE == 'yes'):
            try:
                validation.add_plddt()
            except Exception as e:
                logger.warning("Failed to add pLDDT scores; pLDDT filter will be skipped. Error: %s", e)

        elif args.conf_file: #replace with a check to see if plddts can be taken from conf_file in future
            # TODO: implement reading pLDDT from external confidence file and passing to add_plddt()
            logger.warning("External confidence file supplied but reading pLDDT from external files is not yet implemented; pLDDT filter will be skipped.")
            
        if args.gesamt_exe and (args.distformat in ['pdb', 'mmcif']):
            try:
                validation.Run_gesamt_filter(usable_model, args.distfile, args.gesamt_exe, moltype=args.moltype, experimentfiletype=args.pdbformat)
                logger.info(os.linesep + "Added Q-scores.")
            except Exception as e:
                logger.warning("gesamt Q-score filter failed; Q-score filter will be skipped. Error: %s", e)

        if  {'PLDDT', 'CONTACTS', 'Q_IN_ERROR'}.issubset(validation.data.columns):
            # run the trained combination filters if all filter features calculated
            if args.RUN_MAP_ALIGN=='yes' and args.moltype=='RNA': 
                validation.Run_combined_filter(filter_type = 'CMO', filter_th = 0.54)

            if args.RUN_SVM=='yes' and args.moltype=='RNA': 
                validation.Run_combined_filter(filter_type = 'RF', filter_th = 0.76)
   
    logger.info(os.linesep + "Creating Figure.")
    validation.draw(RUN_SVM=(args.RUN_SVM=='yes'), RUN_MAP_ALIGN=(args.RUN_MAP_ALIGN=='yes'), RUN_FILTERS=(args.RUN_FILTERS=='yes'), svm_threshold=args.score_threshold, moltype=args.moltype, numbering_anomalies=numbering_anomalies)

    validation.savefig(args.output, overwrite=args.overwrite)
    logger.info(os.linesep + "Validation plot written to %s", args.output)

    residue_info = validation.data.loc[:, ['RESNUM', 'SCORE', 'MISALIGNED']]
    for filter_name in ['CMO_FILTER', 'RF_FILTER', 'PLDDT', 'CONTACTS', 'Q_IN_ERROR']:
        if filter_name in validation.data.columns:
            residue_info[filter_name] = validation.data.loc[:, filter_name]
        else:
            residue_info[filter_name] = ''

    residue_info['NEW_REGISTER'] = ''

    table = PrettyTable()
    table.field_names = ["Residue", "Predicted score", "Suggested register", "map align filter", "classifier filter","plddt", "predicted contacts", "Q in error"]

    _error_score_template = '*** {0:.2f} ***'
    _correct_score_template = '    {0:.2f}    '
    _register_template = '*** {} ({}) ***'
    _empty_register = '               '

    def _resnum_display(new_seq_id):
        """Format residue number, showing original→new when they differ."""
        orig_seq_id, orig_icode, _ = original_map.get((new_seq_id, ' '), (new_seq_id, ' ', ''))
        orig_str = f"{orig_seq_id}{orig_icode.strip()}"
        new_str = str(new_seq_id)
        return orig_str if orig_str == new_str else f"{orig_str}→{new_str}"

    def _is_numbering_anomaly(new_seq_id):
        """True if the residue at new_seq_id carried an insertion code in the input file."""
        _, orig_icode, _ = original_map.get((new_seq_id, ' '), (new_seq_id, ' ', ''))
        return orig_icode.strip() != ''

    # Collect rows as (sort_key, row_list) so EXTRA_CANONICAL annotation rows
    # can be merged in sequence order rather than appended at the end.
    table_rows = []

    for residue in residue_info.values:
        resnum, score, misalignment, cmo_filter, rf_filter, plddt, contacts, Qs, register = residue
        num_display = _resnum_display(resnum)
        residue_letter = sequence.seq[resnum - 1]
        if _is_numbering_anomaly(resnum):
            current_residue = f'*** {residue_letter} ({num_display}) ***'
        else:
            current_residue = f'{residue_letter} ({num_display})'

        score = _error_score_template.format(score) if score > args.score_threshold else _correct_score_template.format(score)
        if type(cmo_filter) in [int, float]:
            cmo_filter = _error_score_template.format(cmo_filter) if cmo_filter > args.cmo_filter_threshold else _correct_score_template.format(cmo_filter)
        if type(rf_filter) in [int, float]:
            rf_filter = _error_score_template.format(rf_filter) if rf_filter > args.rf_filter_threshold else _correct_score_template.format(rf_filter)

        if misalignment and resnum in validation.alignment.keys():
            register = _register_template.format(sequence.seq[validation.alignment[resnum] - 1], validation.alignment[resnum])
            residue_info.loc[residue_info['RESNUM'] == resnum, 'NEW_REGISTER'] = register
        else:
            register = _empty_register
            residue_info.loc[residue_info['RESNUM'] == resnum, 'NEW_REGISTER'] = register

        table_rows.append(((resnum, ' '), [current_residue, score, register, cmo_filter, rf_filter, plddt, contacts, Qs]))

    # Inject annotation-only rows for EXTRA_CANONICAL residues (no FASTA slot,
    # skipped by pdb.py — no metrics available, but user should see they exist).
    for new_seq_id, new_icode, orig_seq_id, orig_icode, resname, atype in numbering_anomalies:
        if atype == 'EXTRA_CANONICAL':
            orig_str = f"{orig_seq_id}{orig_icode.strip()}"
            new_str = f"{new_seq_id}{new_icode.strip()}"
            current_residue = f'*** {resname} ({orig_str}→{new_str}) ***'
            table_rows.append(((new_seq_id, new_icode), [current_residue, 'N/A', '', 'N/A', 'N/A', 'N/A', 'N/A', 'N/A']))

    table_rows.sort(key=lambda x: x[0])
    for _, row in table_rows:
        table.add_row(row)

    ### add json format report ###

    if args.output_json:
        residue_info_json = residue_info.to_dict(orient='list')
        residue_info_json['numbering_anomalies'] = [
            {
                'orig_seq_id': orig_seq_id,
                'orig_icode': orig_icode.strip(),
                'new_seq_id': new_seq_id,
                'new_icode': new_icode.strip(),
                'resname': resname,
                'type': atype,
            }
            for new_seq_id, new_icode, orig_seq_id, orig_icode, resname, atype in numbering_anomalies
        ]
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
