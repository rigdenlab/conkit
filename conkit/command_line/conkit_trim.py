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
"""This script provides a command-line interface to trim structures using a single contact probability map.

"""

__author__ = "Aderik Voorspoels"
__date__ = "11 february 2026"
__version__ = "0.14.1"

import argparse
import inspect
import glob
import fnmatch
import re
import os

import conkit.command_line
import conkit.io
from conkit.io import DISTANCE_FILE_PARSERS
from conkit.io.tools import set_contact_definition
import conkit.plot
import conkit.plot.tools
from conkit.core import Contact, ContactMap, ContactFile


logger = None

def create_argument_parser():
    """Create a parser for the command line arguments used in conkit-trim"""

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("contact_map_summary", type=str, help="contact map to use")
    parser.add_argument("contact_map_format", type=str, help="format of that map")
    parser.add_argument("struct_files", type=str, help="expression to find the needed structure files")
    parser.add_argument("struct_format", type=str, help="structure format")


    parser.add_argument("--score_threshold", dest="threshold", default=0.9, type=float,
                        help="score threshold above which contact will be kept")

    parser.add_argument("--output_prefix", dest="output", default="trimmed", type=str,
                        help="prefix to attach to trimmed structure file")
    parser.add_argument("--overwrite", dest="overwrite", default=False, action="store_true",
                        help="overwrite output file if it already exists")

    parser.add_argument("--moltype", dest="moltype", default="Protein", type=str,
                        help="Type of molecule [Protein,RNA].")
    parser.add_argument("--contact_dist", dest="contact_distance_cutoff", default=None, type=float,
                        help="distance cutoff for contacts when using Custom moltype.")
    parser.add_argument("--rep_atom", dest="rep_atom", default=None, type=str,
                        help="representative atom for contacts when using Custom moltype.")  

    return parser

def main():
    """The main routine for conkit-trim functionality"""
    parser = create_argument_parser()
    args = parser.parse_args()

    global logger
    logger = conkit.command_line.setup_logging(level="info")

    if os.path.isfile(args.output) and not args.overwrite:
        raise FileExistsError('The output file {} already exists!'.format(args.output))

    logger.info(os.linesep + "Working directory:                           %s", os.getcwd())

    logger.info(os.linesep + "Reading the template contact map from:                           %s", args.contact_map_summary)

    rep_atom, cutoff = set_contact_definition(args.moltype,rep_atom=args.rep_atom,cutoff=args.contact_distance_cutoff)

    if args.contact_map_format in ['pdb', 'mmcif']:
        contact_file = conkit.io.read(args.contact_map_summary, args.contact_map_format, distance_cutoff=cutoff, atom_type=rep_atom)
    if args.contact_map_format in DISTANCE_FILE_PARSERS:
        contact_file = conkit.io.read(args.contact_map_summary, args.contact_map_format)
        if args.rep_atom:
            logger.info(os.linesep + "It looks looks like you are manually adjusting the representative atom for contact definition.")
            logger.info(os.linesep + "However the contactmap you passed did'nt retain more than one atom's information.")
            logger.info(os.linesep + "This change to the representative atom is being ignored.")
            logger.info(os.linesep + "A chage to the distance cuttoff done explicitly or through the molecule type will however be respected")
    else:
        contact_file = conkit.io.read(args.contact_map_summary, args.contact_map_format)
        if args.rep_atom or args.contact_distance_cutoff or (args.moltype != 'Protein'):
            logger.info(os.linesep + "It looks looks like you are manually adjusting the contact definition or specifying the moltype.")
            logger.info(os.linesep + "However the contactmap you passed did'nt retain the information needed to adjust those.")
            logger.info(os.linesep + "This change is being ignored, willl use the info provided in %s.", args.contact_map_summary)
            logger.info(os.linesep + "Consider changing the contents of that file if you want to change the contact definition", args.contact_map_summary)

    if args.contact_map_format in DISTANCE_FILE_PARSERS:
        contact_map = (contact_file.top).as_contactmap( distance_cutoff=cutoff )
    else:
        contact_map = (contact_file.top)

    contact_map = contact_map.remove_neighbors(min_distance=5)
    contact_set = contact_map.as_set()

    logger.info(os.linesep + "Building list of residues to keep.")

    residues_to_keep = set()
    conserved_contacts = 0
    for contact in contact_map:
        if contact.raw_score >= args.threshold:
            conserved_contacts +=1
            residues_to_keep.add(contact.res1_seq)
            residues_to_keep.add(contact.res2_seq)

    
    logger.info(os.linesep + "Found %d residues to keep from %d conserved contacts.", len(residues_to_keep),conserved_contacts)

    pattern = args.struct_files
    struct_fn_list = glob.glob(pattern)
    number_of_structures = len(struct_fn_list)

    # make regex version of structurefile search string to extract output names
    regex_pattern = fnmatch.translate(pattern)
    regex_pattern = regex_pattern.replace('.*', '(.*)')
    regex = re.compile(regex_pattern)

    logger.info(os.linesep + f"Found {number_of_structures} structure files to trim.")

    if args.struct_format == 'pdb':
        from Bio.PDB.PDBParser import PDBParser as STRUCTParser
        from Bio.PDB.PDBIO import PDBIO as STRUCTIO
        from Bio.PDB.PDBIO import Select
    elif args.struct_format == 'mmcif':
        from Bio.PDB.MMCIFParser import MMCIFParser
        from Bio.PDB.mmcifio  import MMCIFIO as STRUCTIO
        from Bio.PDB.mmcifio import Select
    else:
        raise ValueError('Unrecognized structure file type {}'.format(args.struct_format))

    for fn in struct_fn_list:

        # read in the structure
        parser = STRUCTParser()
        s = parser.get_structure('structure',fn)
        io = STRUCTIO()
        io.set_structure(s)
        
        class Select_from_set(Select):
            def accept_residue(self, residue):
                if residue.get_id()[1] in residues_to_keep:
                    return 1
                else:
                    return 0

        # build output names for the trimmed structures
        # currently this is done by using a coomn prefix (including file location) provided by the user
        # to this we append the wildcards that were detected when looking for input files
        # this is slightly messy and might lead to not the best output naming scheme 
        # but it allows one to match trimed structures to the input structure 
        # and for the user to specify at least some of the output location/name

        match = regex.match(fn)
        out_fn = args.output
        for m in match.groups():
            out_fn += '_'+m
        out_fn += '.'+args.struct_format

        io.save(out_fn, Select_from_set())

    logger.info(os.linesep + f"Trimmed {number_of_structures} structure files and saved them under {args.output}.")
        

if __name__ == "__main__":
    main()
