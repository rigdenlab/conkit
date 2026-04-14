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
"""This script provides a command-line interface to summmarise multple structures into a single contact probability map.

"""

__author__ = "Aderik Voorspoels"
__date__ = "03 december 2025"
__version__ = "0.13.3"

import argparse
import glob
import os


import conkit.command_line
import conkit.io
from conkit.io import DISTANCE_FILE_PARSERS
from conkit.io.tools import set_contact_definition
import conkit.plot
import conkit.plot.tools
from conkit.core import Contact, ContactMap, ContactFile
from conkit.core import Distance, Distogram, DistanceFile

logger = None

def create_argument_parser():
    """Create a parser for the command line arguments used in conkit-summarise"""

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("struct_files", type=str, help="expression to find the needed structure files")
    parser.add_argument("struct_format", type=str, help="structure format")


    parser.add_argument("--overlay_structure", dest="o_struct_fn", default=None, type=str,
                        help="Optional additional structure file to overlay on the summary contact map.")
    parser.add_argument("--overlay_struct_type", dest="o_struct_type", default=None, type=str,
                        help="File type of the addtional structure to be overlayed.")

    parser.add_argument("--other_files", dest="other_struct_files", default=None, type=str,
                        help="Optional additional structure files to summarise below diagonal.")
    parser.add_argument("--other_struct_type", dest="other_struct_type", default=None, type=str,
                        help="file type of additional structure files")

    parser.add_argument("--output", dest="output", default="conkit_summary.png", type=str,
                        help="File to save the summary figure to.")
    parser.add_argument("--output_contact_map", dest="contact_map_out", default="conkit_summary.mat", type=str,
                        help="File to save the summary contact matrix to.")
    parser.add_argument("--output_contact_map_format", dest="contact_map_out_format", default="mapalign", type=str,
                        help="Desired output format of the contact map")
    parser.add_argument("--color_map", dest="cmap", default="Greys", type=str,
                        help="matplotlib cmap for output figure.")    
    parser.add_argument("--overwrite", dest="overwrite", default=False, action="store_true",
                        help="overwrite output figure png file if it already exists")

    parser.add_argument("--moltype", dest="moltype", default="Protein", type=str,
                        help="Type of molecule [Protein,RNA].")
    parser.add_argument("--contact_dist", dest="contact_distance_cutoff", default=None, type=float,
                        help="distance cutoff for contacts when using Custom moltype.")
    parser.add_argument("--rep_atom", dest="rep_atom", default=None, type=str,
                        help="representative atom for contacts when using Custom moltype.")  

    return parser

def main():
    """The main routine for conkit-summarise functionality"""
    parser = create_argument_parser()
    args = parser.parse_args()

    global logger
    logger = conkit.command_line.setup_logging(level="info")

    if os.path.isfile(args.output) and not args.overwrite:
        raise FileExistsError('The output file {} already exists!'.format(args.output))

    logger.info(os.linesep + "Working directory:                           %s", os.getcwd())

    rep_atom, cutoff = set_contact_definition(args.moltype,rep_atom=args.rep_atom,cutoff=args.contact_distance_cutoff)

    struct_fn_list = glob.glob(args.struct_files)
    number_of_structures = len(struct_fn_list)

    logger.info(os.linesep + f"Found {number_of_structures} structure files to compile.")

    #summary_contact_map = ContactMap("summary")
    summary_contact_dict = {}

    for fn in struct_fn_list:
        logger.info(os.linesep + f"extracting contacts from {fn}.")
    
        file = conkit.io.read(fn, args.struct_format, distance_cutoff=cutoff, atom_type=rep_atom)
        print(file)
        struct = file.top
        contact_map = struct.as_contactmap( distance_cutoff=cutoff )
        contact_set = contact_map.as_set()

        for contact in contact_set:
            if contact in summary_contact_dict.keys():
                summary_contact_dict[contact] += 1
            else:
                summary_contact_dict[contact] =1

    summary_contact_map = ContactMap("summary")

    logger.info(os.linesep + f"Creating summary contactmap.")

    for contact in summary_contact_dict.keys():
        summary_contact_dict[contact] = summary_contact_dict[contact]/number_of_structures
        summary_contact_map.add(Contact(contact[0], contact[1], summary_contact_dict[contact]) )

    if args.other_struct_files and args.other_struct_type:
        other_struct_fn_list = glob.glob(args.other_struct_files)
        number_of_structures = len(other_struct_fn_list)
        logger.info(os.linesep + f"Found {number_of_structures} structure files to compile below the diagonal.")
        other_summary_contact_map = summarise_files_to_single_map(other_struct_fn_list,args.other_struct_type,moltype=args.moltype,rep_atom=args.rep_atom,cutoff=args.contact_distance_cutoff)
    else: 
        other_summary_contact_map = None

    logger.info(os.linesep + "Creating Figure.")
    if args.o_struct_fn and args.o_struct_type:
        overlay_file = conkit.io.read(args.o_struct_fn, args.o_struct_type, distance_cutoff=cutoff, atom_type=rep_atom)
        overlay_map = (overlay_file.top).as_contactmap(distance_cutoff=cutoff)

    else: 
        overlay_map = None

    logger.info(os.linesep + "Creating Figure.")
    summary_plot = conkit.plot.ContactMapMatrixFigure( summary_contact_map, other= other_summary_contact_map, overlay=overlay_map, cmap=args.cmap )

    summary_plot.savefig(args.output, overwrite=args.overwrite)
    logger.info(os.linesep + "Validation plot written to %s", args.output)

    if args.contact_map_out:
        contact_mat_out_file = ContactFile("out_matrix_file")
        contact_mat_out_file.add(summary_contact_map)
        conkit.io.write(args.contact_map_out, args.contact_map_out_format, contact_mat_out_file)

def summarise_files_to_single_map(fns,f_type,moltype='Protein',rep_atom=None,cutoff=None):

    number_of_structures = len(fns)
    rep_atom, cutoff = set_contact_definition(moltype,rep_atom=rep_atom,cutoff=cutoff)
    summary_contact_dict = {}

    for fn in fns:
        logger.info(os.linesep + f"extracting contacts from {fn}.")
    
        file = conkit.io.read(fn, f_type, distance_cutoff=cutoff, atom_type=rep_atom)
        contact_map = (file.top).as_contactmap( distance_cutoff=cutoff )
        contact_set = contact_map.as_set()

        for contact in contact_set:
            if contact in summary_contact_dict.keys():
                summary_contact_dict[contact] += 1
            else:
                summary_contact_dict[contact] =1

        summary_contact_map = ContactMap("summary")

    logger.info(os.linesep + f"Creating summary contactmap.")

    for contact in summary_contact_dict.keys():
        summary_contact_dict[contact] = summary_contact_dict[contact]/number_of_structures
        summary_contact_map.add(Contact(contact[0], contact[1], summary_contact_dict[contact]) )

    return summary_contact_map

if __name__ == "__main__":
    main()
