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
from collections import namedtuple
import os
from operator import attrgetter
from prettytable import PrettyTable
import shutil

import conkit.applications
import conkit.command_line
import conkit.io
import conkit.plot

logger = None

Outlier = namedtuple('Outlier', ("resnum", "wrmsd", "fn_count", "misalignment"))
Alignment = namedtuple('Alignment', ("residue_range", "residue_pairs"))


def create_argument_parser():
    """Create a parser for the command line arguments used in conkit-validate"""

    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("seqfile", type=check_file_exists, help="Path to sequence file")
    parser.add_argument("seqformat", type=str, help="Sequence format")
    parser.add_argument("distfile", type=check_file_exists, help="Path to distance prediction file")
    parser.add_argument("distformat", type=str, help="Format of distance prediction file")
    parser.add_argument("pdbfile", type=check_file_exists, help="Path to structure file")
    parser.add_argument("pdbformat", type=str, help="Format of structure file")
    parser.add_argument("--outdir", dest="outdir", type=str, help="Output directory",
                        default=os.path.join(os.getcwd(), 'conkit-validate'))
    parser.add_argument("--overwrite", dest="overwrite", default=False, action="store_true",
                        help="overwrite output directory if exists")
    parser.add_argument("--skip_alignment", dest="skip_alignment", default=False, action="store_true",
                        help="skip contact map alignment step")
    parser.add_argument("--map_align_exe", dest="map_align_exe", default=None,
                        type=check_file_exists, help="Path to the map_align executable")
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


def prepare_output_directory(outdir, overwrite=False):
    """Prepare the output directory for conkit-validate.

    Parameters
    ----------
    outdir : str
       Path to the output directory
    overwrite : bool
       Whether the output directory should be overwritten or not [default: False]

    Raises
    ------
    :exc:`ValueError`
        The output directory already exists and overwrite is False
    """

    if os.path.isdir(outdir):
        if not overwrite:
            raise ValueError('Output directory already exists')
        else:
            shutil.rmtree(outdir)
    os.mkdir(outdir)


def get_residue_ranges(resnums):
    nums = sorted(set(resnums))
    gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s + 3 < e]
    edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
    return list(zip(edges, edges))


def parse_map_align_stdout(stdout):
    """Parse the stdout of map_align and extract the alignment of residues.

    Parameters
    ----------
    stdout : str
       Standard output created with map_align

    Returns
    ------
    alignments: list
        A list of :obj:`~conkit.command_line.conkit_validate.Alignment` that described regions of at
        least 5 misaligned residues
    """

    residues_map_a = []
    residues_map_b = []
    resnum = 0
    misalinged_residues = []

    for line in stdout.split('\n'):
        if line and line.split()[0] == "MAX":
            line = line.rstrip().lstrip().split()
            for residue_pair in line[8:]:
                resnum += 1
                residue_pair = residue_pair.split(":")
                residues_map_a.append(int(residue_pair[0]))
                residues_map_b.append(int(residue_pair[1]))
                if residue_pair[0] != residue_pair[1]:
                    misalinged_residues.append(resnum)

    alignments = []

    for start, stop in get_residue_ranges(misalinged_residues):
        alingment_slice = slice(start - 1, stop)
        alingment_range = range(residues_map_a[start - 1], residues_map_a[stop - 1] + 1)
        residue_pairs = []

        if len(alingment_range) < 5:
            continue
        for res_a, res_b in zip(residues_map_a[alingment_slice], residues_map_b[alingment_slice]):
            residue_pairs.append((res_a, res_b))

        new_alignment = Alignment(alingment_range, residue_pairs)
        alignments.append(new_alignment)

    return alignments


def get_outliers(figure_outliers, rmsd_profile, fn_profile, alignment_list=None):
    outliers = []
    if alignment_list is not None:
        _copy_alignment_list = alignment_list.copy()
    else:
        _copy_alignment_list = alignment_list = []

    for idx, resnum in enumerate(figure_outliers, 1):
        outlier_alignment = []
        for alignment in alignment_list:
            if any([x in alignment.residue_range for x in (resnum, resnum + 10, resnum - 10)]):
                outlier_alignment += alignment.residue_pairs
                if alignment in _copy_alignment_list:
                    _copy_alignment_list.remove(alignment)
        outlier = Outlier(resnum, rmsd_profile[resnum], fn_profile[resnum], outlier_alignment)
        outliers.append(outlier)

    for alignment in _copy_alignment_list:
        resnum = alignment.residue_range[0]
        outlier = Outlier(resnum, rmsd_profile[resnum], fn_profile[resnum], alignment.residue_pairs)
        outliers.append(outlier)

    return sorted(outliers, key=attrgetter('resnum'))


def main():
    """The main routine for conkit-validate functionality"""
    parser = create_argument_parser()
    args = parser.parse_args()

    global logger
    logger = conkit.command_line.setup_logging(level="info")

    logger.info("Output directory:                           %s", args.outdir)

    prepare_output_directory(args.outdir, args.overwrite)

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

    if args.skip_alignment or args.map_align_exe is None:
        logger.info("Skipping contact map alignment, no fixes will be suggested.")
        outliers = get_outliers(figure.outliers, figure.rmsd_profile, figure.fn_profile)
    else:
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

        logger.info(os.linesep + "Executing: %s", map_align_cline)
        stdout, stderr = map_align_cline()
        map_align_log = os.path.join(args.outdir, 'map_align.log')
        with open(map_align_log, 'w') as fhandle:
            fhandle.write(stdout)

        alignments = parse_map_align_stdout(stdout)
        outliers = get_outliers(figure.outliers, figure.rmsd_profile, figure.fn_profile, alignments)

    if any(outliers):
        table = PrettyTable()
        table.field_names = ["Outlier no.", "Residue no.", "wRMSD", "FN Count", "Misalignment"]
        for idx, outlier in enumerate(outliers, 1):
            table.add_row([idx, outlier.resnum, '{0:.2f}'.format(outlier.wrmsd),
                           '{0:.2f}'.format(outlier.fn_count), 'Yes' if outlier.misalignment else 'No'])
        logger.info(os.linesep + "List of detected outliers:")
        logger.info(table)

        for idx, outlier in enumerate(outliers, 1):
            if not outlier.misalignment:
                logger.info(os.linesep + "Cannot find optimal re-alignment for outlier no. {}.".format(idx))
                continue
            table = PrettyTable()
            table.field_names = ["Current Residue", "New Residue"]
            for residue_pairs in outlier.misalignment:
                table.add_row(['{} ({})'.format(sequence.seq[residue_pairs[1] - 1], residue_pairs[1]),
                               '{} ({})'.format(sequence.seq[residue_pairs[0] - 1], residue_pairs[0])])
            logger.info(os.linesep + "List of proposed changes to fix outlier no. {}:".format(idx))
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
