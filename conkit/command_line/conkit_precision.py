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
"""This script calculates the precision score for a contact prediction
when compared against contacts extracted from a protein structure.

"""

__author__ = "Felix Simkovic"
__date__ = "21 Nov 2016"
__version__ = "0.1"

import argparse

import conkit.command_line
import conkit.io

logger = None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument(
        '-c',
        dest='pdbchain',
        default=None,
        help='PDB chain to use [default: first in file]. '
        'Inter-molecular predictions use two letter '
        'convention, i.e AD for contacts between A and D.')
    parser.add_argument('-d', dest='dtn', default=5, type=int, help='Minimum sequence separation [default: 5]')
    parser.add_argument(
        '-f',
        dest='dfactor',
        default=1.0,
        type=float,
        help='number of contacts to include relative to sequence length [default: 1.0]')
    parser.add_argument('pdbfile')
    parser.add_argument('pdbformat')
    parser.add_argument('seqfile')
    parser.add_argument('seqformat')
    parser.add_argument('confile')
    parser.add_argument('conformat')
    args = parser.parse_args()

    global logger
    logger = conkit.command_line.setup_logging(level='info')

    if args.pdbchain:
        pdb = conkit.io.read(args.pdbfile, args.pdbformat)[args.pdbchain]
    else:
        pdb = conkit.io.read(args.pdbfile, args.pdbformat)[0]
    seq = conkit.io.read(args.seqfile, args.seqformat)[0]
    con = conkit.io.read(args.confile, args.conformat)[0]

    con.sequence = seq
    con.set_sequence_register()

    logger.info('Min sequence separation for contacting residues: %d', args.dtn)
    logger.info('Contact list cutoff factor: %f * L', args.dfactor)

    con.remove_neighbors(min_distance=args.dtn, inplace=True)
    ncontacts = int(seq.seq_len * args.dfactor)
    con.sort('raw_score', reverse=True, inplace=True)
    con_sliced = con[:ncontacts]

    con_matched = con_sliced.match(pdb)
    precision = con_matched.precision

    logger.info('Precision score: %f', precision)


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
