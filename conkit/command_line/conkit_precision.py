#!/usr/bin/env python
"""This script calculates the precision score for a contact prediction
when compared against contacts extracted from a protein structure.

"""

__author__ = "Felix Simkovic"
__date__ = "21 Nov 2016"
__version__ = "0.1"

import argparse
import sys

import conkit.command_line
import conkit.io

logger = conkit.command_line.get_logger('precision', level='info')


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-c', dest='pdbchain', default=None,
                        help='PDB chain to use [default: first in file]. '
                             'Inter-molecular predictions use two letter '
                             'convention, i.e AD for contacts between A and D.')
    parser.add_argument('-d', dest='dtn', default=5, type=int,
                        help='Minimum sequence separation [default: 5]')
    parser.add_argument('-f', dest='dfactor', default=1.0, type=float,
                        help='number of contacts to include relative to sequence length [default: 1.0]')
    parser.add_argument('pdbfile')
    parser.add_argument('pdbformat')
    parser.add_argument('seqfile')
    parser.add_argument('seqformat')
    parser.add_argument('confile')
    parser.add_argument('conformat')
    args = parser.parse_args()
    
    # Compute all the data
    if args.pdbchain:
        pdb = conkit.io.read(args.pdbfile, args.pdbformat)[args.pdbchain]
    else:
        pdb = conkit.io.read(args.pdbfile, args.pdbformat)[0]
    seq = conkit.io.read(args.seqfile, args.seqformat)[0]
    con = conkit.io.read(args.confile, args.conformat)[0]
    
    con.sequence = seq
    con.assign_sequence_register()

    logger.info('Min sequence separation for contacting residues: %d', args.dtn)
    logger.info('Contact list cutoff factor: %f * L', args.dfactor)

    con.remove_neighbors(min_distance=args.dtn, inplace=True)
    ncontacts = int(seq.seq_len * args.dfactor)
    con.sort('raw_score', reverse=True, inplace=True)
    con_sliced = con[:ncontacts]


    con_matched = con_sliced.match(pdb)
    precision = con_matched.precision
    
    logger.info('Precision score: %f', precision)

    return 

if __name__ == "__main__":
    sys.exit(main())
