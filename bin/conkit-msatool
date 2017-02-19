#!/usr/bin/env python
"""This script analyses a Multiple Sequence Alignment.

It provides you estimates about the quality, which will allow
you to determine its usefulness for covariance-based contact
prediction.

"""

__author__ = "Felix Simkovic"
__date__ = "21 Nov 2016"
__version__ = 0.1

import argparse
import conkit
import logging
import sys

logging.basicConfig(format='%(message)s', level=logging.INFO)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('-id', default=0.7, type=float, help='sequence identity [default: 0.7]')
    parser.add_argument('msafile', help='Multiple Sequence Alignment file')
    parser.add_argument('msaformat', help='Multiple Sequence Alignment format')
    args = parser.parse_args()
    
    # Compute all the data
    hierarchy = conkit.io.read(args.msafile, args.msaformat)
    seq_len = hierarchy.top_sequence.seq_len
    nseqs = hierarchy.nseqs
    meff = hierarchy.calculate_meff(identity=args.id)
    plot = args.msafile.rsplit('.', 1)[0] + '.png'
    conkit.plot.SequenceCoverageFigure(hierarchy, file_name=plot)

    logging.info('Input MSA File:                            {0}'.format(args.msafile))
    logging.info('Input MSA Format:                          {0}'.format(args.msaformat))
    logging.info('Pairwise Sequence Identity Threshold:      {0}'.format(args.id))
    logging.info('Length of the Target Sequence:             {0}'.format(seq_len))
    logging.info('Total Number of Sequences:                 {0}'.format(nseqs))
    logging.info('Number of Effective Sequences:             {0}'.format(meff))
    logging.info('Sequence Coverage Plot:                    {0}'.format(plot))

    return 0


if __name__ == "__main__":
    sys.exit(main())
