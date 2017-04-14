#!/usr/bin/env python
"""This script analyses a Multiple Sequence Alignment.

It provides you estimates about the quality, which will allow
you to determine its usefulness for covariance-based contact
prediction.

"""

__author__ = "Felix Simkovic"
__date__ = "21 Nov 2016"
__version__ = "0.1"

import argparse
import sys

import conkit.command_line
import conkit.io
import conkit.plot

logger = conkit.command_line.get_logger('msatool', level='info')


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

    logger.info('Input MSA File:                            %s', args.msafile)
    logger.info('Input MSA Format:                          %s', args.msaformat)
    logger.info('Pairwise Sequence Identity Threshold:      %f', args.id)
    logger.info('Length of the Target Sequence:             %d', seq_len)
    logger.info('Total Number of Sequences:                 %d', nseqs)
    logger.info('Number of Effective Sequences:             %d', meff)
    logger.info('Sequence Coverage Plot:                    %s', plot)

    return


if __name__ == "__main__":
    sys.exit(main())
