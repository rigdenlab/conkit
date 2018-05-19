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
"""This script analyses a Multiple Sequence Alignment.

It provides you estimates about the quality, which will allow
you to determine its usefulness for covariance-based contact
prediction.

"""

__author__ = "Felix Simkovic"
__date__ = "21 Nov 2016"
__version__ = "0.1"

import argparse

import conkit.command_line
import conkit.io
import conkit.plot
import conkit.plot.tools

logger = None


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('msafile', help='Multiple Sequence Alignment file')
    parser.add_argument('msaformat', help='Multiple Sequence Alignment format')
    args = parser.parse_args()

    global logger
    logger = conkit.command_line.setup_logging(level='info')

    msa = conkit.io.read(args.msafile, args.msaformat)

    plot = args.msafile.rsplit('.', 1)[0] + '.png'
    figure = conkit.plot.SequenceCoverageFigure(msa, legend=True)
    figure.ax.set_aspect(conkit.plot.tools.get_adjusted_aspect(figure.ax, 0.3))
    figure.savefig(plot)

    logger.info('Input MSA File:                            %s', args.msafile)
    logger.info('Input MSA Format:                          %s', args.msaformat)
    logger.info('Length of the Target Sequence:             %d', msa.top_sequence.seq_len)
    logger.info('Total Number of Sequences:                 %d', msa.nseq)
    logger.info('Number of Effective Sequences:             %d', msa.meff)
    logger.info('Sequence Coverage Plot:                    %s', plot)


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
