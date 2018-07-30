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
"""Command line facility for ConKit scripts"""

__author__ = "Felix Simkovic"
__date__ = "14 Apr 2017"
__version__ = "0.1"

import logging
import os
import sys


def setup_logging(level='info', logfile=None):
    """Set up logging to the console for the root logger.

    Parameters
    ----------
    level : str, optional
       The console logging level to be used [default: info]
       To change, use one of
           [ notset | info | debug | warning | error | critical ]
    logfile : str, optional
       The path to a full file log

    Returns
    -------
    :obj:`~logging.Logger`
       Instance of a :obj:`~logging.Logger`

    """

    class ColorFormatter(logging.Formatter):
        """Formatter to color console logging output"""

        # ANSI foreground color codes
        colors = {
            logging.DEBUG: 34,  # blue
            logging.WARNING: 33,  # yellow
            logging.ERROR: 31,  # red
            logging.CRITICAL: 31,  # red
        }

        def format(self, record):
            if record.levelno in self.colors:
                prefix = '\033[1;{}m'.format(ColorFormatter.colors[record.levelno])
                postfix = '\033[0m'
                record.msg = os.linesep.join([prefix + l + postfix for l in str(record.msg).splitlines()])
            return logging.Formatter.format(self, record)

    # Reset any Handlers or Filters already in the logger to start from scratch
    # https://stackoverflow.com/a/16966965
    map(logging.getLogger().removeHandler, logging.getLogger().handlers[:])
    map(logging.getLogger().removeFilter, logging.getLogger().filters[:])

    logging_levels = {
        'notset': logging.NOTSET,
        'info': logging.INFO,
        'debug': logging.DEBUG,
        'warning': logging.WARNING,
        'error': logging.ERROR,
        'critical': logging.CRITICAL
    }

    # Create logger and default settings
    logging.getLogger().setLevel(logging.DEBUG)

    # create console handler with a higher log level
    ch = logging.StreamHandler(stream=sys.stdout)
    ch.setLevel(logging_levels.get(level, logging.INFO))
    ch.setFormatter(ColorFormatter('%(message)s'))
    logging.getLogger().addHandler(ch)

    # create file handler which logs even debug messages
    if logfile:
        fh = logging.FileHandler(logfile)
        fh.setLevel(logging.NOTSET)
        fh.setFormatter(logging.Formatter('%(asctime)s\t%(name)s [%(lineno)d]\t%(levelname)s\t%(message)s'))
        logging.getLogger().addHandler(fh)

    logging.getLogger().debug('Console logger level: %s', logging_levels.get(level, logging.INFO))
    logging.getLogger().debug('File logger level: %s', logging.NOTSET)

    return logging.getLogger()
