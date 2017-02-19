#!/usr/bin/env python
"""This script provides direct access to ConKit's conversion algorithms.

This script can convert either contact prediction files or sequence files.
In case of the latter, a file with a single or multiple sequences can be
converted.

!!! IMPORTANT
=============
Do not attempt to mix formats, i.e. convert from a contact file format
to a sequence file format.

"""

__author__ = "Felix Simkovic"
__date__ = "01 Oct 2016"
__version__ = 0.1

import argparse
import conkit
import sys

_OPTIONS = sorted(
    conkit.io.CONTACT_FILE_PARSERS.keys() 
    + conkit.io.SEQUENCE_FILE_PARSERS.keys()
)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('infile')
    parser.add_argument('informat')
    parser.add_argument('outfile')
    parser.add_argument('outformat')
    args = parser.parse_args()
    
    # Perform the conversion
    conkit.io.convert(args.infile, args.informat, args.outfile, args.outformat)

    return 0


if __name__ == "__main__":
    sys.exit(main())
