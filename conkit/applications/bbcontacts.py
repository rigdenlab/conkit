# coding=utf-8
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
"""
Command line object for bbcontacts contact filtering application
"""

__author__ = "Felix Simkovic"
__date__ = "10 Aug 2016"
__version__ = "0.1"

from Bio.Application import _Argument
from Bio.Application import _Option
from Bio.Application import _Switch
from Bio.Application import AbstractCommandline


class BbcontactsCommandline(AbstractCommandline):
    """
    Command line object for bbcontacts [#]_ contact filtering application

    https://github.com/soedinglab/bbcontacts

    The bbcontacts program is a Python program predicting residue-level
    contacts between beta-strands by detecting patterns in matrices of
    predicted couplings. bbcontacts can make use of a secondary structure
    assignment or a secondary structure prediction.
    
    .. [#] Andreani J., Söding J. (2015). bbcontacts: prediction of beta-strand
       pairing from direct coupling patterns. Bioinformatics 31(11), 1729-1737.


    Examples
    --------
    To filter a contact map using a Multiple Sequence Alignment in
    CCMpred format, use:

    >>> from conkit.applications import BbcontactsCommandline
    >>> bbcontacts_cline = BbcontactsCommandline(
    ...     matfile='test.mat', diversity_score=0.482, prefix='test'
    ... )
    >>> print(bbcontacts_cline)
    bbcontacts

    You would typically run the command line with :func:`bbcontacts_cline` or via
    the :mod:`~subprocess` module.

    Note
    ----
    Installation instructions are available via the `GitHub repository 
    <https://github.com/soedinglab/bbcontacts>`_.

    """

    def __init__(self, cmd="bbcontacts", **kwargs):

        # TODO: figure a way to group CL arguments as in `mutually_exclusive_group`
        if 'dssp_file' in list(kwargs.keys()) and 'psipred_file' in list(kwargs.keys()):
            msg = 'Provide only one of [dssp_file|psipred_file]!'
            raise RuntimeError(msg)
        elif not ('dssp_file' in list(kwargs.keys()) or 'psipred_file' in list(kwargs.keys())):
            msg = 'Provide one of [dssp_file|psipred_file]!'
            raise RuntimeError(msg)

        self.parameters = [
            _Option(['-c', 'config_file'], 'bbcontacts configuration file', filename=True, equate=False),
            _Option(
                ['-s', 'smoothing_size'],
                'Perform local background correction of the coupling matrix '
                'before decoding: from each coupling, subtract the average '
                'coupling (smoothed background) over an area extending by '
                'SMOOTHINGSIZE in each direction [default=10, use 0 for no '
                'local background correction]',
                equate=False),
            _Switch(['-l', 'long_predictions'], 'Turn off (slow) prediction-shortening mode (this mode is on '
                    'by default but will only get triggered when long predictions occur)'),
            _Option(
                ['-n', 'pdb_name'],
                'Provide a PDB identifier (when also using -e, this will be the '
                'PDB name to look for in EVALUATIONFILE)',
                equate=False),
            _Option(
                ['-e', 'evaluation_file'],
                'Provide a file containing the true contacts (BetaSheet916.dat, '
                'BetaSheet1452.dat or same format) for evaluation',
                filename=True,
                equate=False),
            _Argument(['matfile'], 'CCMpred-like coupling matrix', filename=True, is_required=True),
            _Argument(['diversity_score'], 'sequence-dependent diversity score', is_required=True),
            _Argument(['prefix'], 'output prefix', is_required=True),
            _Option(['-d', 'dssp_file'], 'DSSP secondary structure prediction file', filename=True, equate=False),
            _Option(['-p', 'psipred_file'], 'PSIPRED secondary structure prediction file', filename=True, equate=False),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
