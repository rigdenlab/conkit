# BSD 3-Clause License
#
# Copyright (c) 2016-19, University of Liverpool
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
"""Energy function templates for restraint generation"""

__author__ = "Felix Simkovic"
__date__ = "13 Aug 2018"
__version__ = "1.0"


class RosettaFunctionConstructs(object):
    """Storage for string formats of different Rosetta energy function constructs

    For more information on the different energy functions, please refer to the
    corresponding references or the official `RosettaCommons documentation
    <https://www.rosettacommons.org/docs/latest/rosetta_basics/file_types/constraint-file>`_

    """

    _ATOMPAIR = "AtomPair {atom1: >2} {res1_seq: >4} {atom2: >2} {res2_seq: >4} "
    _SCALARWEIGHTED = "SCALARWEIGHTEDFUNC {scalar_score: .3f} "

    @property
    def BOUNDED_default(self):
        """Simple bounded energy function"""
        sstream = RosettaFunctionConstructs._ATOMPAIR
        sstream += "BOUNDED {lower_bound: >.3f} {upper_bound: >.3f} 1 0.5 #"
        return sstream

    @property
    def BOUNDED_gremlin(self):
        """Energy function according to [#]_

        References
        ----------
        .. [#] Ovchinnekov et al. (2015). Large-scale determination of previously unsolved
           protein structures using evolutionary information. Elife 3(4), e09248.

        """
        sstream = RosettaFunctionConstructs._ATOMPAIR
        sstream += RosettaFunctionConstructs._SCALARWEIGHTED
        sstream += "BOUNDED 0 {lower_bound: >.3f} 1 0.5"
        return sstream

    @property
    def FADE(self):
        """Energy function according to [#]_ and [#]_

        References
        ----------
        .. [#] Simkovic et al. (2016). Residue contacts predicted by evolutionary covariance
           extend the application of ab initio molecular replacement to larger and more
           challenging protein folds. IUCrJ 3(Pt 4), 259-270.
        .. [#] Michel et al. (2014). PconsFold: improved contact predictions improve protein
           models. Bioinformatics 30(17), i482-i488

        """
        sstream = RosettaFunctionConstructs._ATOMPAIR
        sstream += "FADE -10 19 10 {energy_bonus: >5.2f} 0"
        return sstream

    @property
    def FADE_default(self):
        """Energy function according to [#]_

        References
        ----------
        .. [#] Michel et al. (2014). PconsFold: improved contact predictions improve protein
           models. Bioinformatics 30(17), i482-i488

        """
        sstream = RosettaFunctionConstructs._ATOMPAIR
        sstream += "FADE -10 19 10 -15.00 0"
        return sstream

    @property
    def SIGMOID_default(self):
        """Simple sigmoidal energy function"""
        sstream = RosettaFunctionConstructs._ATOMPAIR
        sstream += "SIGMOID 8.00 1.00 #ContactMap: {raw_score: >4.3f}"
        return sstream

    @property
    def SIGMOID_gremlin(self):
        """Energy function according to [#]_

        References
        ----------
        .. [#] Ovchinnekov et al. (2015). Large-scale determination of previously unsolved
           protein structures using evolutionary information. Elife 4, e09248.

        """
        sstream = RosettaFunctionConstructs._ATOMPAIR
        sstream += RosettaFunctionConstructs._SCALARWEIGHTED
        sstream += "SUMFUNC 2 SIGMOID {sigmoid_cutoff: >6.3f} {sigmoid_slope: >6.3f} CONSTANTFUNC -0.5"
        return sstream
