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

from multiprocessing import Pool

from conkit.io import read
from conkit.misc.selectalg import SUBSELECTION_ALGORITHMS
from conkit.misc.selectalg import SubselectionAlgorithm


class StructureSelector(object):
    """Structure selection class for assessment by short-, medium- and long-range contact satisfaction"""

    def __init__(self, contactmap, nprocesses=1):
        """Instantiate a new :obj:`~conkit.misc.selector.StructureSelector` object

        Parameters
        ----------
        contactmap : :obj:`~conkit.core.contactmap.ContactMap`
           An instance of a :obj:`~conkit.core.contactmap.ContactMap`
        nprocesses : int, optional
           The number of processes

        """
        self.contactmap = contactmap
        self.nprocesses = nprocesses

    def assess(self, decoys, decoy_format, mode="linear"):
        """Subselect decoys excluding those not satisfying long-distance restraints

        Parameters
        ----------
        decoys : list, tuple
           A list containing paths to decoy files
        decoy_format : str
           The file format of `decoys`
        mode : str, optional
           The :obj:`~conkit.misc.selectalg.SubselectionAlgorithm` mode to use

        Returns
        -------
        list
           A :obj:`bool` list to categorize decoy as keeper

        Raises
        ------
        :exc:`ValueError`
           Unknown subselection mode

        """
        if mode not in SUBSELECTION_ALGORITHMS:
            raise ValueError("Unknown subselection mode: {}".format(mode))
        scores = self.compute_precision_by_range(decoys, decoy_format)
        _, _, longrange = zip(*scores)
        f = getattr(SubselectionAlgorithm, mode)
        keep, throw = f(longrange)
        keep = set(keep)
        return [True if i in keep else False for i in range(len(decoys))]

    def compute_precision_by_range(self, decoys, decoy_format):
        """Compute restraint precision score by sequence separation range

        Parameters
        ----------
        decoys : list, tuple
           A list containing paths to decoy files
        decoy_format : str
           The file format of ``decoys``

        Returns
        -------
        tuple
           A 2-D tuple containing short-range, medium-range and long-range scores for all decoys

        """
        formats = [decoy_format for _ in range(len(decoys))]
        cmap_copies = [self.contactmap.deepcopy() for _ in range(len(decoys))]
        args = zip(decoys, formats, cmap_copies)
        return Pool(self.nprocesses).map(_compute_single, args)


# This needs to be outside for the function to be pickleable by Pool
def _compute_single(args):
    decoy, decoy_format, cmap = args
    dmap = read(decoy, decoy_format).top_map
    matched = cmap.match(dmap)
    shortrange = matched.short_range_contacts
    mediumrange = matched.medium_range_contacts
    longrange = matched.long_range_contacts
    sprec, mprec, lprec = float("NaN"), float("NaN"), float("NaN")
    if shortrange.ncontacts > 0:
        sprec = shortrange.precision
    if mediumrange.ncontacts > 0:
        mprec = mediumrange.precision
    if longrange.ncontacts > 0:
        lprec = longrange.precision
    return sprec, mprec, lprec
