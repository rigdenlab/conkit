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

import inspect
import numpy as np


class SubselectionAlgorithm(object):
    """A class to collect all subselection algorithms"""

    @classmethod
    def _numpify(cls, data):
        """Convert a Python array to a :obj:`~numpy.ndarray`"""
        return np.asarray(data)

    @classmethod
    def cutoff(cls, data, cutoff=0.287):
        """A cutoff-defined subselection algorithm

        This algorithm removes a decoy, if its score is less
        than the cutoff.

        Parameters
        ----------
        data : list, tuple
           A 1D array of scores
        cutoff : float, optional
           The cutoff of keeping decoys

        Returns
        -------
        list
           The decoy indices to keep
        list
           The decoy indices to throw

        """
        data = cls._numpify(data)
        keep = np.where(data >= cutoff)[0]
        throw = np.where(data < cutoff)[0]
        return keep.tolist(), throw.tolist()

    @classmethod
    def linear(cls, data, cutoff=0.5):
        """A linearly-defined subselection algorithm

        This algorithm removes the worst ``x``% decoys, where ``x``
        is defined by ``cutoff``.

        Parameters
        ----------
        data : list, tuple
           A 1D array of scores
        cutoff : float, optional
           The porportion of the total number of decoys to keep

        Returns
        -------
        list
           The decoy indices to keep
        list
           The decoy indices to throw

        """
        sorted_indices = cls._numpify(data).argsort()[::-1]
        pivot = np.ceil(sorted_indices.shape[0] * cutoff).astype(np.int)
        keep = sorted_indices[:pivot]
        throw = sorted_indices[pivot:]
        return keep.tolist(), throw.tolist()

    @classmethod
    def scaled(cls, data, cutoff=0.5):
        """A scaling-defined subselection algorithm

        This algorithm removes a decoy, if its scaled score
        is less than 0.5. The scaled score is calculated by
        dividing the precision score by the average of the
        set.

        Parameters
        ----------
        data : list, tuple
           A 1D array of scores
        cutoff : float, optional
           The cutoff of keeping decoys

        Returns
        -------
        list
           The decoy indices to keep
        list
           The decoy indices to throw

        """
        data = cls._numpify(data)
        data_scaled = data / np.mean(data)
        return SubselectionAlgorithm.cutoff(data_scaled, cutoff=cutoff)

    @classmethod
    def ignore(cls, data):
        """"A subselection algorithm to keep all

        This algorithm doesn't do anything except mimic others.

        It will not discard any decoys and keep all!!

        Parameters
        ----------
        data : list, tuple
           A 1D array of scores

        Returns
        -------
        list
           The decoy indices to keep
        list
           The decoy indices to throw

        """
        data = cls._numpify(data)
        return cls.cutoff(data, cutoff=0)


SUBSELECTION_ALGORITHMS = [
    func_name for func_name, _ in inspect.getmembers(SubselectionAlgorithm) if not func_name.startswith("_")
]
