# coding=utf-8
#
# BSD 3-Clause License
#
# Copyright (c) 2016-21, University of Liverpool
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
"""ConfidenceFile container used throughout ConKit"""

from __future__ import division
from __future__ import print_function

__author__ = "Aderik Voorspoels"
__date__ = "07 Apr 2025"
__version__ = "0.13.3"

from conkit.core.entity import Entity


class ConfidenceFile(Entity):
    """A confidence file object representing seperate confidences of a single prediction file

    The confidence file class represents a data structure to hold confidences associated to predictions. 
    It contains functions to store,
    manipulate and organise confidence functions.

    Examples
    --------
    >>> from conkit.core import ConfidenceMap, ConfidenceFile
    >>> confidence_file = ConfidenceFile("example")
    >>> confidence_file.add(ConfidenceMap("foo"))
    >>> confidence_file.add(ConfidenceMap("bar"))
    >>> print(confidence_file)
    ConfidenceFile(id="example" nmaps=2)

    Attributes
    ----------
    author : str
       The author of the :obj:`~conkit.core.confidencefile.ConfidenceFile`
    method : list, str
       The :obj:`~conkit.core.confidencefile.ConfidenceFile`-specific method
    remark : list, str
       The :obj:`~conkit.core.confidencefile.ConfidenceFile`-specific remarks
    target : str
       The target name
    top_map : :obj:`~conkit.core.confidencemap.ConfidenceMap`
       The first :obj:`~conkit.core.confidencemap.ConfidenceMap` entry in :obj:`~conkit.core.confidencefile.ConfidenceFile`

    """

    __slots__ = ["author", "target", "_method", "_remark"]
    
    def __init__(self, id):
        """Initialise a new confidence map

        Parameters
        ----------
        id : str
           A unique identifier for this :obj:`~conkit.core.confidencefile.ConfidenceFile`

        """
        self.author = None
        self.target = None
        self._method = []
        self._remark = []
        super(ConfidenceFile, self).__init__(id)

    def __repr__(self):
        return '{}(id="{}" nmaps={})'.format(self.__class__.__name__, self.id, len(self))

    @property
    def method(self):
        """The :obj:`~conkit.core.confidencefile.ConfidenceFile`-specific method"""
        return self._method

    @method.setter
    def method(self, method):
        """Set the :obj:`~conkit.core.confidencefile.ConfidenceFile` method

        Parameters
        ----------
        method : str, list
           The method will be added to the list of methods

        """
        self._method += Entity.listify(method)

    @property
    def remark(self):
        """The :obj:`~conkit.core.confidencefile.ConfidenceFile`-specific remarks"""
        return self._remark

    @remark.setter
    def remark(self, remark):
        """Set the :obj:`~conkit.core.confidencefile.ConfidenceFile` remark

        Parameters
        ----------
        remark : str, list
           The remark will be added to the list of remarks

        """
        self._remark += Entity.listify(remark)

    @property
    def top_map(self):
        """The first :obj:`~conkit.core.confidencemap.ConfidenceMap` entry"""
        return self.top

    def sort(self, kword, reverse=False, inplace=False):
        """Sort the :obj:`~conkit.core.confidencefile.ConfidenceFile`

        Parameters
        ----------
        kword : str
           The dictionary key to sort confidences by
        reverse : bool, optional
           Sort the confidence pairs in descending order [default: False]
        inplace : bool, optional
           Replace the saved order of confidences [default: False]

        Returns
        -------
        :obj:`~conkit.core.confidencemap.ConfidenceMap`
           The reference to the :obj:`~conkit.core.confidencemap.ConfidenceMap`, regardless of `inplace`

        Raises
        ------
        :exc:`ValueError`
           `kword` not in :obj:`~conkit.core.confidencefile.ConfidenceFile`

        """
        confidence_file = self._inplace(inplace)
        confidence_file._sort(kword, reverse)
        return confidence_file
