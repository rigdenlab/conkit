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
"""ContactFile container used throughout ConKit"""

from __future__ import division
from __future__ import print_function

__author__ = "Felix Simkovic"
__date__ = "03 Aug 2016"
__version__ = "1.0"

from conkit.core.entity import Entity


class ContactFile(Entity):
    """A contact file object representing a single prediction file

    The contact file class represents a data structure to hold all predictions
    with a single contact map file. It contains functions to store,
    manipulate and organise contact maps.

    Examples
    --------
    >>> from conkit.core import ContactMap, ContactFile
    >>> contact_file = ContactFile("example")
    >>> contact_file.add(ContactMap("foo"))
    >>> contact_file.add(ContactMap("bar"))
    >>> print(contact_file)
    ContactFile(id="example" nseq=2)

    Attributes
    ----------
    author : str
       The author of the :obj:`~conkit.core.contactfile.ContactFile`
    method : list, str
       The :obj:`~conkit.core.contactfile.ContactFile`-specific method
    remark : list, str
       The :obj:`~conkit.core.contactfile.ContactFile`-specific remarks
    target : str
       The target name
    top_map : :obj:`~conkit.core.contactmap.ContactMap`
       The first :obj:`~conkit.core.contactmap.ContactMap` entry in :obj:`~conkit.core.contactfile.ContactFile`

    """
    __slots__ = ['author', 'target', '_method', '_remark']

    def __init__(self, id):
        """Initialise a new contact map

        Parameters
        ----------
        id : str
           A unique identifier for this :obj:`~conkit.core.contactfile.ContactFile`

        """
        self.author = None
        self.target = None
        self._method = []
        self._remark = []
        super(ContactFile, self).__init__(id)

    def __repr__(self):
        return '{}(id="{}" nmaps={})'.format(self.__class__.__name__, self.id, len(self))

    @property
    def method(self):
        """The :obj:`~conkit.core.contactfile.ContactFile`-specific method"""
        return self._method

    @method.setter
    def method(self, method):
        """Set the :obj:`~conkit.core.contactfile.ContactFile` method

        Parameters
        ----------
        method : str, list
           The method will be added to the list of methods

        """
        self._method += Entity.listify(method)

    @property
    def remark(self):
        """The :obj:`~conkit.core.contactfile.ContactFile`-specific remarks"""
        return self._remark

    @remark.setter
    def remark(self, remark):
        """Set the :obj:`~conkit.core.contactfile.ContactFile` remark

        Parameters
        ----------
        remark : str, list
           The remark will be added to the list of remarks

        """
        self._remark += Entity.listify(remark)

    @property
    def top_map(self):
        """The first :obj:`~conkit.core.contactmap.ContactMap` entry"""
        return self.top

    def sort(self, kword, reverse=False, inplace=False):
        """Sort the :obj:`~conkit.core.contactfile.ContactFile`

        Parameters
        ----------
        kword : str
           The dictionary key to sort contacts by
        reverse : bool, optional
           Sort the contact pairs in descending order [default: False]
        inplace : bool, optional
           Replace the saved order of contacts [default: False]

        Returns
        -------
        :obj:`~conkit.core.contactmap.ContactMap` 
           The reference to the :obj:`~conkit.core.contactmap.ContactMap`, regardless of `inplace`

        Raises
        ------
        :exc:`ValueError`
           `kword` not in :obj:`~conkit.core.contactfile.ContactFile`

        """
        contact_file = self._inplace(inplace)
        contact_file._sort(kword, reverse)
        return contact_file
