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
"""Distance prediction file container used throughout ConKit"""

from conkit.core.contactfile import ContactFile


class DistanceFile(ContactFile):
    """A distance file object representing a single distance prediction file. This class inherits methods and attributes
    from :obj:`~conkit.core.contactfile.ContactFile`


    Examples
    --------
    >>> from conkit.core import Distogram, DistanceFile
    >>> distance_file = DistanceFile("example")
    >>> distance_file.add(Distogram("foo"))
    >>> distance_file.add(Distogram("bar"))
    >>> print(distance_file)
    DistanceFile(id="example" ndistograms=2)

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
    original_file_format : str
       The original file format used to create the :obj:`~conkit.core.distogram.Distogram` instance

    """

    __slots__ = ["author", "target", "_method", "_remark", "_original_file_format"]

    def __init__(self, id):
        """Initialise a new distance file"""
        self._original_file_format = None
        super(DistanceFile, self).__init__(id)

    def __repr__(self):
        return '{}(id="{}" ndistograms={})'.format(self.__class__.__name__, self.id, len(self))

    @property
    def original_file_format(self):
        """The original file format used to create the :obj:`~conkit.core.distancefile.DistanceFile` instance"""
        return self._original_file_format

    @original_file_format.setter
    def original_file_format(self, value):
        self._original_file_format = value

    def add(self, entity):
        entity.original_file_format = self.original_file_format
        super(ContactFile, self).add(entity)
