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
"""Entity container used throughout ConKit"""

from __future__ import division
from __future__ import print_function

__author__ = "Felix Simkovic"
__date__ = "03 Aug 2016"
__version__ = "1.0"

import copy
import operator


class Entity(object):
    """Base class for all entities used in this interface.

    It handles the storage of data. It also provides a high-efficiency
    methods to allow fast lookup and iterations of each entity. It also
    provides a hierarchical structure to remember parent and child
    entities.

    Warning 
    -------
    It is strongly advised against the use of the :obj:`~conkit.core.entity.Entity` class directly.
    Instead, use one or more of the the remaining data models.

    Attributes
    ----------
    id : str, list, tuple
       The ID of the selected entity
    full_id : tuple
       A traceback id including all parent classes
    parent : :obj:`~conkit.core.entity.Entity`
       An attribute to store the reference to the parent :obj:`~conkit.core.entity.Entity`
    child_list : list
       A list storing the child entities
    child_dict : dict
       A dictionary storing the child entities

    """
    __slots__ = ['parent', '_id', 'child_list', 'child_dict']

    def __init__(self, id):
        """Initialise a generic :obj:`~conkit.core.entity.Entity`

        Parameters
        ----------
        id : str, list, tuple
           The ID of the selected entity

        """
        self._id = None
        self.parent = None
        self.child_list = []
        self.child_dict = {}
        self.id = id

    def __contains__(self, id):
        """True if there is a child element with the given id"""
        return id in self.child_dict

    def __delitem__(self, id):
        """Remove a child with given id"""
        child = self[id]
        child.parent = None
        self.child_dict.pop(id)
        self.child_list.remove(child)

    def __getitem__(self, id):
        """Return the child with the given id"""
        if isinstance(id, slice):
            indexes_to_keep = set(range(*id.indices(len(self))))
            copy_to_return = self.copy()
            for i, child in enumerate(self):
                if i not in indexes_to_keep:
                    copy_to_return.remove(child.id)
            return copy_to_return
        elif isinstance(id, int):
            return self.child_list[id]
        else:
            return self.child_dict[id]

    def __iter__(self):
        """Iterate over children"""
        for child in self.child_list:
            yield child

    def __len__(self):
        """Return the number of children"""
        return len(self.child_list)

    def __reversed__(self):
        """Reversed list of the children"""
        for child in reversed(self.child_list):
            yield child

    @property
    def full_id(self):
        """A traceback id including all parent classes

        The full id is a tuple containing all id's starting from
        the top object (:obj:`~conkit.core.contactfile.ContactFile`) down to the current object.
        A full id for a :obj:`~conkit.core.contact.Contact` e.g. is something like:
        ('1aa', 1, (1, 10))

        This corresponds to:

        :obj:`~conkit.core.contactfile.ContactFile` identifier => 1aaa
        :obj:`~conkit.core.contactmap.ContactMap` identifier => 1
        :obj:`~conkit.core.contact.Contact` identifier => (1, 10)

        """
        traceback = [self.id]
        mother = self.parent
        while mother is not None:
            traceback.append(mother.id)
            mother = mother.parent
        return tuple(reversed(traceback))

    @property
    def id(self):
        """The ID of the selected entity"""
        return self._id

    @id.setter
    def id(self, id):
        """Set the ID of the selected entity

        Parameters
        ----------
        id : str, list, tuple
           The unique ID for an :obj:`~conkit.core.entity.Entity`

        Warning
        -------
        You cannot provide an :obj:`int` or :obj:`float` as ID.

        Raises
        ------
        :obj:`TypeError`
           Please provide data type of str, list, or tuple

        """
        if isinstance(id, int) or isinstance(id, float):
            raise TypeError('Please provide data type of str, list, or tuple')
        elif isinstance(id, list):
            id = tuple(id)
        self._id = id

    @property
    def top(self):
        """The first child in the :obj:`~conkit.core.entity.Entity`"""
        if len(self) > 0:
            return self.child_list[0]
        else:
            return None

    def _inplace(self, inplace):
        """Modify the current version using a copy

        Parameters
        ----------
        inplace : bool

        """
        if inplace:
            return self
        else:
            return self.deepcopy()

    def _sort(self, kword, reverse):
        """Sort the :obj:`~conkit.core.entity.Entity`"""
        if any(not hasattr(e, kword) for e in self.child_list):
            raise ValueError('Attribute not defined')
        self.child_list.sort(key=operator.attrgetter(kword), reverse=reverse)

    def add(self, entity):
        """Add a child to the :obj:`~conkit.core.entity.Entity`

        Parameters
        ----------
        entity : :obj:`~conkit.core.entity.Entity`

        """
        if entity.id in self:
            raise ValueError("%s defined twice" % str(entity.id))
        entity.parent = self
        self.child_list.append(entity)
        self.child_dict[entity.id] = entity

    def copy(self):
        """Create a shallow copy of :obj:`~conkit.core.entity.Entity`"""
        shallow = copy.copy(self)

        shallow.child_list = []
        shallow.child_dict = {}
        shallow.parent = None

        for child in self:
            shallow.add(child.copy())
        return shallow

    def deepcopy(self):
        """Create a deep copy of :obj:`~conkit.core.entity.Entity`"""
        deep = copy.deepcopy(self)

        deep.child_list = []
        deep.child_dict = {}
        deep.parent = None

        for child in self:
            deep.add(child.copy())
        return deep

    def remove(self, id):
        """Remove a child

        Parameters
        ----------
        id : str, int, list, tuple

        Warning
        -------
        If `id` is of type :obj:`int`, then the :obj:`~conkit.core.entity.Entity`
        in the :attr:`~conkit.core.entity.Entity.child_list`` at index `id` will 
        be deleted

        """
        del self[id]

    @staticmethod
    def listify(s):
        """Convert unknown input to a list

        Parameters
        ----------
        s : str, int, float, list, tuple
        
        Returns
        -------
        list
           The input as list

        """
        if isinstance(s, list):
            return s
        elif isinstance(s, tuple):
            return list(s)
        else:
            return [s]
