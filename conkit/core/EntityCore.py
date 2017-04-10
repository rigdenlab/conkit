"""
Base class for all entities used in this interface.
"""

__author__ = "Felix Simkovic"
__date__ = "17 Aug 2016"
__version__ = 0.1

import copy
import operator


class Entity(object):
    """Base class for all entities used in this interface.

    It handles the storage of data. It also provides a high-efficiency
    methods to allow fast lookup and iterations of each entity. It also
    provides a hierarchical structure to remember parent and child
    entities.

    Attributes
    ----------
    id : str, list, tuple
       The ID of the selected entity
    full_id : tuple
       A traceback id including all parent classes
    parent : :obj:`Entity <conkit.core.Entity>`
       An attribute to store the reference to the parent :obj:`Entity <conkit.core.Entity>`
    child_list : list
       A list storing the child entities
    child_dict : dict
       A dictionary storing the child entities

    Notes
    -----
    It is strongly advised against the use of the :obj:`Entity <conkit.core.Entity>` class directly.
    Instead, use one or more of the the remaining data models.

    """
    __slots__ = ['_id', '_parent', '_child_list', '_child_dict']

    def __init__(self, id):
        """Initialise a generic :obj:`Entity <conkit.core.Entity>`

        Parameters
        ----------
        id : str, list, tuple
           The ID of the selected entity

        """
        self._id = None
        self._parent = None
        self._child_list = []
        self._child_dict = {}

        # Assign values post creation to use setter/getter methods
        # Possibly very bad practice but no better alternative for now
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
            indexes_to_keep = list(range(*id.indices(len(self))))
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
    def child_dict(self):
        """A dictionary storing the child entities"""
        return self._child_dict

    @child_dict.setter
    def child_dict(self, child_dict):
        """Define a dictionary storing the child entities

        Parameters
        ----------
        child_dict : dict

        """
        self._child_dict = child_dict

    @property
    def child_list(self):
        """A list storing the child entities"""
        return self._child_list

    @child_list.setter
    def child_list(self, child_list):
        """Define a list storing the child entities

        Parameters
        ----------
        child_list : dict

        """
        self._child_list = child_list

    @property
    def full_id(self):
        """A traceback id including all parent classes

        The full id is a tuple containing all id's starting from
        the top object (:obj:`ContactFile <conkit.core.ContactFile>`) down to the current object.
        A full id for a :obj:`Contact <conkit.core.Contact>` e.g. is something like:
        ('1aa', 1, (1, 10))

        This corresponds to:

        :obj:`ContactFile <conkit.core.ContactFile>` identifier => 1aaa
        :obj:`ContactMap <conkit.core.ContactMap>` identifier => 1
        :obj:`Contact <conkit.core.Contact>` identifier => (1, 10)

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
           The unique ID for an :obj:`Entity <conkit.core.Entity>`

        Warnings
        --------
        You cannot provide an :obj:`int` or :obj:`float` as ID.

        """
        if isinstance(id, int) or isinstance(id, float):
            raise TypeError('Please provide data type of str, list, or tuple')
        elif isinstance(id, list):
            id = tuple(id)
        self._id = id

    @property
    def parent(self):
        """An attribute to store the reference to the parent :obj:`Entity <conkit.core.Entity>`"""
        return self._parent

    @parent.setter
    def parent(self, parent):
        """Define the reference to the parent :obj:`Entity <conkit.core.Entity>`

        Parameters
        ----------
        parent : :obj:`Entity <conkit.core.Entity>`

        """
        self._parent = parent

    def _inplace(self, inplace):
        """Modify the current version using a copy

        Parameters
        ----------
        inplace : bool

        """
        if inplace:
            return self
        else:
            return self.copy()

    def _sort(self, kword, reverse):
        """Sort the :obj:`Entity <conkit.core.Entity>`"""
        if any(not hasattr(e, kword) for e in self._child_list):
            raise ValueError('Attribute not defined')
        self.child_list.sort(key=operator.attrgetter(kword), reverse=reverse)

    def add(self, entity):
        """Add a child to the :obj:`Entity <conkit.core.Entity>`

        Parameters
        ----------
        entity : :obj:`Entity <conkit.core.Entity>`

        """
        if entity.id in self:
            raise ValueError("%s defined twice" % str(entity.id))
        entity.parent = self
        self.child_list.append(entity)
        self.child_dict[entity.id] = entity

    def copy(self):
        """Create a shallow copy of :obj:`Entity <conkit.core.Entity>`"""
        shallow = copy.copy(self)

        shallow.child_list = []
        shallow.child_dict = {}
        shallow.parent = None

        for child in self:
            shallow.add(child.copy())
        return shallow

    def deepcopy(self):
        """Create a deep copy of :obj:`Entity <conkit.core.Entity>`"""
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

        Warnings
        --------
        If ``id`` is of type :obj:`int`, then the :obj:`Entity <conkit.core.Entity>`
        in the ``child_list`` at index ``id`` will be deleted

        """
        del self[id]
