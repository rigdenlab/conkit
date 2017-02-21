"""
Storage space for a contact file
"""

__author__ = "Felix Simkovic"
__date__ = "03 Aug 2016"
__version__ = 0.1

from conkit.core.Entity import Entity


class ContactFile(Entity):
    """A contact file object representing a single prediction file

    The contact file class represents a data structure to hold all predictions
    with a single contact map file. It contains functions to store,
    manipulate and organise contact maps.

    Attributes
    ----------
    author : str
       The author of the :obj:`ContactFile <conkit.core.ContactFile>`
    method : list, str
       The :obj:`ContactFile <conkit.core.ContactFile>`-specific method
    remark : list, str
       The :obj:`ContactFile <conkit.core.ContactFile>`-specific remarks
    target : str
       The target name
    top_map : :obj:`ContactMap <conkit.core.ContactMap>`
       The first :obj:`ContactMap <conkit.core.ContactMap>` entry in :obj:`ContactFile <conkit.core.ContactFile>`

    Examples
    --------
    >>> from conkit.core import ContactMap, ContactFile
    >>> contact_file = ContactFile("example")
    >>> contact_file.add(ContactMap("foo"))
    >>> contact_file.add(ContactMap("bar"))
    >>> print(contact_file)
    ContactFile(id="example" nseqs=2)

    """
    __slots__ = ['_author', '_method', '_remark', '_target']

    def __init__(self, id):
        """Initialise a new contact map

        Parameters
        ----------
        id : str
           A unique identifier for the contact file

        """
        self._author = None
        self._method = []
        self._remark = []
        self._target = None
        super(ContactFile, self).__init__(id)

    def __repr__(self):
        return "ContactFile(id=\"{0}\" nmaps={1})".format(self.id, len(self))

    @property
    def author(self):
        """The author of the :obj:`ContactFile <conkit.core.ContactFile>`"""
        return self._author

    @author.setter
    def author(self, author):
        """Define the author of the :obj:`ContactFile <conkit.core.ContactFile>`

        Parameters
        ----------
        author : str

        """
        self._author = author

    @property
    def method(self):
        """The :obj:`ContactFile <conkit.core.ContactFile>`-specific method"""
        return self._method

    @method.setter
    def method(self, method):
        """Set the :obj:`ContactFile <conkit.core.ContactFile>` method

        Parameters
        ----------
        method : str, list
           The method will be added to the list of methods

        """
        if isinstance(method, list):
            self._method += method
        elif isinstance(method, tuple):
            self._method += list(method)
        else:
            self._method += [method]

    @property
    def remark(self):
        """The :obj:`ContactFile <conkit.core.ContactFile>`-specific remarks"""
        return self._remark

    @remark.setter
    def remark(self, remark):
        """Set the :obj:`ContactFile <conkit.core.ContactFile>` remark

        Parameters
        ----------
        remark : str, list
           The remark will be added to the list of remarks

        """
        if isinstance(remark, list):
            self._remark += remark
        elif isinstance(remark, tuple):
            self._remark += list(remark)
        else:
            self._remark += [remark]

    @property
    def target(self):
        """The target name"""
        return self._target

    @target.setter
    def target(self, target):
        """Define the target name

        Parameters
        ----------
        target : str

        """
        self._target = target

    @property
    def top_map(self):
        """The first :obj:`ContactMap <conkit.core.ContactMap>` entry in :obj:`ContactFile <conkit.core.ContactFile>`

        Returns
        -------
        top_map : :obj:`ContactMap <conkit.core.ContactMap>`, None
           The first :obj:`ContactMap <conkit.core.ContactMap>` entry in :obj:`ContactFile <conkit.core.ContactFile>`

        """
        if len(self) > 0:
            return self[0]
        else:
            return None

    def sort(self, kword, reverse=False, inplace=False):
        """Sort the :obj:`ContactFile <conkit.core.ContactFile>`

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
        contact_map : :obj:`ContactMap <conkit.core.ContactMap>`
           The reference to the :obj:`ContactMap <conkit.core.ContactMap>`, regardless of inplace

        Raises
        ------
        ValueError
           ``kword`` not in :obj:`ContactFile <conkit.core.ContactFile>`

        """
        contact_file = self._inplace(inplace)
        contact_file._sort(kword, reverse)
        return contact_file
