"""
Parser module specific to bbcontacts predictions
"""

__author__ = "Felix Simkovic"
__date__ = "26 Oct 2016"
__version__ = 0.1

from conkit.core import Contact
from conkit.core import ContactMap
from conkit.core import ContactFile
from conkit.io._ParserIO import _ContactFileParser

import re

RE_COMMENT = re.compile(r'^#+.*$')


class BbcontactsParser(_ContactFileParser):
    """Class to parse a Bbcontacts contact file
    """
    def __init__(self):
        super(BbcontactsParser, self).__init__()

    def read(self, f_handle, f_id="bbcontacts"):
        """Read a contact file

        Parameters
        ----------
        f_handle
           Open file handle [read permissions]
        f_id : str, optional
           Unique contact file identifier

        Returns
        -------
        :obj:`ContactFile <conkit.core.ContactFile>`

        """

        contact_file = ContactFile(f_id)
        contact_map = ContactMap("map_1")
        contact_file.add(contact_map)

        for line in f_handle:
            line = line.rstrip()

            if not line:
                continue

            elif RE_COMMENT.match(line):
                continue

            else:
                # bbcontacts reverse residue numbering so swap
                _, _, _, raw_score, _, _, res2_seq, res1_seq = line.split()
                contact = Contact(int(res1_seq), int(res2_seq), float(raw_score))
                contact_map.add(contact)

        contact_file.method = 'Contact map predicted using Bbcontacts'

        return contact_file

    def write(self, f_handle, hierarchy):
        """Write a contact file instance to to file

        Parameters
        ----------
        f_handle
           Open file handle [write permissions]
        hierarchy : :obj:`ContactFile <conkit.core.ContactFile>`, :obj:`ContactMap <conkit.core.ContactMap>`
                    or :obj:`Contact <conkit.core.Contact>`

        Notes
        -----
        Creating a :func`write` function for the Bbcontacts parser
        would come with a lot of issues, such as the parallel/antiparallel
        direction, scoring etc.

        Thus, no :func:`write` method is available.

        Raises
        ------
        RuntimeError
           Not available

        """
        raise RuntimeError('Not available')
