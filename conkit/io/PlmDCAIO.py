"""
Parser module specific to plmDCA predictions
"""

__author__ = "Felix Simkovic"
__date__ = "03 Aug 2016"
__version__ = 0.1

from conkit.core import Contact
from conkit.core import ContactMap
from conkit.core import ContactFile
from conkit.io._ParserIO import _ContactFileParser

import os

class PlmDCAParser(_ContactFileParser):
    """Class to parse a plmDCA contact prediction
    """
    def __init__(self):
        super(PlmDCAParser, self).__init__()

    def read(self, f_handle, f_id="plmdca"):
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
            line = line.strip()

            if not line or line[0].isalpha():
                continue

            elif line[0].isdigit():
                res1_seq, res2_seq, raw_score = line.split(',')
                contact = Contact(int(res1_seq), int(res2_seq), float(raw_score))
                contact_map.add(contact)

        contact_file.method = 'Contact map predicted using plmDCA'

        return contact_file

    def write(self, f_handle, hierarchy):
        """Write a contact file instance to to file

        Parameters
        ----------
        f_handle
           Open file handle [write permissions]
        hierarchy : :obj:`ContactFile <conkit.core.ContactFile>`, :obj:`ContactMap <conkit.core.ContactMap>`
                    or :obj:`Contact <conkit.core.Contact>`

        Raises
        ------
        RuntimeError
           More than one contact map in the hierarchy

        """
        # Double check the type of hierarchy and reconstruct if necessary
        contact_file = self._reconstruct(hierarchy)

        if len(contact_file) > 1:
            raise RuntimeError('More than one contact map provided')

        for contact_map in contact_file:
            for contact in contact_map:
                line = "{res1_seq},{res2_seq},{raw_score:.6f}"
                line = line.format(res1_seq=contact.res1_seq, res2_seq=contact.res2_seq, raw_score=contact.raw_score)
                f_handle.write(line + os.linesep)

        return
