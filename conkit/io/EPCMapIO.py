"""
Parser module specific to EPC-Map predictions
"""

__author__ = "Felix Simkovic"
__date__ = "12 Dec 2016"
__version__ = "0.1"

from conkit.core import Contact
from conkit.core import ContactMap
from conkit.core import ContactFile
from conkit.io._ParserIO import _ContactFileParser

import os


class EPCMapParser(_ContactFileParser):
    """Class to parse a EPC-Map contact prediction
    """
    def read(self, f_handle, f_id="epcmap"):
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

        hierarchy = ContactFile(f_id)
        _map = ContactMap("map_1")
        hierarchy.add(_map)

        for line in f_handle:
            line = line.strip().split()

            if not line or line[0].isalpha():
                continue

            elif line[0].isdigit():
                _contact = Contact(int(line[0]), int(line[1]), float(line[4]),
                                   distance_bound=(int(line[2]), int(line[3])))
                _map.add(_contact)

        hierarchy.method = 'Contact map predicted using EPC-Map'

        return hierarchy

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
                line = "{res1_seq} {res2_seq} {lb} {ub} {raw_score:.6f}"
                line = line.format(res1_seq=contact.res1_seq, res2_seq=contact.res2_seq, raw_score=contact.raw_score,
                                   lb=contact.distance_bound[0], ub=contact.distance_bound[1])
                f_handle.write(line + os.linesep)

        return
