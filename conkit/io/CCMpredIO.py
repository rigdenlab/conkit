"""
Parser module specific to CCMpred predictions
"""

from __future__ import division

__author__ = "Felix Simkovic"
__date__ = "03 Aug 2016"
__version__ = 0.1

from conkit.core import Contact
from conkit.core import ContactMap
from conkit.core import ContactFile
from conkit.io._ParserIO import _ContactFileParser

import numpy
import sys


class CCMpredParser(_ContactFileParser):
    """
    Class to parse a CCMpred contact matrix

    """
    def __init__(self):
        super(CCMpredParser, self).__init__()
    
    def read(self, f_handle, f_id="ccmpred"):
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

        # Bits ripped from Stefan Seemayer's script shipped with CCMpred
        mat = numpy.loadtxt(f_handle)
        raw_contacts = self._get_contact_pairs(mat)

        for res1_seq, res2_seq, raw_score in zip(raw_contacts[0], raw_contacts[1], mat[raw_contacts]):
            if res1_seq > res2_seq:
                continue
            # Matrix starts count at 0 so increment numbers by one straight away
            contact = Contact(
                int(res1_seq+1),
                int(res2_seq+1),
                float(raw_score)
            )
            contact_map.add(contact)
        contact_file.method = 'Contact map predicted using CCMpred'

        return contact_file
        
    def _get_contact_pairs(self, mat):
        """Get all contact pairs in the matrix

        Parameters
        ----------
        mat : numpy.ndarray
           A numpy arranged matrix

        Returns
        -------
        contacts : list
           A list of contact pairs

        """
        contacts = mat.argsort(axis=None)[::-1]
        contacts = (contacts % mat.shape[0]).astype(numpy.uint16), \
                    numpy.floor(contacts / mat.shape[0]).astype(numpy.uint16)
        return contacts

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
        TypeError
           Python3 requires f_handle to be in 'wb' or 'ab' mode

        """
        # Python3 support requires bytes mode
        if sys.version_info.major == 3 and not (f_handle.mode == 'wb' or f_handle.mode == 'ab'):
            raise TypeError("Python3 requires f_handle to be in 'wb' or 'ab' mode")

        # Double check the type of hierarchy and reconstruct if necessary
        contact_file = self._reconstruct(hierarchy)

        if len(contact_file) > 1:
            raise RuntimeError('More than one contact map provided')

        for contact_map in contact_file:
            len_mat = max([c.res1_seq for c in contact_map] + [c.res2_seq for c in contact_map])
            mat = numpy.zeros((len_mat, len_mat), numpy.float64)

            for contact in contact_map:
                mat[contact.res1_seq - 1][contact.res2_seq - 1] = contact.raw_score
                mat[contact.res2_seq - 1][contact.res1_seq - 1] = contact.raw_score
            
            numpy.savetxt(f_handle, mat)

        return
