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
"""
Parser module specific to CCMpred predictions
"""

from __future__ import division

__author__ = "Felix Simkovic"
__date__ = "03 Aug 2016"
__version__ = "0.1"

import numpy as np
import sys

from conkit.io._parser import ContactFileParser
from conkit.core.contact import Contact
from conkit.core.contactmap import ContactMap
from conkit.core.contactfile import ContactFile


class CCMpredParser(ContactFileParser):
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
        :obj:`~conkit.core.contactfile.ContactFile`

        """
        contact_file = ContactFile(f_id)
        contact_file.method = 'Contact map predicted using CCMpred'
        contact_map = ContactMap("map_1")
        contact_file.add(contact_map)

        # Bits ripped from Stefan Seemayer's script shipped with CCMpred
        mat = np.loadtxt(f_handle)
        if mat.size > 0:
            raw_contacts = self._get_contact_pairs(mat)
            for res1_seq, res2_seq, raw_score in zip(raw_contacts[0], raw_contacts[1], mat[raw_contacts]):
                if res1_seq > res2_seq:
                    continue
                # Matrix starts count at 0 so increment numbers by one straight away
                contact = Contact(int(res1_seq + 1), int(res2_seq + 1), float(raw_score))
                contact_map.add(contact)

        return contact_file

    def _get_contact_pairs(self, mat):
        """Get all contact pairs in the matrix

        Parameters
        ----------
        mat : :obj:`~numpy.ndarray`
           A :mod:`numpy` matrix

        Returns
        -------
        list
           A list of contact pairs

        """
        contacts = mat.argsort(axis=None)[::-1]
        contacts = (contacts % mat.shape[0]).astype(np.uint16), \
                    np.floor(contacts / mat.shape[0]).astype(np.uint16)
        return contacts

    def write(self, f_handle, hierarchy):
        """Write a contact file instance to to file

        Parameters
        ----------
        f_handle
           Open file handle [write permissions]
        hierarchy : :obj:`~conkit.core.contactfile.ContactFile`, :obj:`~conkit.core.contactmap.ContactMap`
                    or :obj:`~conkit.core.contact.Contact`

        Raises
        ------
        :exc:`RuntimeError`
           More than one contact map in the hierarchy
        :exc:`TypeError`
           Python3 requires f_handle to be in `wb` or `ab` mode

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
            mat = np.zeros((len_mat, len_mat), np.float64)

            for contact in contact_map:
                mat[contact.res1_seq - 1][contact.res2_seq - 1] = contact.raw_score
                mat[contact.res2_seq - 1][contact.res1_seq - 1] = contact.raw_score

            np.savetxt(f_handle, mat, delimiter="\t")

        return
