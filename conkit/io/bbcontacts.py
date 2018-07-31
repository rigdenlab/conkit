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
Parser module specific to bbcontacts predictions
"""

__author__ = "Felix Simkovic"
__date__ = "23 Jul 2018"
__version__ = "0.2"

import re

from conkit.io._parser import ContactFileParser
from conkit.core.contact import Contact
from conkit.core.contactmap import ContactMap
from conkit.core.contactfile import ContactFile


class BbcontactsParser(ContactFileParser):
    """Class to parse a Bbcontacts contact file
    """

    def __init__(self):
        super(BbcontactsParser, self).__init__()

    def read(self, f_handle, f_id="bbcontacts", del_one_two=False):
        """Read a contact file

        Parameters
        ----------
        f_handle
           Open file handle [read permissions]
        f_id : str, optional
           Unique contact file identifier
        del_one_two : bool
           Remove one- & two-strand sheets

        Returns
        -------
        :obj:`~conkit.core.contactfile.ContactFile`

        """

        contact_file = ContactFile(f_id)
        contact_map = ContactMap("map_1")
        contact_file.add(contact_map)

        previous = 'first'
        for line in f_handle:
            line = line.strip()

            if line and not line.startswith('#'):
                _, _, _, raw_score, _, current, res2_seq, res1_seq = line.split()
                if del_one_two and previous == 'first' and current == 'last':
                    contact_map.child_list.pop()
                elif any(value == "NA" for value in [raw_score, res2_seq, res1_seq]):
                    pass
                else:
                    contact = Contact(int(res1_seq), int(res2_seq), float(raw_score))
                    contact_map.add(contact)
                previous = current

        if del_one_two and previous == 'first' and len(contact_map) > 0:
            contact_map.child_list.pop()

        contact_file.method = 'Contact map predicted using Bbcontacts'

        return contact_file

    def write(self, f_handle, hierarchy):
        """Write a contact file instance to to file

        Parameters
        ----------
        f_handle
           Open file handle [write permissions]
        hierarchy : :obj:`~conkit.core.contactfile.ContactFile`, :obj:`~conkit.core.contactmap.ContactMap`
                    or :obj:`~conkit.core.contact.Contact`

        Note
        ----
        Creating a :meth:`~conkit.io.bbcontacts.BbcontactsParser.write` method
        would come with a lot of issues, such as the parallel/antiparallel
        direction, scoring etc ... thus, this function is unvailable.

        Raises
        ------
        :exc:`NotImplementedError`
           Write function not available

        """
        raise NotImplementedError("Write function not available")
