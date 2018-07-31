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
Parser module specific to ncont contact comparison
"""

__author__ = "Felix Simkovic"
__date__ = "16 Nov 2017"
__version__ = "0.1"

import re
import warnings

from conkit.io._parser import ContactFileParser
from conkit.core.contact import Contact
from conkit.core.contactmap import ContactMap
from conkit.core.contactfile import ContactFile

_re_hit_part = "\/\d+\/([A-Za-z])\/\s*(\d+)\(([A-Z]{3})\)\.\s+\/\s*[A-Z]+\s+\[\s*[A-Z]\]:\s+"
RE_CONTACT = re.compile(_re_hit_part + _re_hit_part + "(\d+\.\d+)\s*")


class NcontParser(ContactFileParser):
    """Class to parse a Ncont contact file
    """

    def __init__(self):
        super(NcontParser, self).__init__()

    def read(self, f_handle, f_id="ncont"):
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
        contact_map = ContactMap("map_1")
        contact_file.add(contact_map)

        for line in f_handle:
            line = line.strip()

            if RE_CONTACT.match(line):
                matches = RE_CONTACT.match(line)
                res1_seq = int(matches.group(2))
                res2_seq = int(matches.group(5))
                lb = ub = float(matches.group(7))

                if (res1_seq, res2_seq) in contact_map:
                    msg = "This parser cannot handle multiple atoms of the same residue. " \
                          "If your contact map contains such entries, only the first will be stored!"
                    warnings.warn(msg, Warning)
                    continue

                contact = Contact(res1_seq, res2_seq, 1.0, distance_bound=(lb, ub))
                contact.res1_chain = matches.group(1)
                contact.res2_chain = matches.group(4)
                contact.res1 = matches.group(3)
                contact.res2 = matches.group(6)
                contact_map.add(contact)

        contact_file.method = 'Contact map generated using Ncont'
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
        Creating a :meth:`~conkit.io.ncont.NcontParser.write` function 
        would come with a lot of issues, such as the parallel/antiparallel
        direction, scoring etc ... thus, this method is unavailable.

        Raises
        ------
        :exc:`NotImplementedError`
           Write function not available

        """
        raise NotImplementedError("Write function not available")
