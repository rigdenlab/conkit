# BSD 3-Clause License
#
# Copyright (c) 2016-21, University of Liverpool
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
Parser module specific to COMSAT predictions
"""

__author__ = "Felix Simkovic"
__date__ = "03 Aug 2016"
__version__ = "0.13.3"

import re

from conkit.io._parser import ContactFileParser
from conkit.core.contact import Contact
from conkit.core.contactmap import ContactMap
from conkit.core.contactfile import ContactFile

RE_SPLIT = re.compile(r"\s+")


class ComsatParser(ContactFileParser):
    """Class to parse a COMSAT contact file
    """

    def __init__(self):
        super(ComsatParser, self).__init__()

    def read(self, f_handle, f_id="comsat"):
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
            line = line.rstrip()

            if not line:
                continue

            else:
                res1_seq, res1, res2_seq, res2, _ = RE_SPLIT.split(line)
                contact = Contact(int(res1_seq), int(res2_seq), 0.0)
                contact.res1 = res1
                contact.res2 = res2

                contact_map.add(contact)

        contact_file.method = "Contact map predicted using COMSAT"

        return contact_file

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

        """
        contact_file = self._reconstruct(hierarchy)
        if len(contact_file) > 1:
            raise RuntimeError("More than one contact map provided")
        content = ""
        for contact_map in contact_file:
            for contact in contact_map:
                line = "{res1_seq}{sep}{res1}{sep}{res2_seq}{sep}{res2}{sep}Hx-Hx\n"
                line = line.format(
                    res1_seq=contact.res1_seq, res2_seq=contact.res2_seq, res1=contact.res1, res2=contact.res2, sep="\t"
                )
                content += line
        f_handle.write(content)
