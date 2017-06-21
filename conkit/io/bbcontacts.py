# BSD 3-Clause License
#
# Copyright (c) 2016-17, University of Liverpool
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
__date__ = "26 Oct 2016"
__version__ = "0.1"

import re

from conkit.io._parser import _ContactFileParser
from conkit.core.contact import Contact
from conkit.core.contactmap import ContactMap
from conkit.core.contactfile import ContactFile

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
                if any(value == "NA" for value in [raw_score, res2_seq, res1_seq]):
                    continue
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
