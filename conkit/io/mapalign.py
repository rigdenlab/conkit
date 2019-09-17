# BSD 3-Clause License
#
# Copyright (c) 2016-19, University of Liverpool
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
Parser module specific to map_align map files
"""

from conkit.io._parser import ContactFileParser
from conkit.core.contact import Contact
from conkit.core.contactmap import ContactMap
from conkit.core.contactfile import ContactFile


class MapAlignParser(ContactFileParser):
    """Class to parse a map_align map file
    """

    def read(self, f_handle, f_id="map_align"):
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

        hierarchy = ContactFile(f_id)
        _map = ContactMap("map_1")
        hierarchy.add(_map)

        for line in f_handle:
            line = line.strip().split()

            if line[0] == "CON" and line[1].isdigit() and line[2].isdigit():
                _contact = Contact(int(line[1]), int(line[2]), float(line[3]))
                _map.add(_contact)

        hierarchy.method = "Contact map compatible with map_algin"

        return hierarchy

    def write(self, f_handle, hierarchy):
        """Write a contact file instance to a file

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
        cmap = contact_file.top_map
        content = "LEN {}\n".format(cmap.highest_residue_number)
        line_template = "CON {} {} {:.6f}\n"
        for contact in cmap:
            content += line_template.format(contact.res1_seq, contact.res2_seq, contact.raw_score)
        f_handle.write(content)
