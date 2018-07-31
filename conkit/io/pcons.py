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
Parser module specific to Pcons predictions

This module can be used to parse all versions
of the Pcons programs, i.e. PconsC, PconsC2,
and PconsC3.
"""

__author__ = "Felix Simkovic"
__date__ = "26 Oct 2016"
__version__ = "0.1"

import os
import re

from conkit.io._parser import ContactFileParser
from conkit.core.contact import Contact
from conkit.core.contactmap import ContactMap
from conkit.core.contactfile import ContactFile
from conkit.core.sequence import Sequence

RE_COMMENT = re.compile(r"^#+\s*$")
RE_JUNK = re.compile(r"^(PconsC3|Total|Sequence number|Sequence length).*\s*$")
RE_GENERATED = re.compile(r'^Generated.*$')
RE_SEQUENCE_NAME = re.compile(r"^Sequence name:\s+(.*)\s*$")
RE_SEQUENCE = re.compile(r"^Sequence:\s*$")
RE_PRED_CONTACTS = re.compile(r"^Predicted\s+contacts:\s*$")
RE_CONTACT_HEADER = re.compile(r"^Res1\s+Res2\s+Score\s*$")
RE_CONTACT = re.compile(r"^\s*(\d+)\s+(\d+)\s+(-?\d*\.\d+|\d+)\s*$")


class PconsParser(ContactFileParser):
    """Class to parse a Pcons output

    This module can be used to parse all versions of the
    Pcons programs, i.e. PconsC, PconsC2, and PconsC3.

    """

    def __init__(self):
        super(PconsParser, self).__init__()

    def read(self, f_handle, f_id="pcons"):
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
        contact_map = ContactMap("1")
        contact_file.add(contact_map)

        lines = iter([l.rstrip() for l in f_handle if l.rstrip()])
        done = object()
        line = next(lines, done)

        seq = ''
        seq_id = 'seq_1'

        while line is not done:

            if not line:
                pass

            elif RE_GENERATED.match(line):
                contact_file.remark = line

            elif RE_SEQUENCE_NAME.match(line):
                seq_id = RE_SEQUENCE_NAME.match(line).group(1)

            elif RE_SEQUENCE.match(line):
                line = next(lines, done)
                while line is not done:
                    if not line:
                        break
                    elif RE_CONTACT_HEADER.match(line):
                        break
                    elif RE_PRED_CONTACTS.match(line):
                        break
                    elif RE_CONTACT.match(line):
                        break
                    else:
                        seq += line
                    line = next(lines, done)

            if RE_CONTACT.match(line):
                res1_seq, res2_seq, raw_score = line.split()
                contact = Contact(int(res1_seq), int(res2_seq), float(raw_score))
                contact_map.add(contact)

            line = next(lines, done)

        if seq:
            contact_map.sequence = Sequence(seq_id, seq)

        contact_file.method = 'Contact map predicted using Pcons'

        return contact_file

    def write(self, f_handle, hierarchy, write_header_footer=True):
        """Write a contact file instance to to file

        Default format is ``PconsC3`` style, including comments
        and sequence information (if provided).

        Parameters
        ----------
        f_handle
           Open file handle [write permissions]
        hierarchy : :obj:`~conkit.core.contactfile.ContactFile`, :obj:`~conkit.core.contactmap.ContactMap`
                    or :obj:`~conkit.core.contact.Contact`
        write_header_footer : bool
           Write a PconsC3-typical header

        Raises
        ------
        :exc:`RuntimeError`
           More than one contact map in the hierarchy

        """
        # Double check the type of hierarchy and reconstruct if necessary
        contact_file = self._reconstruct(hierarchy)

        if len(contact_file) > 1:
            raise RuntimeError('More than one contact map provided')

        comment_line = "##############################################################################"

        content = ""

        for contact_map in contact_file:
            if write_header_footer:
                content += comment_line + os.linesep
                content += "PconsC3 result file" + os.linesep
                content += "Generated using ConKit" + os.linesep
                content += comment_line + os.linesep

                if contact_map.sequence is not None:
                    content += "Sequence number: 1" + os.linesep
                    content += "Sequence name: {0}".format(contact_map.sequence.id) + os.linesep
                    content += "Sequence length: {0} aa.".format(contact_map.sequence.seq_len) + os.linesep
                    content += "Sequence:" + os.linesep
                    content += contact_map.sequence.seq + os.linesep * 3

                content += "Predicted contacts:" + os.linesep
                content += "Res1 Res2 Score" + os.linesep

            for contact in contact_map:
                res1_seq = contact.res1_seq
                res2_seq = contact.res2_seq
                raw_score = contact.raw_score
                l = "{res1_seq:>4} {res2_seq:>4} {raw_score:>.6f}".format(
                    res1_seq=res1_seq, res2_seq=res2_seq, raw_score=raw_score)
                content += l + os.linesep

            if write_header_footer:
                content += os.linesep + comment_line + os.linesep

        f_handle.write(content)
