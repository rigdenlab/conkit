"""
Parser module specific to Pcons predictions

This module can be used to parse all versions
of the Pcons programs, i.e. PconsC, PconsC2,
and PconsC3.
"""

__author__ = "Felix Simkovic"
__date__ = "26 Oct 2016"
__version__ = 0.1

from conkit.core import Contact
from conkit.core import ContactMap
from conkit.core import ContactFile
from conkit.core import Sequence
from conkit.io._ParserIO import _ContactFileParser

import os
import re

RE_COMMENT = re.compile(r"^#+\s*$")
RE_JUNK = re.compile(r"^(PconsC3|Total|Sequence number|Sequence length).*\s*$")
RE_GENERATED = re.compile(r'^Generated.*$')
RE_SEQUENCE_NAME = re.compile(r"^Sequence name:\s+(.*)\s*$")
RE_SEQUENCE = re.compile(r"^Sequence:\s*$")
RE_PRED_CONTACTS = re.compile(r"^Predicted\s+contacts:\s*$")
RE_CONTACT_HEADER = re.compile(r"^Res1\s+Res2\s+Score\s*$")
RE_CONTACT = re.compile(r"^\s*(\d+)\s+(\d+)\s+(\d*\.\d+|\d+)\s*$")
RE_SPLIT = re.compile(r'\s+')


class PconsParser(_ContactFileParser):
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
        :obj:`ContactFile <conkit.core.ContactFile>`

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
                res1_seq, res2_seq, raw_score = RE_SPLIT.split(line)
                contact = Contact(int(res1_seq), int(res2_seq), float(raw_score))
                contact_map.add(contact)

            line = next(lines, done)

        if seq:
            contact_map.sequence = Sequence(seq_id, seq)

        contact_file.method = 'Contact map predicted using Pcons'

        return contact_file

    def write(self, f_handle, hierarchy):
        """Write a contact file instance to to file

        Default format is ``PconsC3`` style, including comments
        and sequence information (if provided).

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

        comment_line = "##############################################################################"

        for contact_map in contact_file:
            f_handle.write(comment_line + os.linesep)
            f_handle.write("PconsC3 result file" + os.linesep)
            f_handle.write("Generated from ConKit" + os.linesep)
            f_handle.write(comment_line + os.linesep)

            if contact_map.sequence is not None:
                f_handle.write("Sequence number: 1" + os.linesep)
                f_handle.write("Sequence name: {0}".format(contact_map.sequence.id) + os.linesep)
                f_handle.write("Sequence length: {0} aa.".format(contact_map.sequence.seq_len) + os.linesep)
                f_handle.write("Sequence:" + os.linesep)
                f_handle.write(contact_map.sequence.seq + os.linesep * 3)

            f_handle.write("Predicted contacts:" + os.linesep)
            f_handle.write("Res1 Res2 Score" + os.linesep)

            for contact in contact_map:
                res1_seq = contact.res1_seq
                res2_seq = contact.res2_seq
                raw_score = contact.raw_score
                l = "{res1_seq:>4} {res2_seq:>4} {raw_score:>.6f}".format(res1_seq=res1_seq, res2_seq=res2_seq,
                                                                          raw_score=raw_score)
                f_handle.write(l + os.linesep)

            f_handle.write(os.linesep + comment_line + os.linesep)

        return
