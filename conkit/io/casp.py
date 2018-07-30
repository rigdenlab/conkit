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
Parser module specific to CASP-RR predictions
"""

__author__ = "Felix Simkovic"
__date__ = "03 Aug 2016"
__version__ = "1.0"

import collections
import os
import re

from conkit.io._parser import ContactFileParser
from conkit.core.contact import Contact
from conkit.core.contactmap import ContactMap
from conkit.core.contactfile import ContactFile
from conkit.core.sequence import Sequence

# Credits to Stefan Seemayer - bits taken from PyMOL-RR
RE_PRFMAT = re.compile(r'PFRMAT\s+(RR)\s*$')
RE_TARGET = re.compile(r'^TARGET\s+(.*?)\s*$')
RE_AUTHOR = re.compile(r'^AUTHOR\s+(.*?)\s*$')
RE_REMARK = re.compile(r'^REMARK\s+(.*?)\s*$')
RE_METHOD = re.compile(r'^METHOD\s+(.*?)\s*$')
RE_MODEL = re.compile(r'^MODEL\s+(\w+)\s*$')
RE_SEQ = re.compile(r'^([A-Za-z\-]+)$')
RE_SPLIT = re.compile(r'\s+')
RE_RES = re.compile(r'([A-Za-z]+)([0-9]+)')
RE_ENDMDL = re.compile(r'^ENDMDL\s*$')
RE_END = re.compile(r'^END\s*$')

# Intermediate storage structures
ModelTemplate = collections.namedtuple('Model', ['id', 'contacts', 'sequence'])
ContactTemplate = collections.namedtuple(
    'Contact',
    ['res1_seq', 'res2_seq', 'lb', 'ub', 'raw_score', 'res1_chain', 'res2_chain', 'res1_altseq', 'res2_altseq'])


class CaspParser(ContactFileParser):
    """Parser class for CASP RR contact prediction file"""

    def __init__(self):
        super(CaspParser, self).__init__()

    def read(self, f_handle, f_id="casp"):
        """Read a contact file into a :obj:`~conkit.core.contactfile.ContactFile` instance

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
        lines = [l.strip() for l in f_handle.readlines()]

        contact_file = ContactFile(f_id)

        it = iter(lines)
        while True:

            try:
                line = next(it)
            except StopIteration:
                break

            if RE_PRFMAT.match(line):
                continue

            elif RE_TARGET.match(line):
                contact_file.remark = RE_TARGET.match(line).group(1)

            elif RE_AUTHOR.match(line):
                contact_file.author = RE_AUTHOR.match(line).group(1)

            elif RE_REMARK.match(line):
                contact_file.remark = RE_REMARK.match(line).group(1)

            elif RE_METHOD.match(line):
                contact_file.method = RE_METHOD.match(line).group(1)

            elif RE_MODEL.match(line):
                contact_map = ContactMap(RE_MODEL.match(line).group(1))

                seq_chunks = []

                while True:

                    try:
                        line = next(it)
                    except StopIteration:
                        break

                    if not line:
                        break

                    if RE_ENDMDL.match(line):
                        break

                    elif RE_END.match(line):
                        break

                    elif RE_SEQ.match(line):
                        seq_chunks.append(line)

                    else:
                        res1_entry, res2_entry, lb, ub, raw_score = RE_SPLIT.split(line)

                        # Split in case we have chain in inter-molecular scenarios
                        res1_split = RE_RES.split(res1_entry)
                        if len(res1_split) == 1:
                            res1_chain, res1_seq = '', res1_split[0]
                        elif len(res1_split) == 4:
                            res1_chain, res1_seq = res1_split[1], res1_split[2]
                        res2_split = RE_RES.split(res2_entry)

                        if len(res2_split) == 1:
                            res2_chain, res2_seq = '', res2_split[0]
                        elif len(res2_split) == 4:
                            res2_chain, res2_seq = res2_split[1], res2_split[2]

                        contact = Contact(
                            int(res1_seq), int(res2_seq), float(raw_score), distance_bound=(float(lb), float(ub)))
                        contact.res1_chain = res1_chain
                        contact.res2_chain = res2_chain
                        contact.res1_altseq = int(res1_seq)
                        contact.res2_altseq = int(res2_seq)
                        contact_map.add(contact)

                if seq_chunks:
                    seq = "".join(seq_chunks)
                    sequence = Sequence('seq_{0}'.format(contact_map.id), seq)
                    contact_map.sequence = sequence
                    contact_map.set_sequence_register()
                contact_file.add(contact_map)

            elif RE_END.match(line):
                break

            else:
                raise ValueError('Unrecognized line type. Please report this issue')

        return contact_file

    def write(self, f_handle, hierarchy):
        """Write a contact file instance to to file

        Parameters
        ----------
        f_handle
           Open file handle [write permissions]
        hierarchy : :obj:`~conkit.core.contactfile.ContactFile`, :obj:`~conkit.core.contactmap.ContactMap`
                    or :obj:`~conkit.core.contact.Contact`

        """

        # Double check the type of hierarchy and reconstruct if necessary
        contact_file = self._reconstruct(hierarchy)

        content = "PFRMAT RR" + os.linesep
        if contact_file.target:
            content += "TARGET {0}".format(contact_file.target) + os.linesep
        if contact_file.author:
            content += "AUTHOR {0}".format(contact_file.author) + os.linesep
        if contact_file.remark:
            for remark in contact_file.remark:
                content += "REMARK {0}".format(remark) + os.linesep
        if contact_file.method:
            for method in contact_file.method:
                content += "METHOD {0}".format(method) + os.linesep

        for contact_map in contact_file:

            content += "MODEL  {0}".format(contact_map.id) + os.linesep
            if isinstance(contact_map.sequence, Sequence):
                sequence = contact_map.sequence
                for i in range(0, sequence.seq_len, 50):
                    content += sequence.seq[i:i + 50] + os.linesep
            # Casp Roll format specifies raw scores to be in [0, 1]
            if any(c.raw_score > 1.0 or c.raw_score < 0.0 for c in contact_map):
                contact_map.rescale(inplace=True)
            for contact in contact_map:
                s = '{res1_chain: <}{res1_seq: <4} {res2_chain: <}{res2_seq:<4} {lb: <3} {ub: <3} {raw_score: <.6f}'
                if contact.res1_chain == contact.res2_chain:
                    res1_chain = res2_chain = ""
                else:
                    res1_chain = contact.res1_chain
                    res2_chain = contact.res2_chain
                lb = int(contact.lower_bound) if float(contact.lower_bound).is_integer() else contact.lower_bound
                ub = int(contact.upper_bound) if float(contact.upper_bound).is_integer() else contact.upper_bound
                s = s.format(
                    res1_chain=res1_chain,
                    res1_seq=contact.res1_seq,
                    res2_chain=res2_chain,
                    res2_seq=contact.res2_seq,
                    lb=lb,
                    ub=ub,
                    raw_score=contact.raw_score)
                content += s + os.linesep
            content += "ENDMDL" + os.linesep

        content += "END" + os.linesep

        f_handle.write(content)
