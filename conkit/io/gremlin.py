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
Parser module specific to Gremlin predictions
"""

__author__ = "Felix Simkovic"
__date__ = "04 Oct 2016"
__version__ = "0.1"

import os
import re

from conkit.io._parser import ContactFileParser
from conkit.core.contact import Contact
from conkit.core.contactmap import ContactMap
from conkit.core.contactfile import ContactFile

RE_HEADER_INTRA = re.compile(r'^i\s+j\s+i_id\s+j_id\s+r_sco\s+s_sco\s+prob$')
RE_HEADER_INTER = re.compile(r'^i\s+j\s+gene\s+i_id\s+j_id\s+r_sco\s+s_sco\s+prob\s+I_prob$')
RE_COMMENT = re.compile(r'^#+(.*)$')
RE_SPLIT = re.compile(r'\s+')


class GremlinParser(ContactFileParser):
    """Parser class for GREMLIN contact prediction file
    """

    def __init__(self):
        super(GremlinParser, self).__init__()

    def read(self, f_handle, f_id="gremlin"):
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

        lines = iter([l.rstrip() for l in f_handle if l.rstrip()])
        done = object()
        line = next(lines, done)

        inter = False
        chain_list = set()
        contact_list = []
        while line is not done:

            if RE_COMMENT.match(line):
                hierarchy.remark = RE_COMMENT.match(line).group(1)

            elif RE_HEADER_INTRA.match(line):
                inter = False
            elif RE_HEADER_INTER.match(line):
                inter = True
            else:
                if inter:
                    res1_seq, res2_seq, chain, _, _, raw_score, scalar_score, _, _ = RE_SPLIT.split(line)
                else:
                    res1_seq, res2_seq, _, _, raw_score, scalar_score, _ = RE_SPLIT.split(line)
                    chain = 'UNK'

                c = Contact(int(res1_seq), int(res2_seq), float(raw_score))
                c.scalar_score = float(scalar_score)

                if chain == 'UNK':
                    chain_list.add('UNK')
                elif len(chain) == 1:
                    c.res1_chain = chain[0]
                    c.res2_chain = chain[0]
                    chain_list.add((c.res1_chain, c.res2_chain))
                elif len(chain) == 2:
                    c.res1_chain = chain[0]
                    c.res2_chain = chain[1]
                    chain_list.add((c.res1_chain, c.res2_chain))
                elif len(chain) > 2:
                    raise ValueError('Cannot distinguish between chains')

                contact_list.append(c)

            line = next(lines, done)

        chain_list = list(chain_list)
        if len(chain_list) == 1 and chain_list[0] == 'UNK':
            contact_map = ContactMap('1')
            for c in contact_list:
                contact_map.add(c)
            hierarchy.add(contact_map)
        elif len(chain_list) == 1:
            chain = chain_list[0]
            map_id = chain[0] if chain[0] == chain[1] else "".join(chain)
            contact_map = ContactMap(map_id)
            for c in contact_list:
                contact_map.add(c)
            hierarchy.add(contact_map)
        else:
            for chain in chain_list:
                map_id = chain[0] if chain[0] == chain[1] else "".join(chain)
                contact_map = ContactMap(map_id)
                for c in contact_list:
                    if c.res1_chain == chain[0] and c.res2_chain == chain[1]:
                        contact_map.add(c)
                hierarchy.add(contact_map)

        hierarchy.sort('id', inplace=True)
        return hierarchy

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

        content = ""

        if contact_file.top_map.top_contact.res1_chain and contact_file.top_map.top_contact.res2_chain:
            header_args = ['i', 'j', 'gene', 'i_id', 'j_id', 'r_sco', 's_sco', 'prob', 'I_prob']
            content += '\t'.join(header_args) + os.linesep

            out_kwargs = [
                '{res1_seq}', '{res2_seq}', '{chains}', '{res1_code}', '{res2_code}', '{raw_score}', '{scalar_score}',
                '1.0', 'N/A'
            ]

        else:
            header_args = ['i', 'j', 'i_id', 'j_id', 'r_sco', 's_sco', 'prob']
            content += '\t'.join(header_args) + os.linesep

            out_kwargs = [
                '{res1_seq}', '{res2_seq}', '{res1_code}', '{res2_code}', '{raw_score}', '{scalar_score}', '1.0'
            ]

        for contact_map in contact_file:
            contact_map.set_scalar_score()
            for c in contact_map:
                res1_code = str(c.res1_seq) + '_' + c.res1
                res2_code = str(c.res2_seq) + '_' + c.res2

                if c.res1_chain == c.res2_chain:
                    chains = c.res1_chain
                else:
                    chains = "{0}{1}".format(c.res1_chain, c.res2_chain)

                out_line = '\t'.join(out_kwargs)
                out_line = out_line.format(
                    res1_seq=c.res1_seq,
                    res2_seq=c.res2_seq,
                    res1_code=res1_code,
                    res2_code=res2_code,
                    chains=chains,
                    raw_score=c.raw_score,
                    scalar_score=round(c.scalar_score, 1))

                content += out_line + os.linesep

        f_handle.write(content)
