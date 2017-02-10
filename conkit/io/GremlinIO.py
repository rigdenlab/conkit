"""
Parser module specific to Gremlin predictions
"""

__author__ = "Felix Simkovic"
__date__ = "04 Oct 2016"
__version__ = 0.1

from conkit.core import Contact
from conkit.core import ContactMap
from conkit.core import ContactFile
from conkit.io._ParserIO import _ContactFileParser

import os
import re

RE_HEADER_INTRA = re.compile(r'^i\s+j\s+i_id\s+j_id\s+r_sco\s+s_sco\s+prob$')
RE_HEADER_INTER = re.compile(r'^i\s+j\s+gene\s+i_id\s+j_id\s+r_sco\s+s_sco\s+prob\s+I_prob$')
RE_COMMENT = re.compile(r'^#+(.*)$')
RE_SPLIT = re.compile(r'\s+')


class GremlinParser(_ContactFileParser):
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
        :obj:`ContactFile <conkit.core.ContactFile>`

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
        hierarchy : :obj:`ContactFile <conkit.core.ContactFile>`, :obj:`ContactMap <conkit.core.ContactMap>`
                    or :obj:`Contact <conkit.core.Contact>`

        """
        # Double check the type of hierarchy and reconstruct if necessary
        contact_file = self._reconstruct(hierarchy)

        if contact_file.top_map.top_contact.res1_chain and contact_file.top_map.top_contact.res2_chain:
            header_args = ['i', 'j', 'gene', 'i_id', 'j_id', 'r_sco', 's_sco', 'prob', 'I_prob']
            f_handle.write('\t'.join(header_args) + os.linesep)

            out_kwargs = ['{res1_seq}', '{res2_seq}', '{chains}', '{res1_code}', '{res2_code}',
                          '{raw_score}', '{scalar_score}', '1.0', 'N/A']

        else:
            header_args = ['i', 'j', 'i_id', 'j_id', 'r_sco', 's_sco', 'prob']
            f_handle.write('\t'.join(header_args) + os.linesep)

            out_kwargs = ['{res1_seq}', '{res2_seq}', '{res1_code}', '{res2_code}',
                          '{raw_score}', '{scalar_score}', '1.0']

        for contact_map in contact_file:
            contact_map.calculate_scalar_score()
            for c in contact_map:
                res1_code = str(c.res1_seq) + '_' + c.res1
                res2_code = str(c.res2_seq) + '_' + c.res2

                if c.res1_chain == c.res2_chain:
                    chains = c.res1_chain
                else:
                    chains = "{0}{1}".format(c.res1_chain, c.res2_chain)

                out_line = '\t'.join(out_kwargs)
                out_line = out_line.format(res1_seq=c.res1_seq, res2_seq=c.res2_seq, res1_code=res1_code,
                                           res2_code=res2_code, chains=chains,
                                           raw_score=c.raw_score, scalar_score=round(c.scalar_score, 1))

                f_handle.write(out_line + os.linesep)

        return
