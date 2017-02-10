"""
Parser module specific to CASP-RR predictions
"""

__author__ = "Felix Simkovic"
__date__ = "03 Aug 2016"
__version__ = "1.0"

from conkit.core import Contact
from conkit.core import ContactFile
from conkit.core import ContactMap
from conkit.core import Sequence
from conkit.io._ParserIO import _ContactFileParser

import collections
import os
import re


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
ContactTemplate = collections.namedtuple('Contact', ['res1_seq', 'res2_seq', 'lb', 'ub', 'raw_score',
                                                     'res1_chain', 'res2_chain', 'res1_altseq', 'res2_altseq'])

class CaspParser(_ContactFileParser):
    """Parser class for CASP RR contact prediction file

    """

    def __init__(self):
        super(CaspParser, self).__init__()

    def read(self, f_handle, f_id="casp"):
        """Read a contact file into a :obj:`conkit.core.ContactFile` instance

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

                        contact = Contact(int(res1_seq), int(res2_seq), float(raw_score),
                                          distance_bound=(int(lb), int(ub)))
                        contact.res1_chain = res1_chain
                        contact.res2_chain = res2_chain
                        contact.res1_altseq = int(res1_seq)
                        contact.res2_altseq = int(res2_seq)
                        contact_map.add(contact)

                if seq_chunks:
                    seq = "".join(seq_chunks)
                    sequence = Sequence('seq_{0}'.format(contact_map.id), seq)
                    contact_map.sequence = sequence
                    contact_map.assign_sequence_register()
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
        hierarchy : :obj:`ContactFile <conkit.core.ContactFile>`, :obj:`ContactMap <conkit.core.ContactMap>`
                    or :obj:`Contact <conkit.core.Contact>`

        """

        # Double check the type of hierarchy and reconstruct if necessary
        contact_file = self._reconstruct(hierarchy)

        f_handle.write("PFRMAT RR" + os.linesep)
        if contact_file.target:
            f_handle.write("TARGET {0}".format(contact_file.target) + os.linesep)
        if contact_file.author:
            f_handle.write("AUTHOR {0}".format(contact_file.author) + os.linesep)
        if contact_file.remark:
            for remark in contact_file.remark:
                f_handle.write("REMARK {0}".format(remark) + os.linesep)
        if contact_file.method:
            for method in contact_file.method:
                f_handle.write("METHOD {0}".format(method) + os.linesep)

        for contact_map in contact_file:

            f_handle.write("MODEL  {0}".format(contact_map.id) + os.linesep)
            if isinstance(contact_map.sequence, Sequence):
                sequence = contact_map.sequence
                for i in range(0, sequence.seq_len, 50):
                    f_handle.write(sequence.seq[i:i+50] + os.linesep)
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
                s = s.format(res1_chain=res1_chain, res1_seq=contact.res1_seq, res2_chain=res2_chain,
                             res2_seq=contact.res2_seq, lb=contact.distance_bound[0], ub=contact.distance_bound[1],
                             raw_score=contact.raw_score)
                f_handle.write(s + os.linesep)
            f_handle.write("ENDMDL" + os.linesep)

        f_handle.write("END" + os.linesep)

        return
