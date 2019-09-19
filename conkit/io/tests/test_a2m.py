"""Testing facility for conkit.io.a2m"""

__author__ = "Felix Simkovic"
__date__ = "30 Jul 2018"

import unittest

from conkit.io.a2m import A2mParser
from conkit.io.tests.helpers import ParserTestCase


class TestA2mParser(ParserTestCase):

    def test_read_1(self):
        msa = """GSMFTPKPPQDSAVI--GYCVKQGAVMKNWKRRY--LDENTIGYF
EVHK--ECKQSDIMMRD--FEIVTTSRTFYVQADSPEEMHSWIKA
EVHKVQECK--DIMMRDNLFEI--TSRTFWKRRY--LDENTIGYF
EVHKVQECK--DIMMRDNLFEI--TSRTF--RRY--LDENTIGYF
"""
        fname = self._create_iotools_fname(content=msa)
        with open(fname, "r") as f_in:
            sequence_file = A2mParser().read(f_in)
        for i, sequence_entry in enumerate(sequence_file):
            if i == 0:
                self.assertEqual("seq_0", sequence_entry.id)
                self.assertEqual("GSMFTPKPPQDSAVI--GYCVKQGAVMKNWKRRY--LDENTIGYF", sequence_entry.seq)
            elif i == 1:
                self.assertEqual("seq_1", sequence_entry.id)
                self.assertEqual("EVHK--ECKQSDIMMRD--FEIVTTSRTFYVQADSPEEMHSWIKA", sequence_entry.seq)
            elif i == 2:
                self.assertEqual("seq_2", sequence_entry.id)
                self.assertEqual("EVHKVQECK--DIMMRDNLFEI--TSRTFWKRRY--LDENTIGYF", sequence_entry.seq)
            elif i == 3:
                self.assertEqual("seq_3", sequence_entry.id)
                self.assertEqual("EVHKVQECK--DIMMRDNLFEI--TSRTF--RRY--LDENTIGYF", sequence_entry.seq)

    def test_read_2(self):
        msa = """>header1
GSMFTPKPPQDSAVI--GYCVKQGAVMKNWKRRY--LDENTIGYF
>header2
EVHK--ECKQSDIMMRD--FEIVTTSRTFYVQADSPEEMHSWIKA
>header3
EVHKVQECK--DIMMRDNLFEI--TSRTFWKRRY--LDENTIGYF
>header4
EVHKVQECK--DIMMRDNLFEI--TSRTF--RRY--LDENTIGYF
"""
        fname = self._create_iotools_fname(content=msa)
        with open(fname, "r") as f_in:
            with self.assertRaises(ValueError):
                A2mParser().read(f_in)

    def test_write_1(self):
        msa = [
            "GSMFTPKPPQDSAVI--GYCVKQGAVMKNWKRRY--LDENTIGYF",
            "EVHK--ECKQSDIMMRD--FEIVTTSRTFYVQADSPEEMHSWIKA",
            "EVHKVQECK--DIMMRDNLFEI--TSRTFWKRRY--LDENTIGYF",
            "EVHKVQECK--DIMMRDNLFEI--TSRTF--RRY--LDENTIGYF",
        ]
        infname = self._create_iotools_fname(content="\n".join(msa))
        outfname = self._create_iotools_fname()
        with open(infname, "r") as f_in, open(outfname, "w") as f_out:
            sequence_file = parser.read(f_in)
            A2mParser().write(f_out, sequence_file)
        with open(outfname, "r") as f_in:
            output = f_in.read().splitlines()
        self.assertEqual(msa, output)


if __name__ == "__main__":
    unittest.main(verbosity=2)
