"""Testing facility for conkit.io.JonesIO"""

__author__ = "Felix Simkovic"
__date__ = "13 Sep 2016"

from conkit.io.JonesIO import JonesIO
from conkit.io._iotools import create_tmp_f

import os
import unittest


class Test(unittest.TestCase):

    def test_read_1(self):
        msa = """GSMFTPKPPQDSAVI--GYCVKQGAVMKNWKRRY--LDENTIGYF
EVHK--ECKQSDIMMRD--FEIVTTSRTFYVQADSPEEMHSWIKA
EVHKVQECK--DIMMRDNLFEI--TSRTFWKRRY--LDENTIGYF
EVHKVQECK--DIMMRDNLFEI--TSRTF--RRY--LDENTIGYF
"""
        f_name = create_tmp_f(content=msa)
        parser = JonesIO()
        with open(f_name, 'r') as f_in:
            sequence_file = parser.read(f_in)
        for i, sequence_entry in enumerate(sequence_file):
            if i == 0:
                self.assertEqual('seq_0', sequence_entry.id)
                self.assertEqual('GSMFTPKPPQDSAVI--GYCVKQGAVMKNWKRRY--LDENTIGYF', sequence_entry.seq)
            elif i == 1:
                self.assertEqual('seq_1', sequence_entry.id)
                self.assertEqual('EVHK--ECKQSDIMMRD--FEIVTTSRTFYVQADSPEEMHSWIKA', sequence_entry.seq)
            elif i == 2:
                self.assertEqual('seq_2', sequence_entry.id)
                self.assertEqual('EVHKVQECK--DIMMRDNLFEI--TSRTFWKRRY--LDENTIGYF', sequence_entry.seq)
            elif i == 3:
                self.assertEqual('seq_3', sequence_entry.id)
                self.assertEqual('EVHKVQECK--DIMMRDNLFEI--TSRTF--RRY--LDENTIGYF', sequence_entry.seq)
        os.unlink(f_name)

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
        f_name = create_tmp_f(content=msa)
        parser = JonesIO()
        with open(f_name, 'r') as f_in:
            with self.assertRaises(ValueError):
                parser.read(f_in)
        os.unlink(f_name)

    def test_write_1(self):
        msa = [
            "GSMFTPKPPQDSAVI--GYCVKQGAVMKNWKRRY--LDENTIGYF",
            "EVHK--ECKQSDIMMRD--FEIVTTSRTFYVQADSPEEMHSWIKA",
            "EVHKVQECK--DIMMRDNLFEI--TSRTFWKRRY--LDENTIGYF",
            "EVHKVQECK--DIMMRDNLFEI--TSRTF--RRY--LDENTIGYF",
            "",
        ]
        msa = os.linesep.join(msa)
        f_name_in = create_tmp_f(content=msa)
        f_name_out = create_tmp_f()
        parser = JonesIO()
        with open(f_name_in, 'r') as f_in, open(f_name_out, 'w') as f_out:
            sequence_file = parser.read(f_in)
            parser.write(f_out, sequence_file)
        with open(f_name_out, 'r') as f_in:
            output = "".join(f_in.readlines())
        self.assertEqual(msa, output)
        os.unlink(f_name_in)
        os.unlink(f_name_out)


if __name__ == "__main__":
    unittest.main(verbosity=2)
