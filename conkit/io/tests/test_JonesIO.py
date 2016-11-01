"""Testing facility for conkit.io.JonesIO"""

__author__ = "Felix Simkovic"
__date__ = "13 Sep 2016"

from conkit.io.JonesIO import JonesIO

import os
import unittest
import tempfile


def _create_tmp(data=None):
    f_in = tempfile.NamedTemporaryFile(delete=False)
    if data:
        f_in.write(data)
    f_in.close()
    return f_in.name


class Test(unittest.TestCase):

    def test_read(self):
        # ==================================================
        # Multiple sequence alignment
        msa = """GSMFTPKPPQDSAVI--GYCVKQGAVMKNWKRRY--LDENTIGYF
EVHK--ECKQSDIMMRD--FEIVTTSRTFYVQADSPEEMHSWIKA
EVHKVQECK--DIMMRDNLFEI--TSRTFWKRRY--LDENTIGYF
EVHKVQECK--DIMMRDNLFEI--TSRTF--RRY--LDENTIGYF
"""
        f_name = _create_tmp(msa)
        parser = JonesIO()
        sequence_file = parser.read(open(f_name, 'r'))
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
        del parser, sequence_file, sequence_entry
        os.unlink(f_name)

        # ==================================================
        # Multiple sequence alignment - crash
        msa = """>header1
GSMFTPKPPQDSAVI--GYCVKQGAVMKNWKRRY--LDENTIGYF
>header2
EVHK--ECKQSDIMMRD--FEIVTTSRTFYVQADSPEEMHSWIKA
>header3
EVHKVQECK--DIMMRDNLFEI--TSRTFWKRRY--LDENTIGYF
>header4
EVHKVQECK--DIMMRDNLFEI--TSRTF--RRY--LDENTIGYF
"""
        f_name = _create_tmp(msa)
        parser = JonesIO()
        self.assertRaises(ValueError, parser.read, open(f_name, 'r'))
        del parser
        os.unlink(f_name)

    def test_write(self):
        # ==================================================
        # Multiple sequence alignment
        msa = """GSMFTPKPPQDSAVI--GYCVKQGAVMKNWKRRY--LDENTIGYF
EVHK--ECKQSDIMMRD--FEIVTTSRTFYVQADSPEEMHSWIKA
EVHKVQECK--DIMMRDNLFEI--TSRTFWKRRY--LDENTIGYF
EVHKVQECK--DIMMRDNLFEI--TSRTF--RRY--LDENTIGYF
"""
        f_name_in = _create_tmp(msa)
        parser = JonesIO()
        sequence_file = parser.read(open(f_name_in, 'r'))
        f_name_out = _create_tmp()
        parser.write(open(f_name_out, 'w'), sequence_file)
        output = "".join(open(f_name_out, 'r').readlines())
        self.assertEqual(msa, output)
        del parser, sequence_file
        os.unlink(f_name_in)
        os.unlink(f_name_out)

if __name__ == "__main__":
    unittest.main(verbosity=2)
