"""Testing facility for conkit.io.FastaIO"""

__author__ = "Felix Simkovic"
__date__ = "09 Sep 2016"

import os
import unittest

from conkit.io.clustal import ClustalParser
from conkit.io._iotools import create_tmp_f


class TestClustalParser(unittest.TestCase):
    def test_read_1(self):
        seq = """CLUSTAL W

seq_0           MLDLEVVPE-RSLGNEQW-------E-F-TLG-MPLAQAV-AILQKHC--

seq_0           -RIIKNVQV

"""
        f_name = create_tmp_f(content=seq)
        parser = ClustalParser()
        with open(f_name, "r") as f_in:
            sequence_file = parser.read(f_in)
        sequence_entry = sequence_file.top_sequence
        ref_id = "seq_0"
        self.assertEqual(ref_id, sequence_entry.id)
        ref_seq = "MLDLEVVPE-RSLGNEQW-------E-F-TLG-MPLAQAV-AILQKHC---RIIKNVQV"
        self.assertEqual(ref_seq, sequence_entry.seq)
        os.unlink(f_name)

    def test_read_2(self):
        msa = """CLUSTAL W
seq_0           AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
seq_1           BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
seq_2           CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
                **************************************************

seq_0           AAAAAAAAA
seq_1           BBBBBBBBB
seq_2           CCCCCCCCC
                *********

"""
        f_name = create_tmp_f(content=msa)
        parser = ClustalParser()
        with open(f_name, "r") as f_in:
            sequence_file = parser.read(f_in)
        for i, sequence_entry in enumerate(sequence_file):
            if i == 0:
                self.assertEqual("seq_0", sequence_entry.id)
                self.assertEqual("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", sequence_entry.seq)
            elif i == 1:
                self.assertEqual("seq_1", sequence_entry.id)
                self.assertEqual("BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB", sequence_entry.seq)
            elif i == 2:
                self.assertEqual("seq_2", sequence_entry.id)
                self.assertEqual("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", sequence_entry.seq)
        os.unlink(f_name)

    def test_read_3(self):
        msa = """CLUSTAL FORMAT for

seq_0           AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA
seq_1           BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB
seq_2           CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

seq_0           AAAAAAAAA
seq_1           BBBBBBBBB
seq_2           CCCCCCCCC

"""
        f_name = create_tmp_f(content=msa)
        parser = ClustalParser()
        with open(f_name, "r") as f_in:
            sequence_file = parser.read(f_in)
        for i, sequence_entry in enumerate(sequence_file):
            if i == 0:
                self.assertEqual("seq_0", sequence_entry.id)
                self.assertEqual("AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA", sequence_entry.seq)
            elif i == 1:
                self.assertEqual("seq_1", sequence_entry.id)
                self.assertEqual("BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB", sequence_entry.seq)
            elif i == 2:
                self.assertEqual("seq_2", sequence_entry.id)
                self.assertEqual("CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC", sequence_entry.seq)
        os.unlink(f_name)

    def test_write_1(self):
        seq = [
            "CLUSTAL FORMAT written with ConKit",
            "",
            "seq_0\tAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
            "seq_1\tBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBB",
            "" "seq_0\tAAAAAAAAAAAAAAAAAAAAAA",
            "seq_1\tBBBBBBBBBBBBBBBBBBBBBB",
        ]
        joinedseq = "\n".join(seq)
        f_name_in = create_tmp_f(content=joinedseq)
        f_name_out = create_tmp_f()
        parser = ClustalParser()
        with open(f_name_in, "r") as f_in, open(f_name_out, "w") as f_out:
            sequence_file = parser.read(f_in)
            parser.write(f_out, sequence_file)
        with open(f_name_out, "r") as f_in:
            output = f_in.read().splitlines()
        self.assertEqual(seq, output)
        map(os.unlink, [f_name_in, f_name_out])


if __name__ == "__main__":
    unittest.main(verbosity=2)
