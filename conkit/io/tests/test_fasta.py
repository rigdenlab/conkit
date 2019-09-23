"""Testing facility for conkit.io.FastaIO"""

__author__ = "Felix Simkovic"
__date__ = "09 Sep 2016"

import os
import unittest

from conkit.io.fasta import FastaParser
from conkit.io.tests.helpers import ParserTestCase


class TestFastaParser(ParserTestCase):
    def test_read_1(self):
        seq = """>00FAF_A <unknown description>
GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYFKSELEKEPLRVIPLK
EVHKVQECKQSDIMMRDNLFEIVTTSRTFYVQADSPEEMHSWIKAVSGAIVAQRGPGRSA
SSEHP
"""
        f_name = self.tempfile(content=seq)
        parser = FastaParser()
        with open(f_name, "r") as f_in:
            sequence_file = parser.read(f_in)
        sequence_entry = sequence_file.top_sequence
        ref_id = "00FAF_A <unknown description>"
        self.assertEqual(ref_id, sequence_entry.id)
        ref_seq = "GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYFKSELEKEPLRVIPLKEVHKVQECKQSDIMMRDNLFEIVTTSRTFYVQADSPEEMHSWIKAVSGAIVAQRGPGRSASSEHP"
        self.assertEqual(ref_seq, sequence_entry.seq)

    def test_read_2(self):
        seq = """# Hello World
>00FAF_A <unknown description>
GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYFKSELEKEPLRVIPLK
"""
        f_name = self.tempfile(content=seq)
        parser = FastaParser()
        with open(f_name, "r") as f_in:
            sequence_file = parser.read(f_in)
        sequence_entry = sequence_file.top_sequence
        ref_f_remark = [" Hello World"]
        self.assertEqual(ref_f_remark, sequence_file.remark)
        ref_id = "00FAF_A <unknown description>"
        self.assertEqual(ref_id, sequence_entry.id)
        ref_seq = "GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYFKSELEKEPLRVIPLK"
        self.assertEqual(ref_seq, sequence_entry.seq)

    def test_read_3(self):
        msa = """#foo
#bar
>seq1
GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYF
>seq2
EVHKVQECKQSDIMMRDNLFEIVTTSRTFYVQADSPEEMHSWIKA
>seq3
EVHKVQECKQSDIMMRDNLFEIVTTSRTFWKRRYFQLDENTIGYF
"""
        f_name = self.tempfile(content=msa)
        parser = FastaParser()
        with open(f_name, "r") as f_in:
            sequence_file = parser.read(f_in)
        self.assertEqual(["foo", "bar"], sequence_file.remark)
        for i, sequence_entry in enumerate(sequence_file):
            if i == 0:
                self.assertEqual("seq1", sequence_entry.id)
                self.assertEqual("GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYF", sequence_entry.seq)
            elif i == 1:
                self.assertEqual("seq2", sequence_entry.id)
                self.assertEqual("EVHKVQECKQSDIMMRDNLFEIVTTSRTFYVQADSPEEMHSWIKA", sequence_entry.seq)
            elif i == 2:
                self.assertEqual("seq3", sequence_entry.id)
                self.assertEqual("EVHKVQECKQSDIMMRDNLFEIVTTSRTFWKRRYFQLDENTIGYF", sequence_entry.seq)

    def test_write_1(self):
        seq = [
            ">00FAF_A|<unknown description>",
            "GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYFKSELEKEPLRVIPLK",
            "EVHKVQECKQSDIMMRDNLFEIVTTSRTFYVQADSPEEMHSWIKAVSGAIVAQRGPGRSA",
            "SSEHP",
        ]
        f_name_in = self.tempfile(content="\n".join(seq))
        f_name_out = self.tempfile()
        parser = FastaParser()
        with open(f_name_in, "r") as f_in, open(f_name_out, "w") as f_out:
            sequence_file = parser.read(f_in)
            parser.write(f_out, sequence_file)
        with open(f_name_out, "r") as f_in:
            output = f_in.read().splitlines()
        self.assertEqual(seq, output)

    def test_write_2(self):
        seq = [
            "# Hello World",
            ">00FAF_A|<unknown description>",
            "GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYFKSELEKEPLRVIPLK",
        ]
        f_name_in = self.tempfile(content="\n".join(seq))
        f_name_out = self.tempfile()
        parser = FastaParser()
        with open(f_name_in, "r") as f_in, open(f_name_out, "w") as f_out:
            sequence_file = parser.read(f_in)
            parser.write(f_out, sequence_file)
        with open(f_name_out, "r") as f_in:
            output = f_in.read().splitlines()
        self.assertEqual(seq, output)

    def test_write_3(self):
        msa = [
            "#foo",
            "#bar",
            ">seq1",
            "GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYF",
            ">seq2",
            "EVHKVQECKQSDIMMRDNLFEIVTTSRTFYVQADSPEEMHSWIKA",
            ">seq3",
            "EVHKVQECKQSDIMMRDNLFEIVTTSRTFWKRRYFQLDENTIGYF",
        ]
        f_name_in = self.tempfile(content="\n".join(msa))
        f_name_out = self.tempfile()
        parser = FastaParser()
        with open(f_name_in, "r") as f_in, open(f_name_out, "w") as f_out:
            sequence_file = parser.read(f_in)
            parser.write(f_out, sequence_file)
        with open(f_name_out, "r") as f_in:
            output = f_in.read().splitlines()
        self.assertEqual(msa, output)


if __name__ == "__main__":
    unittest.main(verbosity=2)
