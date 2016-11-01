"""Testing facility for conkit.io.FastaIO"""

__author__ = "Felix Simkovic"
__date__ = "09 Sep 2016"

from conkit.io.FastaIO import FastaIO

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
        # Normal sequence only mode
        seq = """>00FAF_A <unknown description>
GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYFKSELEKEPLRVIPLK
EVHKVQECKQSDIMMRDNLFEIVTTSRTFYVQADSPEEMHSWIKAVSGAIVAQRGPGRSA
SSEHP
"""
        f_name = _create_tmp(seq)
        parser = FastaIO()
        sequence_file = parser.read(open(f_name, 'r'))
        sequence_entry = sequence_file.top_sequence
        ref_id = "00FAF_A <unknown description>"
        self.assertEqual(ref_id, sequence_entry.id)
        ref_seq = "GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYFKSELEKEPLRVIPLKEVHKVQECKQSDIMMRDNLFEIVTTSRTFYVQADSPEEMHSWIKAVSGAIVAQRGPGRSASSEHP"
        self.assertEqual(ref_seq, sequence_entry.seq)
        del parser, sequence_file, sequence_entry
        os.unlink(f_name)

        # ==================================================
        # Comments at the beginning of the file
        seq = """# Hello World
>00FAF_A <unknown description>
GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYFKSELEKEPLRVIPLK
"""
        f_name = _create_tmp(seq)
        parser = FastaIO()
        sequence_file = parser.read(open(f_name, 'r'))
        sequence_entry = sequence_file.top_sequence
        ref_f_remark = [' Hello World']
        self.assertItemsEqual(ref_f_remark, sequence_file.remark)
        ref_id = "00FAF_A <unknown description>"
        self.assertEqual(ref_id, sequence_entry.id)
        ref_seq = "GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYFKSELEKEPLRVIPLK"
        self.assertEqual(ref_seq, sequence_entry.seq)
        del parser, sequence_file, sequence_entry
        os.unlink(f_name)

        # ==================================================
        # Multiple sequence alignment
        msa = """#foo
#bar
>seq1
GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYF
>seq2
EVHKVQECKQSDIMMRDNLFEIVTTSRTFYVQADSPEEMHSWIKA
>seq3
EVHKVQECKQSDIMMRDNLFEIVTTSRTFWKRRYFQLDENTIGYF
"""
        f_name = _create_tmp(msa)
        parser = FastaIO()
        sequence_file = parser.read(open(f_name, 'r'))
        self.assertItemsEqual(['foo', 'bar'], sequence_file.remark)
        for i, sequence_entry in enumerate(sequence_file):
            if i == 0:
                self.assertEqual('seq1', sequence_entry.id)
                self.assertEqual('GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYF', sequence_entry.seq)
            elif i == 1:
                self.assertEqual('seq2', sequence_entry.id)
                self.assertEqual('EVHKVQECKQSDIMMRDNLFEIVTTSRTFYVQADSPEEMHSWIKA', sequence_entry.seq)
            elif i == 2:
                self.assertEqual('seq3', sequence_entry.id)
                self.assertEqual('EVHKVQECKQSDIMMRDNLFEIVTTSRTFWKRRYFQLDENTIGYF', sequence_entry.seq)
        del parser, sequence_file, sequence_entry
        os.unlink(f_name)

    def test_write(self):
        # ==================================================
        # Normal sequence mode
        seq = """>00FAF_A|<unknown description>
GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYFKSELEKEPLRVIPLK
EVHKVQECKQSDIMMRDNLFEIVTTSRTFYVQADSPEEMHSWIKAVSGAIVAQRGPGRSA
SSEHP
"""
        f_name_in = _create_tmp(seq)
        parser = FastaIO()
        sequence_file = parser.read(open(f_name_in, 'r'))
        f_name_out = _create_tmp()
        parser.write(open(f_name_out, 'w'), sequence_file)
        output = "".join(open(f_name_out, 'r').readlines())
        self.assertEqual(seq, output)
        del parser, sequence_file
        os.unlink(f_name_in)
        os.unlink(f_name_out)

        # ==================================================
        # Normal sequence mode - with comment
        seq = """# Hello World
>00FAF_A|<unknown description>
GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYFKSELEKEPLRVIPLK
"""
        f_name_in = _create_tmp(seq)
        parser = FastaIO()
        sequence_file = parser.read(open(f_name_in, 'r'))
        f_name_out = _create_tmp()
        parser.write(open(f_name_out, 'w'), sequence_file)
        output = "".join(open(f_name_out, 'r').readlines())
        self.assertEqual(seq, output)
        del parser, sequence_file
        os.unlink(f_name_in)
        os.unlink(f_name_out)

        # ==================================================
        # Multiple sequence alignment
        msa = """#foo
#bar
>seq1
GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYF
>seq2
EVHKVQECKQSDIMMRDNLFEIVTTSRTFYVQADSPEEMHSWIKA
>seq3
EVHKVQECKQSDIMMRDNLFEIVTTSRTFWKRRYFQLDENTIGYF
"""
        f_name_in = _create_tmp(msa)
        parser = FastaIO()
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
