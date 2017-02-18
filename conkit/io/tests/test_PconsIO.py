"""Testing facility for conkit.io.PconsIO"""

__author__ = "Felix Simkovic"
__date__ = "26 Oct 2016"

from conkit.core import Contact
from conkit.core import ContactFile
from conkit.core import ContactMap
from conkit.core import Sequence
from conkit.io.PconsIO import PconsParser
from conkit.io._iotools import create_tmp_f

import os
import unittest


class Test(unittest.TestCase):

    def test_read_1(self):
        content = """1 2 0.93514
1 3 0.67324
1 4 0.23692
1 5 0.13166
1 6 0.09188
1 7 0.07957
1 8 0.06556
1 9 0.05188
1 10 0.04146
1 11 0.03264
1 12 0.02515
1 13 0.02137
1 14 0.01961
1 15 0.01710
1 16 0.01397
1 17 0.01192
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = PconsParser().read(f_in)
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(16, len(contact_map1))
        self.assertEqual([1] * 16, [c.res1_seq for c in contact_map1])
        self.assertEqual(list(range(2, 18)), [c.res2_seq for c in contact_map1])
        self.assertEqual([0.93514, 0.67324, 0.23692, 0.13166, 0.09188], [c.raw_score for c in contact_map1][:5])
        os.unlink(f_name)

    def test_read_2(self):
        content = """# Check one two
Hello WOrld
1 2 0.93514
1 3 0.67324
1 4 0.23692
1 5 0.13166
1 6 0.09188
1 7 0.07957
1 8 0.06556
1 9 0.05188
1 10 0.04146
1 11 0.03264
1 12 0.02515
1 13 0.02137
1 14 0.01961
1 15 0.01710
1 16 0.01397
1 17 0.01192
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = PconsParser().read(f_in)
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(16, len(contact_map1))
        self.assertEqual([1] * 16, [c.res1_seq for c in contact_map1])
        self.assertEqual(list(range(2, 18)), [c.res2_seq for c in contact_map1])
        self.assertEqual([0.93514, 0.67324, 0.23692, 0.13166, 0.09188], [c.raw_score for c in contact_map1][:5])
        os.unlink(f_name)

    def test_read_3(self):
        content = """##############################################################################
PconsC3 result file
Generated from test_remark
Total request time: 83148.0 seconds.
##############################################################################
Sequence number: 1
Sequence name: shorter_test
Sequence length: 33 aa.
Sequence:
HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD


Predicted contacts:
Res1 Res2 Score
1 2 0.93514
1 3 0.67324
1 4 0.23692
1 5 0.13166
1 6 0.09188
1 7 0.07957
1 8 0.06556
1 9 0.05188
1 10 0.04146
1 11 0.03264
1 12 0.02515
1 13 0.02137
1 14 0.01961
1 15 0.01710
1 16 0.01397
1 17 0.01192

##############################################################################
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = PconsParser().read(f_in)
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(['Generated from test_remark'], contact_file.remark)
        self.assertEqual(16, len(contact_map1))
        self.assertEqual([1] * 16, [c.res1_seq for c in contact_map1])
        self.assertEqual(list(range(2, 18)), [c.res2_seq for c in contact_map1])
        self.assertEqual([0.93514, 0.67324, 0.23692, 0.13166, 0.09188], [c.raw_score for c in contact_map1][:5])
        self.assertEqual('shorter_test', contact_map1.sequence.id)
        self.assertEqual('HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD', contact_map1.sequence.seq)
        os.unlink(f_name)

    def test_read_4(self):
        content = """##############################################################################
PconsC3 result file
Generated from test_remark
Total request time: 83148.0 seconds.
##############################################################################
Sequence number: 1
Sequence name: longer_test
Sequence length: 132 aa.
Sequence:
HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSDHLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD
HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSDHLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD


Predicted contacts:
Res1 Res2 Score
1 2 0.93514
1 3 0.67324
1 4 0.23692
1 5 0.13166
1 6 0.09188
1 7 0.07957
1 8 0.06556
1 9 0.05188
1 10 0.04146
1 11 0.03264
1 12 0.02515
1 13 0.02137
1 14 0.01961
1 15 0.01710
1 16 0.01397
1 17 0.01192

##############################################################################
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = PconsParser().read(f_in)
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(['Generated from test_remark'], contact_file.remark)
        self.assertEqual(16, len(contact_map1))
        self.assertEqual([1] * 16, [c.res1_seq for c in contact_map1])
        self.assertEqual(list(range(2, 18)), [c.res2_seq for c in contact_map1])
        self.assertEqual([0.93514, 0.67324, 0.23692, 0.13166, 0.09188], [c.raw_score for c in contact_map1][:5])
        self.assertEqual('longer_test', contact_map1.sequence.id)
        self.assertEqual('HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSDHLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD'
                         'HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSDHLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD',
                         contact_map1.sequence.seq)
        os.unlink(f_name)

    def test_write_1(self):
        contact_file = ContactFile('RR')
        contact_file.target = 'R9999'
        contact_file.author = '1234-5678-9000'
        contact_file.remark = ['Predictor remarks']
        contact_file.method = ['Description of methods used', 'Description of methods used']
        contact_map = ContactMap('1')
        contact_file.add(contact_map)
        for c in [(1, 9, 0, 8, 0.7), (1, 10, 0, 8, 0.7), (2, 8, 0, 8, 0.9), (3, 12, 0, 8, 0.4)]:
            contact = Contact(c[0], c[1], c[4], distance_bound=(c[2], c[3]))
            contact_map.add(contact)
        contact_map.sequence = Sequence('sequence_1', 'HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD')
        contact_map.assign_sequence_register()
        f_name = create_tmp_f()
        with open(f_name, 'w') as f_out:
            PconsParser().write(f_out, contact_file)
        content = [
            "##############################################################################",
            "PconsC3 result file",
            "Generated from ConKit",
            "##############################################################################",
            "Sequence number: 1",
            "Sequence name: sequence_1",
            "Sequence length: 33 aa.",
            "Sequence:",
            "HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD",
            "",
            "",
            "Predicted contacts:",
            "Res1 Res2 Score",
            "   1    9 0.700000",
            "   1   10 0.700000",
            "   2    8 0.900000",
            "   3   12 0.400000",
            "",
            "##############################################################################",
            "",
        ]
        content = os.linesep.join(content)
        with open(f_name, 'r') as f_in:
            data = "".join(f_in.readlines())
        self.assertEqual(content, data)
        os.unlink(f_name)


if __name__ == "__main__":
    unittest.main(verbosity=2)
