"""Testing facility for conkit.io.ComsatIO"""

__author__ = "Felix Simkovic"
__date__ = "14 Sep 2016"

from conkit.core import Contact
from conkit.core import ContactFile
from conkit.core import ContactMap
from conkit.core import Sequence
from conkit.io.ComsatIO import ComsatParser
from conkit.io._iotools import create_tmp_f

import os
import unittest


class Test(unittest.TestCase):

    def test_read_1(self):
        content = """19   A   41   A   H1-H2
19   A   42   C   H1-H2
11   L   47   L   H1-H2
11   L   48   L   H1-H2
12   L   47   L   H1-H2
12   L   48   L   H1-H2
40   I   66   I   H2-H3
41   A   66   I   H2-H3
33   Y   73   H   H2-H3
33   Y   74   A   H2-H3
46   L   62   L   H2-H3
47   L   62   L   H2-H3
69   M   88   V   H3-H4
69   M   89   A   H3-H4
96   A   117  V   H4-H5
96   A   118  A   H4-H5
82   A   129  A   H4-H5
82   A   130  G   H4-H5
82   A   133  F   H4-H5
83   F   133  F   H4-H5
128  I   154  I   H5-H6
129  A   154  I   H5-H6
118  A   163  A   H5-H6
119  V   163  A   H5-H6
20   A   160  A   H1-H6
21   L   160  A   H1-H6
8    N   171  V   H1-H6
9    V   171  V   H1-H6
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = ComsatParser().read(f_in, 'r')
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(28, len(contact_map1))
        self.assertEqual(
            [19, 19, 11, 11, 12, 12, 40, 41, 33, 33, 46, 47, 69, 69, 96, 96, 82, 82, 82, 83, 128, 129, 118, 119, 20, 21, 8, 9],
            [c.res1_seq for c in contact_map1]
        )
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
        contact_map.sequence = Sequence('1', 'HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD')
        contact_map.assign_sequence_register()
        f_name = create_tmp_f()
        with open(f_name, 'w') as f_out:
            ComsatParser().write(f_out, contact_file)
        content = [
            "1	H	9	L	Hx-Hx",
            "1	H	10	L	Hx-Hx",
            "2	L	8	I	Hx-Hx",
            "3	E	12	K	Hx-Hx",
            "",
        ]
        content = os.linesep.join(content)
        with open(f_name, 'r') as f_in:
            data = "".join(f_in.readlines())
        self.assertEqual(content, data)
        os.unlink(f_name)


if __name__ == "__main__":
    unittest.main(verbosity=2)
