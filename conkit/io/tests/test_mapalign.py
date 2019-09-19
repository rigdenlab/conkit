"""Testing facility for conkit.io.map_align"""

import os
import unittest

from conkit.core.contact import Contact
from conkit.core.contactfile import ContactFile
from conkit.core.contactmap import ContactMap
from conkit.core.sequence import Sequence
from conkit.io.mapalign import MapAlignParser
from conkit.io._iotools import create_tmp_f


class TestMapAlignParser(unittest.TestCase):
    def test_read_1(self):
        content = """LEN     77
CON     10      14      1
CON     1       7       0.883
CON     10      13      0.871
CON     6       9       0.847
CON     5       71      0.816
CON     36      42      0.807
CON     1       73      0.806
CON     33      37      0.8
CON     21      68      0.563
CON     38      57      0.561
PRF     0       T       X       0.0321  0.0078  0.0582  0.1021  0.0207  0.0384  0.038   0.0386  0.0697  0.0376  0.0586
PRF     1       M       X       0.0228  0.0052  0.0103  0.0112  0.1236  0.0083  0.0106  0.1182  0.0104  0.3046  0.1922
PRF     2       K       X       0.0322  0.0048  0.0563  0.128   0.0116  0.0481  0.0332  0.0319  0.1099  0.0366  0.0228
PRF     3       I       X       0.0147  0.0062  0.0025  0.0034  0.2786  0.0058  0.0057  0.1876  0.0053  0.2779  0.0302
PRF     4       I       X       0.0404  0.0057  0.1139  0.2398  0.0104  0.028   0.0275  0.0278  0.0491  0.0283  0.0051
"""
        f_name = create_tmp_f(content=content)
        self.addCleanup(os.remove, f_name)
        with open(f_name, "r") as f_in:
            contact_file = MapAlignParser().read(f_in)
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(10, len(contact_map1))
        self.assertEqual([10, 1, 10, 6, 5, 36, 1, 33, 21, 38], [c.res1_seq for c in contact_map1])
        self.assertEqual([14, 7, 13, 9, 71, 42, 73, 37, 68, 57], [c.res2_seq for c in contact_map1])
        self.assertEqual(
            [1.0, 0.883, 0.871, 0.847, 0.816, 0.807, 0.806, 0.8, 0.563, 0.561], [c.raw_score for c in contact_map1]
        )

    def test_write_1(self):
        contact_file = ContactFile("RR")
        contact_file.target = "R9999"
        contact_file.author = "1234-5678-9000"
        contact_file.remark = ["Predictor remarks"]
        contact_file.method = ["Description of methods used", "Description of methods used"]
        contact_map = ContactMap("1")
        contact_file.add(contact_map)
        for c in [(1, 9, 0, 8, 0.7), (1, 10, 0, 8, 0.7), (2, 8, 0, 8, 0.9), (3, 12, 0, 8, 0.4)]:
            contact = Contact(c[0], c[1], c[4], distance_bound=(c[2], c[3]))
            contact_map.add(contact)
        contact_map.sequence = Sequence("1", "HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD")
        contact_map.set_sequence_register()
        f_name = create_tmp_f()
        self.addCleanup(os.remove, f_name)
        with open(f_name, "w") as f_out:
            MapAlignParser().write(f_out, contact_file)
        content = ["LEN 12", "CON 1 9 0.700000", "CON 1 10 0.700000", "CON 2 8 0.900000", "CON 3 12 0.400000"]
        with open(f_name, "r") as f_in:
            output = f_in.read().splitlines()
        self.assertEqual(content, output)

    def tearDown(self):
        self.doCleanups()


if __name__ == "__main__":
    unittest.main(verbosity=2)
