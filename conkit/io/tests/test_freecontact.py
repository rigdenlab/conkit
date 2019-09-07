"""Testing facility for conkit.io.FreeContactIO"""

__author__ = "Felix Simkovic"
__date__ = "26 Oct 2016"

import os
import unittest

from conkit.core.contact import Contact
from conkit.core.contactfile import ContactFile
from conkit.core.contactmap import ContactMap
from conkit.core.sequence import Sequence
from conkit.io.freecontact import FreeContactParser
from conkit.io._iotools import create_tmp_f


class TestFreeContactParser(unittest.TestCase):
    def test_read_1(self):
        content = """1 M 2 V 0.0338619 0
1 M 3 G 0.0307956 0
1 M 4 L 0.0268079 0
1 M 5 T 0.0219783 0
1 M 6 T 0.0222061 0
1 M 7 L 0.0213079 0
1 M 8 F 0.0119054 0
1 M 9 W 0.0275182 0
1 M 10 L 0.0134577 0
1 M 11 G 0.0234555 0
"""
        f_name = create_tmp_f(content=content)
        self.addCleanup(os.remove, f_name)
        with open(f_name, "r") as f_in:
            contact_file = FreeContactParser().read(f_in)
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(10, len(contact_map1))
        self.assertEqual([1, 1, 1, 1, 1, 1, 1, 1, 1, 1], [c.res1_seq for c in contact_map1])
        self.assertEqual([2, 3, 4, 5, 6, 7, 8, 9, 10, 11], [c.res2_seq for c in contact_map1])
        self.assertEqual(
            [
                0.0338619,
                0.0307956,
                0.0268079,
                0.0219783,
                0.0222061,
                0.0213079,
                0.0119054,
                0.0275182,
                0.0134577,
                0.0234555,
            ],
            [c.raw_score for c in contact_map1],
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
            FreeContactParser().write(f_out, contact_file)
        content = ["1 H 9 L 0.7 0", "1 H 10 L 0.7 0", "2 L 8 I 0.9 0", "3 E 12 K 0.4 0"]
        with open(f_name, "r") as f_in:
            output = f_in.read().splitlines()
        self.assertEqual(content, output)

    def tearDown(self):
        self.doCleanups()


if __name__ == "__main__":
    unittest.main(verbosity=2)
