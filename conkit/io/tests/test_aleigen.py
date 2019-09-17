"""Testing facility for conkit.io.aleigen"""

import os
import unittest

from conkit.core.contact import Contact
from conkit.core.contactfile import ContactFile
from conkit.core.contactmap import ContactMap
from conkit.core.sequence import Sequence
from conkit.io.aleigen import AleigenParser
from conkit.io._iotools import create_tmp_f


class TestAleigenParser(unittest.TestCase):
    def test_read_1(self):
        content = """77
10 14
1 7
10 13
6 9
5 71
36 42
1 73
33 37
21 68
38 57
"""
        f_name = create_tmp_f(content=content)
        self.addCleanup(os.remove, f_name)
        with open(f_name, "r") as f_in:
            contact_file = AleigenParser().read(f_in)
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(10, len(contact_map1))
        self.assertEqual([10, 1, 10, 6, 5, 36, 1, 33, 21, 38], [c.res1_seq for c in contact_map1])
        self.assertEqual([14, 7, 13, 9, 71, 42, 73, 37, 68, 57], [c.res2_seq for c in contact_map1])
        self.assertEqual([0.5 for x in range(0, 10)], [c.raw_score for c in contact_map1])

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
            AleigenParser().write(f_out, contact_file)
        content = ["12", "1 9", "1 10", "2 8", "3 12"]
        with open(f_name, "r") as f_in:
            output = f_in.read().splitlines()
        self.assertEqual(content, output)

    def tearDown(self):
        self.doCleanups()


if __name__ == "__main__":
    unittest.main(verbosity=2)
