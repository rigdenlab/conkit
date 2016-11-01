"""Testing facility for conkit.io.PsicovIO"""

__author__ = "Felix Simkovic"
__date__ = "26 Oct 2016"

from conkit.core import Contact
from conkit.core import ContactFile
from conkit.core import ContactMap
from conkit.core import Sequence
from conkit.io.PsicovIO import PsicovParser

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
        # Test Case 1
        content = """46 78 0 8 9.301869
80 105 0 8 8.856009
111 129 0 8 7.252451
75 205 0 8 6.800462
19 44 0 8 6.588349
111 130 0 8 6.184269
23 41 0 8 6.163786
171 205 0 8 5.519271
53 126 0 8 5.440612
100 140 0 8 5.382865
"""
        f_name = _create_tmp(content)
        contact_file = PsicovParser().read(open(f_name, 'r'))
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(10, len(contact_map1))
        self.assertEqual([46, 80, 111, 75, 19, 111, 23, 171, 53, 100], [c.res1_seq for c in contact_map1])
        self.assertEqual([78, 105, 129, 205, 44, 130, 41, 205, 126, 140], [c.res2_seq for c in contact_map1])
        self.assertItemsEqual(
            [9.301869, 8.856009, 7.252451, 6.800462, 6.588349, 6.184269, 6.163786,
             5.519271, 5.440612, 5.382865],
            [c.raw_score for c in contact_map1]
        )
        os.unlink(f_name)

    def test_write(self):
        # ======================================================
        # Test Case 1
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
        f_name = _create_tmp()
        PsicovParser().write(open(f_name, 'w'), contact_file)
        content = """1 9 0 8 0.700000
1 10 0 8 0.700000
2 8 0 8 0.900000
3 12 0 8 0.400000
"""
        data = "".join(open(f_name, 'r').readlines())
        self.assertEqual(content, data)
        os.unlink(f_name)


if __name__ == "__main__":
    unittest.main(verbosity=2)
