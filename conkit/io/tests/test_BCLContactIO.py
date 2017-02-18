"""Testing facility for conkit.io.BCLContactIO"""

__author__ = "Felix Simkovic"
__date__ = "12 Dec 2016"

from conkit.io.BCLContactIO import BCLContactParser
from conkit.io._iotools import create_tmp_f

import os
import unittest


class Test(unittest.TestCase):

    def test_read_1(self):
        content = """5 I    9 Q 0.000 0.286 0.185 0.836 0.875 0.749
5 I   10 R 0.000 0.000 0.105 0.875 0.482 0.634
5 I   11 I 0.000 0.178 0.066 0.730 0.876 0.727
5 I   21 I 0.030 0.021 0.233 0.645 0.733 0.557
5 I   58 G 0.000 0.054 0.010 0.642 0.799 0.535
6 T   62 V 0.000 0.000 0.027 0.485 0.428 0.585
6 T   63 S 0.000 0.004 0.051 0.547 0.387 0.529
6 T   78 L 0.000 0.000 0.039 0.624 0.384 0.581
6 T   79 T 0.000 0.000 0.036 0.657 0.415 0.679
6 T   80 I 0.000 0.076 0.003 0.513 0.386 0.578
6 T   94 Q 0.000 0.068 0.041 0.534 0.489 0.679
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = BCLContactParser().read(f_in)
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(11, len(contact_map1))
        self.assertEqual([5, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6], [c.res1_seq for c in contact_map1])
        self.assertEqual([9, 10, 11, 21, 58, 62, 63, 78, 79, 80, 94], [c.res2_seq for c in contact_map1])
        self.assertEqual(
            [0.749, 0.634, 0.727, 0.557, 0.535, 0.585, 0.529, 0.581, 0.679, 0.578, 0.679],
            [c.raw_score for c in contact_map1]
        )
        os.unlink(f_name)


if __name__ == "__main__":
    unittest.main(verbosity=2)
