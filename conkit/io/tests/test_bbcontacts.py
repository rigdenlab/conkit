"""Testing facility for conkit.io.Bbcontacts"""

__author__ = "Felix Simkovic"
__date__ = "26 Oct 2016"

import os
import unittest

from conkit.io.bbcontacts import BbcontactsParser
from conkit.io._iotools import create_tmp_f


class TestBbcontactsParser(unittest.TestCase):
    def test_read_1(self):
        content = """#identifier diversity     direction viterbiscore indexpred        state  res1  res2
1EAZ      0.65  Antiparallel     9.860725         1        first    29    24
1EAZ      0.65  Antiparallel     9.860725         1     internal    30    23
1EAZ      0.65  Antiparallel     9.860725         1         last    31    22
1EAZ      0.65      Parallel    -6.855870        29        first    87    54
1EAZ      0.65      Parallel    -6.855870        29     internal    88    55
1EAZ      0.65      Parallel    -6.855870        29         last    89    56
"""
        f_name = create_tmp_f(content=content)
        self.addCleanup(os.remove, f_name)
        with open(f_name, "r") as f_in:
            contact_file = BbcontactsParser().read(f_in)
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(6, len(contact_map1))
        self.assertEqual([24, 23, 22, 54, 55, 56], [c.res1_seq for c in contact_map1])
        self.assertEqual([29, 30, 31, 87, 88, 89], [c.res2_seq for c in contact_map1])
        self.assertEqual(
            sorted([9.860725, 9.860725, 9.860725, -6.855870, -6.855870, -6.855870]),
            sorted([c.raw_score for c in contact_map1]),
        )

    def test_read_2(self):
        content = """#identifier diversity     direction viterbiscore indexpred        state  res1  res2
1EAZ      0.65  Antiparallel     9.860725         1        first    29    24
1EAZ      0.65  Antiparallel     9.860725         1         last    30    23
1EAZ      0.65      Parallel    -6.855870        29        first    87    54
"""
        f_name = create_tmp_f(content=content)
        self.addCleanup(os.remove, f_name)
        with open(f_name, "r") as f_in:
            contact_file = BbcontactsParser().read(f_in, del_one_two=True)
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(0, len(contact_map1))

    def test_read_3(self):
        content = """#identifier diversity     direction viterbiscore indexpred        state  res1  res2
1EAZ      0.65  Antiparallel     9.860725         1        first    29    24
1EAZ      0.65  Antiparallel     9.860725         1     internal    30    23
1EAZ      0.65  Antiparallel     9.860725         1         last    31    22
1EAZ      0.65      Parallel    -6.855870        29        first    87    54
1EAZ      0.65      Parallel    -6.855870        29     internal    88    55
1EAZ      0.65      Parallel    -6.855870        29         last    89    56
1EAZ      0.65  Antiparallel     0.000000         1        first   100    24
1EAZ      0.65  Antiparallel     0.000000         1         last   101    23
1EAZ      0.65      Parallel     0.000000        29        first   100    15
"""
        f_name = create_tmp_f(content=content)
        self.addCleanup(os.remove, f_name)
        with open(f_name, "r") as f_in:
            contact_file = BbcontactsParser().read(f_in, del_one_two=False)
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(9, len(contact_map1))
        self.assertEqual([24, 23, 22, 54, 55, 56, 24, 23, 15], [c.res1_seq for c in contact_map1])
        self.assertEqual([29, 30, 31, 87, 88, 89, 100, 101, 100], [c.res2_seq for c in contact_map1])
        self.assertEqual(
            sorted([9.860725, 9.860725, 9.860725, -6.855870, -6.855870, -6.855870, 0.0, 0.0, 0.0]),
            sorted([c.raw_score for c in contact_map1]),
        )

    def tearDown(self):
        self.doCleanups()


if __name__ == "__main__":
    unittest.main(verbosity=2)
