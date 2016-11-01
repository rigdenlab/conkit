"""Testing facility for conkit.io.Bbcontacts"""

__author__ = "Felix Simkovic"
__date__ = "26 Oct 2016"

from conkit.io.BbcontactsIO import BbcontactsParser

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
        content = """#identifier diversity     direction viterbiscore indexpred        state  res1  res2
1EAZ      0.65  Antiparallel     9.860725         1        first    29    24
1EAZ      0.65  Antiparallel     9.860725         1     internal    30    23
1EAZ      0.65  Antiparallel     9.860725         1         last    31    22
1EAZ      0.65      Parallel    -6.855870        29        first    87    54
1EAZ      0.65      Parallel    -6.855870        29     internal    88    55
1EAZ      0.65      Parallel    -6.855870        29         last    89    56
"""
        f_name = _create_tmp(content)
        contact_file = BbcontactsParser().read(open(f_name, 'r'))
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(6, len(contact_map1))
        self.assertEqual([24, 23, 22, 54, 55, 56], [c.res1_seq for c in contact_map1])
        self.assertEqual([29, 30, 31, 87, 88, 89], [c.res2_seq for c in contact_map1])
        self.assertItemsEqual(
            [9.860725, 9.860725, 9.860725, -6.855870, -6.855870, -6.855870], [c.raw_score for c in contact_map1]
        )
        os.unlink(f_name)


if __name__ == "__main__":
    unittest.main(verbosity=2)
