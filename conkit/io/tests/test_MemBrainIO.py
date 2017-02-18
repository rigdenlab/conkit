"""Testing facility for conkit.io.MemBrainIO"""

__author__ = "Felix Simkovic"
__date__ = "26 Oct 2016"

from conkit.core import Contact
from conkit.core import ContactFile
from conkit.core import ContactMap
from conkit.core import Sequence
from conkit.io.MemBrainIO import MemBrainParser
from conkit.io._iotools import create_tmp_f

import os
import unittest


class Test(unittest.TestCase):

    def test_read_1(self):
        content = """Helix   Position        Residue Helix   Position        Residue Probability
H1      30      F       H2      55      F       1.000000
H1      33      L       H2      51      A       0.944091
H1      18      G       H2      65      C       0.942259
H1      30      F       H2      54      G       0.919241
H1      26      I       H2      57      L       0.817638
H1      18      G       H2      58      S       0.797449
H1      33      L       H2      63      L       0.795520
H1      12      A       H2      68      V       0.795462
H1      29      V       H2      55      F       0.791829
H1      24      I       H2      51      A       0.790044
H1      19      L       H2      62      G       0.784613
H1      19      L       H2      55      F       0.782741
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = MemBrainParser().read(f_in)
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(12, len(contact_map1))
        self.assertEqual([30, 33, 18, 30, 26, 18, 33, 12, 29, 24, 19, 19], [c.res1_seq for c in contact_map1])
        self.assertEqual([55, 51, 65, 54, 57, 58, 63, 68, 55, 51, 62, 55], [c.res2_seq for c in contact_map1])
        self.assertEqual(
            [1.000000, 0.944091, 0.942259, 0.919241, 0.817638, 0.797449, 0.795520,
             0.795462, 0.791829, 0.790044, 0.784613, 0.782741],
            [c.raw_score for c in contact_map1]
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
            MemBrainParser().write(f_out, contact_file)
        content = [
            "Helix   Position        Residue Helix   Position        Residue Probability",
            "Hx      1       H       Hx      9       L       0.700000",
            "Hx      1       H       Hx      10      L       0.700000",
            "Hx      2       L       Hx      8       I       0.900000",
            "Hx      3       E       Hx      12      K       0.400000",
            "",
        ]
        content = os.linesep.join(content)
        with open(f_name, 'r') as f_in:
            data = "".join(f_in.readlines())
        self.assertEqual(content, data)
        os.unlink(f_name)


if __name__ == "__main__":
    unittest.main(verbosity=2)
