"""Testing facility for conkit.io.CaspIO"""

__author__ = "Felix Simkovic"
__date__ = "17 Aug 2016"

from conkit.core import Contact
from conkit.core import ContactFile
from conkit.core import ContactMap
from conkit.core import Sequence
from conkit.io.CaspIO import CaspParser
from conkit.io._iotools import create_tmp_f

import os
import unittest


class Test(unittest.TestCase):

    def test_read_1(self):
        content = """PFRMAT RR
TARGET R9999
AUTHOR 1234-5678-9000
REMARK Predictor remarks
METHOD Description of methods used
METHOD Description of methods used
MODEL  1
HLEGSIGILLKKHEIVFDGC
HDFGRTYIWQMSD
1  9  0  8  0.70
1 10  0  8  0.70
1 12  0  8  0.60
2  8  0  8  0.90
3  7  0  8  0.70
3 12  0  8  0.40
4  6  0  8  0.90
7 14  0  8  0.30
9 14  0  8  0.50
END
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = CaspParser().read(f_in)
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(9, len(contact_map1))
        self.assertEqual("HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD", contact_map1.sequence.seq)
        self.assertEqual("HLEG-IGILL-K-E-------------------", contact_map1.repr_sequence.seq)
        os.unlink(f_name)

    def test_read_2(self):
        content = """PFRMAT RR
TARGET R9999
AUTHOR 1234-5678-9000
REMARK Predictor remarks
METHOD Description of methods used
MODEL  1
1  9  0  8  0.70
1 10  0  8  0.70
1 12  0  8  0.60
2  8  0  8  0.90
3  7  0  8  0.70
3 12  0  8  0.40
4  6  0  8  0.90
7 14  0  8  0.30
9 14  0  8  0.50
END
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = CaspParser().read(f_in)
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(9, len(contact_map1))
        self.assertIsNone(contact_map1.sequence)
        with self.assertRaises(TypeError):
            _ = contact_map1.repr_sequence
        os.unlink(f_name)

    def test_read_3(self):
        content = """PFRMAT RR
MODEL  1
1  9  0  8  0.70
1 10  0  8  0.70
1 12  0  8  0.60
2  8  0  8  0.90
3  7  0  8  0.70
3 12  0  8  0.40
4  6  0  8  0.90
7 14  0  8  0.30
9 14  0  8  0.50
END
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = CaspParser().read(f_in)
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(9, len(contact_map1))
        self.assertIsNone(contact_map1.sequence)
        with self.assertRaises(TypeError):
            _ = contact_map1.repr_sequence
        os.unlink(f_name)

    def test_read_4(self):
        content = """PFRMAT RR
TARGET R9999
AUTHOR 1234-5678-9000
REMARK Predictor remarks
METHOD Description of methods used
METHOD Description of methods used
MODEL  1
HLEGSIGILLKKHEIVFDGC
HDFGRTYIWQMSD
1  9  0  8  0.70
1 10  0  8  0.70
1 12  0  8  0.60
2  8  0  8  0.90
3  7  0  8  0.70
3 12  0  8  0.40
4  6  0  8  0.90
7 14  0  8  0.30
9 14  0  8  0.50
ENDMDL
MODEL  2
HLEGSIGILLKKHEIVFDGC
HDFGRTYIWQM
1 10  0  8  0.70
1 12  0  8  0.60
2  8  0  8  0.90
3  7  0  8  0.70
3 12  0  8  0.40
4  6  0  8  0.90
7 14  0  8  0.30
9 14  0  8  0.50
ENDMDL
END
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = CaspParser().read(f_in)
        self.assertEqual(2, len(contact_file))
        contact_map1 = contact_file[0]
        self.assertEqual(9, len(contact_map1))
        self.assertEqual("HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD", contact_map1.sequence.seq)
        self.assertEqual("HLEG-IGILL-K-E-------------------", contact_map1.repr_sequence.seq)
        contact_map2 = contact_file['2']
        self.assertEqual(8, len(contact_map2))
        self.assertEqual("HLEGSIGILLKKHEIVFDGCHDFGRTYIWQM", contact_map2.sequence.seq)
        self.assertEqual("HLEG-IGILL-K-E-----------------", contact_map2.repr_sequence.seq)
        os.unlink(f_name)

    def test_read_5(self):
        content = """PFRMAT RR
TARGET R9999
MODEL  1
1  9  0  8  0.70
1 10  0  8  0.70
1 12  0  8  0.60
2  8  0  8  0.90
3  7  0  8  0.70
3 12  0  8  0.40
4  6  0  8  0.90
7 14  0  8  0.30
9 14  0  8  0.50
ENDMDL
MODEL  2
HLEGSIGILLKKHEIVFDGC
HDFGRTYIWQMSD
1  9  0  8  0.70
1 10  0  8  0.70
1 12  0  8  0.60
2  8  0  8  0.90
3  7  0  8  0.70
3 12  0  8  0.40
4  6  0  8  0.90
7 14  0  8  0.30
9 14  0  8  0.50
END
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = CaspParser().read(f_in)
        self.assertEqual(2, len(contact_file))
        contact_map1 = contact_file["1"]
        self.assertEqual(9, len(contact_map1))
        self.assertIsNone(contact_map1.sequence)
        contact_map2 = contact_file["2"]
        self.assertEqual(9, len(contact_map2))
        self.assertEqual("HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD", contact_map2.sequence.seq)
        self.assertEqual("HLEG-IGILL-K-E-------------------", contact_map2.repr_sequence.seq)
        os.unlink(f_name)

    def test_read_6(self):
        content = """1  9  0  8  0.70
1 10  0  8  0.70
1 12  0  8  0.60
2  8  0  8  0.90
3  7  0  8  0.70
3 12  0  8  0.40
4  6  0  8  0.90
7 14  0  8  0.30
9 14  0  8  0.50
"""
        f_name = create_tmp_f(content=content)
        parser = CaspParser()
        with open(f_name, 'r') as f_in:
            with self.assertRaises(ValueError):
                CaspParser().read(f_in)
        os.unlink(f_name)

    @unittest.skip("Not yet implemented")
    def test_read_7(self):
        content = """PFRMAT RR
TARGET R9999
AUTHOR 1234-5678-9000
REMARK Predictor remarks
METHOD Description of methods used
METHOD Description of methods used
MODEL  1
HLEGSIGILLKKHEIVFDGC
HDFGRTYIWQMSD
A1 B9   0  8  0.70
A1 B10  0  8  0.70
A1 B12  0  8  0.60
A1 B14  0  8  0.20
A1 B15  0  8  0.10
A1 B17  0  8  0.30
A1 B19  0  8  0.50
A2 B8   0  8  0.90
A3 B7   0  8  0.70
A3 B12  0  8  0.40
A3 B14  0  8  0.70
A3 B15  0  8  0.30
A4 B6   0  8  0.90
A7 B14  0  8  0.30
A9 B14  0  8  0.50
END
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = CaspParser().read(f_in)
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(15, len(contact_map1))
        self.assertEqual("HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD", contact_map1.sequence.seq)
        self.assertEqual("HLEG-IGILL-K-E-------------------", contact_map1.repr_sequence.seq)
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
            CaspParser().write(f_out, contact_file)
        content = [
            "PFRMAT RR",
            "TARGET R9999",
            "AUTHOR 1234-5678-9000",
            "REMARK Predictor remarks",
            "METHOD Description of methods used",
            "METHOD Description of methods used",
            "MODEL  1",
            "HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD",
            "1    9    0   8   0.700000",
            "1    10   0   8   0.700000",
            "2    8    0   8   0.900000",
            "3    12   0   8   0.400000",
            "ENDMDL",
            "END",
            "",
        ]
        content = os.linesep.join(content)
        with open(f_name, 'r') as f_in:
            data = "".join(f_in.readlines())
        self.assertEqual(content, data)
        os.unlink(f_name)

    def test_write_2(self):
        contact_file = ContactFile('RR')
        contact_map = ContactMap('1')
        contact_file.add(contact_map)
        for c in [(1, 9, 0, 8, 0.7), (1, 10, 0, 8, 0.7), (2, 8, 0, 8, 0.9), (3, 12, 0, 8, 0.4)]:
            contact = Contact(c[0], c[1], c[4], distance_bound=(c[2], c[3]))
            contact_map.add(contact)
        contact_map.sequence = Sequence('1', 'HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD')
        f_name = create_tmp_f()
        with open(f_name, 'w') as f_out:
            CaspParser().write(f_out, contact_file)
        content = [
            "PFRMAT RR",
            "MODEL  1",
            "HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD",
            "1    9    0   8   0.700000",
            "1    10   0   8   0.700000",
            "2    8    0   8   0.900000",
            "3    12   0   8   0.400000",
            "ENDMDL",
            "END",
            "",
        ]
        content = os.linesep.join(content)
        with open(f_name, 'r') as f_in:
            data = "".join(f_in.readlines())
        self.assertEqual(content, data)
        os.unlink(f_name)

    def test_write_3(self):
        contact_file = ContactFile('RR')
        contact_map = ContactMap('1')
        contact_file.add(contact_map)
        for c in [(1, 9, 0, 8, 0.7), (1, 10, 0, 8, 0.7), (2, 8, 0, 8, 0.9), (3, 12, 0, 8, 0.4)]:
            contact = Contact(c[0], c[1], c[4], distance_bound=(c[2], c[3]))
            contact_map.add(contact)
        f_name = create_tmp_f()
        with open(f_name, 'w') as f_out:
            CaspParser().write(f_out, contact_file)
        content = [
            "PFRMAT RR",
            "MODEL  1",
            "1    9    0   8   0.700000",
            "1    10   0   8   0.700000",
            "2    8    0   8   0.900000",
            "3    12   0   8   0.400000",
            "ENDMDL",
            "END",
            "",
        ]
        content = os.linesep.join(content)
        with open(f_name, 'r') as f_in:
            data = "".join(f_in.readlines())
        self.assertEqual(content, data)
        os.unlink(f_name)

    def test_write_4(self):
        contact_file = ContactFile('RR')
        contact_map = ContactMap('1')
        contact_file.add(contact_map)
        for c in [(1, 9, 0, 8, 0.7), (1, 10, 0, 8, 0.7), (2, 8, 0, 8, 0.9), (3, 12, 0, 8, 0.4)]:
            contact = Contact(c[0], c[1], c[4], distance_bound=(c[2], c[3]))
            contact_map.add(contact)
        contact_map.sequence = Sequence('1', 'HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSDHLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD')
        f_name = create_tmp_f()
        with open(f_name, 'w') as f_out:
            CaspParser().write(f_out, contact_file)
        content = [
            "PFRMAT RR",
            "MODEL  1",
            "HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSDHLEGSIGILLKKHEIVF",
            "DGCHDFGRTYIWQMSD",
            "1    9    0   8   0.700000",
            "1    10   0   8   0.700000",
            "2    8    0   8   0.900000",
            "3    12   0   8   0.400000",
            "ENDMDL",
            "END",
            "",
        ]
        content = os.linesep.join(content)
        with open(f_name, 'r') as f_in:
            data = "".join(f_in.readlines())
        self.assertEqual(content, data)
        os.unlink(f_name)

    def test_write_5(self):
        contact_file = ContactFile('RR')
        contact_map = ContactMap('1')
        contact_file.add(contact_map)
        for c in [(1, 9, 0, 8, 1.5), (1, 10, 0, 8, -0.3), (2, 8, 0, 8, 0.9), (3, 12, 0, 8, 0.4)]:
            contact = Contact(c[0], c[1], c[4], distance_bound=(c[2], c[3]))
            contact_map.add(contact)
        contact_map.sequence = Sequence('1', 'HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSDHLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD')
        f_name = create_tmp_f()
        with open(f_name, 'w') as f_out:
            CaspParser().write(f_out, contact_file)
        content = [
            "PFRMAT RR",
            "MODEL  1",
            "HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSDHLEGSIGILLKKHEIVF",
            "DGCHDFGRTYIWQMSD",
            "1    9    0   8   1.000000",
            "1    10   0   8   0.000000",
            "2    8    0   8   0.666667",
            "3    12   0   8   0.388889",
            "ENDMDL",
            "END",
            "",
        ]
        content = os.linesep.join(content)
        with open(f_name, 'r') as f_in:
            data = "".join(f_in.readlines())
        self.assertEqual(content, data)
        os.unlink(f_name)

    def test_write_6(self):
        contact_file = ContactFile('RR')
        contact_map = ContactMap('1')
        contact_file.add(contact_map)
        for c in [('A', 1, 'B', 9, 0, 8, 0.7),
                  ('A', 1, 'B', 10, 0, 8, 0.7),
                  ('A', 2, 'B', 8, 0, 8, 0.9),
                  ('A', 3, 'B', 12, 0, 8, 0.4)]:
            contact = Contact(c[1], c[3], c[6], distance_bound=(c[4], c[5]))
            contact.res1_chain = c[0]
            contact.res2_chain = c[2]
            contact_map.add(contact)
        contact_map.sequence = Sequence('1', 'HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSDHLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD')
        f_name = create_tmp_f()
        with open(f_name, 'w') as f_out:
            CaspParser().write(f_out, contact_file)
        content = [
            "PFRMAT RR",
            "MODEL  1",
            "HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSDHLEGSIGILLKKHEIVF",
            "DGCHDFGRTYIWQMSD",
            "A1    B9    0   8   0.700000",
            "A1    B10   0   8   0.700000",
            "A2    B8    0   8   0.900000",
            "A3    B12   0   8   0.400000",
            "ENDMDL",
            "END",
            "",
        ]
        content = os.linesep.join(content)
        with open(f_name, 'r') as f_in:
            data = "".join(f_in.readlines())
        self.assertEqual(content, data)
        os.unlink(f_name)


if __name__ == "__main__":
    unittest.main(verbosity=2)

