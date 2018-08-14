__author__ = "Felix Simkovic"
__date__ = "18 Oct 2016"

from conkit.core.contact import Contact
from conkit.core.contactfile import ContactFile
from conkit.core.contactmap import ContactMap
from conkit.core.sequence import Sequence
from conkit.core.sequencefile import SequenceFile
from conkit.io._parser import Parser

import unittest


class Test_Parser(unittest.TestCase):
    def test__reconstruct_1(self):
        hierarchy = Parser._reconstruct(Contact(1, 3, 1.0))
        self.assertTrue(isinstance(hierarchy, ContactFile))
        self.assertTrue(isinstance(hierarchy[0], ContactMap))
        self.assertTrue(isinstance(hierarchy[0][0], Contact))
        self.assertEqual((1, 3), (hierarchy[0][0].res1_seq, hierarchy[0][0].res2_seq))

    def test__reconstruct_2(self):
        hierarchy = Parser._reconstruct(ContactMap('test'))
        self.assertTrue(isinstance(hierarchy, ContactFile))
        self.assertTrue(isinstance(hierarchy[0], ContactMap))
        self.assertEqual('test', hierarchy[0].id)

    def test__reconstruct_3(self):
        hierarchy = Parser._reconstruct(ContactFile('test'))
        self.assertTrue(isinstance(hierarchy, ContactFile))
        self.assertEqual('test', hierarchy.id)

    def test__reconstruct_4(self):
        hierarchy = Parser._reconstruct(SequenceFile('test'))
        self.assertTrue(isinstance(hierarchy, SequenceFile))

    def test__reconstruct_5(self):
        hierarchy = Parser._reconstruct(Sequence('test', 'AAA'))
        self.assertTrue(isinstance(hierarchy, SequenceFile))


if __name__ == "__main__":
    unittest.main(verbosity=2)
