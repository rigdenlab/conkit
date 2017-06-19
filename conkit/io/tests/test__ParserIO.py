"""Testing facility for conkit.io._ParserIO"""

__author__ = "Felix Simkovic"
__date__ = "18 Oct 2016"

from conkit.core.ContactCore import Contact
from conkit.core.ContactFileCore import ContactFile
from conkit.core.ContactMapCore import ContactMap
from conkit.core.SequenceCore import Sequence
from conkit.core.SequenceFileCore import SequenceFile
from conkit.io._ParserIO import _ContactFileParser
from conkit.io._ParserIO import _SequenceFileParser

import unittest


class Test_ContactFileParser(unittest.TestCase):
    def test__reconstruct_1(self):
        hierarchy = _ContactFileParser()._reconstruct(Contact(1, 3, 1.0))
        self.assertTrue(isinstance(hierarchy, ContactFile))
        self.assertTrue(isinstance(hierarchy[0], ContactMap))
        self.assertTrue(isinstance(hierarchy[0][0], Contact))
        self.assertEqual((1, 3), (hierarchy[0][0].res1_seq, hierarchy[0][0].res2_seq))

    def test__reconstruct_2(self):
        hierarchy = _ContactFileParser()._reconstruct(ContactMap('test'))
        self.assertTrue(isinstance(hierarchy, ContactFile))
        self.assertTrue(isinstance(hierarchy[0], ContactMap))
        self.assertEqual('test', hierarchy[0].id)

    def test__reconstruct_3(self):
        hierarchy = _ContactFileParser()._reconstruct(ContactFile('test'))
        self.assertTrue(isinstance(hierarchy, ContactFile))
        self.assertEqual('test', hierarchy.id)


class Test_SequenceFileParser(unittest.TestCase):
    def test__reconstruct_1(self):
        hierarchy = _SequenceFileParser()._reconstruct(SequenceFile('test'))
        self.assertTrue(isinstance(hierarchy, SequenceFile))

    def test__reconstruct_2(self):
        hierarchy = _SequenceFileParser()._reconstruct(Sequence('test', 'AAA'))
        self.assertTrue(isinstance(hierarchy, SequenceFile))

if __name__ == "__main__":
    unittest.main(verbosity=2)
