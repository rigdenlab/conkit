"""Testing facility for conkit.io._ParserIO"""

__author__ = "Felix Simkovic"
__date__ = "18 Oct 2016"

from conkit.core import Contact
from conkit.core import ContactFile
from conkit.core import ContactMap
from conkit.core import Sequence
from conkit.core import SequenceFile
from conkit.io._ParserIO import _ContactFileParser, _SequenceFileParser

import unittest


class Test1(unittest.TestCase):
    """Test for _ContactFileParser"""
    def test__reconstruct(self):
        parser = _ContactFileParser()
        # ======================================================
        # Test Case 1
        hierarchy = parser._reconstruct(Contact(1, 3, 1.0))
        self.assertTrue(isinstance(hierarchy, ContactFile))
        self.assertTrue(isinstance(hierarchy[0], ContactMap))
        self.assertTrue(isinstance(hierarchy[0][0], Contact))
        self.assertEqual((1, 3), (hierarchy[0][0].res1_seq, hierarchy[0][0].res2_seq))
        # ======================================================
        # Test Case 2
        hierarchy = parser._reconstruct(ContactMap('test'))
        self.assertTrue(isinstance(hierarchy, ContactFile))
        self.assertTrue(isinstance(hierarchy[0], ContactMap))
        self.assertEqual('test', hierarchy[0].id)
        # ======================================================
        # Test Case 3
        hierarchy = parser._reconstruct(ContactFile('test'))
        self.assertTrue(isinstance(hierarchy, ContactFile))
        self.assertEqual('test', hierarchy.id)


class Test2(unittest.TestCase):
    """Test for _SequenceFileParser"""
    def test__reconstruct(self):
        parser = _SequenceFileParser()
        # ======================================================
        # Test Case 1
        hierarchy = parser._reconstruct(SequenceFile('test'))
        self.assertTrue(isinstance(hierarchy, SequenceFile))
        # ======================================================
        # Test Case 2
        hierarchy = parser._reconstruct(Sequence('test', 'AAA'))
        self.assertTrue(isinstance(hierarchy, SequenceFile))

if __name__ == "__main__":
    unittest.main(verbosity=2)
