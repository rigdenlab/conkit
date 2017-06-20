"""Testing facility for conkit.io._cache"""

__author__ = "Felix Simkovic"
__date__ = "19 Jun 2017"

import unittest

from conkit.io._cache import _ParserCache


class Test_ParserCache(unittest.TestCase):

    def test_1(self):
        c = _ParserCache()
        self.assertTrue("casprr" in c)

    def test_2(self):
        c = _ParserCache()
        self.assertFalse("foo" in c)

    def test_3(self):
        c = _ParserCache()
        self.assertTrue("casprr" in c.contact_file_parsers)

    def test_4(self):
        c = _ParserCache()
        self.assertFalse("bar" in c.contact_file_parsers)

    def test_5(self):
        c = _ParserCache()
        self.assertTrue("fasta" in c.sequence_file_parsers)

    def test_6(self):
        c = _ParserCache()
        self.assertFalse("baz" in c.sequence_file_parsers)

    def test_7(self):
        c = _ParserCache()
        self.assertFalse("fasta" in c.contact_file_parsers)

    def test_8(self):
        c = _ParserCache()
        self.assertFalse("casprr" in c.sequence_file_parsers)


if __name__ == "__main__":
    unittest.main(verbosity=2)
