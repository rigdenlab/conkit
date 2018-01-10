"""Testing facility for conkit.misc.__init__"""

__author__ = "Felix Simkovic"
__date__ = "10 Jan 2018"

import unittest

from conkit.misc import *


class Test(unittest.TestCase):
    def test_normalize_1(self):
        self.assertListEqual([0.0, 0.5, 1.0], normalize([1, 2, 3]))

    def test_normalize_2(self):
        self.assertListEqual([0.0, 0.5, 1.0], normalize([0.0, 0.5, 1.0]))

    def test_normalize_3(self):
        self.assertListEqual([0.0, 0.5, 1.0], normalize([-3, -2, -1]))

    def test_normalize_4(self):
        self.assertListEqual([0.0, 1.0], normalize([1, 2]))

    def test_normalize_5(self):
        self.assertListEqual([-1.0, 1.0], normalize([1, 2], vmin=-1))

    def test_normalize_6(self):
        self.assertListEqual([0.0, 2.0], normalize([1, 2], vmax=2))

    def test_normalize_7(self):
        self.assertListEqual([0.0, -1.0], normalize([1, 2], vmax=-1))

    def test_normalize_8(self):
        self.assertListEqual([0.2, 0.8], normalize([1, 2], vmin=0.2, vmax=0.8))

    def test_normalize_9(self):
        self.assertListEqual([0.2, 0.5, 0.8], normalize([1, 2, 3], vmin=0.2, vmax=0.8))


if __name__ == "__main__":
    unittest.main(verbosity=2)
