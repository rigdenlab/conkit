__author__ = "Felix Simkovic"

import unittest

from conkit.misc.selectalg import SubselectionAlgorithm


class TestSubselectionAlgorithm(unittest.TestCase):
    def test_cutoff_1(self):
        data = [1.0, 0.6, 0.5, 0.4, 0.4, 0.3, 0.2, 0.1]
        keep, throw = SubselectionAlgorithm.cutoff(data)
        self.assertEqual([0, 1, 2, 3, 4, 5], keep)
        self.assertEqual([6, 7], throw)

    def test_cutoff_2(self):
        data = [1.0, 0.3, 0.2, 0.1, 0.6, 0.5, 0.4, 0.4]
        keep, throw = SubselectionAlgorithm.cutoff(data)
        self.assertEqual([0, 1, 4, 5, 6, 7], keep)
        self.assertEqual([2, 3], throw)

    def test_cutoff_3(self):
        data = [0.2, 0.1]
        keep, throw = SubselectionAlgorithm.cutoff(data)
        self.assertEqual([], keep)
        self.assertEqual([0, 1], throw)

    def test_cutoff_4(self):
        data = [0.286, 0.287, 0.288]
        keep, throw = SubselectionAlgorithm.cutoff(data)
        self.assertEqual([1, 2], keep)
        self.assertEqual([0], throw)

    def test_linear_1(self):
        data = [1.0, 0.6, 0.5, 0.45, 0.4, 0.3, 0.2, 0.1]
        keep, throw = SubselectionAlgorithm.linear(data)
        self.assertEqual([0, 1, 2, 3], keep)
        self.assertEqual([4, 5, 6, 7], throw)

    def test_linear_2(self):
        data = [0.1, 0.2, 0.3, 0.45, 0.4, 0.5, 0.6, 1.0]
        keep, throw = SubselectionAlgorithm.linear(data)
        self.assertEqual([7, 6, 5, 3], keep)
        self.assertEqual([4, 2, 1, 0], throw)

    def test_linear_3(self):
        data = [1.0, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1]
        keep, throw = SubselectionAlgorithm.linear(data)
        self.assertEqual([0, 1, 2, 3], keep)
        self.assertEqual([4, 5, 6], throw)

    def test_linear_4(self):
        data = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 1.0]
        keep, throw = SubselectionAlgorithm.linear(data)
        self.assertEqual([6, 5, 4, 3], keep)
        self.assertEqual([2, 1, 0], throw)

    def test_linear_5(self):
        data = [1.0, 0.21, 0.4, 0.6, 0.5, 0.3, 0.1, 0.2]
        keep, throw = SubselectionAlgorithm.linear(data)
        self.assertEqual([0, 3, 4, 2], keep)
        self.assertEqual([5, 1, 7, 6], throw)

    def test_scaled_1(self):
        data = [1.0, 0.6, 0.5, 0.45, 0.4, 0.3, 0.2, 0.1]
        keep, throw = SubselectionAlgorithm.scaled(data)
        self.assertEqual([0, 1, 2, 3, 4, 5], keep)
        self.assertEqual([6, 7], throw)

    def test_scaled_2(self):
        data = [1.0, 1.0, 1.0, 1.0]
        keep, throw = SubselectionAlgorithm.scaled(data)
        self.assertEqual([0, 1, 2, 3], keep)
        self.assertEqual([], throw)

    def test_scaled_3(self):
        data = [100.0, 1.0, 1.0, 1.0]
        keep, throw = SubselectionAlgorithm.scaled(data)
        self.assertEqual([0], keep)
        self.assertEqual([1, 2, 3], throw)

    def test_ignore_1(self):
        data = [1.0, 0.6, 0.5, 0.45, 0.4, 0.3, 0.2, 0.1]
        keep, throw = SubselectionAlgorithm.ignore(data)
        self.assertEqual([0, 1, 2, 3, 4, 5, 6, 7], keep)
        self.assertEqual([], throw)

    def test_ignore_2(self):
        data = [1.0, 1.0, 1.0, 1.0]
        keep, throw = SubselectionAlgorithm.ignore(data)
        self.assertEqual([0, 1, 2, 3], keep)
        self.assertEqual([], throw)

    def test_ignore_3(self):
        data = [100.0, 1.0, 1.0, 1.0]
        keep, throw = SubselectionAlgorithm.ignore(data)
        self.assertEqual([0, 1, 2, 3], keep)
        self.assertEqual([], throw)


if __name__ == "__main__":
    unittest.main(verbosity=2)
