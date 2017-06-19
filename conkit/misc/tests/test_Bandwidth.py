"""Testing facility for conkit.core._Entity"""

__author__ = "Felix Simkovic"
__date__ = "19 Jun 2017"

import numpy as np
import unittest

from conkit.misc import Bandwidth


class TestAmiseBW(unittest.TestCase):

    def test_1(self):
        xy = np.array([(1, 5), (3, 3), (2, 4)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(Bandwidth.AmiseBW(x).bw, 7), 1.1455243)

    def test_2(self):
        xy = np.array([(1, 5), (3, 3), (2, 4), (1, 10), (4, 9)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(Bandwidth.AmiseBW(x).bw, 7), 1.5310027)

    def test_3(self):
        xy = np.array([(3, 5), (2, 4), (3, 4)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(Bandwidth.AmiseBW(x).bw, 7), 0.3758801)


class TestBowmanBW(unittest.TestCase):

    def test_1(self):
        xy = np.array([(1, 5), (3, 3), (2, 4)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(Bandwidth.BowmanBW(x).bw, 7), 0.7881495)

    def test_2(self):
        xy = np.array([(1, 5), (3, 3), (2, 4), (1, 10), (4, 9)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(Bandwidth.BowmanBW(x).bw, 7), 1.4223373)

    def test_3(self):
        xy = np.array([(3, 5), (2, 4), (3, 4)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(Bandwidth.BowmanBW(x).bw, 7), 0.6052020)


class TestScottBW(unittest.TestCase):

    def test_1(self):
        xy = np.array([(1, 5), (3, 3), (2, 4)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(Bandwidth.ScottBW(x).bw, 7), 0.8357821)

    def test_2(self):
        xy = np.array([(1, 5), (3, 3), (2, 4), (1, 10), (4, 9)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(Bandwidth.ScottBW(x).bw, 7), 1.4513602)

    def test_3(self):
        xy = np.array([(3, 5), (2, 4), (3, 4)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(Bandwidth.ScottBW(x).bw, 7), 0.5179240)


class TestSilvermanBW(unittest.TestCase):

    def test_1(self):
        xy = np.array([(1, 5), (3, 3), (2, 4)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(Bandwidth.SilvermanBW(x).bw, 7), 0.7523629)

    def test_2(self):
        xy = np.array([(1, 5), (3, 3), (2, 4), (1, 10), (4, 9)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(Bandwidth.SilvermanBW(x).bw, 7), 1.3065002)

    def test_3(self):
        xy = np.array([(3, 5), (2, 4), (3, 4)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(Bandwidth.SilvermanBW(x).bw, 7), 0.4662301)


if __name__ == "__main__":
    unittest.main(verbosity=2)
