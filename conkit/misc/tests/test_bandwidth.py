"""Testing facility for conkit.misc.bandwidth"""

__author__ = "Felix Simkovic"
__date__ = "19 Jun 2017"

import numpy as np
import unittest

from conkit.misc import bandwidth


class TestAmiseBW(unittest.TestCase):
    def test_extended_range_1(self):
        min_, max_ = bandwidth.AmiseBW.extended_range(1, 5, 1.5, 3)
        self.assertEqual(min_, -3.5)
        self.assertEqual(max_, 9.5)

    def test_gauss_curvature_1(self):
        data = np.array([1, 2, 3, 4, 5, 3, 2, 3, 4], dtype=np.int64)[:, np.newaxis]
        curvature = bandwidth.AmiseBW.gauss_curvature(data, -1.5, 0.5)
        self.assertEqual(round(curvature, 7), 3.17e-5)

    def test_stiffness_integral_1(self):
        data = np.array([1, 2, 3, 4, 5, 3, 2, 3, 4], dtype=np.int64)[:, np.newaxis]
        integral = bandwidth.AmiseBW.stiffness_integral(data, 2.0, 1e-4)
        self.assertEqual(round(integral, 7), 0.0031009)

    def test_optimal_bandwidth_1(self):
        data = np.array([1, 2, 3, 4, 5, 3, 2, 3, 4], dtype=np.int64)[:, np.newaxis]
        optimal = bandwidth.AmiseBW.optimal_bandwidth(data, 2.0)
        self.assertEqual(round(optimal, 7), 0.4116948)

    def test_bandwidth_1(self):
        xy = np.array([(1, 5), (3, 3), (2, 4)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(bandwidth.AmiseBW(x).bandwidth, 7), 1.1455243)

    def test_bandwidth_2(self):
        xy = np.array([(1, 5), (3, 3), (2, 4), (1, 10), (4, 9)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(bandwidth.AmiseBW(x).bandwidth, 7), 1.5310027)

    def test_bandwidth_3(self):
        xy = np.array([(3, 5), (2, 4), (3, 4)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(bandwidth.AmiseBW(x).bandwidth, 7), 0.3758801)


class TestBowmanBW(unittest.TestCase):
    def test_bandwidth_1(self):
        xy = np.array([(1, 5), (3, 3), (2, 4)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(bandwidth.BowmanBW(x).bandwidth, 7), 0.7881495)

    def test_bandwidth_2(self):
        xy = np.array([(1, 5), (3, 3), (2, 4), (1, 10), (4, 9)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(bandwidth.BowmanBW(x).bandwidth, 7), 1.4223373)

    def test_bandwidth_3(self):
        xy = np.array([(3, 5), (2, 4), (3, 4)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(bandwidth.BowmanBW(x).bandwidth, 7), 0.6052020)


class TestLinearBW(unittest.TestCase):
    def test_bandwidth_1(self):
        xy = np.array([(1, 5), (3, 3), (2, 4)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(bandwidth.LinearBW(x, threshold=8).bandwidth, 0.625)

    def test_bandwidth_2(self):
        xy = np.array([(1, 5), (3, 3), (2, 4), (1, 10), (4, 9)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(bandwidth.LinearBW(x, threshold=10).bandwidth, 1.0)

    def test_bandwidth_3(self):
        xy = np.array([(3, 5), (2, 4), (3, 4)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(bandwidth.LinearBW(x, threshold=15).bandwidth, 7), 0.3333333)


class TestScottBW(unittest.TestCase):
    def test_bandwidth_1(self):
        xy = np.array([(1, 5), (3, 3), (2, 4)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(bandwidth.ScottBW(x).bandwidth, 7), 0.8357821)

    def test_bandwidth_2(self):
        xy = np.array([(1, 5), (3, 3), (2, 4), (1, 10), (4, 9)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(bandwidth.ScottBW(x).bandwidth, 7), 1.4513602)

    def test_bandwidth_3(self):
        xy = np.array([(3, 5), (2, 4), (3, 4)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(bandwidth.ScottBW(x).bandwidth, 7), 0.5179240)


class TestSilvermanBW(unittest.TestCase):
    def test_bandwidth_1(self):
        xy = np.array([(1, 5), (3, 3), (2, 4)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(bandwidth.SilvermanBW(x).bandwidth, 7), 0.7523629)

    def test_bandwidth_2(self):
        xy = np.array([(1, 5), (3, 3), (2, 4), (1, 10), (4, 9)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(bandwidth.SilvermanBW(x).bandwidth, 7), 1.3065002)

    def test_bandwidth_3(self):
        xy = np.array([(3, 5), (2, 4), (3, 4)])
        x = np.asarray([i for (x, y) in xy for i in np.arange(x, y + 1)])[:, np.newaxis]
        self.assertEqual(round(bandwidth.SilvermanBW(x).bandwidth, 7), 0.4662301)


if __name__ == "__main__":
    unittest.main(verbosity=2)
