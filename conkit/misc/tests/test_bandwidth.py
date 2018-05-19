"""Testing facility for conkit.misc.bandwidth"""

__author__ = "Felix Simkovic"
__date__ = "19 Jun 2017"

import numpy as np
import unittest

from conkit.misc import bandwidth
from conkit.misc.ext import c_bandwidth


class TestAmiseBW(unittest.TestCase):
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


class Test(unittest.TestCase):
    def test_bandwidth_factory_1(self):
        obj = bandwidth.bandwidth_factory("amise")
        self.assertEqual(str(obj), "<class 'conkit.misc.bandwidth.AmiseBW'>")

    def test_bandwidth_factory_2(self):
        obj = bandwidth.bandwidth_factory("bowman")
        self.assertEqual(str(obj), "<class 'conkit.misc.bandwidth.BowmanBW'>")

    def test_bandwidth_factory_3(self):
        obj = bandwidth.bandwidth_factory("linear")
        self.assertEqual(str(obj), "<class 'conkit.misc.bandwidth.LinearBW'>")

    def test_bandwidth_factory_4(self):
        obj = bandwidth.bandwidth_factory("scott")
        self.assertEqual(str(obj), "<class 'conkit.misc.bandwidth.ScottBW'>")

    def test_bandwidth_factory_5(self):
        obj = bandwidth.bandwidth_factory("silverman")
        self.assertEqual(str(obj), "<class 'conkit.misc.bandwidth.SilvermanBW'>")

    def test_bandwidth_factory_6(self):
        with self.assertRaises(ValueError):
            bandwidth.bandwidth_factory("SILVERMAN")

    def test_bandwidth_factory_7(self):
        with self.assertRaises(ValueError):
            bandwidth.bandwidth_factory("Silverman")

    def test_bandwidth_factory_8(self):
        with self.assertRaises(ValueError):
            bandwidth.bandwidth_factory("silvermn")

    def test_bandwidth_factory_9(self):
        with self.assertRaises(ValueError):
            bandwidth.bandwidth_factory("garbage")


class TestExt(unittest.TestCase):
    def test_gauss_curvature_1(self):
        A = np.array(
            [[1], [2], [3], [4], [5], [3], [2], [3], [4]], dtype=np.int64)
        curvature = c_bandwidth.c_get_gauss_curvature(A, -1.5, 0.5)
        self.assertAlmostEqual(3.171746247735917e-05, curvature)

    def test_gauss_curvature_2(self):
        A = np.array([[1]], dtype=np.int64)
        curvature = c_bandwidth.c_get_gauss_curvature(A, -1.5, 0.5)
        self.assertAlmostEqual(0.0002854501468289852, curvature)

    def test_gauss_curvature_3(self):
        A = np.array([[1]], dtype=np.int64)
        curvature = c_bandwidth.c_get_gauss_curvature(A, 0, 0.5)
        self.assertAlmostEqual(1.2957831963165134, curvature)

    def test_stiffness_integral_1(self):
        A = np.array([[1], [2], [3], [4], [5], [3], [2], [3], [4]], dtype=np.int64)
        stiff_integ = c_bandwidth.c_get_stiffness_integral(A, 2.0, 0.0001)
        self.assertAlmostEqual(0.003100864697366348, stiff_integ)

    def test_stiffness_integral_2(self):
        A = np.array([[1], [2], [3], [4], [5], [3], [2], [3], [4]], dtype=np.int64)
        stiff_integ = c_bandwidth.c_get_stiffness_integral(A, 2.0, 0.1)
        self.assertAlmostEqual(0.003100864697366348, stiff_integ)

    def test_stiffness_integral_3(self):
        A = np.array([[1], [2], [3], [4], [5], [3], [2], [3], [4]], dtype=np.int64)
        stiff_integ = c_bandwidth.c_get_stiffness_integral(A, 1000.0, 0.0001)
        self.assertAlmostEqual(2.1106164693083826e-16, stiff_integ)

    def test_optimize_bandwidth_1(self):
        A = np.array([[1], [2], [3], [4], [5], [3], [2], [3], [4]], dtype=np.int64)
        optimized = c_bandwidth.c_optimize_bandwidth(A, 2.0)
        self.assertAlmostEqual(0.4116948343202962, optimized)

    def test_optimize_bandwidth_2(self):
        A = np.array([[1], [2], [3], [4], [5], [3], [2], [3], [4]], dtype=np.int64)
        optimized = c_bandwidth.c_optimize_bandwidth(A, 1000.0)
        self.assertAlmostEqual(317.11331138268406, optimized)


if __name__ == "__main__":
    unittest.main(verbosity=2)
