"""Testing facility for conkit.plot._iotools"""

__author__ = "Felix Simkovic"
__date__ = "16 Feb 2017"

from conkit.core import *
from conkit.plot import tools

import unittest


class Test(unittest.TestCase):
    def test_find_minima_1(self):
        with self.assertRaises(ValueError):
            tools.find_minima([])

    def test_find_minima_2(self):
        with self.assertRaises(ValueError):
            tools.find_minima([5])

    def test_find_minima_3(self):
        self.assertEqual([1], tools.find_minima([5, 4, 5], order=1))

    def test_find_minima_4(self):
        self.assertEqual([2, 4], tools.find_minima([4, 5, 3, 5, 2, 6], order=1))

    def test_find_minima_5(self):
        self.assertEqual([4], tools.find_minima([4, 5, 3, 5, 2, 6], order=3))

    def test_find_minima_6(self):
        self.assertEqual([4], tools.find_minima([4, 5, 3, 5, 2, 6], order=10))

    def test_find_minima_7(self):
        with self.assertRaises(ValueError):
            tools.find_minima([4, 5, 3, 5, 2, 6], order=0)

    def test_find_minima_8(self):
        with self.assertRaises(ValueError):
            tools.find_minima([4, 5, 3, 5, 2, 6], order=-1)

    def test_get_points_on_circle_1(self):
        coords = tools.get_points_on_circle(0)
        self.assertEqual([[]], coords)

    def test_get_points_on_circle_2(self):
        coords = tools.get_points_on_circle(1)
        self.assertEqual([[1.0, 0.0]], coords)

    def test_get_points_on_circle_3(self):
        coords = tools.get_points_on_circle(2)
        self.assertEqual([2.0, 0.0], coords[0])
        self.assertEqual([-2.0, 0], coords[1])

    def test_get_points_on_circle_4(self):
        coords = tools.get_points_on_circle(3)
        self.assertEqual([3.0, 0.0], coords[0])
        self.assertEqual([-1.5, 2.598076], coords[1])
        self.assertEqual([-1.5, -2.598076], coords[2])

    def test_get_points_on_circle_5(self):
        coords = tools.get_points_on_circle(4)
        self.assertEqual([4.0, 0.0], coords[0])
        self.assertEqual([0.0, 4.0], coords[1])
        self.assertEqual([-4.0, 0.0], coords[2])
        self.assertEqual([0, -4], coords[3])

    def test_get_points_on_circle_6(self):
        coords = tools.get_points_on_circle(4, h=2)
        self.assertEqual([6.0, 0.0], coords[0])
        self.assertEqual([2.0, 4.0], coords[1])
        self.assertEqual([-2.0, 0.0], coords[2])
        self.assertEqual([2, -4], coords[3])

    def test_get_points_on_circle_7(self):
        coords = tools.get_points_on_circle(4, k=-3)
        self.assertEqual([4.0, -3.0], coords[0])
        self.assertEqual([0.0, 1.0], coords[1])
        self.assertEqual([-4.0, -3.0], coords[2])
        self.assertEqual([0, -7], coords[3])

    def test_get_radius_around_circle_1(self):
        p1, p2 = [0.0, 0.0], [0.0, 0]
        self.assertEqual(0.0, tools.get_radius_around_circle(p1, p2))

    def test_get_radius_around_circle_2(self):
        p1, p2 = [2.0, 0.0], [-2.0, 0]
        self.assertEqual(1.6, tools.get_radius_around_circle(p1, p2))

    def test_get_radius_around_circle_3(self):
        p1, p2 = [2.0, 1.0], [-2.0, 4]
        self.assertAlmostEqual(2.0, tools.get_radius_around_circle(p1, p2))

    def test__isinstance_1(self):
        h = Contact(1, 2, 2)
        self.assertTrue(tools._isinstance(h, "Contact"))

    def test__isinstance_2(self):
        h = ContactMap("test")
        self.assertTrue(tools._isinstance(h, "ContactMap"))

    def test__isinstance_3(self):
        h = ContactFile("test")
        self.assertTrue(tools._isinstance(h, "ContactFile"))

    def test__isinstance_4(self):
        h = Sequence("test", "AAAA")
        self.assertTrue(tools._isinstance(h, "Sequence"))

    def test__isinstance_5(self):
        h = SequenceFile("test")
        self.assertTrue(tools._isinstance(h, "SequenceFile"))

    def test__isinstance_6(self):
        self.assertTrue(tools._isinstance(5, int))

    def test__isinstance_7(self):
        self.assertFalse(tools._isinstance(5, str))


if __name__ == "__main__":
    unittest.main(verbosity=2)
