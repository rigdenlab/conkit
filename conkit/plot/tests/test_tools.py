"""Testing facility for conkit.plot._iotools"""

__author__ = "Felix Simkovic"
__date__ = "16 Feb 2017"

from conkit.core import *
from conkit.plot import tools

import unittest


class Test(unittest.TestCase):
    def test_points_on_circle_1(self):
        coords = tools.points_on_circle(2)
        self.assertEqual([2., 0.], coords[0])
        self.assertEqual([-2., 0], coords[1])

    def test_points_on_circle_2(self):
        coords = tools.points_on_circle(3)
        self.assertEqual([3., 0.], coords[0])
        self.assertEqual([-1.5, 2.598076], coords[1])
        self.assertEqual([-1.5, -2.598076], coords[2])

    def test_points_on_circle_3(self):
        coords = tools.points_on_circle(4)
        self.assertEqual([4., 0.], coords[0])
        self.assertEqual([0., 4.], coords[1])
        self.assertEqual([-4., 0.], coords[2])
        self.assertEqual([0, -4], coords[3])

    def test_points_on_circle_4(self):
        coords = tools.points_on_circle(4, h=2)
        self.assertEqual([6., 0.], coords[0])
        self.assertEqual([2., 4.], coords[1])
        self.assertEqual([-2., 0.], coords[2])
        self.assertEqual([2, -4], coords[3])

    def test_points_on_circle_5(self):
        coords = tools.points_on_circle(4, k=-3)
        self.assertEqual([4., -3.], coords[0])
        self.assertEqual([0., 1.], coords[1])
        self.assertEqual([-4., -3.], coords[2])
        self.assertEqual([0, -7], coords[3])

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


if __name__ == "__main__":
    unittest.main(verbosity=2)
