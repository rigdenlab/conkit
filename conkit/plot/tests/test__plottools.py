"""Testing facility for conkit.plot._iotools"""

__author__ = "Felix Simkovic"
__date__ = "16 Feb 2017"

from conkit.plot import _plottools

import unittest


class Test(unittest.TestCase):

    def test_points_on_circle1(self):
        coords = _plottools.points_on_circle(2)
        self.assertEqual([2., 0.], coords[0])
        self.assertEqual([-2., 0], coords[1])

    def test_points_on_circle2(self):
        coords = _plottools.points_on_circle(3)
        self.assertEqual([3., 0.], coords[0])
        self.assertEqual([-1.5, 2.598076], coords[1]) 
        self.assertEqual([-1.5, -2.598076], coords[2]) 

    def test_points_on_circle3(self):
        coords = _plottools.points_on_circle(4)
        self.assertEqual([4., 0.], coords[0])
        self.assertEqual([0., 4.], coords[1]) 
        self.assertEqual([-4., 0.], coords[2]) 
        self.assertEqual([0, -4], coords[3])

    def test_points_on_circle4(self):
        coords = _plottools.points_on_circle(4, h=2)
        self.assertEqual([6., 0.], coords[0])
        self.assertEqual([2., 4.], coords[1]) 
        self.assertEqual([-2., 0.], coords[2]) 
        self.assertEqual([2, -4], coords[3])

    def test_points_on_circle5(self):
        coords = _plottools.points_on_circle(4, k=-3)
        self.assertEqual([4., -3.], coords[0])
        self.assertEqual([0., 1.], coords[1]) 
        self.assertEqual([-4., -3.], coords[2]) 
        self.assertEqual([0, -7], coords[3])


if __name__ == "__main__":
    unittest.main(verbosity=2)
