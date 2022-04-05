"""Testing facility for conkit.core.distance.Distance"""

import unittest
import numpy as np
from conkit.core.distance import Distance
from conkit.core.distogram import Distogram


class TestDistance(unittest.TestCase):

    def test_predicted_distance_1(self):
        distance = Distance(1, 25, (0.15, 0.45, 0.25, 0.05, 0.1), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf)))
        self.assertEqual(distance.max_score, 0.45)
        self.assertTupleEqual(distance.predicted_distance_bin, (4, 6))
        self.assertEqual(distance.predicted_distance, 5)

    def test_predicted_distance_2(self):
        distogram = Distogram('test')
        distogram.original_file_format = 'pdb'
        distance = Distance(36, 86, (1,), ((6.589181, 6.589181),), 0.934108)
        distogram.add(distance)
        self.assertEqual(6.589181, distance.predicted_distance)

    def test_predicted_distance_3(self):
        distance = Distance(2, 3, (0.2, 0.3, 0.3, 0.2), ((0, 4), (4, 6), (6, 8), (8, np.inf)))
        self.assertEqual(distance.max_score, 0.3)
        self.assertTupleEqual(distance.predicted_distance_bin, (4, 6))
        self.assertEqual(distance.predicted_distance, 5)

    def test_get_probability_within_distance_1(self):
        distance = Distance(1, 25, (0.15, 0.45, 0.25, 0.05, 0.1), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf)))
        self.assertEqual(distance.raw_score, 0.85)
        self.assertEqual(distance.get_probability_within_distance(5), 0.375)
        self.assertEqual(distance.get_probability_within_distance(8), 0.85)
        self.assertEqual(distance.get_probability_within_distance(10), 0.9)
        self.assertEqual(distance.get_probability_within_distance(25), 0.999999969409768)
        self.assertEqual(distance.get_probability_within_distance(np.inf), 1)
        self.assertEqual(distance.get_probability_within_distance(0), 0)
        with self.assertRaises(ValueError):
            distance.get_probability_within_distance(-5)

    def test_get_probability_within_distance_2(self):
        distogram = Distogram('test')
        distogram.original_file_format = 'pdb'
        distance = Distance(36, 86, (1,), ((6.589181, 6.589181),), 0.934108)
        distogram.add(distance)
        self.assertEqual(1, distance.get_probability_within_distance(8))
        self.assertEqual(0, distance.get_probability_within_distance(5))

    def test_reshape_bins_1(self):
        distance = Distance(1, 25, (0.15, 0.45, 0.25, 0.05, 0.1), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf)))
        new_bins = ((0, 2), (2, 8), (8, np.inf))
        distance.reshape_bins(new_bins)
        self.assertEqual(distance.raw_score, 0.85)
        self.assertEqual(round(distance.get_probability_within_distance(8), 2), 0.85)
        self.assertTupleEqual(new_bins, distance.distance_bins)
        self.assertTupleEqual((0.075, 0.775, 0.15000000000000002), distance.distance_scores)

    def test_reshape_bins_2(self):
        distogram = Distogram('test')
        distogram.original_file_format = 'pdb'
        distance = Distance(36, 86, (1,), ((6.589181, 6.589181),), 0.934108)
        distogram.add(distance)
        with self.assertRaises(ValueError):
            distance.reshape_bins(((0, 1), (1, 10), (10, np.inf)))

    def test__assert_valid_bins_1(self):
        distance = Distance(1, 25, (0.15, 0.45, 0.25, 0.05, 0.1), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf)))
        with self.assertRaises(ValueError):
            distance._assert_valid_bins(((0, 1), (np.inf, 10), (10, np.inf)))

    def test__assert_valid_bins_2(self):
        distance = Distance(1, 25, (0.15, 0.45, 0.25, 0.05, 0.1), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf)))
        with self.assertRaises(ValueError):
            distance._assert_valid_bins(((0, 1), (1, 10), (10, 20)))

    def test__assert_valid_bins_3(self):
        distance = Distance(1, 25, (0.15, 0.45, 0.25, 0.05, 0.1), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf)))
        with self.assertRaises(ValueError):
            distance._assert_valid_bins(((0, 1), (5, 10, 15), (15, np.inf)))

    def test__assert_valid_bins_4(self):
        distance = Distance(1, 25, (0.15, 0.45, 0.25, 0.05, 0.1), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf)))
        with self.assertRaises(ValueError):
            distance._assert_valid_bins(((0, 1), (5, 10), (10, np.inf)))

    def test__assert_valid_bins_5(self):
        distance = Distance(1, 25, (0.15, 0.45, 0.25, 0.05, 0.1), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf)))
        with self.assertRaises(ValueError):
            distance._assert_valid_bins(((0, 1), (1, 10), (10, 20), (25, 30), (30, 31), (45, np.inf)))

    def test__assert_valid_bins_6(self):
        distance = Distance(1, 25, (0.15, 0.45, 0.25, 0.05, 0.1), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf)))
        with self.assertRaises(ValueError):
            distance._assert_valid_bins(((0, 1), (10, 10), (10, np.inf)))

    def test__assert_valid_bins_7(self):
        distance = Distance(1, 25, (0.15, 0.45, 0.25, 0.05, 0.1), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf)))
        with self.assertRaises(ValueError):
            distance._assert_valid_bins(((0, np.inf),))


if __name__ == "__main__":
    unittest.main(verbosity=2)
