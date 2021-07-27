"""Testing facility for conkit.core.distance.Distance"""


import unittest
import numpy as np
from conkit.core.distance import Distance


class TestDistance(unittest.TestCase):

    def test_predicted_distance(self):
        distance = Distance(1, 25, (0.15, 0.45, 0.25, 0.05, 0.1), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf)))
        self.assertEqual(distance.max_score, 0.45)
        self.assertTupleEqual(distance.predicted_distance_bin, (4, 6))
        self.assertEqual(distance.predicted_distance, 5)

    def test_get_probability_within_distance(self):
        distance = Distance(1, 25, (0.15, 0.45, 0.25, 0.05, 0.1), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf)))
        self.assertEqual(distance.raw_score, 0.85)
        self.assertEqual(distance.get_probability_within_distance(5), 0.375)
        self.assertEqual(distance.get_probability_within_distance(8), 0.85)
        self.assertEqual(distance.get_probability_within_distance(10), 0.9)
        self.assertEqual(distance.get_probability_within_distance(25), 0.999999969409768)
        self.assertEqual(distance.get_probability_within_distance(np.inf), 1)
        with self.assertRaises(ValueError):
            distance.get_probability_within_distance(0)


if __name__ == "__main__":
    unittest.main(verbosity=2)
