"""Testing facility for conkit.core.Distogram"""

import unittest
import numpy as np
from conkit.core.distance import Distance
from conkit.core.distogram import Distogram
from conkit.core.sequence import Sequence


class TestDistogram(unittest.TestCase):
    def test_ndistances_1(self):
        distogram = Distogram("test")
        distogram.add(Distance(1, 25, (0.25, 0.45, 0.25, 0.05), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(7, 19, (0.15, 0.15, 0.60, 0.1), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        self.assertEqual(2, distogram.ndistances)

    def test_get_unique_distances_1(self):
        distogram = Distogram("test")
        print(Distance(1, 25, (0.25, 0.45, 0.25, 0.05), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(1, 25, (0.25, 0.45, 0.25, 0.05), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(25, 1, (0.25, 0.45, 0.25, 0.05), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(7, 19, (0.15, 0.15, 0.60, 0.1), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(19, 7, (0.15, 0.15, 0.60, 0.1), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(1, 7, (0.1, 0.2, 0.55, 0.15), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.get_unique_distances(inplace=True)
        self.assertListEqual([[25, 1], [19, 7], [1, 7]], distogram.as_list())

    def test_get_unique_distances_2(self):
        distogram = Distogram("test")
        distogram.add(Distance(1, 25, (0.25, 0.45, 0.25, 0.05), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(25, 1, (0.25, 0.45, 0.25, 0.05), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(7, 19, (0.15, 0.15, 0.60, 0.1), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(19, 7, (0.15, 0.15, 0.60, 0.1), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(1, 7, (0.1, 0.2, 0.55, 0.15), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        new_distogram = distogram.get_unique_distances(inplace=False)
        self.assertListEqual([[1, 25], [25, 1], [7, 19], [19, 7], [1, 7]], distogram.as_list())
        self.assertListEqual([[25, 1], [19, 7], [1, 7]], new_distogram.as_list())

    def test_get_absent_residues_1(self):
        distogram = Distogram("test")
        distogram.add(Distance(1, 5, (0.25, 0.45, 0.25, 0.05), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(2, 3, (0.15, 0.15, 0.60, 0.1), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        self.assertListEqual([4, 6, 7], distogram.get_absent_residues(7))

    def test_get_absent_residues_2(self):
        distogram = Distogram("test")
        distogram.add(Distance(1, 5, (0.25, 0.45, 0.25, 0.05), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(2, 3, (0.15, 0.15, 0.60, 0.1), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.sequence = Sequence('test', 'AAAAAAA')
        self.assertListEqual([4, 6, 7], distogram.get_absent_residues())

    def test_as_array_1(self):
        distogram = Distogram("test")
        distogram.add(Distance(1, 5, (0.25, 0.45, 0.25, 0.05), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(2, 3, (0.15, 0.15, 0.60, 0.1), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(1, 4, (0.05, 0.25, 0.70, 0.0), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(3, 5, (0.5, 0.1, 0.35, 0.05), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        output = np.nan_to_num(distogram.as_array(seq_len=5)).tolist()

        expected = [
            [0., 0., 0., 7.0, 5.0],
            [0., 0., 7.0, 0., 0.],
            [0., 7.0, 0., 0., 2.0],
            [7.0, 0., 0., 0., 0.],
            [5.0, 0., 2.0, 0., 0.]
        ]

        self.assertListEqual(output, expected)

    def test_as_array_2(self):
        distogram = Distogram("test")
        distogram.add(Distance(1, 5, (0.25, 0.45, 0.25, 0.05), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(2, 3, (0.15, 0.15, 0.60, 0.1), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(1, 4, (0.05, 0.25, 0.70, 0.0), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(3, 5, (0.5, 0.1, 0.35, 0.05), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.sequence = Sequence('test', 'AAAAA')
        output = np.nan_to_num(distogram.as_array()).tolist()

        expected = [
            [0., 0., 0., 7.0, 5.0],
            [0., 0., 7.0, 0., 0.],
            [0., 7.0, 0., 0., 2.0],
            [7.0, 0., 0., 0., 0.],
            [5.0, 0., 2.0, 0., 0.]
        ]

        self.assertListEqual(output, expected)

    def test_as_array_3(self):
        distogram = Distogram("test")
        distogram.add(Distance(1, 5, (0.25, 0.45, 0.25, 0.05), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(2, 3, (0.15, 0.15, 0.60, 0.1), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(1, 4, (0.05, 0.25, 0.70, 0.0), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram.add(Distance(3, 5, (0.5, 0.1, 0.35, 0.05), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        output = np.nan_to_num(distogram.as_array(seq_len=5, get_weigths=True)).tolist()

        expected = [
            [0., 0., 0., 0.7, 0.45],
            [0., 0., 0.6, 0., 0.],
            [0., 0.6, 0., 0., 0.5],
            [0.7, 0., 0., 0., 0.],
            [0.45, 0., 0.5, 0., 0.]
        ]

        self.assertListEqual(output, expected)
