"""Testing facility for conkit.plot._iotools"""

__author__ = "Felix Simkovic"
__date__ = "16 Feb 2017"

import os.path
import sys

from conkit.core import *
from conkit.plot import tools

import numpy as np
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

    def test_convolution_smooth_values(self):
        output = tools.convolution_smooth_values([3, 10, 8, 8, 10, 5]).tolist()
        expected = [4.2, 5.8, 7.8, 8.2, 6.2, 4.6]
        self.assertListEqual(expected, [round(x, 2) for x in output])

    def test_parse_map_align_stdout_1(self):
        stdout_contents = """OPT -------------------------------------------------------------------
    OPT                           MAP_ALIGN                                
    OPT -------------------------------------------------------------------
    OPT   -a          /home/filo/opt/map_align_v1/map_align/3u97_A.gremlin.map
    OPT   -b          /home/filo/opt/map_align_v1/map_align/2pd0_A.pdb.map
    OPT   -gap_o      -1
    OPT   -gap_e      -0.01
    OPT   -sep_cut    3
    OPT   -iter       20
    OPT   -silent     0
    OPT -------------------------------------------------------------------
    OPT   -use_gap_ss  0
    OPT   -use_prf     0
    OPT -------------------------------------------------------------------
    TMP	0_1_0	21.7832	-1.06	20.7233
    TMP	0_1_1	21.7832	-1.06	20.7233
    TMP	0_1_2	24.8138	-2.525	22.2888
    TMP	2_32_2	31.8905	-4.39	27.5005
    TMP	2_32_3	31.8905	-4.39	27.5005
    MAX 2_2_2	/home/filo/opt/map_align_v1/map_align/3u97_A.gremlin.map	/home/filo/opt/map_align_v1/map_align/2pd0_A.pdb.map	44.5272	-4.38	40.1472	74	1:1 2:2 3:3 4:4 5:6 6:7 7:8 8:9 10:10 11:11 12:13 13:14 14:15 16:16 17:17 18:18 19:19 20:20 21:23 22:24 23:25 24:26 25:27
"""

        expected = {6: 5, 7: 6, 8: 7, 9: 8, 13: 12, 14: 13, 15: 14, 23: 21, 24: 22, 25: 23, 26: 24, 27: 25}
        output = tools.parse_map_align_stdout(stdout_contents)
        self.assertDictEqual(output, expected)

    def test_get_residue_ranges_1(self):
        residues = [0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 19, 25, 26, 30, 31, 32, 33, 34, 35]
        expected = [(0, 12), (19, 19), (25, 26), (30, 35)]

        output = tools.get_residue_ranges(residues)
        self.assertListEqual(output, expected)

    def test_calculate_zscore_1(self):
        observed_score = 1
        population = [0, 0, 5, 0, 0, 2, 3, 4, 5, 2, 0, 0, 1, 9, 3, 0, 0, 5]
        expected = -0.46

        output = tools.calculate_zscore(observed_score, population)
        self.assertEqual(round(output, 2), expected)

    def test_calculate_zscore_2(self):
        observed_score = 1
        population = [1, ]
        expected = 0

        output = tools.calculate_zscore(observed_score, population)
        self.assertEqual(output, expected)

    def test_calculate_zscore_3(self):
        observed_score = 1
        population = [0, 0, 0, 0, 0, 0, 0]
        expected = 0

        output = tools.calculate_zscore(observed_score, population)
        self.assertEqual(output, expected)

    def test_get_rmsd_1(self):
        distogram_1 = Distogram("test_1")
        distogram_1.add(Distance(1, 5, (0.25, 0.45, 0.05, 0.05, 0.2), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf))))
        distogram_1.add(Distance(2, 3, (0.15, 0.15, 0.60, 0.1, 0.0), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf))))
        distogram_1.add(Distance(1, 4, (0.05, 0.2, 0.0, 0.6, 0.15), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf))))
        distogram_1.add(Distance(3, 5, (0.4, 0.1, 0.35, 0.05, 0.1), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf))))
        distogram_1.sequence = Sequence("test_seq", "AAAAA")

        distogram_2 = Distogram("test_2")
        distogram_2.add(Distance(1, 5, (0.45, 0.05, 0.25, 0.25), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram_2.add(Distance(2, 3, (0.1, 0.15, 0.15, 0.6), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram_2.add(Distance(1, 4, (0.75, 0.20, 0.05, 0.0), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram_2.add(Distance(3, 5, (0.05, 0.1, 0.35, 0.5), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram_2.sequence = Sequence("test_seq", "AAAAA")

        output = tools.get_rmsd(distogram_1, distogram_2, calculate_wrmsd=False)
        expected_0 = [5.39, 3.0, 6.04, 7.0, 6.04]
        expected_1 = [1.44, 2.14, 2.75, 2.75, 2.75, 2.75, 2.75, 2.75, 2.21, 1.91]
        self.assertListEqual(expected_0, [round(x, 2) for x in output[0]])
        self.assertListEqual(expected_1, [round(x, 2) for x in output[1]])

    def test_get_rmsd_2(self):
        distogram_1 = Distogram("test_1")
        distogram_1.add(Distance(1, 5, (0.25, 0.45, 0.05, 0.05, 0.2), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf))))
        distogram_1.add(Distance(2, 3, (0.15, 0.15, 0.60, 0.1, 0.0), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf))))
        distogram_1.add(Distance(1, 4, (0.05, 0.2, 0.0, 0.6, 0.15), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf))))
        distogram_1.add(Distance(3, 5, (0.4, 0.1, 0.35, 0.05, 0.1), ((0, 4), (4, 6), (6, 8), (8, 10), (10, np.inf))))
        distogram_1.sequence = Sequence("test_seq", "AAAAA")

        distogram_2 = Distogram("test_2")
        distogram_2.add(Distance(1, 5, (0.45, 0.05, 0.25, 0.25), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram_2.add(Distance(2, 3, (0.1, 0.15, 0.15, 0.6), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram_2.add(Distance(1, 4, (0.75, 0.20, 0.05, 0.0), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram_2.add(Distance(3, 5, (0.05, 0.1, 0.35, 0.5), ((0, 4), (4, 6), (6, 8), (8, np.inf))))
        distogram_2.sequence = Sequence("test_seq", "AAAAA")

        output = tools.get_rmsd(distogram_1, distogram_2, calculate_wrmsd=True)
        expected_0 = [4.09, 2.32, 3.94, 5.42, 3.85]
        expected_1 = [1.04, 1.58, 1.96, 1.96, 1.96, 1.96, 1.96, 1.96, 1.55, 1.32]
        self.assertListEqual(expected_0, [round(x, 2) for x in output[0]])
        self.assertListEqual(expected_1, [round(x, 2) for x in output[1]])

    def test_get_cmap_validation_metrics_1(self):
        cmap_1 = ContactMap("cmap_1")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2), Contact(1, 1, 0)]:
            cmap_1.add(c)

        cmap_2 = ContactMap("cmap_2")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(3, 4, 0.1), Contact(2, 1, 0.2), Contact(1, 1, 0)]:
            cmap_2.add(c)

        cmap_1.add(Contact(6, 1, 0))
        sequence = Sequence('seq', 'AAAAAAA')
        absent_residues = {6, 7}

        output = tools.get_cmap_validation_metrics(cmap_1.as_dict(), cmap_2.as_dict(), sequence, absent_residues)
        expected_accuracy = [0.6, 0.6, 0.8, 0.6, 0.8]
        expected_fn = [1, 1, 1, 1, 0]
        expected_fnr = [0.5, 0.25, 0.25, 0.25, 0.0]
        expected_fp = [1, 1, 0, 1, 1]
        expected_fpr = [0.33, 1.0, 0.0, 1.0, 0.5]
        expected_sensitivity = [0.67, 0.0, 0.5, 0.0, 1.0]
        expected_specificity = [0.5, 0.75, 1.0, 0.75, 0.75]
        expected = [expected_accuracy, expected_fn, expected_fnr, expected_fp,
                    expected_fpr, expected_sensitivity, expected_specificity]
        self.assertListEqual([round(x, 2) for l in expected for x in l], [round(x, 2) for l in output[0] for x in l])

    def test_get_cmap_validation_metrics_2(self):
        cmap_1 = ContactMap("cmap_1")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2), Contact(1, 1, 0)]:
            cmap_1.add(c)

        cmap_2 = ContactMap("cmap_2")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(3, 4, 0.1), Contact(2, 1, 0.2), Contact(1, 1, 0)]:
            cmap_2.add(c)

        cmap_1.add(Contact(6, 1, 0))
        sequence = Sequence('seq', 'AAAAAAA')
        absent_residues = {6, 7}

        output = tools.get_cmap_validation_metrics(cmap_1.as_dict(), cmap_2.as_dict(), sequence, absent_residues)
        expected_accuracy = [0.4, 0.52, 0.68, 0.56, 0.44]
        expected_fn = [0.6, 0.8, 0.8, 0.6, 0.4]
        expected_fnr = [0.2, 0.25, 0.25, 0.15, 0.1]
        expected_fp = [0.4, 0.6, 0.8, 0.6, 0.4]
        expected_fpr = [0.27, 0.47, 0.57, 0.5, 0.3]
        expected_sensitivity = [0.23, 0.23, 0.43, 0.3, 0.3]
        expected_specificity = [0.45, 0.6, 0.75, 0.65, 0.5]
        expected = [expected_accuracy, expected_fn, expected_fnr, expected_fp,
                    expected_fpr, expected_sensitivity, expected_specificity]

        self.assertListEqual([round(x, 2) for l in expected for x in l], [round(x, 2) for l in output[1] for x in l])

    def test_get_cmap_validation_metrics_3(self):
        cmap_1 = ContactMap("cmap_1")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2), Contact(1, 1, 0)]:
            cmap_1.add(c)

        cmap_2 = ContactMap("cmap_2")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(3, 4, 0.1), Contact(2, 1, 0.2), Contact(1, 1, 0)]:
            cmap_2.add(c)

        sequence = Sequence('seq', 'AAAAAAA')
        absent_residues = {2}

        output = tools.get_cmap_validation_metrics(cmap_1.as_dict(), cmap_2.as_dict(), sequence, absent_residues)

        expected_accuracy = [0.8333333333333334, np.nan, 0.8333333333333334, 0.8333333333333334, 0.8333333333333334]
        expected_fn = [0, np.nan, 1, 1, 0]
        expected_fnr = [0.0, np.nan, 0.2, 0.16666666666666666, 0.0]
        expected_fp = [1, np.nan, 0, 0, 1]
        expected_fpr = [0.3333333333333333, np.nan, 0.0, 0, 0.5]
        expected_sensitivity = [1.0, np.nan, 0.5, 0.0, 1.0]
        expected_specificity = [0.75, np.nan, 1.0, 1.0, 0.8]
        expected = np.array([expected_accuracy, expected_fn, expected_fnr, expected_fp,
                             expected_fpr, expected_sensitivity, expected_specificity])

        self.assertIsNone(np.testing.assert_allclose(expected, np.array(output[0])))

    def test_is_executable_1(self):
        self.assertEqual(sys.executable, tools.is_executable(sys.executable))

    def test_is_executable_2(self):
        self.assertEqual(sys.executable, tools.is_executable(os.path.basename(sys.executable)))

    def test_is_executable_3(self):
        with self.assertRaises(ValueError):
            tools.is_executable('qweasdzxcpoilkjmnb')


if __name__ == "__main__":
    unittest.main(verbosity=2)
