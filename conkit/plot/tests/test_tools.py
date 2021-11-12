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

        expected = {
            range(5, 15): [(5, 6), (6, 7), (7, 8), (8, 9), (10, 10), (11, 11), (12, 13), (13, 14), (14, 15)],
            range(21, 26): [(21, 23), (22, 24), (23, 25), (24, 26), (25, 27)]
        }

        output = {i.residue_range: i.residue_pairs for i in tools.parse_map_align_stdout(stdout_contents)}
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


if __name__ == "__main__":
    unittest.main(verbosity=2)
