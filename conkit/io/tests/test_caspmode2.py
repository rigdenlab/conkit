"""Testing facility for conkit.io.CaspMode2"""

import numpy as np
from conkit.core.distance import Distance
from conkit.core.distancefile import DistanceFile
from conkit.core.distogram import Distogram
from conkit.io.caspmode2 import CaspMode2Parser
from conkit.io.tests.helpers import ParserTestCase


class TestCaspMode2Parser(ParserTestCase):

    def test_read_1(self):
        content = """PFRMAT RR
TARGET T0999
AUTHOR 1234-5678-9000
REMARK Predictor remarks
METHOD Description of methods used
METHOD Description of methods used
RMODE 2
MODEL 1
1 8 .72 .345 .225 .15 .15 .1 .03 0 0 0 0
1 10 .715 .34 .225 .15 .15 .1 .035 0 0 0 0
31 38 .71 0 .56 .15 .12 .13 .04 0 0 0 0
10 20 .69 .34 .22 .13 .11 .1 .04 .06 0 0 0    
30 37 .67 .33 .21 .13 .1 .15 .04 .04 0 0 0
11 29 .65 .33 .21 .11 .12 .15 .04 .04 0 0 0
1 9 .63 0 .51 .12 .15 .15 .04 .03 0 0 0
21 37 .505 0 .2 .305 .25 .15 .07 .025 0 0 0
8 15 .4 0 .1 .3 .2 .1 .15 .05 .05 .05 0
3 14 .4 0 .2 .2 .1 .1 .15 .05 .05 .05 .1
5 15 .35 0 .15 .2 .1 .1 .15 .05 .05 .1 .1
7 14 .2 0 0 .2 .3 .1 .1 .05 .05 .1 .1
END"""
        expected_res1 = [1, 1, 31, 10, 30, 11, 1, 21, 8, 3, 5, 7]
        expected_res2 = [8, 10, 38, 20, 37, 29, 9, 37, 15, 14, 15, 14]
        expected_raw_score = [0.72, 0.715, 0.71, 0.69, 0.67, 0.65, 0.63, 0.505, 0.4, 0.4, 0.35, 0.2]
        expected_bin_distance = [(0, 4), (0, 4), (4, 6), (0, 4), (0, 4), (0, 4),
                                 (4, 6), (6, 8), (6, 8), (4, 6), (6, 8), (8, 10)]
        expected_bin_score = [0.345, 0.34, 0.56, 0.34, 0.33, 0.33, 0.51, 0.305, 0.3, 0.2, 0.2, 0.3]

        f_name = self.tempfile(content=content)
        with open(f_name, "r") as f_in:
            distancefile = CaspMode2Parser().read(f_in)

        self.assertIsInstance(distancefile, DistanceFile)
        self.assertEqual(1, len(distancefile))
        distogram = distancefile.top
        self.assertEqual('caspmode2', distogram.original_file_format)
        self.assertIsInstance(distogram, Distogram)
        self.assertEqual(12, distogram.ndistances)
        self.assertListEqual(expected_res1, [distance.res1_seq for distance in distogram])
        self.assertListEqual(expected_res2, [distance.res2_seq for distance in distogram])
        self.assertListEqual(expected_raw_score, [contact.raw_score for contact in distogram])
        self.assertListEqual(expected_bin_score, [distance.max_score for distance in distogram])
        self.assertListEqual(expected_bin_distance, [distance.predicted_distance_bin for distance in distogram])

    def test_write_1(self):
        expected_output = """PFRMAT RR
RMODE 2
1 6 0.199696 0.043889 0.085795 0.070011 0.071518 0.054028 0.213284 0.069087 0.097959 0.090083 0.204345
1 7 0.233644 0.049411 0.075135 0.109098 0.150810 0.096584 0.092398 0.096662 0.093350 0.123176 0.113375
1 8 0.246451 0.106886 0.039024 0.100540 0.082028 0.108344 0.078788 0.105980 0.130109 0.113708 0.134592
1 9 0.267139 0.072002 0.083053 0.112084 0.124356 0.128044 0.097491 0.132106 0.047198 0.110915 0.092751
1 10 0.351914 0.081445 0.069721 0.200748 0.099755 0.090368 0.117449 0.127677 0.050879 0.101965 0.059993
2 7 0.228459 0.085973 0.091366 0.051120 0.085890 0.070657 0.119253 0.082744 0.180051 0.097734 0.135213
2 8 0.256177 0.081094 0.077748 0.097335 0.060811 0.138077 0.130496 0.106911 0.101101 0.121346 0.085081
2 9 0.216631 0.046454 0.053018 0.117160 0.196036 0.144154 0.125199 0.090720 0.052621 0.098583 0.076055
2 10 0.284653 0.087567 0.125308 0.071778 0.071988 0.095966 0.099270 0.174715 0.109563 0.062611 0.101233
3 8 0.345583 0.117500 0.110134 0.117950 0.085312 0.098812 0.072826 0.079326 0.196758 0.059058 0.062325
3 9 0.203586 0.036574 0.050725 0.116287 0.174339 0.070881 0.116388 0.083683 0.060738 0.160257 0.130128
3 10 0.293849 0.059364 0.135117 0.099368 0.113124 0.135930 0.066876 0.075962 0.114771 0.127034 0.072454
4 9 0.234649 0.077170 0.048841 0.108638 0.107559 0.119732 0.116349 0.077063 0.111788 0.119497 0.113362
4 10 0.322930 0.090789 0.133412 0.098729 0.099123 0.084633 0.107534 0.137072 0.096560 0.042234 0.109913
5 10 0.279782 0.054314 0.114427 0.111042 0.069073 0.083048 0.105829 0.073806 0.119769 0.088666 0.180028"""

        distancefile = DistanceFile("test")
        distancefile.original_file_format = 'alphafold2'
        distogram = Distogram("1")
        distancefile.add(distogram)

        list_res1 = [1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5]
        list_res2 = [6, 7, 8, 9, 10, 7, 8, 9, 10, 8, 9, 10, 9, 10, 10]
        bin_edges = (2.3125, 2.625, 2.9375, 3.25, 3.5625, 3.875, 4.1875, 4.5, 4.8125, 5.125, 5.4375, 5.75,
                     6.0625, 6.375, 6.6875, 6.9999995, 7.3125, 7.625, 7.9375, 8.25, 8.5625, 8.875, 9.1875,
                     9.5, 9.812499, 10.124999, 10.4375, 10.75, 11.0625, 11.375, 11.687499, 12., 12.3125,
                     12.625, 12.9375, 13.25, 13.5625, 13.874999, 14.187501, 14.499999, 14.812499,
                     15.124999, 15.437499, 15.75, 16.0625, 16.375, 16.687502, 16.999998, 17.312498,
                     17.624998, 17.937498, 18.25, 18.5625, 18.875, 19.1875, 19.5, 19.8125, 20.125,
                     20.437498, 20.75, 21.062498, 21.374998, 21.6875)
        distance_bins = [(0, bin_edges[0])]
        distance_bins += [(bin_edges[idx], bin_edges[idx + 1]) for idx in range(len(bin_edges) - 1)]
        distance_bins.append((bin_edges[-1], np.inf))
        distance_bins = tuple(distance_bins)

        np.random.seed(41)
        for res_1, res_2 in zip(list_res1, list_res2):
            distance_scores = np.random.dirichlet(np.ones(64)).tolist()
            distance = Distance(res_1, res_2, distance_scores, distance_bins)
            distogram.add(distance)

        f_name = self.tempfile()
        with open(f_name, "w") as f_out:
            CaspMode2Parser().write(f_out, distogram)

        with open(f_name, "r") as f_in:
            output = f_in.read().splitlines()

        self.assertListEqual(expected_output.split('\n'), output)

    def test_read_and_write_1(self):
        content = """PFRMAT RR
RMODE 2
1 8 .72 .345 .225 .15 .15 .1 .03 0 0 0 0
1 10 .715 .34 .225 .15 .15 .1 .035 0 0 0 0
31 38 .71 0 .56 .15 .12 .13 .04 0 0 0 0
10 20 .69 .34 .22 .13 .11 .1 .04 .06 0 0 0    
30 37 .67 .33 .21 .13 .1 .15 .04 .04 0 0 0
11 29 .65 .33 .21 .11 .12 .15 .04 .04 0 0 0
1 9 .63 0 .51 .12 .15 .15 .04 .03 0 0 0
21 37 .505 0 .2 .305 .25 .15 .07 .025 0 0 0
8 15 .4 0 .1 .3 .2 .1 .15 .05 .05 .05 0
3 14 .4 0 .2 .2 .1 .1 .15 .05 .05 .05 .1
5 15 .35 0 .15 .2 .1 .1 .15 .05 .05 .1 .1
7 14 .2 0 0 .2 .3 .1 .1 .05 .05 .1 .1
END"""
        f_name_1 = self.tempfile(content=content)
        with open(f_name_1, "r") as f_in:
            distancefile = CaspMode2Parser().read(f_in)

        f_name_2 = self.tempfile()
        with open(f_name_2, "w") as f_out:
            CaspMode2Parser().write(f_out, distancefile)

        with open(f_name_2, "r") as f_in:
            output = CaspMode2Parser().read(f_in)

        expected_res1 = [1, 1, 31, 10, 30, 11, 1, 21, 8, 3, 5, 7]
        expected_res2 = [8, 10, 38, 20, 37, 29, 9, 37, 15, 14, 15, 14]
        expected_raw_score = [0.72, 0.715, 0.71, 0.69, 0.67, 0.65, 0.63, 0.505, 0.4, 0.4, 0.35, 0.2]
        expected_bin_distance = [(0, 4), (0, 4), (4, 6), (0, 4), (0, 4), (0, 4),
                                 (4, 6), (6, 8), (6, 8), (4, 6), (6, 8), (8, 10)]
        expected_bin_score = [0.345, 0.34, 0.56, 0.34, 0.33, 0.33, 0.51, 0.305, 0.3, 0.2, 0.2, 0.3]
        self.assertIsInstance(output, DistanceFile)
        self.assertEqual(1, len(output))
        distogram = output.top
        self.assertEqual('caspmode2', distogram.original_file_format)
        self.assertIsInstance(distogram, Distogram)
        self.assertEqual(12, distogram.ndistances)
        self.assertListEqual(expected_res1, [distance.res1_seq for distance in distogram])
        self.assertListEqual(expected_res2, [distance.res2_seq for distance in distogram])
        self.assertListEqual(expected_raw_score, [contact.raw_score for contact in distogram])
        self.assertListEqual(expected_bin_score, [round(distance.max_score, 3) for distance in distogram])
        self.assertListEqual(expected_bin_distance, [distance.predicted_distance_bin for distance in distogram])
