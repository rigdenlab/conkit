"""Testing facility for conkit.io.CaspMode2"""

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
        self.assertEqual('CASPRR_MODE_2', distogram.original_file_format)
        self.assertIsInstance(distogram, Distogram)
        self.assertEqual(12, distogram.ndistances)
        self.assertListEqual(expected_res1, [distance.res1_seq for distance in distogram])
        self.assertListEqual(expected_res2, [distance.res2_seq for distance in distogram])
        self.assertListEqual(expected_raw_score, [contact.raw_score for contact in distogram])
        self.assertListEqual(expected_bin_score, [distance.max_score for distance in distogram])
        self.assertListEqual(expected_bin_distance, [distance.predicted_distance_bin for distance in distogram])
