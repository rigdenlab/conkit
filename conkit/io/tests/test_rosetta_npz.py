"""Testing facility for conkit.io.RosettaNpzParser"""

import numpy as np
import pickle
from conkit.core.distancefile import DistanceFile
from conkit.core.distogram import Distogram
from conkit.io.rosetta_npz import RosettaNpzParser
from conkit.io.tests.helpers import ParserTestCase


class TestRosettaNpzParser(ParserTestCase):

    def test_read_1(self):
        np.random.seed(41)
        prediction = {
            'dist': np.array([[np.random.dirichlet(np.ones(37)).tolist() for x in range(5)] for x in range(5)])
        }

        f_name = self.tempfile(content=None)
        with open(f_name, 'wb') as fhandle:
            pickle.dump(prediction, fhandle)
        with open(f_name, "rb") as fhandle:
            distancefile = RosettaNpzParser().read(fhandle)

        self.assertIsInstance(distancefile, DistanceFile)
        self.assertEqual(1, len(distancefile))
        distogram = distancefile.top
        self.assertEqual('ROSETTA_NPZ', distogram.original_file_format)
        self.assertIsInstance(distogram, Distogram)
        self.assertEqual(15, len(distogram))

        expected_res1 = [1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5]
        expected_res2 = [1, 2, 3, 4, 5, 2, 3, 4, 5, 3, 4, 5, 4, 5, 5]
        expected_raw_score = [0.2542024725712152, 0.2680150396842064, 0.31951115321526485, 0.2872068644639303,
                              0.3270142179786386, 0.4065752412086529, 0.19140953877084063, 0.325080041345469,
                              0.2342294456365818, 0.15548141805069934, 0.41178206457702593, 0.26026142628351123,
                              0.2798115875666083, 0.29054612017825876, 0.3008523890315092]

        expected_bin_distance = [(18.5, 19), (12, 12.5), (9, 9.5), (6, 6.5), (14, 14.5), (17, 17.5), (8, 8.5),
                                 (8.5, 9), (17.5, 18), (13, 13.5), (4.5, 5), (18.5, 19), (16, 16.5), (6, 6.5), (4, 4.5)]

        expected_bin_score = [0.25350480891411237, 0.08088144636134563, 0.08190930855105415, 0.1383295289038701,
                              0.09824995972452338, 0.1079242231548592, 0.12906386318101604, 0.11786428170511123,
                              0.12079995199505335, 0.15691065654436132, 0.13879490895313662, 0.10757779064007214,
                              0.10870485113910544, 0.0702085454511652, 0.0968560988412983]


        self.assertListEqual(expected_res1, [distance.res1_seq for distance in distogram])
        self.assertListEqual(expected_res2, [distance.res2_seq for distance in distogram])
        self.assertListEqual(expected_raw_score, [contact.raw_score for contact in distogram])
        self.assertListEqual(expected_bin_score, [distance.max_score for distance in distogram])
        self.assertListEqual(expected_bin_distance, [distance.predicted_distance_bin for distance in distogram])
