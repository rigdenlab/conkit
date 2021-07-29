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
            'dist': np.array([[np.random.dirichlet(np.ones(37)).tolist() for x in range(10)] for x in range(10)])
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
        expected_res2 = [6, 7, 8, 9, 10, 7, 8, 9, 10, 8, 9, 10, 9, 10, 10]
        expected_raw_score = [0.2777339274998595, 0.4065752412086529, 0.19140953877084063, 0.325080041345469,
                              0.2342294456365818, 0.2958748042140037, 0.20564014562328256, 0.2798115875666083,
                              0.29054612017825876, 0.3015477616596743, 0.2110736170498582, 0.3986234230845984,
                              0.2336248477198487, 0.28146943296360666, 0.25153867636261357]
        expected_bin_distance = [(14.5, 15), (17, 17.5), (8, 8.5), (8.5, 9), (17.5, 18), (17, 17.5), (19.5, 20),
                                 (16, 16.5), (6, 6.5), (2.5, 3), (17, 17.5), (10, 10.5), (19.5, 20), (16.5, 17),
                                 (14.5, 15)]
        expected_bin_score = [0.07978289246628485, 0.1079242231548592, 0.12906386318101604, 0.11786428170511123,
                              0.12079995199505335, 0.13050922159897016, 0.10735698279949522, 0.10870485113910544,
                              0.0702085454511652, 0.14528258648177836, 0.07600243873216728, 0.10630498169797178,
                              0.11173478543658984, 0.13264723820879454, 0.1152472852119953]

        self.assertListEqual(expected_res1, [distance.res1_seq for distance in distogram])
        self.assertListEqual(expected_res2, [distance.res2_seq for distance in distogram])
        self.assertListEqual(expected_raw_score, [contact.raw_score for contact in distogram])
        self.assertListEqual(expected_bin_score, [distance.max_score for distance in distogram])
        self.assertListEqual(expected_bin_distance, [distance.predicted_distance_bin for distance in distogram])
