"""Testing facility for conkit.io.AlphaFold2Parser"""

import numpy as np
import pickle
from conkit.core.distancefile import DistanceFile
from conkit.core.distogram import Distogram
from conkit.io.alphafold import AlphaFold2Parser
from conkit.io.tests.helpers import ParserTestCase


class TestAlphaFold2Parser(ParserTestCase):

    def test_read_1(self):
        np.random.seed(41)
        prediction = {'distogram': {
            'bin_edges': np.array([2.3125, 2.625, 2.9375, 3.25, 3.5625, 3.875, 4.1875, 4.5, 4.8125, 5.125, 5.4375, 5.75,
                                   6.0625, 6.375, 6.6875, 6.9999995, 7.3125, 7.625, 7.9375, 8.25, 8.5625, 8.875, 9.1875,
                                   9.5, 9.812499, 10.124999, 10.4375, 10.75, 11.0625, 11.375, 11.687499, 12., 12.3125,
                                   12.625, 12.9375, 13.25, 13.5625, 13.874999, 14.187501, 14.499999, 14.812499,
                                   15.124999, 15.437499, 15.75, 16.0625, 16.375, 16.687502, 16.999998, 17.312498,
                                   17.624998, 17.937498, 18.25, 18.5625, 18.875, 19.1875, 19.5, 19.8125, 20.125,
                                   20.437498, 20.75, 21.062498, 21.374998, 21.6875]),
            'logits': np.array([[np.random.dirichlet(np.ones(64)).tolist() for x in range(5)] for x in range(5)])
        }}

        f_name = self.tempfile(content=None)
        with open(f_name, 'wb') as fhandle:
            pickle.dump(prediction, fhandle)
        with open(f_name, "rb") as fhandle:
            distancefile = AlphaFold2Parser().read(fhandle)

        self.assertIsInstance(distancefile, DistanceFile)
        self.assertEqual(1, len(distancefile))
        distogram = distancefile.top
        self.assertEqual('alphafold2', distogram.original_file_format)
        self.assertIsInstance(distogram, Distogram)
        self.assertEqual(15, len(distogram))

        expected_res1 = [1, 1, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 5]
        expected_res2 = [1, 2, 3, 4, 5, 2, 3, 4, 5, 3, 4, 5, 4, 5, 5]
        expected_raw_score = [0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3]
        expected_bin_distance = [(12.625, 12.9375), (19.5, 19.8125), (16.0625, 16.375), (21.062498, 21.374998),
                                 (6.9999995, 7.3125), (10.124999, 10.4375), (7.9375, 8.25), (14.499999, 14.812499),
                                 (16.687502, 16.999998), (8.25, 8.5625), (5.125, 5.4375), (20.75, 21.062498),
                                 (6.375, 6.6875), (18.5625, 18.875), (7.9375, 8.25)]
        expected_bin_score = [0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02, 0.02]

        self.assertListEqual(expected_res1, [distance.res1_seq for distance in distogram])
        self.assertListEqual(expected_res2, [distance.res2_seq for distance in distogram])
        self.assertListEqual(expected_raw_score, [round(contact.raw_score, 2) for contact in distogram])
        self.assertListEqual(expected_bin_score, [round(distance.max_score, 2) for distance in distogram])
        self.assertListEqual(expected_bin_distance, [distance.predicted_distance_bin for distance in distogram])
