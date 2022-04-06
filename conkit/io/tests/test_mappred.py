"""Testing facility for conkit.io.MapPred"""

import numpy as np
from conkit.core.distance import Distance
from conkit.core.distancefile import DistanceFile
from conkit.core.distogram import Distogram
from conkit.io.mappred import MapPredParser
from conkit.io.tests.helpers import ParserTestCase


class TestMapPredParser(ParserTestCase):

    def test_read_1(self):
        content = """#REMARK MapPred 1.1
#REMARK idx_i, idx_j, distance distribution of 34 bins
#REMARK 34 bins consist of 32 normal bins (4-20A with a step of 0.5A) and two boundary bins ( [0,4) and [20, inf) ), as follows: [0,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,inf]
5 10 0.0129 0.0358 0.1879 0.2688 0.4436 0.0364 0.0016 0.0011 0.0005 0.0004 0.0004 0.0005 0.0005 0.0004 0.0004 0.0003 0.0002 0.0002 0.0002 0.0001 0.0002 0.0002 0.0002 0.0003 0.0002 0.0002 0.0002 0.0001 0.0001 0.0001 0.0001 0.0000 0.0000 0.0058
1 35 0.0002 0.0017 0.0046 0.0158 0.0374 0.0660 0.1207 0.2161 0.2476 0.1666 0.0729 0.0354 0.0049 0.0013 0.0008 0.0004 0.0002 0.0002 0.0002 0.0002 0.0002 0.0001 0.0002 0.0002 0.0002 0.0001 0.0001 0.0000 0.0000 0.0001 0.0000 0.0000 0.0000 0.0058
43 85 0.0009 0.0027 0.0044 0.0115 0.0209 0.0149 0.0120 0.0190 0.0210 0.0383 0.0718 0.1105 0.1392 0.1948 0.1630 0.1048 0.0407 0.0159 0.0033 0.0006 0.0004 0.0003 0.0002 0.0002 0.0002 0.0002 0.0001 0.0002 0.0002 0.0001 0.0001 0.0001 0.0001 0.0073
85 43 0.0009 0.0027 0.0044 0.0115 0.0209 0.0149 0.0120 0.0190 0.0210 0.0383 0.0718 0.1105 0.1392 0.1948 0.1630 0.1048 0.0407 0.0159 0.0033 0.0006 0.0004 0.0003 0.0002 0.0002 0.0002 0.0002 0.0001 0.0002 0.0002 0.0001 0.0001 0.0001 0.0001 0.0073
50 50 0.0009 0.0027 0.0044 0.0115 0.0209 0.0149 0.0120 0.0190 0.0210 0.0383 0.0718 0.1105 0.1392 0.1948 0.1630 0.1048 0.0407 0.0159 0.0033 0.0006 0.0004 0.0003 0.0002 0.0002 0.0002 0.0002 0.0001 0.0002 0.0002 0.0001 0.0001 0.0001 0.0001 0.0073
18 50 0.0006 0.0009 0.0010 0.0018 0.0024 0.0027 0.0027 0.0032 0.0043 0.0052 0.0068 0.0105 0.0156 0.0222 0.0298 0.0526 0.0895 0.1389 0.1769 0.1865 0.1278 0.0709 0.0316 0.0096 0.0021 0.0005 0.0002 0.0001 0.0001 0.0001 0.0001 0.0001 0.0001 0.0028
"""
        expected_res1 = [5, 1, 43, 85, 50, 18]
        expected_res2 = [10, 35, 85, 43, 50, 50]
        expected_raw_score = [0.9885999999999999, 0.7101, 0.1073, 0.1073, 0.1073, 0.0196]
        expected_bin_distance = [(5.5, 6), (7.5, 8), (10, 10.5), (10, 10.5), (10, 10.5), (13, 13.5)]
        expected_bin_score = [0.4436, 0.2476, 0.1948, 0.1948, 0.1948, 0.1865]

        f_name = self.tempfile(content=content)
        with open(f_name, "r") as f_in:
            distancefile = MapPredParser().read(f_in)

        self.assertIsInstance(distancefile, DistanceFile)
        self.assertEqual(1, len(distancefile))
        distogram = distancefile.top
        self.assertEqual('mappred', distogram.original_file_format)
        self.assertIsInstance(distogram, Distogram)
        self.assertEqual(6, distogram.ndistances)
        self.assertListEqual(expected_res1, [distance.res1_seq for distance in distogram])
        self.assertListEqual(expected_res2, [distance.res2_seq for distance in distogram])
        self.assertListEqual(expected_raw_score, [contact.raw_score for contact in distogram])
        self.assertListEqual(expected_bin_score, [distance.max_score for distance in distogram])
        self.assertListEqual(expected_bin_distance, [distance.predicted_distance_bin for distance in distogram])

    def test_write_1(self):
        expected_output = """#REMARK MapPred 1.1
#REMARK idx_i, idx_j, distance distribution of 34 bins
#REMARK 34 bins consist of 32 normal bins (4-20A with a step of 0.5A) and two boundary bins ( [0,4) and [20, inf) ), as follows: [0,4,4.5,5,5.5,6,6.5,7,7.5,8,8.5,9,9.5,10,10.5,11,11.5,12,12.5,13,13.5,14,14.5,15,15.5,16,16.5,17,17.5,18,18.5,19,19.5,20,inf]
5 10 0.013746 0.002245 0.053742 0.002115 0.005889 0.044058 0.010081 0.052535 0.118677 0.025818 0.019215 0.015831 0.009808 0.018148 0.031220 0.003428 0.058081 0.017978 0.065069 0.024163 0.044585 0.062025 0.026062 0.023824 0.012573 0.027729 0.022212 0.041685 0.005015 0.064340 0.004133 0.006420 0.018552 0.048998
1 35 0.187103 0.008180 0.021642 0.051089 0.038619 0.006100 0.010553 0.031697 0.010831 0.015310 0.006949 0.008237 0.043400 0.051436 0.003820 0.008148 0.018467 0.057307 0.022873 0.029184 0.008235 0.008025 0.004214 0.027027 0.070948 0.028355 0.049284 0.060124 0.041885 0.043900 0.000681 0.006836 0.007679 0.011862
43 85 0.024968 0.014838 0.021987 0.031265 0.019144 0.033038 0.018177 0.008716 0.017331 0.046459 0.051147 0.043912 0.004041 0.007990 0.027690 0.073997 0.001269 0.008161 0.067709 0.055700 0.028615 0.091884 0.021842 0.025949 0.025295 0.006136 0.031655 0.028990 0.082802 0.005069 0.002322 0.015611 0.039637 0.016654
85 43 0.015871 0.013765 0.006593 0.014670 0.029273 0.042705 0.058513 0.014858 0.050493 0.014216 0.010146 0.037020 0.018679 0.003142 0.031215 0.011736 0.008920 0.007325 0.144325 0.003512 0.018591 0.005043 0.001607 0.043659 0.068744 0.052532 0.050643 0.039295 0.003413 0.035119 0.102032 0.004150 0.005737 0.032456
50 50 0.000490 0.027392 0.001090 0.009625 0.011421 0.002011 0.015100 0.018622 0.008785 0.114531 0.044962 0.019562 0.022973 0.008111 0.042691 0.061367 0.001060 0.032753 0.073944 0.006790 0.002509 0.073759 0.025060 0.031361 0.039123 0.043318 0.032752 0.004280 0.044655 0.000556 0.000111 0.095043 0.028036 0.056157
18 50 0.002704 0.015000 0.024442 0.105520 0.014259 0.027628 0.002832 0.035063 0.038354 0.055931 0.039683 0.035546 0.004621 0.019932 0.012316 0.087781 0.006637 0.043857 0.008459 0.053482 0.016937 0.083507 0.031733 0.000793 0.004304 0.066937 0.009968 0.006859 0.038950 0.064003 0.003185 0.008042 0.007331 0.023401"""
        distancefile = DistanceFile("test")
        distancefile.original_file_format = 'mappred'
        distogram = Distogram("1")
        distancefile.add(distogram)

        list_res1 = [5, 1, 43, 85, 50, 18]
        list_res2 = [10, 35, 85, 43, 50, 50]
        distance_bins = ((0, 4), (4, 4.5), (4.5, 5), (5, 5.5), (5.5, 6), (6, 6.5), (6.5, 7), (7, 7.5), (7.5, 8),
                         (8, 8.5), (8.5, 9), (9, 9.5), (9.5, 10), (10, 10.5), (10.5, 11), (11, 11.5), (11.5, 12),
                         (12, 12.5), (12.5, 13), (13, 13.5), (13.5, 14), (14, 14.5), (14.5, 15), (15, 15.5),
                         (15.5, 16), (16, 16.5), (16.5, 17), (17, 17.5), (17.5, 18), (18, 18.5), (18.5, 19),
                         (19, 19.5), (19.5, 20), (20, np.inf))

        np.random.seed(41)
        for res_1, res_2 in zip(list_res1, list_res2):
            distance_scores = np.random.dirichlet(np.ones(34)).tolist()
            distance = Distance(res_1, res_2, distance_scores, distance_bins)
            distogram.add(distance)

        f_name = self.tempfile()
        with open(f_name, "w") as f_out:
            MapPredParser().write(f_out, distogram)

        with open(f_name, "r") as f_in:
            output = f_in.read().splitlines()

        self.assertListEqual(expected_output.split("\n"), output)
