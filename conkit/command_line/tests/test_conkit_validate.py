"""Testing facility for conkit.conkit_validate"""

from conkit.command_line import conkit_validate

import unittest


class Test(unittest.TestCase):
    def test_parse_map_align_stdout(self):
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
            range(5, 15): [(5, 6), (6, 7), (7, 8), (8, 9), (9, 10), (11, 11), (12, 13), (13, 14), (14, 15)],
            range(21, 26): [(21, 23), (22, 24), (23, 25), (24, 26), (25, 27)]
        }

        output = {i.residue_range: i.residue_pairs for i in conkit_validate.parse_map_align_stdout(stdout_contents)}
        self.assertDictEqual(output, expected)

    def test_get_residue_ranges(self):
        residues = [0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 19, 25, 26, 30, 31, 32, 33, 34, 35]
        expected = [(0, 12), (19, 19), (25, 26), (30, 35)]

        output = conkit_validate.get_residue_ranges(residues)
        self.assertListEqual(output, expected)
