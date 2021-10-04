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
MAX 2_2_2	/home/filo/opt/map_align_v1/map_align/3u97_A.gremlin.map	/home/filo/opt/map_align_v1/map_align/2pd0_A.pdb.map	44.5272	-4.38	40.1472	74	1:1	2:2	3:3	4:4	5:5	6:6	7:7	8:8	9:9	10:10	11:11	12:12	13:13	14:14	15:20	16:21	17:22	18:23	19:24 20:25	21:26	22:27	23:30	24:31	25:32	26:33	27:34	28:35	29:36	30:37	31:42	32:43	33:44	34:45	35:46	36:47	37:48	38:49	39:52	40:53	42:54	43:55	44:56	45:57	46:58	47:60	48:61	49:62	50:63	51:64	52:65	53:113	54:114	55:115	56:116	57:117	58:118	59:119	60:120 61:143	62:144	63:145	64:146	65:147	66:148	67:149	68:150	69:151	70:152	71:153	72:154	73:155	74:156	75:157
"""
        expected = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9,
                    10: 10, 11: 11, 12: 12, 13: 13, 14: 14, 20: 15, 21: 16, 22: 17,
                    23: 18, 24: 19, 25: 20, 26: 21, 27: 22, 30: 23, 31: 24, 32: 25,
                    33: 26, 34: 27, 35: 28, 36: 29, 37: 30, 42: 31, 43: 32, 44: 33,
                    45: 34, 46: 35, 47: 36, 48: 37, 49: 38, 52: 39, 53: 40, 54: 42,
                    55: 43, 56: 44, 57: 45, 58: 46, 60: 47, 61: 48, 62: 49, 63: 50,
                    64: 51, 65: 52, 113: 53, 114: 54, 115: 55, 116: 56, 117: 57,
                    118: 58, 119: 59, 120: 60, 143: 61, 144: 62, 145: 63, 146: 64,
                    147: 65, 148: 66, 149: 67, 150: 68, 151: 69, 152: 70, 153: 71,
                    154: 72, 155: 73, 156: 74, 157: 75}

        output = conkit_validate.parse_map_align_stdout(stdout_contents)
        self.assertDictEqual(output, expected)

    def test_get_residue_ranges(self):
        residues = [0, 1, 2, 3, 4, 5, 6, 9, 10, 11, 12, 19, 25, 26, 30, 31, 32, 33, 34, 35]
        expected = [(0, 12), (19, 19), (25, 26), (30, 35)]

        output = conkit_validate.get_residue_ranges(residues)
        self.assertListEqual(output, expected)
