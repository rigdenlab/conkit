"""Testing facility for conkit.command_line.tools"""

from conkit.command_line import cli_tools

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
MAX 2_2_2	/home/filo/opt/map_align_v1/map_align/3u97_A.gremlin.map	/home/filo/opt/map_align_v1/map_align/2pd0_A.pdb.map	44.5272	-4.38	40.1472	74	1:1	2:2	3:3	4:4	5:5	6:6	7:7	8:8	9:9	10:10	11:11	12:12	13:13	14:14	15:20	16:21	17:22	18:23	19:2420:25	21:26	22:27	23:30	24:31	25:32	26:33	27:34	28:35	29:36	30:37	31:42	32:43	33:44	34:45	35:46	36:47	37:48	38:49	39:52	40:53	42:54	43:55	44:56	45:57	46:58	47:60	48:61	49:62	50:63	51:64	52:65	53:113	54:114	55:115	56:116	57:117	58:118	59:119	60:12061:143	62:144	63:145	64:146	65:147	66:148	67:149	68:150	69:151	70:152	71:153	72:154	73:155	74:156	75:157
"""
        expected = {1: 1, 2: 2, 3: 3, 4: 4, 5: 5, 6: 6, 7: 7, 8: 8, 9: 9,
                    10: 10, 11: 11, 12: 12, 13: 13, 14: 14, 15: 20, 16: 21, 17: 22,
                    18: 23, 19: 24, 20: 25, 21: 26, 22: 27, 23: 30, 24: 31, 25: 32,
                    26: 33, 27: 34, 28: 35, 29: 36, 30: 37, 31: 42, 32: 43, 33: 44,
                    34: 45, 35: 46, 36: 47, 37: 48, 38: 49, 39: 52, 40: 53, 42: 54,
                    43: 55, 44: 56, 45: 57, 46: 58, 47: 60, 48: 61, 49: 62, 50: 63,
                    51: 64, 52: 65, 53: 113, 54: 114, 55: 115, 56: 116, 57: 117,
                    58: 118, 59: 119, 60: 120, 61: 143, 62: 144, 63: 145, 64: 146,
                    65: 147, 66: 148, 67: 149, 68: 150, 69: 151, 70: 152, 71: 153,
                    72: 154, 73: 155, 74: 156,
                    75: 157}

        output = cli_tools.parse_map_align_stdout(stdout_contents)
        self.assertDictEqual(output, expected)
