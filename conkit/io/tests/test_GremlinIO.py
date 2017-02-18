"""Testing facility for conkit.io.GremlinIO"""

__author__ = "Felix Simkovic"
__date__ = "04 Oct 2016"

from conkit.core import Contact
from conkit.core import ContactFile
from conkit.core import ContactMap
from conkit.core import Sequence
from conkit.io.GremlinIO import GremlinParser
from conkit.io._iotools import create_tmp_f

import os
import unittest


class Test(unittest.TestCase):

    def test_read_1(self):
        content = """i	j	i_id	j_id	r_sco	s_sco	prob
179	246	179_C	246_L	0.2019	4.740	1.000
262	305	262_G	305_Y	0.1742	4.090	1.000
428	448	428_A	448_N	0.1638	3.846	1.000
214	231	214_F	231_V	0.1342	3.150	1.000
457	488	457_L	488_Y	0.1254	2.945	1.000
220	223	220_A	223_A	0.1187	2.786	0.999
143	209	143_I	209_D	0.1139	2.674	0.999
79	365	79_M	365_I	0.1114	2.615	0.998
215	268	215_V	268_A	0.1109	2.604	0.998
262	266	262_G	266_K	0.1040	2.442	0.997
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = GremlinParser().read(f_in)
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(10, len(contact_map1))
        self.assertEqual(
            [179, 262, 428, 214, 457, 220, 143, 79, 215, 262],
            [c.res1_seq for c in contact_map1]
        )
        self.assertEqual(
            [0.2019, 0.1742, 0.1638, 0.1342, 0.1254, 0.1187, 0.1139, 0.1114, 0.1109, 0.1040],
            [c.raw_score for c in contact_map1]
        )
        os.unlink(f_name)

    def test_read_2(self):
        content = """# Some comments
# That are here for whatever reason
i	j	i_id	j_id	r_sco	s_sco	prob
179	246	179_C	246_L	0.2019	4.740	1.000
262	305	262_G	305_Y	0.1742	4.090	1.000
428	448	428_A	448_N	0.1638	3.846	1.000
214	231	214_F	231_V	0.1342	3.150	1.000
457	488	457_L	488_Y	0.1254	2.945	1.000
220	223	220_A	223_A	0.1187	2.786	0.999
143	209	143_I	209_D	0.1139	2.674	0.999
79	365	79_M	365_I	0.1114	2.615	0.998
215	268	215_V	268_A	0.1109	2.604	0.998
262	266	262_G	266_K	0.1040	2.442	0.997
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = GremlinParser().read(f_in)
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(10, len(contact_map1))
        self.assertEqual(
            [179, 262, 428, 214, 457, 220, 143, 79, 215, 262],
            [c.res1_seq for c in contact_map1]
        )
        self.assertEqual(
            [0.2019, 0.1742, 0.1638, 0.1342, 0.1254, 0.1187, 0.1139, 0.1114, 0.1109, 0.1040],
            [c.raw_score for c in contact_map1]
        )
        os.unlink(f_name)

    def test_read_3(self):
        content = """i	j	gene	i_id	j_id	r_sco	s_sco	prob	I_prob
127	187	A	127_V	187_I	0.183	3.635	1.000	N/A
83	87	A	83_E	87_Q	0.183	3.633	1.000	N/A
108	111	A	108_P	111_P	0.105	2.095	0.989	N/A
431	435	B	241_L	245_L	0.104	2.076	0.988	N/A
63	83	A	63_T	83_E	0.098	1.952	0.980	N/A
23	434	AB	23_T	244_L	0.082	1.624	0.924	0.519
20	438	AB	20_Y	248_T	0.059	1.178	0.647	0.181
265	275	B	75_E	85_V	0.059	1.175	0.644	N/A
263	267	B	73_A	77_G	0.059	1.172	0.641	N/A
19	438	AB	19_L	248_T	0.059	1.17	0.640	0.176
211	215	B	21_D	25_A	0.054	1.069	0.536	N/A
30	65	A	30_A	65_T	0.054	1.065	0.532	N/A
24	434	AB	24_A	244_L	0.054	1.064	0.531	0.123
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = GremlinParser().read(f_in)
        self.assertEqual(3, len(contact_file))
        chain_a_res1seq = [127, 83, 108, 63, 30]
        chain_a_rawscore = [0.183, 0.183, 0.105, 0.098, 0.054]
        chain_b_res1seq = [431, 265, 263, 211]
        chain_b_rawscore = [0.104, 0.059, 0.059, 0.054]
        chain_ab_res1seq = [23, 20, 19, 24]
        chain_ab_rawscore = [0.082, 0.059, 0.059, 0.054]
        for count, res1_seqs, raw_scores, cmap in zip([5, 4, 4],
                                                      [chain_a_res1seq, chain_ab_res1seq, chain_b_res1seq],
                                                      [chain_a_rawscore, chain_ab_rawscore, chain_b_rawscore],
                                                      contact_file):
            self.assertEqual(count, len(cmap))
            self.assertEqual(res1_seqs, [c.res1_seq for c in cmap])
            self.assertEqual(raw_scores, [c.raw_score for c in cmap])
        os.unlink(f_name)

    def test_write_1(self):
        contact_file = ContactFile('test')
        contact_map = ContactMap('A')
        contact_file.add(contact_map)
        for c in [(1, 9, 0, 8, 0.7), (1, 10, 0, 8, 0.7), (2, 8, 0, 8, 0.9), (3, 12, 0, 8, 0.4)]:
            contact = Contact(c[0], c[1], c[4], distance_bound=(c[2], c[3]))
            contact_map.add(contact)
        contact_map.sequence = Sequence('1', 'HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSDHLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD')
        contact_map.assign_sequence_register()
        f_name = create_tmp_f()
        with open(f_name, 'w') as f_out:
            GremlinParser().write(f_out, contact_file)
        content = [
            "i	j	i_id	j_id	r_sco	s_sco	prob",
            "1	9	1_H	9_L	0.7	1.0	1.0",
            "1	10	1_H	10_L	0.7	1.0	1.0",
            "2	8	2_L	8_I	0.9	1.3	1.0",
            "3	12	3_E	12_K	0.4	0.6	1.0",
            "",
        ]
        content = os.linesep.join(content)
        with open(f_name, 'r') as f_in:
            data = "".join(f_in.readlines())
        self.assertEqual(content, data)
        os.unlink(f_name)

    def test_write_2(self):
        contact_file = ContactFile('TEST')
        contact_map = ContactMap('1')
        contact_file.add(contact_map)
        for c in [(1, 9, 0, 8, 0.7), (1, 10, 0, 8, 0.7), (2, 8, 0, 8, 0.9), (3, 12, 0, 8, 0.4)]:
            contact = Contact(c[0], c[1], c[4], distance_bound=(c[2], c[3]))
            contact_map.add(contact)
        f_name = create_tmp_f()
        with open(f_name, 'w') as f_out:
            GremlinParser().write(f_out, contact_file)
        content = [
            "i	j	i_id	j_id	r_sco	s_sco	prob",
            "1	9	1_X	9_X	0.7	1.0	1.0",
            "1	10	1_X	10_X	0.7	1.0	1.0",
            "2	8	2_X	8_X	0.9	1.3	1.0",
            "3	12	3_X	12_X	0.4	0.6	1.0",
            "",
        ]
        content = os.linesep.join(content)
        with open(f_name, 'r') as f_in:
            data = "".join(f_in.readlines())
        self.assertEqual(content, data)
        os.unlink(f_name)

    def test_write_3(self):
        contact_file = ContactFile('TEST')
        contact_maps = [ContactMap('A'), ContactMap('AB'), ContactMap('B')]
        contacts = [(Contact(1, 9, 0.7), Contact(1, 10, 0.7), Contact(2, 8, 0.9), Contact(3, 12, 0.4)),
                    (Contact(1, 9, 0.7), Contact(1, 10, 0.7), Contact(2, 8, 0.9), Contact(3, 12, 0.4)),
                    (Contact(1, 9, 0.7), Contact(1, 10, 0.7), Contact(2, 8, 0.9), Contact(3, 12, 0.4))]
        chains = [('A', 'A'), ('A', 'B'), ('B', 'B')]
        for contact_map, contacts, chain in zip(contact_maps, contacts, chains):
            contact_file.add(contact_map)
            for c in contacts:
                c.res1_chain = chain[0]
                c.res2_chain = chain[1]
                contact_map.add(c)
            contact_map.sequence = Sequence('1', 'HLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSDHLEGSIGILLKKHEIVFDGCHDFGRTYIWQMSD')
            contact_map.assign_sequence_register()
        f_name = create_tmp_f()
        with open(f_name, 'w') as f_out:
            GremlinParser().write(f_out, contact_file)
        content = [
            "i	j	gene	i_id	j_id	r_sco	s_sco	prob	I_prob",
            "1	9	A	1_H	9_L	0.7	1.0	1.0	N/A",
            "1	10	A	1_H	10_L	0.7	1.0	1.0	N/A",
            "2	8	A	2_L	8_I	0.9	1.3	1.0	N/A",
            "3	12	A	3_E	12_K	0.4	0.6	1.0	N/A",
            "1	9	AB	1_H	9_L	0.7	1.0	1.0	N/A",
            "1	10	AB	1_H	10_L	0.7	1.0	1.0	N/A",
            "2	8	AB	2_L	8_I	0.9	1.3	1.0	N/A",
            "3	12	AB	3_E	12_K	0.4	0.6	1.0	N/A",
            "1	9	B	1_H	9_L	0.7	1.0	1.0	N/A",
            "1	10	B	1_H	10_L	0.7	1.0	1.0	N/A",
            "2	8	B	2_L	8_I	0.9	1.3	1.0	N/A",
            "3	12	B	3_E	12_K	0.4	0.6	1.0	N/A",
            "",
        ]
        content = os.linesep.join(content)
        with open(f_name, 'r') as f_in:
            data = "".join(f_in.readlines())
        self.assertEqual(content, data)
        os.unlink(f_name)


if __name__ == "__main__":
    unittest.main(verbosity=2)
