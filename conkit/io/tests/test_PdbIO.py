"""Testing facility for conkit.io.PdbIO"""

__author__ = "Felix Simkovic"
__date__ = "26 Oct 2016"

from conkit.io.PdbIO import PdbParser
from conkit.io._iotools import create_tmp_f

import os
import unittest


class Test(unittest.TestCase):

    def test_read_1(self):
        content = """ATOM      1  N   TYR A  36      39.107  51.628   3.103  0.50 43.13           N
ATOM      2  CA  TYR A  36      38.300  50.814   2.204  0.50 41.80           C
ATOM      3  O   TYR A  36      38.712  48.587   1.405  0.50 41.03           O
ATOM      4  CB  TYR A  36      37.586  51.694   1.175  0.50 41.61           C
ATOM      5  N   PHE A  86      32.465  47.498   5.487  0.50 25.81           N
ATOM      6  CA  PHE A  86      32.670  48.303   4.288  0.50 26.45           C
ATOM      7  O   PHE A  86      31.469  50.326   3.758  0.50 28.47           O
ATOM      8  CB  PHE A  86      32.977  47.392   3.090  0.50 25.35           C
ATOM      9  N   TRP A 171      23.397  37.507  -1.161  0.50 18.04           N
ATOM     10  CA  TRP A 171      23.458  36.846   0.143  0.50 20.46           C
ATOM     11  O   TRP A 171      22.235  34.954   0.951  0.50 22.45           O
ATOM     12  CB  TRP A 171      23.647  37.866   1.275  0.50 18.83           C
ATOM     13  N   PHE A 208      32.221  42.624  -5.829  0.50 19.96           N
ATOM     14  CA  PHE A 208      31.905  43.710  -4.909  0.50 20.31           C
ATOM     15  O   PHE A 208      32.852  45.936  -5.051  0.50 17.69           O
ATOM     16  CB  PHE A 208      31.726  43.102  -3.518  0.50 19.90           C
END
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = PdbParser().read(f_in, distance_cutoff=8, atom_type='CB')
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(6, len(contact_map1))
        self.assertEqual([36, 86], [c.res1_seq for c in contact_map1 if c.is_true_positive])
        self.assertEqual([86, 208], [c.res2_seq for c in contact_map1 if c.is_true_positive])
        self.assertEqual([0.934108, 0.920229], [c.raw_score for c in contact_map1 if c.is_true_positive])
        os.unlink(f_name)

    def test_read_2(self):
        content = """ATOM      1  N   TYR A  36      39.107  51.628   3.103  0.50 43.13           N
ATOM      2  CA  TYR A  36      38.300  50.814   2.204  0.50 41.80           C
ATOM      3  O   TYR A  36      38.712  48.587   1.405  0.50 41.03           O
ATOM      4  CB  TYR A  36      37.586  51.694   1.175  0.50 41.61           C
ATOM      5  N   PHE A  86      32.465  47.498   5.487  0.50 25.81           N
ATOM      6  CA  PHE A  86      32.670  48.303   4.288  0.50 26.45           C
ATOM      7  O   PHE A  86      31.469  50.326   3.758  0.50 28.47           O
ATOM      8  CB  PHE A  86      32.977  47.392   3.090  0.50 25.35           C
ATOM      9  N   TRP A 171      23.397  37.507  -1.161  0.50 18.04           N
ATOM     10  CA  TRP A 171      23.458  36.846   0.143  0.50 20.46           C
ATOM     11  O   TRP A 171      22.235  34.954   0.951  0.50 22.45           O
ATOM     12  CB  TRP A 171      23.647  37.866   1.275  0.50 18.83           C
ATOM     13  N   PHE A 208      32.221  42.624  -5.829  0.50 19.96           N
ATOM     14  CA  PHE A 208      31.905  43.710  -4.909  0.50 20.31           C
ATOM     15  O   PHE A 208      32.852  45.936  -5.051  0.50 17.69           O
ATOM     16  CB  PHE A 208      31.726  43.102  -3.518  0.50 19.90           C
END
        """
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = PdbParser().read(f_in, distance_cutoff=8, atom_type='CA')
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(6, len(contact_map1))
        self.assertEqual([36], [c.res1_seq for c in contact_map1 if c.is_true_positive])
        self.assertEqual([86], [c.res2_seq for c in contact_map1 if c.is_true_positive])
        self.assertEqual([0.934927], [c.raw_score for c in contact_map1 if c.is_true_positive])
        os.unlink(f_name)

    def test_read_3(self):
        content = """ATOM      1  N   TYR A  36      39.107  51.628   3.103  0.50 43.13           N
ATOM      2  CA  TYR A  36      38.300  50.814   2.204  0.50 41.80           C
ATOM      3  O   TYR A  36      38.712  48.587   1.405  0.50 41.03           O
ATOM      4  CB  TYR A  36      37.586  51.694   1.175  0.50 41.61           C
ATOM      5  N   PHE A  86      32.465  47.498   5.487  0.50 25.81           N
ATOM      6  CA  PHE A  86      32.670  48.303   4.288  0.50 26.45           C
ATOM      7  O   PHE A  86      31.469  50.326   3.758  0.50 28.47           O
ATOM      8  CB  PHE A  86      32.977  47.392   3.090  0.50 25.35           C
ATOM      9  N   TRP A 171      23.397  37.507  -1.161  0.50 18.04           N
ATOM     10  CA  TRP A 171      23.458  36.846   0.143  0.50 20.46           C
ATOM     11  O   TRP A 171      22.235  34.954   0.951  0.50 22.45           O
ATOM     12  CB  TRP A 171      23.647  37.866   1.275  0.50 18.83           C
ATOM     13  N   PHE A 208      32.221  42.624  -5.829  0.50 19.96           N
ATOM     14  CA  PHE A 208      31.905  43.710  -4.909  0.50 20.31           C
ATOM     15  O   PHE A 208      32.852  45.936  -5.051  0.50 17.69           O
ATOM     16  CB  PHE A 208      31.726  43.102  -3.518  0.50 19.90           C
END
        """
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = PdbParser().read(f_in, distance_cutoff=7, atom_type='CB')
        contact_map1 = contact_file.top_map
        self.assertEqual(1, len(contact_file))
        self.assertEqual(6, len(contact_map1))
        self.assertEqual([36], [c.res1_seq for c in contact_map1 if c.is_true_positive])
        self.assertEqual([86], [c.res2_seq for c in contact_map1 if c.is_true_positive])
        self.assertEqual([0.934108], [c.raw_score for c in contact_map1 if c.is_true_positive])
        os.unlink(f_name)

    def test_read_4(self):
        content = """ATOM      1  N   TYR A  36      39.107  51.628   3.103  0.50 43.13           N
ATOM      2  CA  TYR A  36      38.300  50.814   2.204  0.50 41.80           C
ATOM      3  O   TYR A  36      38.712  48.587   1.405  0.50 41.03           O
ATOM      4  CB  TYR A  36      37.586  51.694   1.175  0.50 41.61           C
ATOM      5  N   PHE A  86      32.465  47.498   5.487  0.50 25.81           N
ATOM      6  CA  PHE A  86      32.670  48.303   4.288  0.50 26.45           C
ATOM      7  O   PHE A  86      31.469  50.326   3.758  0.50 28.47           O
ATOM      8  CB  PHE A  86      32.977  47.392   3.090  0.50 25.35           C
TER
ATOM      9  N   TRP B 171      23.397  37.507  -1.161  0.50 18.04           N
ATOM     10  CA  TRP B 171      23.458  36.846   0.143  0.50 20.46           C
ATOM     11  O   TRP B 171      22.235  34.954   0.951  0.50 22.45           O
ATOM     12  CB  TRP B 171      23.647  37.866   1.275  0.50 18.83           C
ATOM     13  N   PHE B 208      32.221  42.624  -5.829  0.50 19.96           N
ATOM     14  CA  PHE B 208      31.905  43.710  -4.909  0.50 20.31           C
ATOM     15  O   PHE B 208      32.852  45.936  -5.051  0.50 17.69           O
ATOM     16  CB  PHE B 208      31.726  43.102  -3.518  0.50 19.90           C
END
"""
        f_name = create_tmp_f(content=content)
        with open(f_name, 'r') as f_in:
            contact_file = PdbParser().read(f_in, distance_cutoff=8, atom_type='CB')
        # Two maps because no contacts in B
        contact_map1 = contact_file['A']     # chain A
        contact_map2 = contact_file['AB']    # chain AB
        contact_map3 = contact_file['B']     # chain B
        self.assertEqual(3, len(contact_file))
        self.assertEqual(1, len(contact_map1))
        self.assertEqual(['A', 'A'], [contact_map1.top_contact.res1_chain, contact_map1.top_contact.res2_chain])
        self.assertEqual([36, 86], [contact_map1.top_contact.res1_seq, contact_map1.top_contact.res2_seq])
        self.assertEqual(4, len(contact_map2))
        self.assertEqual(['A', 'B'], [contact_map2.top_contact.res1_chain, contact_map2.top_contact.res2_chain])
        self.assertEqual([36, 171], [contact_map2.top_contact.res1_seq, contact_map2.top_contact.res2_seq])
        self.assertEqual(1, len(contact_map3))
        self.assertEqual(['B', 'B'], [contact_map3.top_contact.res1_chain, contact_map3.top_contact.res2_chain])
        self.assertEqual([171, 208], [contact_map3.top_contact.res1_seq, contact_map3.top_contact.res2_seq])
        os.unlink(f_name)


if __name__ == "__main__":
    unittest.main(verbosity=2)
