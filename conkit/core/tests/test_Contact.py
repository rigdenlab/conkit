"""Testing facility for conkit.core.contact"""

__author__ = "Felix Simkovic"
__date__ = "12 Aug 2016"

from conkit.core.Contact import Contact

import unittest


class Test(unittest.TestCase):

    def test_distance_bound_1(self):
        contact = Contact(1, 2, 1.0)
        contact.distance_bound = (0, 8)
        self.assertEqual((0, 8), contact.distance_bound)

    def test_distance_bound_2(self):
        contact = Contact(1, 2, 1.0)
        contact.distance_bound = (0, 8)
        self.assertEqual((0, 8), contact.distance_bound)
        self.assertTrue(isinstance(contact.distance_bound, tuple))

    def test_distance_bound_3(self):
        contact = Contact(1, 2, 1.0)
        try:
            contact.distance_bound = (1, 100)
        except TypeError:
            self.fail('contact.distance_bound raised TypeError unexpectedly')
        with self.assertRaises(TypeError):
            contact.distance_bound = 8
        with self.assertRaises(TypeError):
            contact.distance_bound = 'test'

    def test_is_false_positive_1(self):
        contact = Contact(1, 2, 1.0)
        contact.define_false_positive()
        self.assertTrue(contact.is_false_positive)
        self.assertFalse(contact.is_true_positive)

    def test_is_true_positive_1(self):
        contact = Contact(1, 2, 1.0)
        contact.define_true_positive()
        self.assertTrue(contact.is_true_positive)
        self.assertFalse(contact.is_false_positive)

    def test_lower_bound_1(self):
        contact = Contact(1, 2, 1.0)
        self.assertEqual(0, contact.lower_bound)

    def test_lower_bound_2(self):
        contact = Contact(1, 2, 1.0)
        with self.assertRaises(ValueError):
            contact.lower_bound = -1

    def test_lower_bound_3(self):
        contact = Contact(1, 2, 1.0)
        contact.lower_bound = 1
        self.assertEqual(1, contact.lower_bound)

    def test_lower_bound_4(self):
        contact = Contact(1, 2, 1.0)
        contact.lower_bound = 7
        self.assertEqual(7, contact.lower_bound)

    def test_lower_bound_5(self):
        contact = Contact(1, 2, 1.0)
        contact.lower_bound = 1
        contact.upper_bound = 8
        self.assertEqual(8, contact.upper_bound)
        with self.assertRaises(ValueError):
            contact.lower_bound = 8
        with self.assertRaises(ValueError):
            contact.lower_bound = 10
        self.assertEqual(1, contact.lower_bound)

    def test_upper_bound_1(self):
        contact = Contact(1, 2, 1.0)
        self.assertEqual(8, contact.upper_bound)

    def test_upper_bound_2(self):
        contact = Contact(1, 2, 1.0)
        with self.assertRaises(ValueError):
            contact.upper_bound = -1

    def test_upper_bound_3(self):
        contact = Contact(1, 2, 1.0)
        contact.upper_bound = 10
        self.assertEqual(10, contact.upper_bound)

    def test_upper_bound_4(self):
        contact = Contact(1, 2, 1.0)
        contact.upper_bound = 7
        self.assertEqual(7, contact.upper_bound)

    def test_upper_bound_5(self):
        contact = Contact(1, 2, 1.0)
        contact.lower_bound = 4
        contact.upper_bound = 8
        self.assertEqual(4, contact.lower_bound)
        with self.assertRaises(ValueError):
            contact.upper_bound = 4
        with self.assertRaises(ValueError):
            contact.upper_bound = 3
        self.assertEqual(8, contact.upper_bound)

    def test_raw_score_1(self):
        contact = Contact(1, 2, 1.0)
        self.assertEqual(1.0, contact.raw_score)
        contact.raw_score = 0.3
        self.assertEqual(0.3, contact.raw_score)
        contact.raw_score = -0.1
        self.assertEqual(-0.1, contact.raw_score)
        contact.raw_score = 1.0
        self.assertEqual(1.0, contact.raw_score)

    def test_res1_1(self):
        contact = Contact(1, 2, 1.0)
        contact.res1 = 'Ala'
        self.assertEqual('A', contact.res1)

    def test_res1_2(self):
        contact = Contact(1, 2, 1.0)
        contact.res1 = 'T'
        self.assertEqual('T', contact.res1)

    def test_res1_3(self):
        contact = Contact(1, 2, 1.0)
        with self.assertRaises(ValueError):
            contact.res1 = '8'

    def test_res2_1(self):
        contact = Contact(1, 2, 1.0)
        contact.res2 = 'Met'
        self.assertEqual('M', contact.res2)

    def test_res2_2(self):
        contact = Contact(1, 2, 1.0)
        contact.res2 = 'P'
        self.assertEqual('P', contact.res2)

    def test_res2_3(self):
        contact = Contact(1, 2, 1.0)
        with self.assertRaises(ValueError):
            contact.res2 = '?'

    def test_res1_chain_1(self):
        contact = Contact(1, 2, 1.0)
        self.assertEqual('', contact.res1_chain)
        contact.res1_chain = 'A'
        self.assertEqual('A', contact.res1_chain)
        contact.res1_chain = 'd'
        self.assertEqual('d', contact.res1_chain)

    def test_res2_chain_1(self):
        contact = Contact(1, 2, 1.0)
        self.assertEqual('', contact.res2_chain)
        contact.res2_chain = 'b'
        self.assertEqual('b', contact.res2_chain)
        contact.res2_chain = 'X'
        self.assertEqual('X', contact.res2_chain)

    def test_res1_seq_1(self):
        contact = Contact(1, 2, 1.0)
        self.assertEqual(1, contact.res1_seq)
        contact.res1_seq = 2
        self.assertEqual(2, contact.res1_seq)
        contact.res1_seq = 5
        self.assertEqual(5, contact.res1_seq)
        contact.res1_seq = 1
        self.assertEqual(1, contact.res1_seq)

    def test_res2_seq_1(self):
        contact = Contact(1, 2000000, 1.0)
        self.assertEqual(2000000, contact.res2_seq)
        contact.res2_seq = 1
        self.assertEqual(1, contact.res2_seq)
        contact.res2_seq = 2000
        self.assertEqual(2000, contact.res2_seq)
        contact.res2_seq = 5
        self.assertEqual(5, contact.res2_seq)

    def test_res1_altseq(self):
        # ======================================================
        # Test Case 1
        contact = Contact(1, 2000000, 1.0)
        contact.res1_altseq = 1000
        self.assertEqual(1000, contact.res1_altseq)
        self.assertNotEqual(10, contact.res2_altseq)

    def test_res2_altseq_1(self):
        contact = Contact(1, 2000000, 1.0)
        contact.res2_altseq = 10
        self.assertEqual(10, contact.res2_altseq)
        self.assertNotEqual(5, contact.res2_altseq)

    def test_scalar_score_1(self):
        contact = Contact(1, 2000000, 1.0)
        contact.scalar_score = 5.432
        self.assertEqual(5.432, contact.scalar_score)

    def test_scalar_score_2(self):
        contact = Contact(1, 2000000, 1.0)
        contact.scalar_score = 3
        self.assertNotEqual(3, contact.res2_altseq)

    def test_status_1(self):
        contact = Contact(1, 2000000, 1.0)
        self.assertEqual(0, contact.status)
        contact.define_false_positive()
        self.assertEqual(-1, contact.status)
        contact.define_true_positive()
        self.assertEqual(1, contact.status)

    def test_weight_1(self):
        contact = Contact(1, 2000000, 1.0)
        self.assertEqual(1.0, contact.weight)
        contact.weight = 2.5
        self.assertEqual(2.5, contact.weight)

    def test_define_false_positive_1(self):
        contact = Contact(1, 2, 1.0)
        contact.define_false_positive()
        self.assertTrue(contact.is_false_positive)
        self.assertFalse(contact.is_true_positive)

    def test_define_true_positive_1(self):
        contact = Contact(1, 2, 1.0)
        contact.define_true_positive()
        self.assertTrue(contact.is_true_positive)
        self.assertFalse(contact.is_false_positive)

    def test__to_dict_1(self):
        contact = Contact(1, 2, 1.0)
        dict = {
            'id': (1, 2), 'is_false_positive': False, 'is_true_positive': False, 'distance_bound': (0, 8),
            'lower_bound': 0, 'upper_bound': 8, 'raw_score': 1.0, 'res1': 'X', 'res2': 'X', 'res1_chain': '',
            'res2_chain': '', 'res1_seq': 1, 'res2_seq': 2, 'res1_altseq': 0, 'res2_altseq': 0, 'scalar_score': 0.0,
            'status': 0, 'weight': 1.0,
        }
        self.assertEqual(dict, contact._to_dict())

    def test__to_dict_2(self):
        contact = Contact(1, 2, 1.0)
        contact.define_true_positive()
        contact.lower_bound = 4
        dict = {
            'id': (1, 2), 'is_false_positive': False, 'is_true_positive': True, 'distance_bound': (4, 8),
            'lower_bound': 4, 'upper_bound': 8, 'raw_score': 1.0, 'res1': 'X', 'res2': 'X', 'res1_chain': '',
            'res2_chain': '', 'res1_seq': 1, 'res2_seq': 2, 'res1_altseq': 0, 'res2_altseq': 0, 'scalar_score': 0.0,
            'status': 1, 'weight': 1.0,
        }
        self.assertEqual(dict, contact._to_dict())

    def test__set_residue_1(self):
        self.assertEqual("A", Contact._set_residue("ALA"))
        self.assertEqual("A", Contact._set_residue("Ala"))
        self.assertEqual("A", Contact._set_residue("ala"))
        self.assertEqual("A", Contact._set_residue("A"))
        self.assertEqual("A", Contact._set_residue("a"))
        self.assertRaises(ValueError, Contact._set_residue, 'AL')
        self.assertRaises(ValueError, Contact._set_residue, '-')
        self.assertRaises(ValueError, Contact._set_residue, 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
