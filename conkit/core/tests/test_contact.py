"""Testing facility for conkit.core.Contact"""

__author__ = "Felix Simkovic"
__date__ = "12 Aug 2016"

import unittest

from conkit.core.contact import Contact
from conkit.core.mappings import ContactMatchState


class TestContact(unittest.TestCase):
    def test_distance_bound_1(self):
        contact = Contact(1, 2, 1.0)
        contact.distance_bound = (0, 8)
        self.assertEqual((0.0, 8.0), contact.distance_bound)

    def test_distance_bound_2(self):
        contact = Contact(1, 2, 1.0)
        contact.distance_bound = (0, 8)
        self.assertEqual((0.0, 8.0), contact.distance_bound)
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

    def test_is_match_1(self):
        contact = Contact(1, 2, 1.0)
        contact.define_match()
        self.assertTrue(contact.is_match)
        self.assertFalse(contact.is_mismatch)
        self.assertFalse(contact.is_unknown)

    def test_is_mismatch_1(self):
        contact = Contact(1, 2, 1.0)
        contact.define_mismatch()
        self.assertFalse(contact.is_match)
        self.assertTrue(contact.is_mismatch)
        self.assertFalse(contact.is_unknown)

    def test_is_unknown_1(self):
        contact = Contact(1, 2, 1.0)
        contact.define_unknown()
        self.assertFalse(contact.is_match)
        self.assertFalse(contact.is_mismatch)
        self.assertTrue(contact.is_unknown)

    def test_lower_bound_1(self):
        contact = Contact(1, 2, 1.0)
        self.assertEqual(0.0, contact.lower_bound)

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
        contact.res1 = 'ALA'
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
        contact.res2 = 'MET'
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

    def test_res1_altseq_1(self):
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
        self.assertEqual(ContactMatchState.unknown.value, contact.status)
        contact.define_match()
        self.assertEqual(ContactMatchState.matched.value, contact.status)

    def test_status_2(self):
        contact = Contact(1, 2000000, 1.0)
        self.assertEqual(ContactMatchState.unknown.value, contact.status)
        contact.define_mismatch()
        self.assertEqual(ContactMatchState.mismatched.value, contact.status)

    def test_status_3(self):
        contact = Contact(1, 2000000, 1.0)
        self.assertEqual(ContactMatchState.unknown.value, contact.status)
        contact.define_unknown()
        self.assertEqual(ContactMatchState.unknown.value, contact.status)

    def test_weight_1(self):
        contact = Contact(1, 2000000, 1.0)
        self.assertEqual(1.0, contact.weight)
        contact.weight = 2.5
        self.assertEqual(2.5, contact.weight)

    def test_define_match_1(self):
        contact = Contact(1, 2, 1.0)
        contact.define_match()
        self.assertTrue(contact.is_match)
        self.assertFalse(contact.is_mismatch)
        self.assertFalse(contact.is_unknown)

    def test_define_mismatch_1(self):
        contact = Contact(1, 2, 1.0)
        contact.define_mismatch()
        self.assertFalse(contact.is_match)
        self.assertTrue(contact.is_mismatch)
        self.assertFalse(contact.is_unknown)

    def test_define_unknown_1(self):
        contact = Contact(1, 2, 1.0)
        contact.define_unknown()
        self.assertFalse(contact.is_match)
        self.assertFalse(contact.is_mismatch)
        self.assertTrue(contact.is_unknown)

    def test__to_dict_1(self):
        contact = Contact(1, 2, 1.0)
        answer_dict = {
            'id': (1, 2),
            'is_match': False,
            'is_mismatch': False,
            'is_unknown': True,
            'distance_bound': (0.0, 8.0),
            'lower_bound': 0.0,
            'upper_bound': 8.0,
            'raw_score': 1.0,
            'res1': 'X',
            'res2': 'X',
            'res1_chain': '',
            'res2_chain': '',
            'res1_seq': 1,
            'res2_seq': 2,
            'res1_altseq': 0,
            'res2_altseq': 0,
            'scalar_score': 0.0,
            'status': 0,
            'weight': 1.0,
        }
        contact_dict = contact._to_dict()
        for k in answer_dict.keys():
            self.assertEqual(answer_dict[k], contact_dict[k], "Key %s differs" % k)

    def test__to_dict_2(self):
        contact = Contact(1, 2, 1.0)
        contact.define_match()
        contact.lower_bound = 4
        answer_dict = {
            'id': (1, 2),
            'is_match': True,
            'is_mismatch': False,
            'is_unknown': False,
            'distance_bound': (4.0, 8.0),
            'lower_bound': 4.0,
            'upper_bound': 8.0,
            'raw_score': 1.0,
            'res1': 'X',
            'res2': 'X',
            'res1_chain': '',
            'res2_chain': '',
            'res1_seq': 1,
            'res2_seq': 2,
            'res1_altseq': 0,
            'res2_altseq': 0,
            'scalar_score': 0.0,
            'status': 1,
            'weight': 1.0,
        }
        contact_dict = contact._to_dict()
        for k in answer_dict.keys():
            self.assertEqual(answer_dict[k], contact_dict[k], "Key %s differs" % k)

    def test__set_residue_1(self):
        self.assertEqual("A", Contact._set_residue("ALA"))

    def test__set_residue_2(self):
        self.assertEqual("A", Contact._set_residue("A"))

    def test__set_residue_3(self):
        self.assertEqual("A", Contact._set_residue("Ala"))

    def test__set_residue_4(self):
        self.assertEqual("A", Contact._set_residue("ala"))

    def test__set_residue_5(self):
        self.assertEqual("A", Contact._set_residue("a"))

    def test__set_residue_6(self):
        self.assertRaises(ValueError, Contact._set_residue, 'AL')

    def test__set_residue_7(self):
        self.assertRaises(ValueError, Contact._set_residue, '-')

    def test__set_residue_8(self):
        self.assertRaises(AttributeError, Contact._set_residue, 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
