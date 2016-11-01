"""Testing facility for conkit.core.contact"""

__author__ = "Felix Simkovic"
__date__ = "12 Aug 2016"

from conkit.core.Contact import Contact

import unittest


class Test(unittest.TestCase):

    def test_distance_bound(self):
        # ======================================================
        # Test Case 1
        contact = Contact(1, 2, 1.0)
        contact.distance_bound = (0, 8)
        self.assertEqual((0, 8), contact.distance_bound)
        # ======================================================
        # Test Case 2
        contact = Contact(1, 2, 1.0)
        contact.distance_bound = [0, 8]
        self.assertEqual((0, 8), contact.distance_bound)
        self.assertTrue(isinstance(contact.distance_bound, tuple))
        # ======================================================
        # Test Case 3
        contact = Contact(1, 2, 1.0)
        try:
            contact.distance_bound = 8
            self.assertTrue(False)
        except TypeError:
            self.assertTrue(True)
        try:
            contact.distance_bound = 'test'
            self.assertTrue(False)
        except TypeError:
            self.assertTrue(True)
        try:
            contact.distance_bound = (1, 100)
            self.assertTrue(True)
        except TypeError:
            self.assertTrue(False)

    def test_is_false_positive(self):
        # ======================================================
        # Test Case 1
        contact = Contact(1, 2, 1.0)
        contact.define_false_positive()
        self.assertTrue(contact.is_false_positive)
        self.assertFalse(contact.is_true_positive)

    def test_is_true_positive(self):
        # ======================================================
        # Test Case 1
        contact = Contact(1, 2, 1.0)
        contact.define_true_positive()
        self.assertTrue(contact.is_true_positive)
        self.assertFalse(contact.is_false_positive)

    def test_raw_score(self):
        # ======================================================
        # Test Case 1
        contact = Contact(1, 2, 1.0)
        self.assertEqual(1.0, contact.raw_score)
        contact.raw_score = 0.3
        self.assertEqual(0.3, contact.raw_score)
        contact.raw_score = -0.1
        self.assertEqual(-0.1, contact.raw_score)
        contact.raw_score = 1.0
        self.assertEqual(1.0, contact.raw_score)

    def test_res1(self):
        # ======================================================
        # Test Case 1
        contact = Contact(1, 2, 1.0)
        contact.res1 = 'Ala'
        self.assertEqual('A', contact.res1)
        # ======================================================
        # Test Case 2
        contact = Contact(1, 2, 1.0)
        contact.res1 = 'T'
        self.assertEqual('T', contact.res1)
        # ======================================================
        # Test Case 3
        contact = Contact(1, 2, 1.0)
        try:
            contact.res1 = '8'
            self.assertFalse(True)
        except ValueError:
            self.assertTrue(True)

    def test_res2(self):
        # ======================================================
        # Test Case 1
        contact = Contact(1, 2, 1.0)
        contact.res2 = 'Met'
        self.assertEqual('M', contact.res2)
        # ======================================================
        # Test Case 2
        contact = Contact(1, 2, 1.0)
        contact.res2 = 'P'
        self.assertEqual('P', contact.res2)
        # ======================================================
        # Test Case 3
        contact = Contact(1, 2, 1.0)
        try:
            contact.res2 = '?'
            self.assertFalse(True)
        except ValueError:
            self.assertTrue(True)

    def test_res1_chain(self):
        # ======================================================
        # Test Case 1
        contact = Contact(1, 2, 1.0)
        self.assertEqual('', contact.res1_chain)
        contact.res1_chain = 'A'
        self.assertEqual('A', contact.res1_chain)
        contact.res1_chain = 'd'
        self.assertEqual('d', contact.res1_chain)

    def test_res2_chain(self):
        # ======================================================
        # Test Case 1
        contact = Contact(1, 2, 1.0)
        self.assertEqual('', contact.res2_chain)
        contact.res2_chain = 'b'
        self.assertEqual('b', contact.res2_chain)
        contact.res2_chain = 'X'
        self.assertEqual('X', contact.res2_chain)

    def test_res1_seq(self):
        # ======================================================
        # Test Case 1
        contact = Contact(1, 2, 1.0)
        self.assertEqual(1, contact.res1_seq)
        contact.res1_seq = 2
        self.assertEqual(2, contact.res1_seq)
        contact.res1_seq = 5
        self.assertEqual(5, contact.res1_seq)
        contact.res1_seq = 1
        self.assertEqual(1, contact.res1_seq)

    def test_res2_seq(self):
        # ======================================================
        # Test Case 1
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

    def test_res2_altseq(self):
        # ======================================================
        # Test Case 1
        contact = Contact(1, 2000000, 1.0)
        contact.res2_altseq = 10
        self.assertEqual(10, contact.res2_altseq)
        self.assertNotEqual(5, contact.res2_altseq)

    def test_scalar_score(self):
        # ======================================================
        # Test Case 1
        contact = Contact(1, 2000000, 1.0)
        contact.scalar_score = 5.432
        self.assertEqual(5.432, contact.scalar_score)
        # ======================================================
        # Test Case 2
        contact = Contact(1, 2000000, 1.0)
        contact.scalar_score = 3
        self.assertNotEqual(3, contact.res2_altseq)

    def test_status(self):
        # ======================================================
        # Test Case 1
        contact = Contact(1, 2000000, 1.0)
        self.assertEqual(0, contact.status)
        contact.define_false_positive()
        self.assertEqual(-1, contact.status)
        contact.define_true_positive()
        self.assertEqual(1, contact.status)

    def test_weight(self):
        # ======================================================
        # Test Case 1
        contact = Contact(1, 2000000, 1.0)
        self.assertEqual(1.0, contact.weight)
        contact.weight = 2.5
        self.assertEqual(2.5, contact.weight)

    def test_define_false_positive(self):
        # ======================================================
        # Test Case 1
        contact = Contact(1, 2, 1.0)
        contact.define_false_positive()
        self.assertTrue(contact.is_false_positive)
        self.assertFalse(contact.is_true_positive)

    def test_define_true_positive(self):
        # ======================================================
        # Test Case 1
        contact = Contact(1, 2, 1.0)
        contact.define_true_positive()
        self.assertTrue(contact.is_true_positive)
        self.assertFalse(contact.is_false_positive)

    def test__set_residue(self):
        contact = Contact(0, 0, 0)
        self.assertEqual("A", contact._set_residue("ALA"))
        self.assertEqual("A", contact._set_residue("Ala"))
        self.assertEqual("A", contact._set_residue("ala"))
        self.assertEqual("A", contact._set_residue("A"))
        self.assertEqual("A", contact._set_residue("a"))
        self.assertRaises(ValueError, contact._set_residue, 'AL')
        self.assertRaises(ValueError, contact._set_residue, '-')
        self.assertRaises(ValueError, contact._set_residue, 1)


if __name__ == "__main__":
    unittest.main()
