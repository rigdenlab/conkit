"""Testing facility for conkit.core.ContactMap"""

__author__ = "Felix Simkovic"
__date__ = "12 Aug 2016"

import unittest

try:
    import sklearn.neighbors

    SKLEARN = True
except ImportError:
    SKLEARN = False

from conkit.core.struct import Gap, Residue
from conkit.core.contact import Contact
from conkit.core.contactmap import ContactMap
from conkit.core.mappings import ContactMatchState
from conkit.core.sequence import Sequence


def skipUnless(condition):
    if condition:
        return lambda x: x
    else:
        return lambda x: None


TP = ContactMatchState.true_positive.value
FP = ContactMatchState.false_positive.value
FN = ContactMatchState.false_negative.value
UNK = ContactMatchState.unknown.value


class TestContactMap(unittest.TestCase):
    def test_coverage_1(self):
        contact_map = ContactMap("test")
        contact_map.add(Contact(1, 4, 1.0))
        contact_map.add(Contact(2, 4, 1.0))
        contact_map.add(Contact(5, 8, 1.0))
        contact_map.add(Contact(3, 6, 1.0))
        contact_map.sequence = Sequence("TEST", "ACDEFGHK")
        self.assertEqual(0.875, contact_map.coverage)

    def test_coverage_2(self):
        contact_map = ContactMap("test")
        contact_map.add(Contact(1, 4, 1.0))
        contact_map.add(Contact(2, 4, 1.0))
        contact_map.add(Contact(5, 8, 1.0))
        contact_map.add(Contact(3, 6, 1.0))
        contact_map.sequence = Sequence("TEST", "ACDEFGHK")
        self.assertEqual(0.875, contact_map.coverage)

    def test_coverage_3(self):
        contact_map = ContactMap("test")
        contact_map.add(Contact(1, 4, 1.0))
        contact_map.add(Contact(2, 4, 1.0))
        contact_map.add(Contact(5, 8, 1.0))
        contact_map.add(Contact(3, 6, 1.0))
        contact_map.sequence = Sequence("TEST", "ACDEFGHK")
        contact_map.remove((5, 8))
        contact_map.remove((3, 6))
        self.assertEqual(0.375, contact_map.coverage)

    def test_empty_1(self):
        contact_map = ContactMap("test")
        self.assertTrue(contact_map.empty)

    def test_empty_2(self):
        contact_map = ContactMap("test")
        contact_map.add(Contact(1, 5, 1.0))
        self.assertFalse(contact_map.empty)

    def test_ncontacts_1(self):
        contact_map = ContactMap("test")
        self.assertEqual(0, contact_map.ncontacts)

    def test_ncontacts_2(self):
        contact_map = ContactMap("test")
        contact_map.add(Contact(1, 5, 1.0))
        self.assertEqual(1, contact_map.ncontacts)

    def test_ncontacts_3(self):
        contact_map = ContactMap("test")
        contact_map.add(Contact(1, 5, 1.0))
        contact_map.add(Contact(2, 10, 1.0))
        self.assertEqual(2, contact_map.ncontacts)

    def test_precision_1(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map.sequence = Sequence("TEST", "AAAAA")
        for i, contact in enumerate(contact_map):
            if i % 2 == 0:
                contact.status = TP
            else:
                contact.status = FP
        self.assertEqual(0.5, contact_map.precision)

    def test_precision_2(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2), Contact(1, 1, 0)]:
            contact_map.add(c)
        contact_map.sequence = Sequence("TEST", "AAAAA")
        for i, contact in enumerate(contact_map):
            if i % 2 == 0:
                contact.status = TP
            else:
                contact.status = FP
        contact_map[(1, 1)].status = FN
        self.assertEqual(0.5, contact_map.precision)

    def test_recall_1(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        for i, contact in enumerate(contact_map):
            if i % 2 == 0:
                contact.status = TP
            else:
                contact.status = FP
        for c in [Contact(2, 5, 1.0), Contact(3, 6, 1.0), Contact(1, 7, 1.0)]:
            c.status = FN
            contact_map.add(c)
        self.assertEqual(0.4, contact_map.recall)

    def test_repr_sequence_1(self):
        contact_map = ContactMap("test")
        for contact in [Contact(1, 5, 1.0), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        contact_map.sequence = Sequence("foo", "ABCDE")
        self.assertEqual("AB-DE", contact_map.repr_sequence.seq)
        contact_map.remove((2, 4))

    def test_repr_sequence_2(self):
        contact_map = ContactMap("test")
        for contact in [Contact(1, 5, 1.0), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        contact_map.sequence = Sequence("foo", "ABCDE")
        self.assertEqual("AB-DE", contact_map.repr_sequence.seq)
        contact_map.remove((2, 4))
        self.assertEqual("A---E", contact_map.repr_sequence.seq)

    def test_repr_sequence_altloc_1(self):
        contact_map = ContactMap("test")
        for contact, altseq in [
            (Contact(1, 5, 1.0), (3, 4)),
            (Contact(2, 4, 0.1), (1, 5)),
            (Contact(5, 1, 0.2), (4, 3)),
        ]:
            contact.res1_altseq = altseq[0]
            contact.res2_altseq = altseq[1]
            contact_map.add(contact)
        contact_map.sequence = Sequence("foo", "ABCDE")
        self.assertEqual("A-CDE", contact_map.repr_sequence_altloc.seq)
        self.assertEqual("AB-DE", contact_map.repr_sequence.seq)

    def test_repr_sequence_altloc_2(self):
        contact_map = ContactMap("test")
        for contact, altseq in [
            (Contact(1, 5, 1.0), (3, 4)),
            (Contact(2, 4, 0.1), (1, 5)),
            (Contact(5, 1, 0.2), (4, 3)),
        ]:
            contact.res1_altseq = altseq[0]
            contact.res2_altseq = altseq[1]
            contact_map.add(contact)
        contact_map.sequence = Sequence("foo", "ABCDE")
        self.assertEqual("A-CDE", contact_map.repr_sequence_altloc.seq)
        self.assertEqual("AB-DE", contact_map.repr_sequence.seq)
        contact_map.remove((2, 4))
        self.assertEqual("--CD-", contact_map.repr_sequence_altloc.seq)
        self.assertEqual("A---E", contact_map.repr_sequence.seq)

    def test_sequence_1(self):
        contact_map = ContactMap("test")
        for contact in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        seq_obj = Sequence("bar", "ABCDE")
        contact_map.sequence = seq_obj
        self.assertEqual(seq_obj, contact_map.sequence)

    def test_sequence_2(self):
        contact_map = ContactMap("test")
        for contact in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        sequence1 = Sequence("foo", "ABC")
        sequence2 = Sequence("bar", "DE")
        with self.assertRaises(TypeError):
            contact_map.sequence([sequence1, sequence2])

    def test_top_contact_1(self):
        contact_map = ContactMap("test")
        self.assertEqual(None, contact_map.top_contact)

    def test_top_contact_2(self):
        contact_map = ContactMap("test")
        contact = Contact(1, 10, 1.0)
        contact_map.add(contact)
        self.assertEqual(contact, contact_map.top_contact)

    def test_top_contact_3(self):
        contact_map = ContactMap("test")
        contact1 = Contact(1, 10, 1.0)
        contact2 = Contact(2, 100, 1.0)
        contact_map.add(contact1)
        contact_map.add(contact2)
        self.assertEqual(contact1, contact_map.top_contact)

    def test_highest_residue_number_1(self):
        contact_map = ContactMap("test")
        self.assertIsNone(contact_map.highest_residue_number)

    def test_highest_residue_number_2(self):
        contact_map = ContactMap("test")
        contact1 = Contact(1, 10, 1.0)
        contact2 = Contact(2, 100, 1.0)
        contact_map.add(contact1)
        contact_map.add(contact2)
        self.assertEqual(100, contact_map.highest_residue_number)

    def test__construct_repr_sequence_1(self):
        contact_map = ContactMap("test")
        for contact in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        seq_obj = Sequence("bar", "ABCDE")
        contact_map.sequence = seq_obj
        self.assertEqual("ABCDE", contact_map._construct_repr_sequence([1, 2, 3, 4, 5]).seq)

    def test__construct_repr_sequence_2(self):
        contact_map = ContactMap("test")
        for contact in [Contact(1, 5, 1.0), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        seq_obj = Sequence("bar", "ABCDE")
        contact_map.sequence = seq_obj
        self.assertEqual("AB-DE", contact_map._construct_repr_sequence([1, 2, 4, 5]).seq)

    def test__construct_repr_sequence_3(self):
        contact_map = ContactMap("test")
        for contact in [Contact(1, 5, 1.0), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        seq_obj = Sequence("bar", "ABCDE")
        contact_map.sequence = seq_obj
        self.assertEqual("A---E", contact_map._construct_repr_sequence([1, 5]).seq)

    def test__construct_repr_sequence_4(self):
        contact_map = ContactMap("test")
        for contact in [Contact(1, 5, 1.0), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        seq_obj = Sequence("bar", "ABCDE")
        contact_map.sequence = seq_obj
        self.assertEqual("----E", contact_map._construct_repr_sequence([-1, 5]).seq)

    def test_set_scalar_score_1(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map.set_scalar_score()
        self.assertListEqual([2.352941, 0.941176, 0.235294, 0.470588], [round(c.scalar_score, 6) for c in contact_map])

    def test_get_jaccard_index_1(self):
        contact_map1 = ContactMap("foo")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map1.add(c)
        contact_map2 = ContactMap("bar")
        for c in [Contact(1, 7, 1.0), Contact(3, 3, 0.4), Contact(2, 5, 0.1), Contact(5, 1, 0.2)]:
            contact_map2.add(c)
        jindex = contact_map1.get_jaccard_index(contact_map2)
        self.assertEqual(0.333333, round(jindex, 6))

    def test_get_jaccard_index_2(self):
        contact_map1 = ContactMap("foo")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map1.add(c)
        contact_map2 = ContactMap("bar")
        for c in [Contact(1, 7, 1.0), Contact(3, 2, 0.4), Contact(2, 5, 0.1), Contact(5, 2, 0.2)]:
            contact_map2.add(c)
        jindex = contact_map1.get_jaccard_index(contact_map2)
        self.assertEqual(0.0, jindex)

    def test_get_jaccard_index_3(self):
        contact_map1 = ContactMap("foo")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map1.add(c)
        contact_map2 = ContactMap("bar")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map2.add(c)
        jindex = contact_map1.get_jaccard_index(contact_map2)
        self.assertEqual(1.0, jindex)

    def test_get_jaccard_index_4(self):
        contact_map1 = ContactMap("foo")
        contact_map2 = ContactMap("bar")
        jindex = contact_map1.get_jaccard_index(contact_map2)
        self.assertEqual(1.0, jindex)

    @skipUnless(SKLEARN)
    def test_get_contact_density_1(self):
        contact_map1 = ContactMap("foo")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1)]:
            contact_map1.add(c)
        density = [round(x, 7) for x in contact_map1.get_contact_density()]
        answer = [0.1194466, 0.2012433, 0.2386849, 0.2012433, 0.1194466]
        self.assertListEqual(answer, density)

    @skipUnless(SKLEARN)
    def test_get_contact_density_2(self):
        contact_map1 = ContactMap("foo")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(3, 4, 0.4)]:
            contact_map1.add(c)
        density = [round(x, 7) for x in contact_map1.get_contact_density()]
        answer = [0.1001259, 0.1983717, 0.2684149, 0.2313211, 0.1217899]
        self.assertEqual(answer, density)

    @skipUnless(SKLEARN)
    def test_get_contact_density_3(self):
        contact_map1 = ContactMap("foo")
        for c in [Contact(3, 5, 0.4), Contact(2, 4, 0.1), Contact(3, 4, 0.4)]:
            contact_map1.add(c)
        density = [round(x, 7) for x in contact_map1.get_contact_density()]
        answer = [0.1442296, 0.4134216, 0.4134216, 0.1442296]
        self.assertEqual(answer, density)

    def test_find_1(self):
        contact_map1 = ContactMap("1")
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        found = contact_map1.find(1)
        self.assertEqual(2, len(found))
        self.assertEqual([1, 1], [c.res1_seq for c in found])
        self.assertEqual([5, 4], [c.res2_seq for c in found])

    def test_find_2(self):
        contact_map1 = ContactMap("1")
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        found = contact_map1.find([10])
        self.assertEqual(0, len(found))
        self.assertEqual([], [c.res1_seq for c in found])
        self.assertEqual([], [c.res2_seq for c in found])

    def test_find_3(self):
        contact_map1 = ContactMap("1")
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        found = contact_map1.find([1, 2, 3])
        self.assertEqual(5, len(found))
        self.assertEqual([1, 2, 1, 3, 2], [c.res1_seq for c in found])
        self.assertEqual([5, 6, 4, 6, 5], [c.res2_seq for c in found])

    def test_find_4(self):
        contact_map1 = ContactMap("1")
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        found = contact_map1.find([1, 6])
        self.assertEqual(4, len(found))
        self.assertEqual([1, 2, 1, 3], [c.res1_seq for c in found])
        self.assertEqual([5, 6, 4, 6], [c.res2_seq for c in found])

    def test_find_5(self):
        contact_map1 = ContactMap("1")
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        found = contact_map1.find([1, 5], strict=True)
        self.assertEqual(1, len(found))
        self.assertEqual([1], [c.res1_seq for c in found])
        self.assertEqual([5], [c.res2_seq for c in found])

    def test_find_6(self):
        contact_map1 = ContactMap("1")
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        found = contact_map1.find([1, 6], strict=True)
        self.assertEqual(0, len(found))
        self.assertEqual([], [c.res1_seq for c in found])
        self.assertEqual([], [c.res2_seq for c in found])

    def test_find_7(self):
        contact_map1 = ContactMap("1")
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        found = contact_map1.find([1, 2, 5], strict=True)
        self.assertEqual(2, len(found))
        self.assertEqual([1, 2], [c.res1_seq for c in found])
        self.assertEqual([5, 5], [c.res2_seq for c in found])

    def test_find_inverse(self):
        contact_map1 = ContactMap("1")
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        found = contact_map1.find([1, 2, 5], strict=False, inverse=True)
        self.assertEqual(1, len(found))
        self.assertEqual([3], [c.res1_seq for c in found])
        self.assertEqual([6], [c.res2_seq for c in found])

    def test_find_strict_inverse(self):
        contact_map1 = ContactMap("1")
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        found = contact_map1.find([1, 2, 5], strict=True, inverse=True)
        self.assertEqual(3, len(found))
        self.assertEqual([2, 1, 3], [c.res1_seq for c in found])
        self.assertEqual([6, 4, 6], [c.res2_seq for c in found])

    def test_set_sequence_register_1(self):
        contact_map1 = ContactMap("1")
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        contact_map1.sequence = Sequence("foo", "ABCDEF")
        contact_map1.set_sequence_register()
        self.assertEqual(
            [("A", "E"), ("B", "F"), ("A", "D"), ("C", "F"), ("B", "E")], [(c.res1, c.res2) for c in contact_map1]
        )

    def test_set_sequence_register_2(self):
        contact_map1 = ContactMap("1")
        for comb, alt in [
            ((101, 105, 1.0), (1, 5)),
            ((102, 106, 1.0), (2, 6)),
            ((101, 104, 1.0), (1, 4)),
            ((103, 106, 1.0), (3, 6)),
            ((102, 105, 1.0), (2, 5)),
        ]:
            c = Contact(*comb)
            c.res1_altseq = alt[0]
            c.res2_altseq = alt[1]
            contact_map1.add(c)
        contact_map1.sequence = Sequence("foo", "ABCDEF")
        contact_map1.set_sequence_register(altloc=True)
        self.assertEqual(
            [("A", "E"), ("B", "F"), ("A", "D"), ("C", "F"), ("B", "E")], [(c.res1, c.res2) for c in contact_map1]
        )

    def test_set_sequence_register_3(self):
        contact_map1 = ContactMap("1")
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        contact_map1.sequence = Sequence("foo", "ABCDE")
        with self.assertRaises(ValueError):
            contact_map1.set_sequence_register()

    def test_set_sequence_register_4(self):
        contact_map1 = ContactMap("1")
        for comb in [(1, 5, 1.0), (6, 3, 1.0), (1, 4, 1.0), (3, 2, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        contact_map1.sequence = Sequence("foo", "ABCDE")
        with self.assertRaises(ValueError):
            contact_map1.set_sequence_register()

    def test_match_1(self):
        contact_map1 = ContactMap("foo")
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence("foo", "AICDEFGH")
        contact_map1.set_sequence_register()

        contact_map2 = ContactMap("bar")
        for i, params in enumerate([(1, 5, 1.0), (1, 7, 1.0), (2, 7, 1.0), (3, 4, 1.0)]):
            contact = Contact(*params)
            contact.res1_altseq = params[0]
            contact.res2_altseq = params[1]
            contact.status = TP
            contact_map2.add(contact)
        contact_map2.sequence = Sequence("bar", "AICDEFG")
        contact_map2.set_sequence_register(altloc=True)

        contact_map1.match(contact_map2, inplace=True)
        self.assertEqual([TP, FP, TP, FP, UNK], [c.status for c in contact_map1])

    def test_match_2(self):
        contact_map1 = ContactMap("foo")
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence("foo", "AICDEFGH")
        contact_map1.set_sequence_register()

        contact_map2 = ContactMap("bar")
        for i, params in enumerate([(1, 5, 1.0), (1, 6, 1.0), (3, 5, 1.0)]):
            contact = Contact(*params)
            contact.res1_altseq = params[0]
            contact.res2_altseq = params[1]
            contact.status = TP
            contact_map2.add(contact)
        contact_map2.sequence = Sequence("bar", "AICDEFG")
        contact_map2.set_sequence_register(altloc=True)

        contact_map1.match(contact_map2, remove_unmatched=True, inplace=True)
        self.assertEqual([TP, TP, FP, TP], [c.status for c in contact_map1])

    def test_match_3(self):
        contact_map1 = ContactMap("foo")
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence("foo", "AICDEFGH")
        contact_map1.set_sequence_register()

        contact_map2 = ContactMap("bar")
        for i, params in enumerate([(2, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0)]):
            contact = Contact(*params)
            contact.res1_altseq = params[0]
            contact.res2_altseq = params[1]
            if i % 2 == 0:
                contact.status = TP
            else:
                contact.status = FP
            contact_map2.add(contact)
        contact_map2.sequence = Sequence("bar", "AICDEFG")
        contact_map2.set_sequence_register(altloc=True)

        contact_map1.match(contact_map2, remove_unmatched=True, inplace=True)
        self.assertEqual([(1, 5), (1, 6), (2, 7), (3, 5)], [c.id for c in contact_map1])
        self.assertEqual([FP, FP, TP, TP], [c.status for c in contact_map1])

    def test_match_4(self):
        contact_map1 = ContactMap("foo")
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence("foo", "AICDEFGH")
        contact_map1.set_sequence_register()

        contact_map2 = ContactMap("bar")
        for i, params in enumerate([(1, 5, 1.0), (1, 6, 1.0), (2, 4, 1.0)]):
            contact = Contact(*params)
            contact.res1_altseq = params[0]
            contact.res2_altseq = params[1]
            contact.status = TP
            contact_map2.add(contact)
        contact_map2.sequence = Sequence("bar", "ICDEFG")
        contact_map2.set_sequence_register(altloc=True)

        contact_map1.match(contact_map2, match_other=True, remove_unmatched=True, inplace=True)
        self.assertEqual([TP, TP], [c.status for c in contact_map1])
        self.assertEqual([2, 2, 3], [c.res1_altseq for c in contact_map2])
        self.assertEqual([6, 7, 5], [c.res2_altseq for c in contact_map2])

    def test_match_5(self):
        contact_map1 = ContactMap("foo")
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence("foo", "AICDEFGH")
        contact_map1.set_sequence_register()

        contact_map2 = ContactMap("bar")

        contact1 = Contact(95, 30, 1.0)
        contact1.res1_chain = "A"
        contact1.res2_chain = "B"
        contact1.res1_altseq = 1
        contact1.res2_altseq = 5
        contact1.status = TP
        contact_map2.add(contact1)

        contact2 = Contact(95, 31, 1.0)
        contact2.res1_chain = "A"
        contact2.res2_chain = "B"
        contact2.res1_altseq = 1
        contact2.res2_altseq = 6
        contact1.status = TP
        contact_map2.add(contact2)

        contact3 = Contact(97, 30, 1.0)
        contact3.res1_chain = "A"
        contact3.res2_chain = "B"
        contact3.res1_altseq = 3
        contact3.res2_altseq = 5
        contact3.status = TP
        contact_map2.add(contact3)

        contact_map2.sequence = Sequence("bar", "AICDEFG")
        contact_map2.set_sequence_register(altloc=True)

        contact_map1.match(contact_map2, inplace=True)
        self.assertEqual([TP, TP, FP, TP, UNK], [c.status for c in contact_map1])

    def test_match_6(self):
        contact_map1 = ContactMap("foo")
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence("foo", "AICDEFGH")
        contact_map1.set_sequence_register()

        contact_map2 = ContactMap("bar")

        contact1 = Contact(95, 30, 1.0)
        contact1.res1_chain = "A"
        contact1.res2_chain = "B"
        contact1.res1_altseq = 1
        contact1.res2_altseq = 5
        contact1.status = TP
        contact_map2.add(contact1)

        contact2 = Contact(95, 31, 1.0)
        contact2.res1_chain = "A"
        contact2.res2_chain = "B"
        contact2.res1_altseq = 1
        contact2.res2_altseq = 6
        contact2.status = TP
        contact_map2.add(contact2)

        contact3 = Contact(97, 30, 1.0)
        contact3.res1_chain = "A"
        contact3.res2_chain = "B"
        contact3.res1_altseq = 3
        contact3.res2_altseq = 5
        contact3.status = TP
        contact_map2.add(contact3)

        contact_map2.sequence = Sequence("bar", "AICDEFG")
        contact_map2.set_sequence_register(altloc=True)

        contact_map1.match(contact_map2, remove_unmatched=True, inplace=True)
        self.assertEqual([TP, TP, FP, TP], [c.status for c in contact_map1])

    def test_match_7(self):
        contact_map1 = ContactMap("foo")
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence("foo", "AICDEFGH")
        contact_map1.set_sequence_register()

        contact_map2 = ContactMap("bar")

        contact1 = Contact(95, 30, 1.0)
        contact1.res1_chain = "A"
        contact1.res2_chain = "B"
        contact1.res1_altseq = 1
        contact1.res2_altseq = 5
        contact1.status = TP
        contact_map2.add(contact1)

        contact2 = Contact(95, 31, 1.0)
        contact2.res1_chain = "A"
        contact2.res2_chain = "B"
        contact2.res1_altseq = 1
        contact2.res2_altseq = 6
        contact2.status = TP
        contact_map2.add(contact2)

        contact3 = Contact(97, 30, 1.0)
        contact3.res1_chain = "A"
        contact3.res2_chain = "B"
        contact3.res1_altseq = 3
        contact3.res2_altseq = 5
        contact3.status = TP
        contact_map2.add(contact3)

        contact_map2.sequence = Sequence("bar", "AICDEFG")
        contact_map2.set_sequence_register(altloc=True)

        contact_map1.match(contact_map2, renumber=True, inplace=True)
        self.assertEqual([TP, TP, FP, TP, UNK], [c.status for c in contact_map1])
        self.assertEqual([95, 95, Gap.IDENTIFIER, 97, Gap.IDENTIFIER], [c.res1_seq for c in contact_map1])
        self.assertEqual(["A", "A", "", "A", ""], [c.res1_chain for c in contact_map1])
        self.assertEqual([30, 31, Gap.IDENTIFIER, 30, Gap.IDENTIFIER], [c.res2_seq for c in contact_map1])
        self.assertEqual(["B", "B", "", "B", ""], [c.res2_chain for c in contact_map1])

    def test_match_8(self):
        contact_map1 = ContactMap("foo")
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence("foo", "AICDEFGH")
        contact_map1.set_sequence_register()

        contact_map2 = ContactMap("bar")

        contact1 = Contact(95, 30, 1.0)
        contact1.res1_chain = "A"
        contact1.res2_chain = "B"
        contact1.res1_altseq = 1
        contact1.res2_altseq = 5
        contact_map2.add(contact1)

        contact2 = Contact(95, 31, 1.0)
        contact2.res1_chain = "A"
        contact2.res2_chain = "B"
        contact2.res1_altseq = 1
        contact2.res2_altseq = 6
        contact_map2.add(contact2)

        contact3 = Contact(97, 30, 1.0)
        contact3.res1_chain = "A"
        contact3.res2_chain = "B"
        contact3.res1_altseq = 3
        contact3.res2_altseq = 5
        contact_map2.add(contact3)

        contact_map2.sequence = Sequence("bar", "AICDEFG")
        contact_map2.set_sequence_register(altloc=True)

        contact_map1.match(contact_map2, remove_unmatched=True, renumber=True, inplace=True)
        self.assertEqual([TP, TP, FP, TP], [c.status for c in contact_map1])
        self.assertEqual([95, 95, Gap.IDENTIFIER, 97], [c.res1_seq for c in contact_map1])
        self.assertEqual(["A", "A", "", "A"], [c.res1_chain for c in contact_map1])
        self.assertEqual([30, 31, Gap.IDENTIFIER, 30], [c.res2_seq for c in contact_map1])
        self.assertEqual(["B", "B", "", "B"], [c.res2_chain for c in contact_map1])

    def test_match_9(self):
        contact_map1 = ContactMap("foo")
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence("foo", "AICDEFGH")
        contact_map1.set_sequence_register()

        contact_map2 = ContactMap("bar")

        contact1 = Contact(95, 30, 1.0)
        contact1.res1_chain = "A"
        contact1.res2_chain = "B"
        contact1.res1_altseq = 1
        contact1.res2_altseq = 5
        contact_map2.add(contact1)

        contact2 = Contact(95, 31, 1.0)
        contact2.res1_chain = "A"
        contact2.res2_chain = "B"
        contact2.res1_altseq = 1
        contact2.res2_altseq = 6
        contact_map2.add(contact2)

        contact3 = Contact(96, 32, 1.0)
        contact3.res1_chain = "A"
        contact3.res2_chain = "B"
        contact3.res1_altseq = 2
        contact3.res2_altseq = 7
        contact_map2.add(contact3)

        contact4 = Contact(97, 30, 1.0)
        contact4.res1_chain = "A"
        contact4.res2_chain = "B"
        contact4.res1_altseq = 3
        contact4.res2_altseq = 5
        contact_map2.add(contact4)

        contact5 = Contact(96, 33, 1.0)
        contact5.res1_chain = "A"
        contact5.res2_chain = "B"
        contact5.res1_altseq = 2
        contact5.res2_altseq = 8
        contact_map2.add(contact5)

        contact_map2.sequence = Sequence("bar", "AICDEFGF")
        contact_map2.set_sequence_register(altloc=True)

        contact_map1.match(contact_map2, renumber=True, inplace=True)
        self.assertEqual([TP, TP, TP, TP, TP], [c.status for c in contact_map1])
        self.assertEqual([95, 95, 96, 97, 96], [c.res1_seq for c in contact_map1])
        self.assertEqual(["A", "A", "A", "A", "A"], [c.res1_chain for c in contact_map1])
        self.assertEqual([30, 31, 32, 30, 33], [c.res2_seq for c in contact_map1])
        self.assertEqual(["B", "B", "B", "B", "B"], [c.res2_chain for c in contact_map1])

    def test_match_10(self):
        contact_map1 = ContactMap("foo")
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence("foo", "AICDEFGH")
        contact_map1.set_sequence_register()

        # # Encryption matrix
        # seq: A B C D E F G H
        # res: 6 7 8 1 2 3 4 5
        # alt: 1 2 3 4 5 6 7 8
        # cha: A A A B B B B B

        contact_map2 = ContactMap("bar")
        contact1 = Contact(6, 2, 1.0)
        contact1.res1_chain = "A"
        contact1.res2_chain = "B"
        contact1.res1_altseq = 1
        contact1.res2_altseq = 5
        contact_map2.add(contact1)

        contact2 = Contact(7, 4, 1.0)
        contact2.res1_chain = "A"
        contact2.res2_chain = "B"
        contact2.res1_altseq = 2
        contact2.res2_altseq = 7
        contact2.status = TP
        contact_map2.add(contact2)

        contact3 = Contact(8, 2, 1.0)
        contact3.res1_chain = "A"
        contact3.res2_chain = "B"
        contact3.res1_altseq = 3
        contact3.res2_altseq = 5
        contact3.status = TP
        contact_map2.add(contact3)

        contact_map2.sequence = Sequence("bar", "AICDEFGH")
        contact_map2.set_sequence_register(altloc=True)

        contact_map1.match(contact_map2, renumber=True, inplace=True)
        self.assertEqual([TP, FP, TP, TP, FP], [c.status for c in contact_map1])
        self.assertEqual(
            [(6, 2), (6, Gap.IDENTIFIER), (7, 4), (8, 2), (7, Gap.IDENTIFIER)],
            [(c.res1_seq, c.res2_seq) for c in contact_map1],
        )
        self.assertEqual(
            [("A", "B"), ("A", ""), ("A", "B"), ("A", "B"), ("A", "")],
            [(c.res1_chain, c.res2_chain) for c in contact_map1],
        )

    def test_match_11(self):
        reference = ContactMap("ref")
        for cparam in [(1, 5, 1.0), (1, 7, 1.0), (2, 7, 1.0), (3, 4, 1.0)]:
            contact = Contact(*cparam)
            contact.res1_altseq = cparam[0]
            contact.res2_altseq = cparam[1]
            contact.status = TP
            reference.add(contact)
        reference.sequence = Sequence("bar", "AIDEFGH")
        reference.set_sequence_register(altloc=True)

        sequence = Sequence("foo", "AICDEFGH")
        for i in range(2):
            target = ContactMap("tar")
            for cparam in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
                target.add(Contact(*cparam))
            target.remove_neighbors(min_distance=1, inplace=True)
            target.sort("raw_score", reverse=True, inplace=True)

            target.sequence = sequence
            target.set_sequence_register()
            target.id = "foobar_%d" % i

            target.match(reference, match_other=False, inplace=True)

    def test_match_12(self):
        contact_map1 = ContactMap("foo")
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence("foo", "AICDEFGH")
        contact_map1.set_sequence_register()

        contact_map2 = ContactMap("bar")
        for i, params in enumerate([(1, 5, 1.0), (1, 7, 1.0), (2, 7, 1.0), (3, 4, 1.0)]):
            contact = Contact(*params)
            contact.res1_altseq = params[0]
            contact.res2_altseq = params[1]
            contact.status = TP
            contact_map2.add(contact)
        contact_map2.sequence = Sequence("bar", "AICDEFG")
        contact_map2.set_sequence_register(altloc=True)

        contact_map1.match(contact_map2, add_false_negatives=True, inplace=True)
        self.assertEqual([TP, FP, TP, FP, UNK, FN, FN], [c.status for c in contact_map1])
        self.assertListEqual([[1, 5], [1, 6], [2, 7], [3, 5], [2, 8], [1, 7], [3, 4]], contact_map1.as_list())

    def test_match_13(self):
        contact_map1 = ContactMap("foo")
        for i, params in enumerate([(1, 5, 1.0), (1, 7, 1.0), (2, 7, 1.0), (3, 4, 1.0)]):
            contact = Contact(*params)
            contact.res1_altseq = params[0]
            contact.res2_altseq = params[1]
            contact_map1.add(contact)
        contact_map1.sequence = Sequence("foo", "AICDEFGH")
        contact_map1.set_sequence_register()

        contact_map2 = ContactMap("bar")
        contact_map2.sequence = Sequence("bar", "AICDEFGH")
        contact_map2.set_sequence_register()

        contact_map1.match(contact_map2, add_false_negatives=True, inplace=True)
        self.assertEqual([FP, FP, FP, FP], [c.status for c in contact_map1])
        self.assertListEqual([[1, 5], [1, 7], [2, 7], [3, 4]], contact_map1.as_list())

    def test_match_14(self):
        contact_map1 = ContactMap("foo")
        contact_map1.sequence = Sequence("foo", "AICDEFGH")
        contact_map1.set_sequence_register()

        contact_map2 = ContactMap("bar")
        for i, params in enumerate([(1, 5, 1.0), (1, 7, 1.0), (2, 7, 1.0), (3, 4, 1.0)]):
            contact = Contact(*params)
            contact.res1_altseq = params[0]
            contact.res2_altseq = params[1]
            contact.status = TP
            contact_map2.add(contact)
        contact_map2.sequence = Sequence("bar", "AICDEFGH")
        contact_map2.set_sequence_register()

        contact_map1.match(contact_map2, add_false_negatives=True, inplace=True, match_other=True)
        self.assertEqual([FN, FN, FN, FN], [c.status for c in contact_map2])
        self.assertListEqual([[1, 5], [1, 7], [2, 7], [3, 4]], contact_map2.as_list())

    def test_remove_neighbors_1(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_mod = contact_map.remove_neighbors(min_distance=2)
        self.assertListEqual([(1, 5), (2, 4), (5, 1)], [c.id for c in contact_map_mod])
        self.assertEqual([(1, 5), (2, 4), (5, 1)], sorted(contact_map_mod.child_dict.keys()))

    def test_remove_neighbors_2(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_mod = contact_map.remove_neighbors(min_distance=5, inplace=True)
        self.assertListEqual([], [c.get_id() for c in contact_map_mod])
        self.assertDictEqual({}, contact_map_mod.child_dict)
        self.assertEqual(contact_map, contact_map_mod)

    def test_remove_neighbors_3(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 30, 1.0), Contact(2, 10, 0.4), Contact(3, 20, 0.1), Contact(1, 5, 0.2)]:
            contact_map.add(c)
        contact_map_mod = contact_map.remove_neighbors(min_distance=6, max_distance=12, inplace=True)
        self.assertEqual([(2, 10)], [c.id for c in contact_map_mod])
        self.assertEqual([(2, 10)], sorted(contact_map_mod.child_dict.keys()))

    def test_remove_neighbors_4(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 30, 1.0), Contact(2, 10, 0.4), Contact(3, 20, 0.1), Contact(1, 5, 0.2)]:
            contact_map.add(c)
        contact_map_mod = contact_map.remove_neighbors(min_distance=12, max_distance=24, inplace=True)
        self.assertEqual([(3, 20)], [c.id for c in contact_map_mod])
        self.assertEqual([(3, 20)], sorted(contact_map_mod.child_dict.keys()))

    def test_remove_neighbors_5(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 30, 1.0), Contact(2, 10, 0.4), Contact(3, 20, 0.1), Contact(1, 5, 0.2)]:
            contact_map.add(c)
        contact_map_mod = contact_map.remove_neighbors(min_distance=24, inplace=True)
        self.assertEqual([(1, 30)], [c.id for c in contact_map_mod])
        self.assertEqual([(1, 30)], sorted(contact_map_mod.child_dict.keys()))

    def test_remove_neighbors_6(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 30, 1.0), Contact(2, 10, 0.4), Contact(3, 20, 0.1), Contact(1, 5, 0.2)]:
            contact_map.add(c)
        contact_map_mod = contact_map.remove_neighbors(min_distance=6, max_distance=24, inplace=True)
        self.assertEqual([(2, 10), (3, 20)], [c.id for c in contact_map_mod])
        self.assertEqual([(2, 10), (3, 20)], sorted(contact_map_mod.child_dict.keys()))

    def test_filter_1(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_mod = contact_map.filter(threshold=0.5)
        self.assertListEqual([(1, 5)], [c.id for c in contact_map_mod])
        self.assertEqual([(1, 5)], sorted(contact_map_mod.child_dict.keys()))

    def test_filter_2(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_mod = contact_map.filter(threshold=0.3, inplace=True)
        self.assertListEqual([(1, 5), (3, 3)], [c.id for c in contact_map])
        self.assertEqual([(1, 5), (3, 3)], sorted(contact_map.child_dict.keys()))

    def test_filter_3(self):
        contact_map = ContactMap("test")
        for c in [
            Contact(1, 5, 1.0, (0, 8)),
            Contact(3, 3, 0.4, (0, 10)),
            Contact(2, 4, 0.1, (5, 7)),
            Contact(5, 1, 0.2, (3, 10)),
        ]:
            contact_map.add(c)
        contact_map_mod = contact_map.filter(threshold=3, inplace=True, filter_by="lower_bound")
        self.assertListEqual([(2, 4), (5, 1)], [c.id for c in contact_map])
        self.assertEqual([(2, 4), (5, 1)], sorted(contact_map.child_dict.keys()))

    def test_rescale_1(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_rescaled = contact_map.rescale()
        self.assertListEqual([1.0, 0.333, 0.0, 0.111], [round(c.raw_score, 3) for c in contact_map_rescaled])
        self.assertNotEqual(contact_map, contact_map_rescaled)

    def test_rescale_2(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_rescaled = contact_map.rescale(inplace=True)
        self.assertListEqual([1.0, 0.333, 0.0, 0.111], [round(c.raw_score, 3) for c in contact_map])
        self.assertEqual(contact_map, contact_map_rescaled)

    def test_sort_1(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_sorted = contact_map.sort("res1_seq", reverse=False, inplace=False)
        self.assertEqual([(1, 5), (2, 4), (3, 3), (5, 1)], [c.id for c in contact_map_sorted])
        self.assertNotEqual(contact_map, contact_map_sorted)

    def test_sort_2(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_sorted = contact_map.sort("res1_seq", reverse=True, inplace=False)
        self.assertEqual([(5, 1), (3, 3), (2, 4), (1, 5)], [c.id for c in contact_map_sorted])
        self.assertNotEqual(contact_map, contact_map_sorted)

    def test_sort_3(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_sorted = contact_map.sort("raw_score", reverse=False, inplace=True)
        self.assertEqual([(2, 4), (5, 1), (3, 3), (1, 5)], [c.id for c in contact_map_sorted])
        self.assertEqual(contact_map, contact_map_sorted)

    def test_sort_4(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_sorted = contact_map.sort("res2_seq", reverse=True, inplace=True)
        self.assertEqual([(1, 5), (2, 4), (3, 3), (5, 1)], [c.id for c in contact_map_sorted])
        self.assertEqual(contact_map, contact_map_sorted)

    def test__adjust_1(self):
        contact_map = ContactMap("test")
        for contact in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(1, 1, 0.2)]:
            contact_map.add(contact)
            contact.res1_altseq = contact.res1_seq
            contact.res2_altseq = contact.res2_seq
        contact_map_keymap = ContactMap._create_keymap(contact_map)
        sequence = [ord(x) for x in "XXXXX"]
        contact_map_keymap = ContactMap._insert_states(sequence, contact_map_keymap)
        contact_map_keymap = ContactMap._reindex_by_keymap(contact_map_keymap)
        adjusted = ContactMap._adjust(contact_map, contact_map_keymap)
        self.assertEqual([1, 3, 2, 1], [c.res1_altseq for c in adjusted])
        self.assertEqual([5, 3, 4, 1], [c.res2_altseq for c in adjusted])

    def test__adjust_2(self):
        contact_map = ContactMap("test")
        for contact in [Contact(1, 5, 1.0), Contact(2, 4, 0.1), Contact(1, 1, 0.2)]:
            contact_map.add(contact)
            contact.res1_altseq = contact.res1_seq + 3
            contact.res2_altseq = contact.res2_seq + 3
        contact_map_keymap = ContactMap._create_keymap(contact_map)
        sequence = [ord(x) for x in "XX-XX"]
        contact_map_keymap = ContactMap._insert_states(sequence, contact_map_keymap)
        contact_map_keymap = ContactMap._reindex_by_keymap(contact_map_keymap)
        adjusted = ContactMap._adjust(contact_map, contact_map_keymap)
        self.assertEqual([1, 2, 1], [c.res1_altseq for c in adjusted])
        self.assertEqual([5, 4, 1], [c.res2_altseq for c in adjusted])

    def test__create_keymap_1(self):
        contact_map = ContactMap("test")
        for contact in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(1, 1, 0.2)]:
            contact_map.add(contact)
        contact_map.sequence = Sequence("foo", "ABCDE")
        contact_map.set_sequence_register()
        contact_map_keymap = ContactMap._create_keymap(contact_map)
        self.assertEqual([1, 2, 3, 4, 5], [c.res_seq for c in contact_map_keymap])
        self.assertEqual(["A", "B", "C", "D", "E"], [c.res_name for c in contact_map_keymap])
        self.assertEqual([0, 0, 0, 0, 0], [c.res_altseq for c in contact_map_keymap])

    def test__create_keymap_2(self):
        contact_map = ContactMap("test")
        for contact, res_altloc in [
            (Contact(1, 5, 1.0), (10, 20)),
            (Contact(3, 3, 0.4), (12, 12)),
            (Contact(2, 4, 0.1), (11, 13)),
            (Contact(1, 1, 0.2), (10, 10)),
        ]:
            contact.res1_altseq = res_altloc[0]
            contact.res2_altseq = res_altloc[1]
            contact_map.add(contact)
        contact_map.sequence = Sequence("foo", "ABCDE")
        contact_map.set_sequence_register()
        contact_map_keymap = ContactMap._create_keymap(contact_map, altloc=True)
        self.assertEqual([1, 2, 3, 4, 5], [c.res_seq for c in contact_map_keymap])
        self.assertEqual(["A", "B", "C", "D", "E"], [c.res_name for c in contact_map_keymap])
        self.assertEqual([10, 11, 12, 13, 20], [c.res_altseq for c in contact_map_keymap])

    def test__create_keymap_3(self):
        contact_map = ContactMap("test")
        for contact, res_altloc in [
            (Contact(1, 5, 1.0), (10, 20)),
            (Contact(2, 4, 0.1), (11, 13)),
            (Contact(1, 1, 0.2), (10, 10)),
        ]:
            contact.res1_altseq = res_altloc[0]
            contact.res2_altseq = res_altloc[1]
            contact_map.add(contact)
        contact_map.sequence = Sequence("foo", "ABCDE")
        contact_map.set_sequence_register()
        contact_map_keymap = ContactMap._create_keymap(contact_map, altloc=True)
        self.assertEqual([1, 2, 4, 5], [c.res_seq for c in contact_map_keymap])
        self.assertEqual(["A", "B", "D", "E"], [c.res_name for c in contact_map_keymap])
        self.assertEqual([10, 11, 13, 20], [c.res_altseq for c in contact_map_keymap])

    def test__find_single_1(self):
        contact_map = ContactMap("test")
        for contact, res_altloc in [
            (Contact(1, 5, 1.0), (10, 20)),
            (Contact(2, 4, 0.1), (11, 13)),
            (Contact(1, 1, 0.2), (10, 10)),
        ]:
            contact.res1_altseq = res_altloc[0]
            contact.res2_altseq = res_altloc[1]
            contact_map.add(contact)
        found_contacts = ContactMap._find_single(contact_map, 1)
        self.assertEqual([(1, 5), (1, 1)], [c.id for c in found_contacts])
        found_contacts = ContactMap._find_single(contact_map, 2)
        self.assertEqual([(2, 4)], [c.id for c in found_contacts])
        found_contacts = ContactMap._find_single(contact_map, 8)
        self.assertEqual([], list(found_contacts))

    def test__find_single_2(self):
        contact_map = ContactMap("test")
        for contact, res_altloc in [
            (Contact(1, 5, 1.0), (10, 20)),
            (Contact(2, 4, 0.1), (11, 13)),
            (Contact(1, 1, 0.2), (10, 10)),
        ]:
            contact.res1_altseq = res_altloc[0]
            contact.res2_altseq = res_altloc[1]
            contact_map.add(contact)
        found_contacts = ContactMap._find_single(contact_map, 1)
        self.assertEqual([(1, 5), (1, 1)], [c.id for c in found_contacts])
        contact_map[0].res1_altseq = 4
        contact_map[0].res2_altseq = 1
        found_contacts = ContactMap._find_single(contact_map, 1)
        self.assertEqual([(1, 5), (1, 1)], [c.id for c in found_contacts])
        found_contacts = ContactMap._find_single(contact_map, 4)
        self.assertEqual([(2, 4)], [c.id for c in found_contacts])

    def test__insert_states_1(self):
        keymap = [Residue(1, 1, "X", ""), Residue(2, 2, "X", ""), Residue(3, 3, "X", "")]
        sequence = [ord(x) for x in "XXX"]
        inserts_added = ContactMap._insert_states(sequence, keymap)
        self.assertEqual(3, len(inserts_added))
        self.assertEqual([True, True, True], [isinstance(r, Residue) for r in inserts_added])

    def test__insert_states_2(self):
        keymap = [Residue(2, 2, "X", ""), Residue(3, 3, "X", "")]
        sequence = [ord(x) for x in "-XX"]
        inserts_added = ContactMap._insert_states(sequence, keymap)
        self.assertEqual(3, len(inserts_added))
        self.assertEqual([False, True, True], [isinstance(r, Residue) for r in inserts_added])
        self.assertEqual([True, False, False], [isinstance(r, Gap) for r in inserts_added])

    def test__insert_states_3(self):
        keymap = [Residue(2, 10, "X", ""), Residue(3, 11, "X", "")]
        sequence = [ord(x) for x in "-X-X--"]
        inserts_added = ContactMap._insert_states(sequence, keymap)
        self.assertEqual(6, len(inserts_added))
        self.assertEqual([False, True, False, True, False, False], [isinstance(r, Residue) for r in inserts_added])
        self.assertEqual([True, False, True, False, True, True], [isinstance(r, Gap) for r in inserts_added])
        self.assertEqual([2, 3], [r.res_seq for r in inserts_added if isinstance(r, Residue)])
        self.assertEqual([10, 11], [r.res_altseq for r in inserts_added if isinstance(r, Residue)])

    def test__reindex_by_keymap_1(self):
        keymap = [Residue(1, 1, "X", ""), Residue(2, 2, "X", ""), Residue(3, 3, "X", "")]
        reindex = ContactMap._reindex_by_keymap(keymap)
        self.assertEqual([1, 2, 3], [c.res_seq for c in reindex])
        self.assertEqual([1, 2, 3], [c.res_altseq for c in reindex])

    def test__reindex_by_keymap_2(self):
        keymap = [Residue(1, -5, "X", ""), Residue(2, -4, "X", ""), Residue(3, -3, "X", "")]
        reindex = ContactMap._reindex_by_keymap(keymap)
        self.assertEqual([1, 2, 3], [c.res_seq for c in reindex])
        self.assertEqual([1, 2, 3], [c.res_altseq for c in reindex])

    def test__reindex_by_keymap_3(self):
        keymap = [Gap(), Residue(2, 1, "X", ""), Residue(3, 2, "X", "")]
        reindex = ContactMap._reindex_by_keymap(keymap)
        self.assertEqual([Gap.IDENTIFIER, 2, 3], [c.res_seq for c in reindex])
        self.assertEqual([1, 2, 3], [c.res_altseq for c in reindex])

    def test__reindex_by_keymap_4(self):
        keymap = [Gap(), Residue(200000, 10000, "X", ""), Gap(), Gap()]
        reindex = ContactMap._reindex_by_keymap(keymap)
        self.assertEqual([Gap.IDENTIFIER, 200000, Gap.IDENTIFIER, Gap.IDENTIFIER], [c.res_seq for c in reindex])
        self.assertEqual([1, 2, 3, 4], [c.res_altseq for c in reindex])

    def test__reindex_by_keymap_5(self):
        keymap = [Gap(), Gap(), Gap(), Gap()]
        reindex = ContactMap._reindex_by_keymap(keymap)
        self.assertEqual([Gap.IDENTIFIER, Gap.IDENTIFIER, Gap.IDENTIFIER, Gap.IDENTIFIER], [c.res_seq for c in reindex])
        self.assertEqual([1, 2, 3, 4], [c.res_altseq for c in reindex])

    def test_as_list_1(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 30, 1.0), Contact(2, 10, 0.4), Contact(3, 20, 0.1), Contact(1, 5, 0.2)]:
            contact_map.add(c)
        self.assertListEqual([[1, 30], [2, 10], [3, 20], [1, 5]], contact_map.as_list())

    def test_as_list_2(self):
        contact_map = ContactMap("test")
        self.assertListEqual([], contact_map.as_list())

    def test_as_list_3(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 30, 1.0), Contact(2, 10, 0.4), Contact(3, 20, 0.1), Contact(1, 5, 0.2)]:
            c.res1_altseq = c.res1_seq + 1
            c.res2_altseq = c.res2_seq + 2
            contact_map.add(c)
        self.assertListEqual([[2, 32], [3, 12], [4, 22], [2, 7]], contact_map.as_list(altloc=True))

    def test_as_list_4(self):
        contact_map = ContactMap("test")
        self.assertListEqual([], contact_map.as_list(altloc=True))

    def test_as_list_5(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 30, 1.0), Contact(2, 10, 0.4), Contact(3, 20, 0.1), Contact(1, 5, 0.2)]:
            contact_map.add(c)
        self.assertListEqual([[0, 0], [0, 0], [0, 0], [0, 0]], contact_map.as_list(altloc=True))

    def test_reindex_1(self):
        contact_map = ContactMap("test")
        for c in [Contact(5, 30, 1.0), Contact(6, 10, 0.4), Contact(7, 20, 0.1), Contact(5, 10, 0.2)]:
            contact_map.add(c)
        reindexed = contact_map.reindex(1)
        self.assertListEqual([[1, 26], [2, 6], [3, 16], [1, 6]], contact_map.reindex(1).as_list())

    def test_reindex_2(self):
        contact_map = ContactMap("test")
        for c in [Contact(5, 30, 1.0), Contact(6, 10, 0.4), Contact(7, 20, 0.1), Contact(5, 10, 0.2)]:
            contact_map.add(c)
        self.assertListEqual([[2, 27], [3, 7], [4, 17], [2, 7]], contact_map.reindex(2).as_list())

    def test_reindex_3(self):
        contact_map = ContactMap("test")
        for c in [Contact(5, 30, 1.0), Contact(6, 10, 0.4), Contact(7, 20, 0.1), Contact(5, 10, 0.2)]:
            contact_map.add(c)
        with self.assertRaises(ValueError):
            contact_map.reindex(-1)

    def test_reindex_4(self):
        contact_map = ContactMap("test")
        for c in [Contact(5, 30, 1.0), Contact(6, 10, 0.4), Contact(7, 20, 0.1), Contact(5, 10, 0.2)]:
            c.res1_altseq = c.res1_seq + 1
            c.res2_altseq = c.res2_seq + 2
            contact_map.add(c)
        reindexed = contact_map.reindex(1, altloc=True)
        self.assertListEqual([[5, 30], [6, 10], [7, 20], [5, 10]], reindexed.as_list())
        self.assertListEqual([[1, 27], [2, 7], [3, 17], [1, 7]], reindexed.as_list(altloc=True))

    def test_singletons_1(self):
        contact_map = ContactMap("test")
        for c in [Contact(5, 5, 0.4), Contact(4, 6, 0.1), Contact(3, 5, 0.2)]:
            contact_map.add(c)
        self.assertListEqual([], contact_map.singletons.as_list())

    def test_singletons_2(self):
        contact_map = ContactMap("test")
        for c in [Contact(5, 5, 0.4), Contact(4, 6, 0.1), Contact(10, 10, 0.2)]:
            contact_map.add(c)
        self.assertListEqual([[10, 10]], contact_map.singletons.as_list())

    def test_singletons_3(self):
        contact_map = ContactMap("test")
        for c in [Contact(5, 5, 0.4), Contact(4, 6, 0.1), Contact(10, 10, 0.2), Contact(10, 11, 0.2)]:
            contact_map.add(c)
        self.assertListEqual([], contact_map.singletons.as_list())

    def test_singletons_4(self):
        contact_map = ContactMap("test")
        for c in [Contact(4, 5, 1.0), Contact(4, 6, 0.4)]:
            contact_map.add(c)
        self.assertListEqual([], contact_map.singletons.as_list())

    def test_singletons_5(self):
        contact_map = ContactMap("test")
        for c in [Contact(3, 4, 1.0), Contact(4, 8, 0.4)]:
            contact_map.add(c)
        self.assertListEqual([[3, 4], [4, 8]], contact_map.singletons.as_list())

    def test_singletons_6(self):
        contact_map = ContactMap("test")
        for c in [Contact(4, 5, 1.0), Contact(6, 7, 0.4)]:
            contact_map.add(c)
        self.assertListEqual([], contact_map.singletons.as_list())

    def test_singletons_7(self):
        contact_map = ContactMap("test")
        for c in [Contact(4, 5, 1.0), Contact(7, 8, 0.4)]:
            contact_map.add(c)
        self.assertListEqual([[4, 5], [7, 8]], contact_map.singletons.as_list())

    def test_remove_false_negatives_1(self):
        contact_map = ContactMap("foo")
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map.add(contact)
        contact_map[(1, 5)].false_negative = True
        contact_map[(3, 5)].false_negative = True
        contact_map.remove_false_negatives(inplace=True)
        self.assertListEqual([[1, 6], [2, 7], [2, 8]], contact_map.as_list())

    def test_remove_false_negatives_2(self):
        contact_map = ContactMap("foo")
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact.false_negative = True
            contact_map.add(contact)
        contact_map.remove_false_negatives(inplace=True)
        self.assertListEqual([], contact_map.as_list())

    def test_remove_false_negatives_3(self):
        contact_map = ContactMap("foo")
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map.add(contact)
        contact_map.remove_false_negatives(inplace=True)
        self.assertListEqual([[1, 5], [1, 6], [2, 7], [3, 5], [2, 8]], contact_map.as_list())

    def test_slice_1(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2), Contact(1, 1, 0)]:
            contact_map.add(c)
        contact_map.sequence = Sequence("TEST", "AAAAAA")
        contact_map.slice_map(0.5, inplace=True)
        self.assertListEqual([[1, 5], [3, 3], [5, 1]], contact_map.as_list())

    def test_slice_2(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2), Contact(1, 1, 0)]:
            contact_map.add(c)
        contact_map.slice_map(0.25, seq_len=12, inplace=True)
        self.assertListEqual([[1, 5], [3, 3], [5, 1]], contact_map.as_list())

    def test_slice_3(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2), Contact(1, 1, 0)]:
            contact_map.add(c)
        with self.assertRaises(ValueError):
            contact_map.slice_map(0.5, inplace=True)

    def test_as_dict(self):
        contact_map = ContactMap("test")
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2), Contact(1, 1, 0)]:
            contact_map.add(c)

        expected = {1: {(1, 1), (5, 1), (1, 5)}, 2: {(2, 4)}, 3: {(3, 3)}, 4: {(2, 4)}, 5: {(5, 1), (1, 5)}}
        output = contact_map.as_dict()
        self.assertDictEqual(expected, output)


if __name__ == "__main__":
    unittest.main(verbosity=2)
