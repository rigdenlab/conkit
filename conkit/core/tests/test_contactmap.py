"""Testing facility for conkit.core.ContactMap"""

__author__ = "Felix Simkovic"
__date__ = "12 Aug 2016"

import unittest

try:
    import sklearn.neighbors
    SKLEARN = True
except ImportError:
    SKLEARN = False

from conkit.core._struct import _Gap, _Residue
from conkit.core.contact import Contact, ContactMatchState
from conkit.core.contactmap import ContactMap
from conkit.core.sequence import Sequence


def skipUnless(condition):
    if condition:
        return lambda x: x
    else:
        return lambda x: None


class TestContactMap(unittest.TestCase):
    def test_coverage_1(self):
        contact_map = ContactMap('test')
        contact_map.add(Contact(1, 4, 1.0))
        contact_map.add(Contact(2, 4, 1.0))
        contact_map.add(Contact(5, 8, 1.0))
        contact_map.add(Contact(3, 6, 1.0))
        contact_map.sequence = Sequence('TEST', 'ABCDEFGH')
        self.assertEqual(0.875, contact_map.coverage)

    def test_coverage_2(self):
        contact_map = ContactMap('test')
        contact_map.add(Contact(1, 4, 1.0))
        contact_map.add(Contact(2, 4, 1.0))
        contact_map.add(Contact(5, 8, 1.0))
        contact_map.add(Contact(3, 6, 1.0))
        contact_map.sequence = Sequence('TEST', 'ABCDEFGH')
        self.assertEqual(0.875, contact_map.coverage)
        contact_map.remove((5, 8))
        contact_map.remove((3, 6))
        self.assertEqual(0.375, contact_map.coverage)

    def test_empty_1(self):
        contact_map = ContactMap("test")
        self.assertTrue(contact_map.empty)

    def test_empty_2(self):
        contact_map = ContactMap('test')
        contact_map.add(Contact(1, 5, 1.0))
        self.assertFalse(contact_map.empty)

    def test_ncontacts_1(self):
        contact_map = ContactMap('test')
        self.assertEqual(0, contact_map.ncontacts)

    def test_ncontacts_2(self):
        contact_map = ContactMap('test')
        contact_map.add(Contact(1, 5, 1.0))
        self.assertEqual(1, contact_map.ncontacts)

    def test_ncontacts_3(self):
        contact_map = ContactMap('test')
        contact_map.add(Contact(1, 5, 1.0))
        contact_map.add(Contact(2, 10, 1.0))
        self.assertEqual(2, contact_map.ncontacts)

    def test_precision_1(self):
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map.sequence = Sequence('TEST', 'AAAAA')
        for i, contact in enumerate(contact_map):
            if i % 2 == 0:
                contact.define_match()
            else:
                contact.define_mismatch()
        self.assertEqual(0.5, contact_map.precision)

    def test_repr_sequence_1(self):
        contact_map = ContactMap('test')
        for contact in [Contact(1, 5, 1.0), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        contact_map.sequence = Sequence('foo', 'ABCDE')
        self.assertEqual('AB-DE', contact_map.repr_sequence.seq)
        contact_map.remove((2, 4))

    def test_repr_sequence_2(self):
        contact_map = ContactMap('test')
        for contact in [Contact(1, 5, 1.0), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        contact_map.sequence = Sequence('foo', 'ABCDE')
        self.assertEqual('AB-DE', contact_map.repr_sequence.seq)
        contact_map.remove((2, 4))
        self.assertEqual('A---E', contact_map.repr_sequence.seq)

    def test_repr_sequence_altloc_1(self):
        contact_map = ContactMap('test')
        for contact, altseq in [(Contact(1, 5, 1.0), (3, 4)), (Contact(2, 4, 0.1), (1, 5)), (Contact(5, 1, 0.2), (4,
                                                                                                                  3))]:
            contact.res1_altseq = altseq[0]
            contact.res2_altseq = altseq[1]
            contact_map.add(contact)
        contact_map.sequence = Sequence('foo', 'ABCDE')
        self.assertEqual('A-CDE', contact_map.repr_sequence_altloc.seq)
        self.assertEqual('AB-DE', contact_map.repr_sequence.seq)

    def test_repr_sequence_altloc_2(self):
        contact_map = ContactMap('test')
        for contact, altseq in [(Contact(1, 5, 1.0), (3, 4)), (Contact(2, 4, 0.1), (1, 5)), (Contact(5, 1, 0.2), (4,
                                                                                                                  3))]:
            contact.res1_altseq = altseq[0]
            contact.res2_altseq = altseq[1]
            contact_map.add(contact)
        contact_map.sequence = Sequence('foo', 'ABCDE')
        self.assertEqual('A-CDE', contact_map.repr_sequence_altloc.seq)
        self.assertEqual('AB-DE', contact_map.repr_sequence.seq)
        contact_map.remove((2, 4))
        self.assertEqual('--CD-', contact_map.repr_sequence_altloc.seq)
        self.assertEqual('A---E', contact_map.repr_sequence.seq)

    def test_sequence_1(self):
        contact_map = ContactMap('test')
        for contact in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        seq_obj = Sequence('bar', 'ABCDE')
        contact_map.sequence = seq_obj
        self.assertEqual(seq_obj, contact_map.sequence)

    def test_sequence_2(self):
        contact_map = ContactMap('test')
        for contact in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        sequence1 = Sequence('foo', 'ABC')
        sequence2 = Sequence('bar', 'DE')
        with self.assertRaises(TypeError):
            contact_map.sequence([sequence1, sequence2])

    def test_top_contact_1(self):
        contact_map = ContactMap('test')
        self.assertEqual(None, contact_map.top_contact)

    def test_top_contact_2(self):
        contact_map = ContactMap('test')
        contact = Contact(1, 10, 1.0)
        contact_map.add(contact)
        self.assertEqual(contact, contact_map.top_contact)

    def test_top_contact_3(self):
        contact_map = ContactMap('test')
        contact1 = Contact(1, 10, 1.0)
        contact2 = Contact(2, 100, 1.0)
        contact_map.add(contact1)
        contact_map.add(contact2)
        self.assertEqual(contact1, contact_map.top_contact)

    def test__construct_repr_sequence_1(self):
        contact_map = ContactMap('test')
        for contact in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        seq_obj = Sequence('bar', 'ABCDE')
        contact_map.sequence = seq_obj
        self.assertEqual('ABCDE', contact_map._construct_repr_sequence([1, 2, 3, 4, 5]).seq)

    def test__construct_repr_sequence_2(self):
        contact_map = ContactMap('test')
        for contact in [Contact(1, 5, 1.0), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        seq_obj = Sequence('bar', 'ABCDE')
        contact_map.sequence = seq_obj
        self.assertEqual('AB-DE', contact_map._construct_repr_sequence([1, 2, 4, 5]).seq)

    def test__construct_repr_sequence_3(self):
        contact_map = ContactMap('test')
        for contact in [Contact(1, 5, 1.0), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        seq_obj = Sequence('bar', 'ABCDE')
        contact_map.sequence = seq_obj
        self.assertEqual('A---E', contact_map._construct_repr_sequence([1, 5]).seq)

    def test__construct_repr_sequence_4(self):
        contact_map = ContactMap('test')
        for contact in [Contact(1, 5, 1.0), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        seq_obj = Sequence('bar', 'ABCDE')
        contact_map.sequence = seq_obj
        self.assertEqual('----E', contact_map._construct_repr_sequence([-1, 5]).seq)

    def test_calculate_scalar_score_1(self):
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map.calculate_scalar_score()
        self.assertListEqual([2.352941, 0.941176, 0.235294, 0.470588], [round(c.scalar_score, 6) for c in contact_map])

    def test_calculate_jaccard_score_1(self):
        contact_map1 = ContactMap('foo')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map1.add(c)
        contact_map2 = ContactMap('bar')
        for c in [Contact(1, 7, 1.0), Contact(3, 3, 0.4), Contact(2, 5, 0.1), Contact(5, 1, 0.2)]:
            contact_map2.add(c)
        jindex = contact_map1.calculate_jaccard_index(contact_map2)
        self.assertEqual(0.333333, round(jindex, 6))

    def test_calculate_jaccard_score_2(self):
        contact_map1 = ContactMap('foo')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map1.add(c)
        contact_map2 = ContactMap('bar')
        for c in [Contact(1, 7, 1.0), Contact(3, 2, 0.4), Contact(2, 5, 0.1), Contact(5, 2, 0.2)]:
            contact_map2.add(c)
        jindex = contact_map1.calculate_jaccard_index(contact_map2)
        self.assertEqual(0.0, jindex)

    def test_calculate_jaccard_score_3(self):
        contact_map1 = ContactMap('foo')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map1.add(c)
        contact_map2 = ContactMap('bar')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map2.add(c)
        jindex = contact_map1.calculate_jaccard_index(contact_map2)
        self.assertEqual(1.0, jindex)

    def test_calculate_jaccard_score_4(self):
        contact_map1 = ContactMap('foo')
        contact_map2 = ContactMap('bar')
        jindex = contact_map1.calculate_jaccard_index(contact_map2)
        self.assertEqual(1.0, jindex)

    # REM: Bug in Sadowski res2(+1) error results in tests to fail.
    # REM: For now, code reverted by @biochunan to match Sadowski's implementation
    # REM: Corresponding results for tests 1-3 are below when code is reverted to correct implementation
    # REM: [0.11944664492151888, 0.20124328187327287, 0.23868489948868032, 0.20124328187327284, 0.11944664492151894]
    # REM: [0.10012592274179978, 0.19837169506444377, 0.26841488628499205, 0.23132109484153407, 0.1217899243718511]
    # REM: [0.14422964774244934, 0.4134215939843202, 0.41342159398432, 0.14422964774244929]

    @skipUnless(SKLEARN)
    def test_calculate_contact_density_1(self):
        contact_map1 = ContactMap('foo')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1)]:
            contact_map1.add(c)
        density = contact_map1.calculate_contact_density()
        self.assertEqual([0.16474084813765252, 0.24964700349773153, 0.24964700349773153, 0.16474084813765252], density)

    @skipUnless(SKLEARN)
    def test_calculate_contact_density_2(self):
        contact_map1 = ContactMap('foo')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(3, 4, 0.4)]:
            contact_map1.add(c)
        density = contact_map1.calculate_contact_density()
        self.assertEqual([0.14936186609839505, 0.26262219502283973, 0.2868154978059642, 0.1768745614827941], density)

    @skipUnless(SKLEARN)
    def test_calculate_contact_density_3(self):
        contact_map1 = ContactMap('foo')
        for c in [Contact(3, 5, 0.4), Contact(2, 4, 0.1), Contact(3, 4, 0.4)]:
            contact_map1.add(c)
        density = contact_map1.calculate_contact_density()
        self.assertEqual([0.2282135747320191, 0.6337372073257503, 0.2282135747320191], density)

    def test_find_1(self):
        contact_map1 = ContactMap('1')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        found = contact_map1.find(1)
        self.assertEqual(2, len(found))
        self.assertEqual([1, 1], [c.res1_seq for c in found])
        self.assertEqual([5, 4], [c.res2_seq for c in found])

    def test_find_2(self):
        contact_map1 = ContactMap('1')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        found = contact_map1.find([10])
        self.assertEqual(0, len(found))
        self.assertEqual([], [c.res1_seq for c in found])
        self.assertEqual([], [c.res2_seq for c in found])

    def test_find_3(self):
        contact_map1 = ContactMap('1')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        found = contact_map1.find([1, 2, 3])
        self.assertEqual(5, len(found))
        self.assertEqual([1, 2, 1, 3, 2], [c.res1_seq for c in found])
        self.assertEqual([5, 6, 4, 6, 5], [c.res2_seq for c in found])

    def test_find_4(self):
        contact_map1 = ContactMap('1')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        found = contact_map1.find([1, 6])
        self.assertEqual(4, len(found))
        self.assertEqual([1, 2, 1, 3], [c.res1_seq for c in found])
        self.assertEqual([5, 6, 4, 6], [c.res2_seq for c in found])

    def test_find_5(self):
        contact_map1 = ContactMap('1')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        found = contact_map1.find([1, 5], strict=True)
        self.assertEqual(1, len(found))
        self.assertEqual([1], [c.res1_seq for c in found])
        self.assertEqual([5], [c.res2_seq for c in found])

    def test_find_6(self):
        contact_map1 = ContactMap('1')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        found = contact_map1.find([1, 6], strict=True)
        self.assertEqual(0, len(found))
        self.assertEqual([], [c.res1_seq for c in found])
        self.assertEqual([], [c.res2_seq for c in found])

    def test_find_7(self):
        contact_map1 = ContactMap('1')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        found = contact_map1.find([1, 2, 5], strict=True)
        self.assertEqual(2, len(found))
        self.assertEqual([1, 2], [c.res1_seq for c in found])
        self.assertEqual([5, 5], [c.res2_seq for c in found])

    def test_assign_sequence_register_1(self):
        contact_map1 = ContactMap('1')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        contact_map1.sequence = Sequence('foo', 'ABCDEF')
        contact_map1.assign_sequence_register()
        self.assertEqual([('A', 'E'), ('B', 'F'), ('A', 'D'), ('C', 'F'), ('B', 'E')],
                         [(c.res1, c.res2) for c in contact_map1])

    def test_assign_sequence_register_2(self):
        contact_map1 = ContactMap('1')
        for comb, alt in [((101, 105, 1.0), (1, 5)), ((102, 106, 1.0), (2, 6)), ((101, 104, 1.0), (1, 4)),
                          ((103, 106, 1.0), (3, 6)), ((102, 105, 1.0), (2, 5))]:
            c = Contact(*comb)
            c.res1_altseq = alt[0]
            c.res2_altseq = alt[1]
            contact_map1.add(c)
        contact_map1.sequence = Sequence('foo', 'ABCDEF')
        contact_map1.assign_sequence_register(altloc=True)
        self.assertEqual([('A', 'E'), ('B', 'F'), ('A', 'D'), ('C', 'F'), ('B', 'E')],
                         [(c.res1, c.res2) for c in contact_map1])

    def test_match_1(self):
        contact_map1 = ContactMap('foo')
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence('foo', 'ABCDEFGH')
        contact_map1.assign_sequence_register()

        contact_map2 = ContactMap('bar')
        for i, params in enumerate([(1, 5, 1.0), (1, 7, 1.0), (2, 7, 1.0), (3, 4, 1.0)]):
            contact = Contact(*params)
            contact.res1_altseq = params[0]
            contact.res2_altseq = params[1]
            contact.define_match()
            contact_map2.add(contact)
        contact_map2.sequence = Sequence('bar', 'ABCDEFG')
        contact_map2.assign_sequence_register(altloc=True)

        contact_map1.match(contact_map2, inplace=True)
        self.assertEqual([
            ContactMatchState.matched.value, ContactMatchState.mismatched.value, ContactMatchState.matched.value,
            ContactMatchState.mismatched.value, ContactMatchState.unknown.value
        ], [c.status for c in contact_map1])

    def test_match_2(self):
        contact_map1 = ContactMap('foo')
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence('foo', 'ABCDEFGH')
        contact_map1.assign_sequence_register()

        contact_map2 = ContactMap('bar')
        for i, params in enumerate([(1, 5, 1.0), (1, 6, 1.0), (3, 5, 1.0)]):
            contact = Contact(*params)
            contact.res1_altseq = params[0]
            contact.res2_altseq = params[1]
            contact.define_match()
            contact_map2.add(contact)
        contact_map2.sequence = Sequence('bar', 'ABCDEFG')
        contact_map2.assign_sequence_register(altloc=True)

        contact_map1.match(contact_map2, remove_unmatched=True, inplace=True)
        self.assertEqual([
            ContactMatchState.matched.value, ContactMatchState.matched.value, ContactMatchState.mismatched.value,
            ContactMatchState.matched.value
        ], [c.status for c in contact_map1])

    def test_match_3(self):
        contact_map1 = ContactMap('foo')
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence('foo', 'ABCDEFGH')
        contact_map1.assign_sequence_register()

        contact_map2 = ContactMap('bar')
        for i, params in enumerate([(2, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0)]):
            contact = Contact(*params)
            contact.res1_altseq = params[0]
            contact.res2_altseq = params[1]
            if i % 2 == 0:
                contact.define_match()
            else:
                contact.define_mismatch()
            contact_map2.add(contact)
        contact_map2.sequence = Sequence('bar', 'ABCDEFG')
        contact_map2.assign_sequence_register(altloc=True)

        contact_map1.match(contact_map2, remove_unmatched=True, inplace=True)
        self.assertEqual([(1, 5), (1, 6), (2, 7), (3, 5)], [c.id for c in contact_map1])
        self.assertEqual([
            ContactMatchState.mismatched.value, ContactMatchState.mismatched.value, ContactMatchState.matched.value,
            ContactMatchState.matched.value
        ], [c.status for c in contact_map1])

    def test_match_4(self):
        contact_map1 = ContactMap('foo')
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence('foo', 'ABCDEFGH')
        contact_map1.assign_sequence_register()

        contact_map2 = ContactMap('bar')
        for i, params in enumerate([(1, 5, 1.0), (1, 6, 1.0), (2, 4, 1.0)]):
            contact = Contact(*params)
            contact.res1_altseq = params[0]
            contact.res2_altseq = params[1]
            contact.define_match()
            contact_map2.add(contact)
        contact_map2.sequence = Sequence('bar', 'BCDEFG')
        contact_map2.assign_sequence_register(altloc=True)

        contact_map1.match(contact_map2, match_other=True, remove_unmatched=True, inplace=True)
        self.assertEqual([ContactMatchState.matched.value, ContactMatchState.matched.value],
                         [c.status for c in contact_map1])
        self.assertEqual([2, 2, 3], [c.res1_altseq for c in contact_map2])
        self.assertEqual([6, 7, 5], [c.res2_altseq for c in contact_map2])

    def test_match_5(self):
        contact_map1 = ContactMap('foo')
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence('foo', 'ABCDEFGH')
        contact_map1.assign_sequence_register()

        contact_map2 = ContactMap('bar')

        contact1 = Contact(95, 30, 1.0)
        contact1.res1_chain = 'A'
        contact1.res2_chain = 'B'
        contact1.res1_altseq = 1
        contact1.res2_altseq = 5
        contact1.define_match()
        contact_map2.add(contact1)

        contact2 = Contact(95, 31, 1.0)
        contact2.res1_chain = 'A'
        contact2.res2_chain = 'B'
        contact2.res1_altseq = 1
        contact2.res2_altseq = 6
        contact2.define_match()
        contact_map2.add(contact2)

        contact3 = Contact(97, 30, 1.0)
        contact3.res1_chain = 'A'
        contact3.res2_chain = 'B'
        contact3.res1_altseq = 3
        contact3.res2_altseq = 5
        contact3.define_match()
        contact_map2.add(contact3)

        contact_map2.sequence = Sequence('bar', 'ABCDEFG')
        contact_map2.assign_sequence_register(altloc=True)

        contact_map1.match(contact_map2, inplace=True)
        self.assertEqual([
            ContactMatchState.matched.value, ContactMatchState.matched.value, ContactMatchState.mismatched.value,
            ContactMatchState.matched.value, ContactMatchState.unknown.value
        ], [c.status for c in contact_map1])

    def test_match_6(self):
        contact_map1 = ContactMap('foo')
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence('foo', 'ABCDEFGH')
        contact_map1.assign_sequence_register()

        contact_map2 = ContactMap('bar')

        contact1 = Contact(95, 30, 1.0)
        contact1.res1_chain = 'A'
        contact1.res2_chain = 'B'
        contact1.res1_altseq = 1
        contact1.res2_altseq = 5
        contact1.define_match()
        contact_map2.add(contact1)

        contact2 = Contact(95, 31, 1.0)
        contact2.res1_chain = 'A'
        contact2.res2_chain = 'B'
        contact2.res1_altseq = 1
        contact2.res2_altseq = 6
        contact2.define_match()
        contact_map2.add(contact2)

        contact3 = Contact(97, 30, 1.0)
        contact3.res1_chain = 'A'
        contact3.res2_chain = 'B'
        contact3.res1_altseq = 3
        contact3.res2_altseq = 5
        contact3.define_match()
        contact_map2.add(contact3)

        contact_map2.sequence = Sequence('bar', 'ABCDEFG')
        contact_map2.assign_sequence_register(altloc=True)

        contact_map1.match(contact_map2, remove_unmatched=True, inplace=True)
        self.assertEqual([
            ContactMatchState.matched.value, ContactMatchState.matched.value, ContactMatchState.mismatched.value,
            ContactMatchState.matched.value
        ], [c.status for c in contact_map1])

    def test_match_7(self):
        contact_map1 = ContactMap('foo')
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence('foo', 'ABCDEFGH')
        contact_map1.assign_sequence_register()

        contact_map2 = ContactMap('bar')

        contact1 = Contact(95, 30, 1.0)
        contact1.res1_chain = 'A'
        contact1.res2_chain = 'B'
        contact1.res1_altseq = 1
        contact1.res2_altseq = 5
        contact1.define_match()
        contact_map2.add(contact1)

        contact2 = Contact(95, 31, 1.0)
        contact2.res1_chain = 'A'
        contact2.res2_chain = 'B'
        contact2.res1_altseq = 1
        contact2.res2_altseq = 6
        contact2.define_match()
        contact_map2.add(contact2)

        contact3 = Contact(97, 30, 1.0)
        contact3.res1_chain = 'A'
        contact3.res2_chain = 'B'
        contact3.res1_altseq = 3
        contact3.res2_altseq = 5
        contact3.define_match()
        contact_map2.add(contact3)

        contact_map2.sequence = Sequence('bar', 'ABCDEFG')
        contact_map2.assign_sequence_register(altloc=True)

        contact_map1.match(contact_map2, renumber=True, inplace=True)
        self.assertEqual([
            ContactMatchState.matched.value, ContactMatchState.matched.value, ContactMatchState.mismatched.value,
            ContactMatchState.matched.value, ContactMatchState.unknown.value
        ], [c.status for c in contact_map1])
        self.assertEqual([95, 95, _Gap.IDENTIFIER, 97, _Gap.IDENTIFIER], [c.res1_seq for c in contact_map1])
        self.assertEqual(['A', 'A', '', 'A', ''], [c.res1_chain for c in contact_map1])
        self.assertEqual([30, 31, _Gap.IDENTIFIER, 30, _Gap.IDENTIFIER], [c.res2_seq for c in contact_map1])
        self.assertEqual(['B', 'B', '', 'B', ''], [c.res2_chain for c in contact_map1])

    def test_match_8(self):
        contact_map1 = ContactMap('foo')
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence('foo', 'ABCDEFGH')
        contact_map1.assign_sequence_register()

        contact_map2 = ContactMap('bar')

        contact1 = Contact(95, 30, 1.0)
        contact1.res1_chain = 'A'
        contact1.res2_chain = 'B'
        contact1.res1_altseq = 1
        contact1.res2_altseq = 5
        contact_map2.add(contact1)

        contact2 = Contact(95, 31, 1.0)
        contact2.res1_chain = 'A'
        contact2.res2_chain = 'B'
        contact2.res1_altseq = 1
        contact2.res2_altseq = 6
        contact_map2.add(contact2)

        contact3 = Contact(97, 30, 1.0)
        contact3.res1_chain = 'A'
        contact3.res2_chain = 'B'
        contact3.res1_altseq = 3
        contact3.res2_altseq = 5
        contact_map2.add(contact3)

        contact_map2.sequence = Sequence('bar', 'ABCDEFG')
        contact_map2.assign_sequence_register(altloc=True)

        contact_map1.match(contact_map2, remove_unmatched=True, renumber=True, inplace=True)
        self.assertEqual([
            ContactMatchState.matched.value, ContactMatchState.matched.value, ContactMatchState.mismatched.value,
            ContactMatchState.matched.value
        ], [c.status for c in contact_map1])
        self.assertEqual([95, 95, _Gap.IDENTIFIER, 97], [c.res1_seq for c in contact_map1])
        self.assertEqual(['A', 'A', '', 'A'], [c.res1_chain for c in contact_map1])
        self.assertEqual([30, 31, _Gap.IDENTIFIER, 30], [c.res2_seq for c in contact_map1])
        self.assertEqual(['B', 'B', '', 'B'], [c.res2_chain for c in contact_map1])

    def test_match_9(self):
        contact_map1 = ContactMap('foo')
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence('foo', 'ABCDEFGH')
        contact_map1.assign_sequence_register()

        contact_map2 = ContactMap('bar')

        contact1 = Contact(95, 30, 1.0)
        contact1.res1_chain = 'A'
        contact1.res2_chain = 'B'
        contact1.res1_altseq = 1
        contact1.res2_altseq = 5
        contact_map2.add(contact1)

        contact2 = Contact(95, 31, 1.0)
        contact2.res1_chain = 'A'
        contact2.res2_chain = 'B'
        contact2.res1_altseq = 1
        contact2.res2_altseq = 6
        contact_map2.add(contact2)

        contact3 = Contact(96, 32, 1.0)
        contact3.res1_chain = 'A'
        contact3.res2_chain = 'B'
        contact3.res1_altseq = 2
        contact3.res2_altseq = 7
        contact_map2.add(contact3)

        contact4 = Contact(97, 30, 1.0)
        contact4.res1_chain = 'A'
        contact4.res2_chain = 'B'
        contact4.res1_altseq = 3
        contact4.res2_altseq = 5
        contact_map2.add(contact4)

        contact5 = Contact(96, 33, 1.0)
        contact5.res1_chain = 'A'
        contact5.res2_chain = 'B'
        contact5.res1_altseq = 2
        contact5.res2_altseq = 8
        contact_map2.add(contact5)

        contact_map2.sequence = Sequence('bar', 'ABCDEFGF')
        contact_map2.assign_sequence_register(altloc=True)

        contact_map1.match(contact_map2, renumber=True, inplace=True)
        self.assertEqual([
            ContactMatchState.matched.value, ContactMatchState.matched.value, ContactMatchState.matched.value,
            ContactMatchState.matched.value, ContactMatchState.matched.value
        ], [c.status for c in contact_map1])
        self.assertEqual([95, 95, 96, 97, 96], [c.res1_seq for c in contact_map1])
        self.assertEqual(['A', 'A', 'A', 'A', 'A'], [c.res1_chain for c in contact_map1])
        self.assertEqual([30, 31, 32, 30, 33], [c.res2_seq for c in contact_map1])
        self.assertEqual(['B', 'B', 'B', 'B', 'B'], [c.res2_chain for c in contact_map1])

    def test_match_10(self):
        contact_map1 = ContactMap('foo')
        for params in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
            contact = Contact(*params)
            contact_map1.add(contact)
        contact_map1.sequence = Sequence('foo', 'ABCDEFGH')
        contact_map1.assign_sequence_register()

        ## Encryption matrix
        # seq: A B C D E F G H
        # res: 6 7 8 1 2 3 4 5
        # alt: 1 2 3 4 5 6 7 8
        # cha: A A A B B B B B

        contact_map2 = ContactMap('bar')
        contact1 = Contact(6, 2, 1.0)
        contact1.res1_chain = 'A'
        contact1.res2_chain = 'B'
        contact1.res1_altseq = 1
        contact1.res2_altseq = 5
        contact_map2.add(contact1)

        contact2 = Contact(7, 4, 1.0)
        contact2.res1_chain = 'A'
        contact2.res2_chain = 'B'
        contact2.res1_altseq = 2
        contact2.res2_altseq = 7
        contact2.define_match()
        contact_map2.add(contact2)

        contact3 = Contact(8, 2, 1.0)
        contact3.res1_chain = 'A'
        contact3.res2_chain = 'B'
        contact3.res1_altseq = 3
        contact3.res2_altseq = 5
        contact3.define_match()
        contact_map2.add(contact3)

        contact_map2.sequence = Sequence('bar', 'ABCDEFGH')
        contact_map2.assign_sequence_register(altloc=True)

        contact_map1.match(contact_map2, renumber=True, inplace=True)
        self.assertEqual([
            ContactMatchState.matched.value, ContactMatchState.mismatched.value, ContactMatchState.matched.value,
            ContactMatchState.matched.value, ContactMatchState.mismatched.value
        ], [c.status for c in contact_map1])
        self.assertEqual([(6, 2), (6, _Gap.IDENTIFIER), (7, 4), (8, 2), (7, _Gap.IDENTIFIER)],
                         [(c.res1_seq, c.res2_seq) for c in contact_map1])
        self.assertEqual([('A', 'B'), ('A', ''), ('A', 'B'), ('A', 'B'), ('A', '')],
                         [(c.res1_chain, c.res2_chain) for c in contact_map1])

    def test_match_11(self):
        reference = ContactMap("ref")
        for cparam in [(1, 5, 1.0), (1, 7, 1.0), (2, 7, 1.0), (3, 4, 1.0)]:
            contact = Contact(*cparam)
            contact.res1_altseq = cparam[0]
            contact.res2_altseq = cparam[1]
            contact.define_match()
            reference.add(contact)
        reference.sequence = Sequence('bar', 'ABDEFGH')
        reference.assign_sequence_register(altloc=True)

        sequence = Sequence("foo", "ABCDEFGH")
        for i in range(2):
            target = ContactMap("tar")
            for cparam in [(1, 5, 1.0), (1, 6, 1.0), (2, 7, 1.0), (3, 5, 1.0), (2, 8, 1.0)]:
                target.add(Contact(*cparam))
            target.remove_neighbors(min_distance=1, inplace=True)
            target.sort("raw_score", reverse=True, inplace=True)

            target.sequence = sequence
            target.assign_sequence_register()
            target.id = "foobar_%d" % i

            target.match(reference, match_other=False, inplace=True)

    def test_remove_neighbors_1(self):
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_mod = contact_map.remove_neighbors(min_distance=2)
        self.assertListEqual([(1, 5), (2, 4), (5, 1)], [c.id for c in contact_map_mod])
        self.assertEqual([(1, 5), (2, 4), (5, 1)], sorted(contact_map_mod.child_dict.keys()))

    def test_remove_neighbors_2(self):
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_mod = contact_map.remove_neighbors(min_distance=5, inplace=True)
        self.assertListEqual([], [c.get_id() for c in contact_map_mod])
        self.assertDictEqual({}, contact_map_mod.child_dict)
        self.assertEqual(contact_map, contact_map_mod)

    def test_remove_neighbors_3(self):
        contact_map = ContactMap('test')
        for c in [Contact(1, 30, 1.0), Contact(2, 10, 0.4), Contact(3, 20, 0.1), Contact(1, 5, 0.2)]:
            contact_map.add(c)
        contact_map_mod = contact_map.remove_neighbors(min_distance=6, max_distance=12, inplace=True)
        self.assertEqual([(2, 10)], [c.id for c in contact_map_mod])
        self.assertEqual([(2, 10)], sorted(contact_map_mod.child_dict.keys()))

    def test_remove_neighbors_4(self):
        contact_map = ContactMap('test')
        for c in [Contact(1, 30, 1.0), Contact(2, 10, 0.4), Contact(3, 20, 0.1), Contact(1, 5, 0.2)]:
            contact_map.add(c)
        contact_map_mod = contact_map.remove_neighbors(min_distance=12, max_distance=24, inplace=True)
        self.assertEqual([(3, 20)], [c.id for c in contact_map_mod])
        self.assertEqual([(3, 20)], sorted(contact_map_mod.child_dict.keys()))

    def test_remove_neighbors_5(self):
        contact_map = ContactMap('test')
        for c in [Contact(1, 30, 1.0), Contact(2, 10, 0.4), Contact(3, 20, 0.1), Contact(1, 5, 0.2)]:
            contact_map.add(c)
        contact_map_mod = contact_map.remove_neighbors(min_distance=24, inplace=True)
        self.assertEqual([(1, 30)], [c.id for c in contact_map_mod])
        self.assertEqual([(1, 30)], sorted(contact_map_mod.child_dict.keys()))

    def test_remove_neighbors_6(self):
        contact_map = ContactMap('test')
        for c in [Contact(1, 30, 1.0), Contact(2, 10, 0.4), Contact(3, 20, 0.1), Contact(1, 5, 0.2)]:
            contact_map.add(c)
        contact_map_mod = contact_map.remove_neighbors(min_distance=6, max_distance=24, inplace=True)
        self.assertEqual([(2, 10), (3, 20)], [c.id for c in contact_map_mod])
        self.assertEqual([(2, 10), (3, 20)], sorted(contact_map_mod.child_dict.keys()))

    def test_rescale_1(self):
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_rescaled = contact_map.rescale()
        self.assertListEqual([1.0, 0.333, 0.0, 0.111], [round(c.raw_score, 3) for c in contact_map_rescaled])
        self.assertNotEqual(contact_map, contact_map_rescaled)

    def test_rescale_2(self):
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_rescaled = contact_map.rescale(inplace=True)
        self.assertListEqual([1.0, 0.333, 0.0, 0.111], [round(c.raw_score, 3) for c in contact_map])
        self.assertEqual(contact_map, contact_map_rescaled)

    def test_sort_1(self):
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_sorted = contact_map.sort('res1_seq', reverse=False, inplace=False)
        self.assertEqual([(1, 5), (2, 4), (3, 3), (5, 1)], [c.id for c in contact_map_sorted])
        self.assertNotEqual(contact_map, contact_map_sorted)

    def test_sort_2(self):
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_sorted = contact_map.sort('res1_seq', reverse=True, inplace=False)
        self.assertEqual([(5, 1), (3, 3), (2, 4), (1, 5)], [c.id for c in contact_map_sorted])
        self.assertNotEqual(contact_map, contact_map_sorted)

    def test_sort_3(self):
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_sorted = contact_map.sort('raw_score', reverse=False, inplace=True)
        self.assertEqual([(2, 4), (5, 1), (3, 3), (1, 5)], [c.id for c in contact_map_sorted])
        self.assertEqual(contact_map, contact_map_sorted)

    def test_sort_4(self):
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_sorted = contact_map.sort('res2_seq', reverse=True, inplace=True)
        self.assertEqual([(1, 5), (2, 4), (3, 3), (5, 1)], [c.id for c in contact_map_sorted])
        self.assertEqual(contact_map, contact_map_sorted)

    def test__adjust_1(self):
        contact_map = ContactMap("test")
        for contact in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(1, 1, 0.2)]:
            contact_map.add(contact)
            contact.res1_altseq = contact.res1_seq
            contact.res2_altseq = contact.res2_seq
        contact_map_keymap = ContactMap._create_keymap(contact_map)
        sequence = [ord(x) for x in 'XXXXX']
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
        sequence = [ord(x) for x in 'XX-XX']
        contact_map_keymap = ContactMap._insert_states(sequence, contact_map_keymap)
        contact_map_keymap = ContactMap._reindex_by_keymap(contact_map_keymap)
        adjusted = ContactMap._adjust(contact_map, contact_map_keymap)
        self.assertEqual([1, 2, 1], [c.res1_altseq for c in adjusted])
        self.assertEqual([5, 4, 1], [c.res2_altseq for c in adjusted])

    def test__create_keymap_1(self):
        contact_map = ContactMap("test")
        for contact in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(1, 1, 0.2)]:
            contact_map.add(contact)
        contact_map.sequence = Sequence('foo', 'ABCDE')
        contact_map.assign_sequence_register()
        contact_map_keymap = ContactMap._create_keymap(contact_map)
        self.assertEqual([1, 2, 3, 4, 5], [c.res_seq for c in contact_map_keymap])
        self.assertEqual(['A', 'B', 'C', 'D', 'E'], [c.res_name for c in contact_map_keymap])
        self.assertEqual([0, 0, 0, 0, 0], [c.res_altseq for c in contact_map_keymap])

    def test__create_keymap_2(self):
        contact_map = ContactMap("test")
        for contact, res_altloc in [(Contact(1, 5, 1.0), (10, 20)), (Contact(3, 3, 0.4), (12, 12)),
                                    (Contact(2, 4, 0.1), (11, 13)), (Contact(1, 1, 0.2), (10, 10))]:
            contact.res1_altseq = res_altloc[0]
            contact.res2_altseq = res_altloc[1]
            contact_map.add(contact)
        contact_map.sequence = Sequence('foo', 'ABCDE')
        contact_map.assign_sequence_register()
        contact_map_keymap = ContactMap._create_keymap(contact_map, altloc=True)
        self.assertEqual([1, 2, 3, 4, 5], [c.res_seq for c in contact_map_keymap])
        self.assertEqual(['A', 'B', 'C', 'D', 'E'], [c.res_name for c in contact_map_keymap])
        self.assertEqual([10, 11, 12, 13, 20], [c.res_altseq for c in contact_map_keymap])

    def test__create_keymap_3(self):
        contact_map = ContactMap("test")
        for contact, res_altloc in [(Contact(1, 5, 1.0), (10, 20)), (Contact(2, 4, 0.1), (11, 13)), (Contact(1, 1, 0.2),
                                                                                                     (10, 10))]:
            contact.res1_altseq = res_altloc[0]
            contact.res2_altseq = res_altloc[1]
            contact_map.add(contact)
        contact_map.sequence = Sequence('foo', 'ABCDE')
        contact_map.assign_sequence_register()
        contact_map_keymap = ContactMap._create_keymap(contact_map, altloc=True)
        self.assertEqual([1, 2, 4, 5], [c.res_seq for c in contact_map_keymap])
        self.assertEqual(['A', 'B', 'D', 'E'], [c.res_name for c in contact_map_keymap])
        self.assertEqual([10, 11, 13, 20], [c.res_altseq for c in contact_map_keymap])

    def test__find_single_1(self):
        contact_map = ContactMap("test")
        for contact, res_altloc in [(Contact(1, 5, 1.0), (10, 20)), (Contact(2, 4, 0.1), (11, 13)), (Contact(1, 1, 0.2),
                                                                                                     (10, 10))]:
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
        for contact, res_altloc in [(Contact(1, 5, 1.0), (10, 20)), (Contact(2, 4, 0.1), (11, 13)), (Contact(1, 1, 0.2),
                                                                                                     (10, 10))]:
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
        keymap = [_Residue(1, 1, 'X', ''), _Residue(2, 2, 'X', ''), _Residue(3, 3, 'X', '')]
        sequence = [ord(x) for x in 'XXX']
        inserts_added = ContactMap._insert_states(sequence, keymap)
        self.assertEqual(3, len(inserts_added))
        self.assertEqual([True, True, True], [isinstance(r, _Residue) for r in inserts_added])

    def test__insert_states_2(self):
        keymap = [_Residue(2, 2, 'X', ''), _Residue(3, 3, 'X', '')]
        sequence = [ord(x) for x in '-XX']
        inserts_added = ContactMap._insert_states(sequence, keymap)
        self.assertEqual(3, len(inserts_added))
        self.assertEqual([False, True, True], [isinstance(r, _Residue) for r in inserts_added])
        self.assertEqual([True, False, False], [isinstance(r, _Gap) for r in inserts_added])

    def test__insert_states_3(self):
        keymap = [_Residue(2, 10, 'X', ''), _Residue(3, 11, 'X', '')]
        sequence = [ord(x) for x in '-X-X--']
        inserts_added = ContactMap._insert_states(sequence, keymap)
        self.assertEqual(6, len(inserts_added))
        self.assertEqual([False, True, False, True, False, False], [isinstance(r, _Residue) for r in inserts_added])
        self.assertEqual([True, False, True, False, True, True], [isinstance(r, _Gap) for r in inserts_added])
        self.assertEqual([2, 3], [r.res_seq for r in inserts_added if isinstance(r, _Residue)])
        self.assertEqual([10, 11], [r.res_altseq for r in inserts_added if isinstance(r, _Residue)])

    def test__reindex_by_keymap_1(self):
        keymap = [_Residue(1, 1, 'X', ''), _Residue(2, 2, 'X', ''), _Residue(3, 3, 'X', '')]
        reindex = ContactMap._reindex_by_keymap(keymap)
        self.assertEqual([1, 2, 3], [c.res_seq for c in reindex])
        self.assertEqual([1, 2, 3], [c.res_altseq for c in reindex])

    def test__reindex_by_keymap_2(self):
        keymap = [_Residue(1, -5, 'X', ''), _Residue(2, -4, 'X', ''), _Residue(3, -3, 'X', '')]
        reindex = ContactMap._reindex_by_keymap(keymap)
        self.assertEqual([1, 2, 3], [c.res_seq for c in reindex])
        self.assertEqual([1, 2, 3], [c.res_altseq for c in reindex])

    def test__reindex_by_keymap_3(self):
        keymap = [_Gap(), _Residue(2, 1, 'X', ''), _Residue(3, 2, 'X', '')]
        reindex = ContactMap._reindex_by_keymap(keymap)
        self.assertEqual([_Gap.IDENTIFIER, 2, 3], [c.res_seq for c in reindex])
        self.assertEqual([1, 2, 3], [c.res_altseq for c in reindex])

    def test__reindex_by_keymap_4(self):
        keymap = [_Gap(), _Residue(200000, 10000, 'X', ''), _Gap(), _Gap()]
        reindex = ContactMap._reindex_by_keymap(keymap)
        self.assertEqual([_Gap.IDENTIFIER, 200000, _Gap.IDENTIFIER, _Gap.IDENTIFIER], [c.res_seq for c in reindex])
        self.assertEqual([1, 2, 3, 4], [c.res_altseq for c in reindex])

    def test__reindex_by_keymap_5(self):
        keymap = [_Gap(), _Gap(), _Gap(), _Gap()]
        reindex = ContactMap._reindex_by_keymap(keymap)
        self.assertEqual([_Gap.IDENTIFIER, _Gap.IDENTIFIER, _Gap.IDENTIFIER, _Gap.IDENTIFIER],
                         [c.res_seq for c in reindex])
        self.assertEqual([1, 2, 3, 4], [c.res_altseq for c in reindex])

    def test_as_list_1(self):
        contact_map = ContactMap('test')
        for c in [Contact(1, 30, 1.0), Contact(2, 10, 0.4), Contact(3, 20, 0.1), Contact(1, 5, 0.2)]:
            contact_map.add(c)
        self.assertListEqual([[1, 30], [2, 10], [3, 20], [1, 5]], contact_map.as_list())

    def test_as_list_2(self):
        contact_map = ContactMap('test')
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
        for c in [Contact(4, 4, 1.0), Contact(5, 5, 0.4), Contact(4, 6, 0.1), Contact(3, 5, 0.2)]:
            contact_map.add(c)
        self.assertListEqual([], contact_map.singletons.as_list())

    def test_singletons_2(self):
        contact_map = ContactMap("test")
        for c in [Contact(4, 4, 1.0), Contact(5, 5, 0.4), Contact(4, 6, 0.1), Contact(10, 10, 0.2)]:
            contact_map.add(c)
        self.assertListEqual([[10, 10]], contact_map.singletons.as_list())

    def test_singletons_3(self):
        contact_map = ContactMap("test")
        for c in [
                Contact(4, 4, 1.0),
                Contact(5, 5, 0.4),
                Contact(4, 6, 0.1),
                Contact(10, 10, 0.2),
                Contact(10, 11, 0.2)
        ]:
            contact_map.add(c)
        self.assertListEqual([], contact_map.singletons.as_list())

    def test_singletons_4(self):
        contact_map = ContactMap("test")
        for c in [Contact(4, 4, 1.0), Contact(4, 6, 0.4)]:
            contact_map.add(c)
        self.assertListEqual([], contact_map.singletons.as_list())

    def test_singletons_5(self):
        contact_map = ContactMap("test")
        for c in [Contact(4, 4, 1.0), Contact(4, 8, 0.4)]:
            contact_map.add(c)
        self.assertListEqual([[4, 4], [4, 8]], contact_map.singletons.as_list())

    def test_singletons_6(self):
        contact_map = ContactMap("test")
        for c in [Contact(4, 4, 1.0), Contact(6, 6, 0.4)]:
            contact_map.add(c)
        self.assertListEqual([], contact_map.singletons.as_list())

    def test_singletons_7(self):
        contact_map = ContactMap("test")
        for c in [Contact(4, 4, 1.0), Contact(7, 7, 0.4)]:
            contact_map.add(c)
        self.assertListEqual([[4, 4], [7, 7]], contact_map.singletons.as_list())


if __name__ == "__main__":
    unittest.main(verbosity=2)
