"""Testing facility for conkit.core.contactmap"""

__author__ = "Felix Simkovic"
__date__ = "12 Aug 2016"

from conkit.core.Contact import Contact
from conkit.core.ContactMap import ContactMap
from conkit.core.Sequence import Sequence

import unittest


class Test(unittest.TestCase):

    def test_coverage(self):
        # ======================================================
        # Test Case 1
        contact_map = ContactMap('test')
        contact_map.add(Contact(1, 4, 1.0))
        contact_map.add(Contact(2, 4, 1.0))
        contact_map.add(Contact(5, 8, 1.0))
        contact_map.add(Contact(3, 6, 1.0))
        contact_map.sequence = Sequence('TEST', 'ABCDEFGH')
        self.assertEqual(0.875, contact_map.coverage)
        # ======================================================
        # Test Case 2
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

    def test_ncontacts(self):
        # ======================================================
        # Test Case 1
        contact_map = ContactMap('test')
        self.assertEqual(0, contact_map.ncontacts)
        # ======================================================
        # Test Case 2
        contact_map = ContactMap('test')
        contact_map.add(Contact(1, 5, 1.0))
        self.assertEqual(1, contact_map.ncontacts)
        # ======================================================
        # Test Case 3
        contact_map = ContactMap('test')
        contact_map.add(Contact(1, 5, 1.0))
        contact_map.add(Contact(2, 10, 1.0))
        self.assertEqual(2, contact_map.ncontacts)

    def test_precision(self):
        # ======================================================
        # Test Case 1
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map.sequence = Sequence('TEST', 'AAAAA')
        for i, contact in enumerate(contact_map):
            if i % 2 == 0:
                contact.define_true_positive()
            else:
                contact.define_false_positive()
        self.assertEqual(0.5, contact_map.precision)

    def test_repr_sequence(self):
        # ======================================================
        # Test Case 1
        contact_map = ContactMap('test')
        for contact in [Contact(1, 5, 1.0), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        contact_map.sequence = Sequence('foo', 'ABCDE')
        self.assertEqual('AB-DE', contact_map.repr_sequence.seq)
        contact_map.remove((2, 4))
        # ======================================================
        # Test Case 2
        contact_map = ContactMap('test')
        for contact in [Contact(1, 5, 1.0), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        contact_map.sequence = Sequence('foo', 'ABCDE')
        self.assertEqual('AB-DE', contact_map.repr_sequence.seq)
        contact_map.remove((2, 4))
        self.assertEqual('A---E', contact_map.repr_sequence.seq)

    def test_repr_sequence_altloc(self):
        # ======================================================
        # Test Case 1
        contact_map = ContactMap('test')
        for contact, altseq in [(Contact(1, 5, 1.0), (3, 4)),
                                (Contact(2, 4, 0.1), (1, 5)),
                                (Contact(5, 1, 0.2), (4, 3))]:
            contact.res1_altseq = altseq[0]
            contact.res2_altseq = altseq[1]
            contact_map.add(contact)
        contact_map.sequence = Sequence('foo', 'ABCDE')
        self.assertEqual('A-CDE', contact_map.repr_sequence_altloc.seq)
        self.assertEqual('AB-DE', contact_map.repr_sequence.seq)
        # ======================================================
        # Test Case 2
        contact_map = ContactMap('test')
        for contact, altseq in [(Contact(1, 5, 1.0), (3, 4)),
                                (Contact(2, 4, 0.1), (1, 5)),
                                (Contact(5, 1, 0.2), (4, 3))]:
            contact.res1_altseq = altseq[0]
            contact.res2_altseq = altseq[1]
            contact_map.add(contact)
        contact_map.sequence = Sequence('foo', 'ABCDE')
        self.assertEqual('A-CDE', contact_map.repr_sequence_altloc.seq)
        self.assertEqual('AB-DE', contact_map.repr_sequence.seq)
        contact_map.remove((2, 4))
        self.assertEqual('--CD-', contact_map.repr_sequence_altloc.seq)
        self.assertEqual('A---E', contact_map.repr_sequence.seq)

    def test_sequence(self):
        # ======================================================
        # Test Case 1
        contact_map = ContactMap('test')
        for contact in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        seq_obj = Sequence('bar', 'ABCDE')
        contact_map.sequence = seq_obj
        self.assertEqual(seq_obj, contact_map.sequence)
        # ======================================================
        # Test Case 2
        contact_map = ContactMap('test')
        for contact in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        sequence1 = Sequence('foo', 'ABC')
        sequence2 = Sequence('bar', 'DE')
        self.assertRaises(TypeError, contact_map.sequence, [sequence1, sequence2])

    def test_top_contact(self):
        # ======================================================
        # Test Case 1
        contact_map = ContactMap('test')
        self.assertEqual(None, contact_map.top_contact)
        # ======================================================
        # Test Case 2
        contact_map = ContactMap('test')
        contact = Contact(1, 10, 1.0)
        contact_map.add(contact)
        self.assertEqual(contact, contact_map.top_contact)
        # ======================================================
        # Test Case 3
        contact_map = ContactMap('test')
        contact1 = Contact(1, 10, 1.0)
        contact2 = Contact(2, 100, 1.0)
        contact_map.add(contact1)
        contact_map.add(contact2)
        self.assertEqual(contact1, contact_map.top_contact)

    def test__construct_repr_sequence(self):
        # ======================================================
        # Test Case 1
        contact_map = ContactMap('test')
        for contact in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        seq_obj = Sequence('bar', 'ABCDE')
        contact_map.sequence = seq_obj
        self.assertEqual('ABCDE', contact_map._construct_repr_sequence([1, 2, 3, 4, 5]).seq)
        # ======================================================
        # Test Case 2
        contact_map = ContactMap('test')
        for contact in [Contact(1, 5, 1.0), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        seq_obj = Sequence('bar', 'ABCDE')
        contact_map.sequence = seq_obj
        self.assertEqual('AB-DE', contact_map._construct_repr_sequence([1, 2, 4, 5]).seq)
        # ======================================================
        # Test Case 3
        contact_map = ContactMap('test')
        for contact in [Contact(1, 5, 1.0), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        seq_obj = Sequence('bar', 'ABCDE')
        contact_map.sequence = seq_obj
        self.assertEqual('A---E', contact_map._construct_repr_sequence([1, 5]).seq)
        # ======================================================
        # Test Case 4
        contact_map = ContactMap('test')
        for contact in [Contact(1, 5, 1.0), Contact(5, 1, 0.2)]:
            contact_map.add(contact)
        seq_obj = Sequence('bar', 'ABCDE')
        contact_map.sequence = seq_obj
        self.assertEqual('----E', contact_map._construct_repr_sequence([-1, 5]).seq)

    def test__create_keymap(self):
        # ======================================================
        # Test Case 1
        contact_map = ContactMap("test")
        for contact in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(1, 1, 0.2)]:
            contact_map.add(contact)
        contact_map.sequence = Sequence('foo', 'ABCDE')
        contact_map.assign_sequence_register()
        keymap = [
            {'res_seq': 1, 'res_name': 'A', 'res_matchseq': 999, 'res_altseq': 0, 'status': -2},
            {'res_seq': 2, 'res_name': 'B', 'res_matchseq': 999, 'res_altseq': 0, 'status': -2},
            {'res_seq': 3, 'res_name': 'C', 'res_matchseq': 999, 'res_altseq': 0, 'status': -2},
            {'res_seq': 4, 'res_name': 'D', 'res_matchseq': 999, 'res_altseq': 0, 'status': -2},
            {'res_seq': 5, 'res_name': 'E', 'res_matchseq': 999, 'res_altseq': 0, 'status': -2},
        ]
        contact_map_keymap = contact_map.create_keymap()
        self.assertEqual(keymap, contact_map_keymap)
        # ======================================================
        # Test Case 2
        contact_map = ContactMap("test")
        for contact, res_altloc in [(Contact(1, 5, 1.0), (10, 20)),
                                    (Contact(3, 3, 0.4), (12, 12)),
                                    (Contact(2, 4, 0.1), (11, 13)),
                                    (Contact(1, 1, 0.2), (10, 10))]:
            contact.res1_altseq = res_altloc[0]
            contact.res2_altseq = res_altloc[1]
            contact_map.add(contact)
        contact_map.sequence = Sequence('foo', 'ABCDE')
        contact_map.assign_sequence_register()
        keymap = [
            {'res_seq': 1, 'res_name': 'A', 'res_matchseq': 999, 'res_altseq': 10, 'status': -2},
            {'res_seq': 2, 'res_name': 'B', 'res_matchseq': 999, 'res_altseq': 11, 'status': -2},
            {'res_seq': 3, 'res_name': 'C', 'res_matchseq': 999, 'res_altseq': 12, 'status': -2},
            {'res_seq': 4, 'res_name': 'D', 'res_matchseq': 999, 'res_altseq': 13, 'status': -2},
            {'res_seq': 5, 'res_name': 'E', 'res_matchseq': 999, 'res_altseq': 20, 'status': -2},
        ]
        contact_map_keymap = contact_map.create_keymap(altloc=True)
        self.assertEqual(keymap, contact_map_keymap)

    def test_calculate_scalar_score(self):
        # ======================================================
        # Test Case 1
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map.calculate_scalar_score()
        self.assertListEqual([2.352941, 0.941176, 0.235294, 0.470588], [round(c.scalar_score, 6) for c in contact_map])

    def test_assign_sequence_register(self):
        # ======================================================
        # Test Case 1
        contact_map1 = ContactMap('1')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        contact_map1.sequence = Sequence('foo', 'ABCDEF')
        contact_map1.assign_sequence_register()
        self.assertEqual(
            [('A', 'E'), ('B', 'F'), ('A', 'D'), ('C', 'F'), ('B', 'E')],
            [(c.res1, c.res2) for c in contact_map1]
        )
        # ======================================================
        # Test Case 2
        contact_map1 = ContactMap('1')
        for comb, alt in [((101, 105, 1.0), (1, 5)),
                          ((102, 106, 1.0), (2, 6)),
                          ((101, 104, 1.0), (1, 4)),
                          ((103, 106, 1.0), (3, 6)),
                          ((102, 105, 1.0), (2, 5))]:
            c = Contact(*comb)
            c.res1_altseq = alt[0]
            c.res2_altseq = alt[1]
            contact_map1.add(c)
        contact_map1.sequence = Sequence('foo', 'ABCDEF')
        contact_map1.assign_sequence_register(altloc=True)
        self.assertEqual(
            [('A', 'E'), ('B', 'F'), ('A', 'D'), ('C', 'F'), ('B', 'E')],
            [(c.res1, c.res2) for c in contact_map1]
        )

    def test_match(self):
        # ======================================================
        # Test Case 1
        contact_map1 = ContactMap('1')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        contact_map1.sequence = Sequence('foo', 'ABCDEF')
        contact_map1.assign_sequence_register()
        contact_map2 = ContactMap('2')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            c = Contact(*comb)
            c.res1_altseq = comb[0]
            c.res2_altseq = comb[1]
            contact_map2.add(c)
        contact_map2.sequence = Sequence('bar', 'ABCDEF')
        contact_map2.assign_sequence_register()
        contact_map1.match(contact_map2, remove_unmatched=False, renumber=False, inplace=True)
        self.assertEqual([1, 2, 1, 3, 2], [c.res1_seq for c in contact_map1])
        self.assertEqual([5, 6, 4, 6, 5], [c.res2_seq for c in contact_map1])

        # ======================================================
        # Test Case 2
        contact_map1 = ContactMap('1')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        contact_map1.sequence = Sequence('foo', 'ABCDEF')
        contact_map1.assign_sequence_register()
        contact_map2 = ContactMap('2')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            c = Contact(*comb)
            c.res1_altseq = comb[0]
            c.res2_altseq = comb[1]
            contact_map2.add(c)
        contact_map2.sequence = Sequence('bar', 'ABCDEF')
        contact_map2.assign_sequence_register()
        contact_map1.match(contact_map2, remove_unmatched=True, renumber=False, inplace=True)
        self.assertEqual([1, 2, 3, 2], [c.res1_seq for c in contact_map1])
        self.assertEqual([5, 6, 6, 5], [c.res2_seq for c in contact_map1])

        # ======================================================
        # Test Case 3
        contact_map1 = ContactMap('1')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        contact_map1.sequence = Sequence('foo', 'ABCDEF')
        contact_map1.assign_sequence_register()
        contact_map2 = ContactMap('2')
        for comb, alt in [((101, 105, 1.0), (1, 5)),
                          ((102, 106, 1.0), (2, 6)),
                          ((101, 104, 1.0), (1, 4)),
                          ((103, 106, 1.0), (3, 6)),
                          ((102, 105, 1.0), (2, 5))]:
            c = Contact(*comb)
            c.res1_altseq = alt[0]
            c.res2_altseq = alt[1]
            contact_map2.add(c)
        contact_map2.sequence = Sequence('bar', 'ABCDEF')
        contact_map2.assign_sequence_register(altloc=True)
        contact_map1.match(contact_map2, remove_unmatched=False, renumber=False, inplace=True)
        self.assertEqual([1, 2, 1, 3, 2], [c.res1_seq for c in contact_map1])
        self.assertEqual([5, 6, 4, 6, 5], [c.res2_seq for c in contact_map1])
        self.assertEqual([True, True, True, True, True], [c.is_true_positive for c in contact_map1])

        # ======================================================
        # Test Case 4 - Do renumber
        contact_map1 = ContactMap('1')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        contact_map1.sequence = Sequence('foo', 'ABCDEF')
        contact_map1.assign_sequence_register()
        contact_map2 = ContactMap('2')
        for comb, alt in [((101, 105, 1.0), (1, 5)),
                          ((102, 106, 1.0), (2, 6)),
                          ((101, 104, 1.0), (1, 4)),
                          ((103, 106, 1.0), (3, 6)),
                          ((102, 105, 1.0), (2, 5))]:
            c = Contact(*comb)
            c.res1_altseq = alt[0]
            c.res2_altseq = alt[1]
            contact_map2.add(c)
        contact_map2.sequence = Sequence('bar', 'ABCDEF')
        contact_map2.assign_sequence_register(altloc=True)
        contact_map1.match(contact_map2, remove_unmatched=False, renumber=True, inplace=True)
        self.assertEqual([101, 102, 101, 103, 102], [c.res1_seq for c in contact_map1])
        self.assertEqual([105, 106, 104, 106, 105], [c.res2_seq for c in contact_map1])
        self.assertEqual([True, True, True, True, True], [c.is_true_positive for c in contact_map1])

        # ======================================================
        # Test Case 5
        contact_map1 = ContactMap('1')
        for comb in [(34, 129, 1.0), (34, 240, 1.0), (36, 129, 1.0), (36, 240, 1.0), (38, 129, 1.0),
                     (38, 240, 1.0), (40, 129, 1.0), (40, 240, 1.0), (42, 129, 1.0), (42, 240, 1.0),
                     (44, 129, 1.0), (44, 240, 1.0), (46, 129, 1.0), (46, 240, 1.0), (48, 129, 1.0),
                     (48, 240, 1.0), (50, 129, 1.0), (50, 240, 1.0), (52, 129, 1.0), (52, 240, 1.0),
                     (54, 129, 1.0), (54, 240, 1.0), (56, 129, 1.0), (56, 240, 1.0), (58, 179, 1.0),
                     (59, 170, 1.0), (60, 179, 1.0), (61, 170, 1.0), (62, 179, 1.0), (63, 170, 1.0),
                     (65, 161, 1.0), (115, 162, 1.0), (116, 161, 1.0), (117, 168, 1.0), (119, 161, 1.0),
                     (120, 162, 1.0)]:
            contact_map1.add(Contact(*comb))
        contact_map1.sequence = Sequence('foo', 'DDLTISSLAKGETTKAAFNQMVQGHKLPAWVMKGGTYTPAQTVTLGDETYQVMSACKPHDCGSQRIAVMW'
                                                'SEKSNQMTGLFSTIDEKTSQEKLTWLNVNDALSIDGKTVLFAALTGSLENHPDGFNFKVFGRCELAAAMK'
                                                'RHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSAL'
                                                'LSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCR')
        contact_map1.assign_sequence_register()
        contact_map2 = ContactMap('pdb')
        for comb, alt in [((3, 43, 1.0), (2, 170)), ((35, 45, 1.0), (34, 172)),
                          ((36, 44, 1.0), (35, 171)), ((37, 35, 1.0), (36, 162)),
                          ((37, 36, 1.0), (36, 163)), ((37, 42, 1.0), (36, 169)),
                          ((37, 43, 1.0), (36, 170)), ((37, 44, 1.0), (36, 171)),
                          ((37, 52, 1.0), (36, 179)), ((37, 57, 1.0), (36, 184)),
                          ((38, 42, 1.0), (37, 169)), ((38, 43, 1.0), (37, 170)),
                          ((38, 44, 1.0), (37, 171)), ((39, 39, 1.0), (38, 166)),
                          ((39, 41, 1.0), (38, 168)), ((39, 42, 1.0), (38, 169)),
                          ((39, 43, 1.0), (38, 170)), ((40, 41, 1.0), (39, 168)),
                          ((57, 36, 1.0), (56, 163)), ((57, 37, 1.0), (56, 164)),
                          ((57, 42, 1.0), (56, 169)), ((59, 44, 1.0), (58, 171)),
                          ((60, 34, 1.0), (59, 161)), ((60, 35, 1.0), (59, 162)),
                          ((60, 57, 1.0), (59, 184)), ((60, 108, 1.0), (59, 235)),
                          ((60, 109, 1.0), (59, 236)), ((60, 110, 1.0), (59, 237)),
                          ((60, 113, 1.0), (59, 240)), ((61, 34, 1.0), (60, 161)),
                          ((61, 110, 1.0), (60, 237)), ((61, 113, 1.0), (60, 240)),
                          ((61, 114, 1.0), (60, 241)), ((62, 33, 1.0), (61, 160)),
                          ((62, 34, 1.0), (61, 161)), ((62, 35, 1.0), (61, 162)),
                          ((62, 36, 1.0), (61, 163)), ((62, 37, 1.0), (61, 164)),
                          ((62, 39, 1.0), (61, 166)), ((62, 42, 1.0), (61, 169)),
                          ((62, 110, 1.0), (61, 237)), ((63, 33, 1.0), (62, 160)),
                          ((63, 34, 1.0), (62, 161)), ((63, 37, 1.0), (62, 164)),
                          ((64, 34, 1.0), (63, 161)), ((64, 114, 1.0), (63, 241)),
                          ((66, 37, 1.0), (65, 164)), ((116, 2, 1.0), (115, 129)),
                          ((117, 2, 1.0), (116, 129)), ((117, 39, 1.0), (116, 166)),
                          ((120, 1, 1.0), (119, 128)), ((120, 2, 1.0), (119, 129)),
                          ((120, 39, 1.0), (119, 166)), ((120, 41, 1.0), (119, 168)),
                          ((121, 1, 1.0), (120, 128)), ((121, 2, 1.0), (120, 129))]:
            c = Contact(*comb)
            c.res1_altseq = alt[0]
            c.res2_altseq = alt[1]
            contact_map2.add(c)
        contact_map2.sequence = Sequence('bar', 'DDLTISSLAKGETTKAAFNQMVQGHKLPAWVMKGGTYTPAQTVTLGDETYQVMSACKPHDCGSQRIAVMW'
                                                'SEKSNQMTGLFSTIDEKTSQEKLTWLNVNDALSIDGKTVLFAALTGSLENHPDGFNFKVFGRCELAAAMK'
                                                'RHGLDNYRGYSLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSAL'
                                                'LSSDITASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCR')
        contact_map2.assign_sequence_register(altloc=True)
        # Test Case 5a
        cmap1_copy1 = contact_map1.copy()
        cmap2_copy1 = contact_map2.copy()
        cmap1_copy1.match(cmap2_copy1, remove_unmatched=False, renumber=True, inplace=True)
        self.assertEqual(
            [35, 35, 37, 37, 39, 39, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 999, 57, 57, 59, 60, 61, 62, 63, 64, 66, 116, 117, 999, 120, 121],
            [c.res1_seq for c in cmap1_copy1]
        )
        self.assertEqual(
            [2, 113, 2, 113, 2, 113, 2, 113, 2, 113, 2, 113, 2, 113, 2, 113, 2, 113, 2, 113, 2, 113, 2, 113, 52, 43, 52,
             43, 52, 43, 34, 35, 34, 41, 34, 35],
            [c.res2_seq for c in cmap1_copy1]
        )

        # Test Case 5b
        cmap1_copy2 = contact_map1.copy()
        cmap2_copy2 = contact_map2.copy()
        cmap1_copy2.match(cmap2_copy2, remove_unmatched=True, renumber=True, inplace=True)
        self.assertItemsEqual(
            [35, 35, 37, 37, 39, 39, 57, 57, 59, 60, 61, 62, 63, 64, 66, 116, 117, 120, 121],
            [c.res1_seq for c in cmap1_copy2]
        )
        self.assertItemsEqual(
            [2, 113, 2, 113, 2, 113, 2, 113, 52, 43, 52, 43, 52, 43, 34, 35, 34, 34, 35],
            [c.res2_seq for c in cmap1_copy2]
        )
        cmap1_copy2_repr_sequence = '---------------------------------G-T-T-----------------C-PHDCGS-R----------------' \
                                    '---------------------------------TG--EN--------V-------------------------------FE' \
                                    '-------T--------D------------------------------------------------------------N---' \
                                    '------------'
        self.assertEqual(cmap1_copy2_repr_sequence, cmap1_copy2.repr_sequence_altloc.seq)

    def test_remove_neighbors(self):
        # ======================================================
        # Test Case 1
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_mod = contact_map.remove_neighbors(min_distance=2)
        self.assertListEqual([(1, 5), (2, 4), (5, 1)], [c.id for c in contact_map_mod])
        self.assertItemsEqual([(1, 5), (2, 4), (5, 1)], contact_map_mod.child_dict.keys())
        # ======================================================
        # Test Case 2
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_mod = contact_map.remove_neighbors(min_distance=5, inplace=True)
        self.assertListEqual([], [c.get_id() for c in contact_map_mod])
        self.assertDictEqual({}, contact_map_mod.child_dict)
        self.assertEqual(contact_map, contact_map_mod)

    def test_rescale(self):
        # ======================================================
        # Test Case 1
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_rescaled = contact_map.rescale()
        self.assertListEqual([1.0, 0.333, 0.0, 0.111], [round(c.raw_score, 3) for c in contact_map_rescaled])
        self.assertNotEqual(contact_map, contact_map_rescaled)
        # ======================================================
        # Test Case 2
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_rescaled = contact_map.rescale(inplace=True)
        self.assertListEqual([1.0, 0.333, 0.0, 0.111], [round(c.raw_score, 3) for c in contact_map])
        self.assertEqual(contact_map, contact_map_rescaled)

    def test_sort(self):
        # ======================================================
        # Test Case 1
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_sorted = contact_map.sort('res1_seq', reverse=False, inplace=False)
        self.assertEqual([(1, 5), (2, 4), (3, 3), (5, 1)], [c.id for c in contact_map_sorted])
        self.assertNotEqual(contact_map, contact_map_sorted)
        # ======================================================
        # Test Case 2
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_sorted = contact_map.sort('res1_seq', reverse=True, inplace=False)
        self.assertEqual([(5, 1), (3, 3), (2, 4), (1, 5)], [c.id for c in contact_map_sorted])
        self.assertNotEqual(contact_map, contact_map_sorted)
        # ======================================================
        # Test Case 3
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_sorted = contact_map.sort('raw_score', reverse=False, inplace=True)
        self.assertEqual([(2, 4), (5, 1), (3, 3), (1, 5)], [c.id for c in contact_map_sorted])
        self.assertEqual(contact_map, contact_map_sorted)
        # ======================================================
        # Test Case 4
        contact_map = ContactMap('test')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map.add(c)
        contact_map_sorted = contact_map.sort('res2_seq', reverse=True, inplace=True)
        self.assertEqual([(1, 5), (2, 4), (3, 3), (5, 1)], [c.id for c in contact_map_sorted])
        self.assertEqual(contact_map, contact_map_sorted)


if __name__ == "__main__":
    unittest.main()
