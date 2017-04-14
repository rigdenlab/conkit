"""Testing facility for conkit.core.contact"""

__author__ = "Felix Simkovic"
__date__ = "12 Aug 2016"

import unittest

from conkit.core import _Entity
from conkit.core import _Gap, _Residue
from conkit.core import SCIPY, SKLEARN
from conkit.core import Contact
from conkit.core import ContactMap
from conkit.core import ContactFile
from conkit.core import Sequence
from conkit.core import SequenceFile


class Test_Entity(unittest.TestCase):

    def test_contains_1(self):
        entity = _Entity('test')
        entity.add(_Entity('foo'))
        self.assertTrue('foo' in entity)

    def test_contains_2(self):
        entity = _Entity('test')
        entity.add(_Entity('foo'))
        entity.add(_Entity('bar'))
        self.assertTrue('foo' in entity)
        self.assertTrue('bar' in entity)
        self.assertFalse('cho' in entity)

    def test_delitem_1(self):
        entity = _Entity('test')
        entity.add(_Entity('foo'))
        del entity['foo']
        self.assertFalse('foo' in entity)

    def test_delitem_2(self):
        entity = _Entity('test')
        entity.add(_Entity('foo'))
        entity.add(_Entity('bar'))
        del entity['foo']
        self.assertFalse('foo' in entity)
        self.assertTrue('bar' in entity)
        del entity['bar']
        self.assertFalse('bar' in entity)

    def test_delitem_3(self):
        entity = _Entity('test')
        with self.assertRaises(KeyError):
            del entity['foo']

    def test_getitem_1(self):
        entity = _Entity('test')
        child_entity = _Entity('foo')
        entity.add(child_entity)
        self.assertEqual(child_entity, entity['foo'])

    def test_getitem_2(self):
        entity = _Entity('test')
        child_entity1 = _Entity('foo')
        child_entity2 = _Entity('bar')
        entity.add(child_entity1)
        entity.add(child_entity2)
        self.assertEqual(child_entity1, entity['foo'])
        self.assertEqual(child_entity2, entity['bar'])

    def test_getitem_3(self):
        entity = _Entity('test')
        child_entity1 = _Entity('foo')
        child_entity2 = _Entity('bar')
        entity.add(child_entity1)
        entity.add(child_entity2)
        self.assertEqual(child_entity1, entity['foo'])
        self.assertEqual(child_entity2, entity['bar'])

    def test_getitem_4(self):
        entity = _Entity('test')
        child_entity1 = _Entity('foo')
        entity.add(child_entity1)
        self.assertEqual(child_entity1, entity[0])

    def test_getitem_5(self):
        entity = _Entity('test')
        child_entity1 = _Entity('foo')
        child_entity2 = _Entity('bar')
        child_entity3 = _Entity('cho')
        entity.add(child_entity1)
        entity.add(child_entity2)
        entity.add(child_entity3)
        self.assertEqual(child_entity1, entity[0])
        self.assertEqual(child_entity2, entity[1])
        self.assertEqual(child_entity2, entity[-2])
        self.assertEqual(child_entity3, entity[-1])

    def test_getitem_6(self):
        entity = _Entity('test')
        for i in range(10):
            entity.add(_Entity('foo_{0}'.format(i)))
        new_entity = entity[:5]
        self.assertEqual(type(entity), type(new_entity))
        self.assertEqual(5, len(new_entity))
        self.assertEqual(['foo_0', 'foo_1', 'foo_2', 'foo_3', 'foo_4'], [e.id for e in new_entity])

    def test_getitem_7(self):
        entity = _Entity('test')
        for i in range(10):
            entity.add(_Entity('foo_{0}'.format(i)))
        new_entity = entity[1::2]
        self.assertEqual(type(entity), type(new_entity))
        self.assertEqual(5, len(new_entity))
        self.assertEqual(['foo_1', 'foo_3', 'foo_5', 'foo_7', 'foo_9'], [e.id for e in new_entity])

    def test_iter_1(self):
        entity = _Entity('test')
        for i in range(10):
            entity.add(_Entity('foo_{0}'.format(i)))
        for i, e in enumerate(entity):
            self.assertEqual('foo_{0}'.format(i), e.id)

    def test_len_1(self):
        entity = _Entity('test')
        for i in range(10):
            entity.add(_Entity('foo_{0}'.format(i)))
        self.assertEqual(10, len(entity))

    def test_len_2(self):
        entity = _Entity('test')
        self.assertEqual(0, len(entity))

    def test_reversed_1(self):
        entity = _Entity('test')
        for i in range(10):
            entity.add(_Entity('foo_{0}'.format(i)))
        rev = list(reversed(list(range(10))))
        for i, e in enumerate(reversed(entity)):
            self.assertNotEqual('foo_{0}'.format(i), e.id)
            self.assertEqual('foo_{0}'.format(rev[i]), e.id)

    def test_child_list_1(self):
        entity = _Entity('test')
        child_entity1 = _Entity('foo')
        entity.add(child_entity1)
        self.assertEqual(1, len(entity.child_list))
        self.assertEqual([child_entity1], entity.child_list)

    def test_child_list_2(self):
        entity = _Entity('test')
        child_entity1 = _Entity('foo')
        child_entity2 = _Entity('bar')
        entity.add(child_entity1)
        entity.add(child_entity2)
        self.assertEqual(2, len(entity.child_list))
        self.assertEqual([child_entity1, child_entity2], entity.child_list)

    def test_child_dict_1(self):
        entity = _Entity('test')
        child_entity1 = _Entity('foo')
        entity.add(child_entity1)
        self.assertDictEqual({'foo': child_entity1}, entity.child_dict)

    def test_child_dict_2(self):
        entity = _Entity('test')
        child_entity1 = _Entity('foo')
        child_entity2 = _Entity('bar')
        entity.add(child_entity1)
        entity.add(child_entity2)
        self.assertDictEqual({'foo': child_entity1, 'bar': child_entity2}, entity.child_dict)

    def test_full_id_1(self):
        entity = _Entity('test')
        entity.add(_Entity('foo'))
        self.assertEqual(('test', 'foo'), entity[0].full_id)

    def test_full_id_2(self):
        entity = _Entity('test')
        self.assertEqual(('test',), entity.full_id)

    def test_full_id_3(self):
        entity = _Entity('test')
        entity.add(_Entity('foo'))
        entity[0].add(_Entity('bar'))
        self.assertEqual(('test', 'foo', 'bar'), entity[0][0].full_id)

    def test_full_id_4(self):
        entity = _Entity('test')
        entity.add(_Entity('foo'))
        entity[0].add(_Entity('bar'))
        entity[0].add(_Entity('cho'))
        self.assertEqual(('test', 'foo', 'bar'), entity[0][0].full_id)
        self.assertEqual(('test', 'foo', 'cho'), entity[0][1].full_id)

    def test_id_1(self):
        entity = _Entity('test')
        self.assertEqual('test', entity.id)

    def test_id_2(self):
        entity = _Entity((1, 2))
        self.assertEqual((1, 2), entity.id)

    def test_id_3(self):
        entity = _Entity([1, 2])
        self.assertEqual((1, 2), entity.id)

    def test_id_4(self):
        entity = _Entity((1., 2.))
        self.assertEqual((1., 2.), entity.id)

    def test_parent_1(self):
        entity = _Entity('test')
        child_entity = _Entity('foo')
        entity.add(child_entity)
        child_child_entity1 = _Entity('bar')
        child_child_entity2 = _Entity('cho')
        child_entity.add(child_child_entity1)
        child_entity.add(child_child_entity2)
        self.assertEqual(entity, entity[0].parent)
        self.assertEqual(child_entity, entity[0][0].parent)
        self.assertEqual(child_entity, entity[0][1].parent)

    def test__inplace_1(self):
        entity = _Entity('foo')
        entity_inplace = entity._inplace(True)
        self.assertEqual(entity, entity_inplace)

    def test__inplace_2(self):
        entity = _Entity('foo')
        entity_inplace = entity._inplace(False)
        self.assertNotEqual(entity, entity_inplace)

    def test__sort_1(self):
        entity = _Entity('test')
        entity.add(_Entity('foo'))
        entity.add(_Entity('bar'))
        entity._sort('id', False)
        self.assertEqual(['bar', 'foo'], [e.id for e in entity])

    def test__sort_2(self):
        entity = _Entity('test')
        entity.add(_Entity('foo'))
        entity.add(_Entity('bar'))
        entity._sort('id', True)
        self.assertEqual(['foo', 'bar'], [e.id for e in entity])

    def test__sort_3(self):
        entity = _Entity('test')
        entity.add(_Entity('foo'))
        entity.add(_Entity('bar'))
        with self.assertRaises(ValueError):
            entity._sort('test', True)

    def test_add_1(self):
        entity = _Entity('test')
        self.assertFalse('foo' in entity.child_dict)
        child_entity = _Entity('foo')
        entity.add(child_entity)
        self.assertTrue(child_entity in entity.child_list)
        self.assertTrue('foo' in entity.child_dict)

    def test_add_2(self):
        entity = _Entity('test')
        child_entity = _Entity('foo')
        entity.add(child_entity)
        self.assertTrue(child_entity in entity.child_list)
        self.assertTrue('foo' in entity.child_dict)
        with self.assertRaises(ValueError):
            entity.add(_Entity('foo'))

    def test_copy_1(self):
        entity = _Entity('test')
        entity.add(_Entity('foo'))
        shallow = entity.copy()
        self.assertNotEqual(entity, shallow)
        self.assertTrue(entity.id, shallow.id)
        self.assertNotEqual(entity[0], shallow[0])
        self.assertEqual(entity[0].id, shallow[0].id)

    def test_copy_2(self):
        entity = _Entity('test')
        entity.add(_Entity('foo'))
        shallow = entity[0].copy()
        self.assertNotEqual(entity[0], shallow)
        self.assertEqual(entity[0].id, shallow.id)
        self.assertIsNone(shallow.parent)
        self.assertEqual('foo', entity[0].id)

    def test_deepcopy_1(self):
        entity = _Entity('test')
        entity.add(_Entity('foo'))
        deep = entity.deepcopy()
        self.assertNotEqual(entity, deep)
        self.assertTrue(entity.id, deep.id)
        self.assertNotEqual(entity[0], deep[0])
        self.assertEqual(entity[0].id, deep[0].id)

    def test_deepcopy_2(self):
        entity = _Entity('test')
        entity.add(_Entity('foo'))
        deep = entity[0].deepcopy()
        self.assertNotEqual(entity[0], deep)
        self.assertEqual(entity[0].id, deep.id)
        self.assertIsNone(deep.parent)
        self.assertEqual('foo', entity[0].id)

    def test_remove_1(self):
        entity = _Entity('test')
        entity.add(_Entity('foo'))
        entity.remove('foo')
        self.assertFalse('foo' in entity)

    def test_remove_2(self):
        entity = _Entity('test')
        entity.add(_Entity('foo'))
        entity.add(_Entity('bar'))
        entity.remove('foo')
        self.assertFalse('foo' in entity)
        self.assertTrue('bar' in entity)
        entity.remove('bar')
        self.assertFalse('bar' in entity)

    def test_remove_3(self):
        entity = _Entity('test')
        with self.assertRaises(KeyError):
            entity.remove('foo')


class TestContact(unittest.TestCase):

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
        self.assertEqual(Contact._UNKNOWN, contact.status)
        contact.define_match()
        self.assertEqual(Contact._MATCH, contact.status)

    def test_status_2(self):
        contact = Contact(1, 2000000, 1.0)
        self.assertEqual(Contact._UNKNOWN, contact.status)
        contact.define_mismatch()
        self.assertEqual(Contact._MISMATCH, contact.status)

    def test_status_3(self):
        contact = Contact(1, 2000000, 1.0)
        self.assertEqual(Contact._UNKNOWN, contact.status)
        contact.define_unknown()
        self.assertEqual(Contact._UNKNOWN, contact.status)

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
        dict = {
            'id': (1, 2), 'is_match': False, 'is_mismatch': False, 'is_unknown': True, 'distance_bound': (0, 8),
            'lower_bound': 0, 'upper_bound': 8, 'raw_score': 1.0, 'res1': 'X', 'res2': 'X', 'res1_chain': '',
            'res2_chain': '', 'res1_seq': 1, 'res2_seq': 2, 'res1_altseq': 0, 'res2_altseq': 0, 'scalar_score': 0.0,
            'status': 0, 'weight': 1.0,
        }
        self.assertEqual(dict, contact._to_dict())

    def test__to_dict_2(self):
        contact = Contact(1, 2, 1.0)
        contact.define_match()
        contact.lower_bound = 4
        dict = {
            'id': (1, 2), 'is_match': True, 'is_mismatch': False, 'is_unknown': False, 'distance_bound': (4, 8),
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
        # ======================================================
        # Test Case 1
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
        for contact, altseq in [(Contact(1, 5, 1.0), (3, 4)),
                                (Contact(2, 4, 0.1), (1, 5)),
                                (Contact(5, 1, 0.2), (4, 3))]:
            contact.res1_altseq = altseq[0]
            contact.res2_altseq = altseq[1]
            contact_map.add(contact)
        contact_map.sequence = Sequence('foo', 'ABCDE')
        self.assertEqual('A-CDE', contact_map.repr_sequence_altloc.seq)
        self.assertEqual('AB-DE', contact_map.repr_sequence.seq)

    def test_repr_sequence_altloc_2(self):
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

    def test_calculate_kernel_density_1(self):
        contact_map1 = ContactMap('foo')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(5, 1, 0.2)]:
            contact_map1.add(c)
        density = contact_map1.calculate_kernel_density()
        self.assertEqual(
            [0.16474084813765252, 0.261112863222513, 0.19866887264998387, 0.06773394884313075], density)

    def test_calculate_kernel_density_2(self):
        contact_map1 = ContactMap('foo')
        for c in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), Contact(3, 4, 0.4)]:
            contact_map1.add(c)
        density = contact_map1.calculate_kernel_density()
        self.assertEqual(
            [0.14936186609839505, 0.2861706643122889, 0.22273876277011645, 0.05660346699219792], density)

    def test_find_1(self):
        contact_map1 = ContactMap('1')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        contact_map1.sequence = Sequence('foo', 'ABCDEF')
        contact_map1.assign_sequence_register()
        found = contact_map1.find([1])
        self.assertEqual(2, len(found))
        self.assertEqual([1, 1], [c.res1_seq for c in found])
        self.assertEqual([5, 4], [c.res2_seq for c in found])

    def test_find_2(self):
        contact_map1 = ContactMap('1')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        contact_map1.sequence = Sequence('foo', 'ABCDEF')
        contact_map1.assign_sequence_register()
        found = contact_map1.find([10])
        self.assertEqual(0, len(found))
        self.assertEqual([], [c.res1_seq for c in found])
        self.assertEqual([], [c.res2_seq for c in found])

    def test_find_3(self):
        contact_map1 = ContactMap('1')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        contact_map1.sequence = Sequence('foo', 'ABCDEF')
        contact_map1.assign_sequence_register()
        found = contact_map1.find([1, 2, 3])
        self.assertEqual(5, len(found))
        self.assertEqual([1, 2, 1, 3, 2], [c.res1_seq for c in found])
        self.assertEqual([5, 6, 4, 6, 5], [c.res2_seq for c in found])

    def test_find_4(self):
        contact_map1 = ContactMap('1')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        contact_map1.sequence = Sequence('foo', 'ABCDEF')
        contact_map1.assign_sequence_register()
        found = contact_map1.find([1, 6])
        self.assertEqual(4, len(found))
        self.assertEqual([1, 2, 1, 3], [c.res1_seq for c in found])
        self.assertEqual([5, 6, 4, 6], [c.res2_seq for c in found])

    def test_assign_sequence_register_1(self):
        contact_map1 = ContactMap('1')
        for comb in [(1, 5, 1.0), (2, 6, 1.0), (1, 4, 1.0), (3, 6, 1.0), (2, 5, 1.0)]:
            contact_map1.add(Contact(*comb))
        contact_map1.sequence = Sequence('foo', 'ABCDEF')
        contact_map1.assign_sequence_register()
        self.assertEqual(
            [('A', 'E'), ('B', 'F'), ('A', 'D'), ('C', 'F'), ('B', 'E')],
            [(c.res1, c.res2) for c in contact_map1]
        )

    def test_assign_sequence_register_2(self):
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
        self.assertEqual(
            [Contact._MATCH, Contact._MISMATCH, Contact._MATCH, Contact._MISMATCH, Contact._UNKNOWN],
            [c.status for c in contact_map1]
        )

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
        self.assertEqual(
            [Contact._MATCH, Contact._MATCH, Contact._MISMATCH, Contact._MATCH],
            [c.status for c in contact_map1]
        )

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
        self.assertEqual(
            [Contact._MISMATCH, Contact._MISMATCH, Contact._MATCH, Contact._MATCH],
            [c.status for c in contact_map1]
        )

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

        contact_map1.match(contact_map2, remove_unmatched=True, inplace=True)
        self.assertEqual(
            [Contact._MATCH, Contact._MATCH],
            [c.status for c in contact_map1]
        )
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
        self.assertEqual(
            [Contact._MATCH, Contact._MATCH, Contact._MISMATCH, Contact._MATCH, Contact._UNKNOWN],
            [c.status for c in contact_map1]
        )

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
        self.assertEqual(
            [Contact._MATCH, Contact._MATCH, Contact._MISMATCH, Contact._MATCH],
            [c.status for c in contact_map1]
        )

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
        self.assertEqual(
            [Contact._MATCH, Contact._MATCH, Contact._MISMATCH, Contact._MATCH, Contact._UNKNOWN],
            [c.status for c in contact_map1]
        )
        self.assertEqual(
            [95, 95, _Gap._IDENTIFIER, 97, _Gap._IDENTIFIER],
            [c.res1_seq for c in contact_map1]
        )
        self.assertEqual(['A', 'A', '', 'A', ''], [c.res1_chain for c in contact_map1])
        self.assertEqual(
            [30, 31, _Gap._IDENTIFIER, 30, _Gap._IDENTIFIER],
            [c.res2_seq for c in contact_map1]
        )
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
        self.assertEqual(
            [Contact._MATCH, Contact._MATCH, Contact._MISMATCH, Contact._MATCH],
            [c.status for c in contact_map1]
        )
        self.assertEqual([95, 95, _Gap._IDENTIFIER, 97], [c.res1_seq for c in contact_map1])
        self.assertEqual(['A', 'A', '', 'A'], [c.res1_chain for c in contact_map1])
        self.assertEqual([30, 31, _Gap._IDENTIFIER, 30], [c.res2_seq for c in contact_map1])
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
        self.assertEqual(
            [Contact._MATCH, Contact._MATCH, Contact._MATCH, Contact._MATCH, Contact._MATCH],
            [c.status for c in contact_map1]
        )
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
        self.assertEqual(
            [Contact._MATCH, Contact._MISMATCH, Contact._MATCH, Contact._MATCH, Contact._MISMATCH],
            [c.status for c in contact_map1]
        )
        self.assertEqual(
                [(6, 2), (6, _Gap._IDENTIFIER), (7, 4), (8, 2), (7, _Gap._IDENTIFIER)],
                [(c.res1_seq, c.res2_seq) for c in contact_map1]
        )
        self.assertEqual(
                [('A', 'B'), ('A', ''), ('A', 'B'), ('A', 'B'), ('A', '')],
                [(c.res1_chain, c.res2_chain) for c in contact_map1]
        )

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
        for contact in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1),
                        Contact(1, 1, 0.2)]:
            contact_map.add(contact)
            contact.res1_altseq = contact.res1_seq
            contact.res2_altseq = contact.res2_seq
        contact_map_keymap = ContactMap._create_keymap(contact_map)
        sequence = [ord(x) for x in 'XXXXX']
        contact_map_keymap = ContactMap._insert_states(sequence, contact_map_keymap)
        contact_map_keymap = ContactMap._reindex(contact_map_keymap)
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
        contact_map_keymap = ContactMap._reindex(contact_map_keymap)
        adjusted = ContactMap._adjust(contact_map, contact_map_keymap)
        self.assertEqual([1, 2, 1], [c.res1_altseq for c in adjusted])
        self.assertEqual([5, 4, 1], [c.res2_altseq for c in adjusted])
    
    def test__create_keymap_1(self):
        contact_map = ContactMap("test")
        for contact in [Contact(1, 5, 1.0), Contact(3, 3, 0.4), Contact(2, 4, 0.1), 
                        Contact(1, 1, 0.2)]:
            contact_map.add(contact)
        contact_map.sequence = Sequence('foo', 'ABCDE')
        contact_map.assign_sequence_register()
        contact_map_keymap = ContactMap._create_keymap(contact_map)
        self.assertEqual([1, 2, 3, 4, 5], [c.res_seq for c in contact_map_keymap])
        self.assertEqual(['A', 'B', 'C', 'D', 'E'], [c.res_name for c in contact_map_keymap])
        self.assertEqual([0, 0, 0, 0, 0], [c.res_altseq for c in contact_map_keymap])

    def test__create_keymap_2(self):
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
        contact_map_keymap = ContactMap._create_keymap(contact_map, altloc=True)
        self.assertEqual([1, 2, 3, 4, 5], [c.res_seq for c in contact_map_keymap])
        self.assertEqual(['A', 'B', 'C', 'D', 'E'], [c.res_name for c in contact_map_keymap])
        self.assertEqual([10, 11, 12, 13, 20], [c.res_altseq for c in contact_map_keymap])

    def test__create_keymap_3(self):
        contact_map = ContactMap("test")
        for contact, res_altloc in [(Contact(1, 5, 1.0), (10, 20)),
                                    (Contact(2, 4, 0.1), (11, 13)),
                                    (Contact(1, 1, 0.2), (10, 10))]:
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
        for contact, res_altloc in [(Contact(1, 5, 1.0), (10, 20)),
                                    (Contact(2, 4, 0.1), (11, 13)),
                                    (Contact(1, 1, 0.2), (10, 10))]:
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
        for contact, res_altloc in [(Contact(1, 5, 1.0), (10, 20)),
                                    (Contact(2, 4, 0.1), (11, 13)),
                                    (Contact(1, 1, 0.2), (10, 10))]:
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

    def test__reindex_1(self):
        keymap = [_Residue(1, 1, 'X', ''), _Residue(2, 2, 'X', ''), _Residue(3, 3, 'X', '')]
        reindex = ContactMap._reindex(keymap)
        self.assertEqual([1, 2, 3], [c.res_seq for c in reindex])
        self.assertEqual([1, 2, 3], [c.res_altseq for c in reindex])

    def test__reindex_2(self):
        keymap = [_Residue(1, -5, 'X', ''), _Residue(2, -4, 'X', ''), _Residue(3, -3, 'X', '')]
        reindex = ContactMap._reindex(keymap)
        self.assertEqual([1, 2, 3], [c.res_seq for c in reindex])
        self.assertEqual([1, 2, 3], [c.res_altseq for c in reindex])

    def test__reindex_3(self):
        keymap = [_Gap(), _Residue(2, 1, 'X', ''), _Residue(3, 2, 'X', '')]
        reindex = ContactMap._reindex(keymap)
        self.assertEqual([_Gap._IDENTIFIER, 2, 3], [c.res_seq for c in reindex])
        self.assertEqual([1, 2, 3], [c.res_altseq for c in reindex])

    def test__reindex_4(self):
        keymap = [_Gap(), _Residue(200000, 10000, 'X', ''), _Gap(), _Gap()]
        reindex = ContactMap._reindex(keymap)
        self.assertEqual(
            [_Gap._IDENTIFIER, 200000, _Gap._IDENTIFIER, _Gap._IDENTIFIER],
            [c.res_seq for c in reindex]
        )
        self.assertEqual([1, 2, 3, 4], [c.res_altseq for c in reindex])

    def test__reindex_5(self):
        keymap = [_Gap(), _Gap(), _Gap(), _Gap()]
        reindex = ContactMap._reindex(keymap)
        self.assertEqual(
            [_Gap._IDENTIFIER, _Gap._IDENTIFIER, _Gap._IDENTIFIER, _Gap._IDENTIFIER],
            [c.res_seq for c in reindex]
        )
        self.assertEqual([1, 2, 3, 4], [c.res_altseq for c in reindex])

    @unittest.skip("Test case missing")
    def test__renumber(self):
        pass


class TestContactFile(unittest.TestCase):

    def test_author_1(self):
        contact_file = ContactFile('test')
        contact_file.author = "John Doe"
        self.assertEqual("John Doe", contact_file.author)

    def test_author_2(self):
        contact_file = ContactFile('test')
        contact_file.author = "John Doe"
        contact_file.author = "Jane Roe"
        self.assertEqual("Jane Roe", contact_file.author)

    def test_method_1(self):
        contact_file = ContactFile('test')
        contact_file.method = 'Hello'
        self.assertEqual(['Hello'], contact_file.method)

    def test_method_2(self):
        contact_file = ContactFile('test')
        contact_file.method = 'Hello'
        contact_file.method = 'World'
        self.assertEqual(['Hello', 'World'], contact_file.method)

    def test_method_3(self):
        contact_file = ContactFile('test')
        contact_file.method = 'Hello'
        contact_file.method = '5'
        contact_file.method = 'World'
        contact_file.method = '!'
        self.assertEqual(['Hello', '5', 'World', '!'], contact_file.method)

    def test_method_4(self):
        contact_file = ContactFile('test')
        self.assertEqual([], contact_file.method)

    def test_method_5(self):
        contact_file = ContactFile('test')
        contact_file.method = 'hello'
        contact_map = ContactMap('foo')
        contact_file.add(contact_map)
        self.assertEqual(['hello'], contact_file.method)

    def test_remark_1(self):
        contact_file = ContactFile('test')
        contact_file.remark = 'Hello'
        self.assertEqual(['Hello'], contact_file.remark)

    def test_remark_2(self):
        contact_file = ContactFile('test')
        contact_file.remark = 'Hello'
        contact_file.remark = 'World'
        self.assertEqual(['Hello', 'World'], contact_file.remark)

    def test_remark_3(self):
        contact_file = ContactFile('test')
        contact_file.remark = 'Hello'
        contact_file.remark = '5'
        contact_file.remark = 'World'
        contact_file.remark = '!'
        self.assertEqual(['Hello', '5', 'World', '!'], contact_file.remark)

    def test_remark_4(self):
        contact_file = ContactFile('test')
        self.assertEqual([], contact_file.remark)

    def test_remark_5(self):
        contact_file = ContactFile('test')
        contact_file.remark = 'hello'
        contact_map = ContactMap('foo')
        contact_file.add(contact_map)
        self.assertEqual(['hello'], contact_file.remark)

    def test_target_1(self):
        contact_file = ContactFile('test')
        contact_file.target = "John Doe"
        self.assertEqual("John Doe", contact_file.target)

    def test_target_2(self):
        contact_file = ContactFile('test')
        contact_file.target = "John Doe"
        contact_file.target = "Jane Roe"
        self.assertEqual("Jane Roe", contact_file.target)

    def test_top_map_1(self):
        contact_file = ContactFile('test')
        self.assertEqual(None, contact_file.top_map)

    def test_top_map_2(self):
        contact_file = ContactFile('test')
        contact_map = ContactMap('foo')
        contact_file.add(contact_map)
        self.assertEqual(contact_map, contact_file.top_map)

    def test_top_map_3(self):
        contact_file = ContactFile('test')
        contact_map1 = ContactMap('foo')
        contact_map2 = ContactMap('bar')
        contact_file.add(contact_map1)
        contact_file.add(contact_map2)
        self.assertEqual(contact_map1, contact_file.top_map)
    
    def test_sort_1(self):
        contact_file = ContactFile('test')
        for map in [ContactMap('foo'), ContactMap('bar'), ContactMap('doe')]:
            contact_file.add(map)
        contact_file_sorted = contact_file.sort('id', reverse=False, inplace=False)
        self.assertEqual(['bar', 'doe', 'foo'], [m.id for m in contact_file_sorted])
        self.assertNotEqual(contact_file, contact_file_sorted)

    def test_sort_2(self):
        contact_file = ContactFile('test')
        for map in [ContactMap('foo'), ContactMap('bar'), ContactMap('doe')]:
            contact_file.add(map)
        contact_file_sorted = contact_file.sort('id', reverse=True, inplace=False)
        self.assertEqual(['foo', 'doe', 'bar'], [m.id for m in contact_file_sorted])
        self.assertNotEqual(contact_file, contact_file_sorted)

    def test_sort_3(self):
        contact_file = ContactFile('test')
        for map in [ContactMap('foo'), ContactMap('bar'), ContactMap('doe')]:
            contact_file.add(map)
        contact_file_sorted = contact_file.sort('id', reverse=False, inplace=True)
        self.assertEqual(['bar', 'doe', 'foo'], [m.id for m in contact_file_sorted])
        self.assertEqual(contact_file, contact_file_sorted)

    def test_sort_4(self):
        contact_file = ContactFile('test')
        for map in [ContactMap('foo'), ContactMap('bar'), ContactMap('doe')]:
            contact_file.add(map)
        contact_file_sorted = contact_file.sort('id', reverse=True, inplace=True)
        self.assertEqual(['foo', 'doe', 'bar'], [m.id for m in contact_file_sorted])
        self.assertEqual(contact_file, contact_file_sorted)


class TestSequence(unittest.TestCase):

    def test_remark_1(self):
        sequence = Sequence('foo', 'GSMFTPK')
        sequence.remark = 'bar'
        self.assertEqual(['bar'], sequence.remark)

    def test_remark_2(self):
        sequence = Sequence('foo', 'GSMFTPK')
        sequence.remark = 'bar'
        sequence.remark = 'baz'
        self.assertEqual(['bar', 'baz'], sequence.remark)

    def test_seq_1(self):
        sequence = Sequence('foo', 'GSMFTPK')
        self.assertEqual('foo', sequence.id)
        self.assertEqual('GSMFTPK', sequence.seq)

    def test_seq_2(self):
        sequence = Sequence('foo', 'GSMFTPK')
        sequence.seq = 'AAAAAA'
        self.assertEqual('foo', sequence.id)
        self.assertEqual('AAAAAA', sequence.seq)

    def test_seq_3(self):
        sequence = Sequence('foo', 'GSMFTPK')
        with self.assertRaises(ValueError):
            sequence.seq = 'A2A'

    def test_seq_4(self):
        sequence = Sequence('foo', 'GSMFTPK')
        sequence.seq = '-------'

    def test_seq_len_1(self):
        sequence = Sequence('foo', 'GSMFTPK')
        self.assertEqual('foo', sequence.id)
        self.assertEqual('GSMFTPK', sequence.seq)
        self.assertEqual(7, sequence.seq_len)

    def test_seq_len_2(self):
        sequence = Sequence('foo', 'GSMFTPK')
        self.assertEqual('foo', sequence.id)
        self.assertEqual('GSMFTPK', sequence.seq)
        self.assertEqual(7, sequence.seq_len)
        sequence.seq = 'AAAAAAAAAA'
        self.assertEqual('foo', sequence.id)
        self.assertEqual('AAAAAAAAAA', sequence.seq)
        self.assertEqual(10, sequence.seq_len)

    def test_align_local_1(self):
        sequence1 = Sequence('foo', 'GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTI'
                                    'GYFKSELEKEPLRVIPLKEVHKVQECKQSDIMMRDNLFEIVT'
                                    'TSRTFYVQADSPEEMHSWIKAVSGAIVAQRGPGRSASSEHP')
        sequence2 = Sequence('bar', 'Q-------YF-------P------------------------'
                                    '--F----------VQADSPEEMHSWIKAVSGAIVAQR')
        sequence1.align_local(sequence2, id_chars=2, nonid_chars=1, gap_open_pen=-0.5, gap_ext_pen=-0.2, inplace=True)
        aligned1 = "GSMFTPKPPQDSAVIKAGYCVKQGAVMKNWKRRYFQLDENTIGYFKSELEKEPLRVIPLKEVHKVQECKQSDIMM" \
                   "RDNLFEIVTTSRTFYVQADSPEEMHSWIKAVSGAIVAQRGPGRSASSEHP"
        aligned2 = "-----------------------------------Q-------YF-------P----------------------" \
                   "----F----------VQADSPEEMHSWIKAVSGAIVAQR-----------"
        self.assertEqual(aligned1, sequence1.seq)
        self.assertEqual(aligned2, sequence2.seq)

    def test_align_local_2(self):
        sequence1 = Sequence('foo', 'DDLTISSLAKGETTKAAFNQMVQGHKLPAWVMKGGTYTPAQTV'
                                    'TLGDETYQVMSACKPHDCGSQRIAVMWSEKSNQMTGLFSTIDE'
                                    'KTSQEKLTWLNVNDALSIDGKTVLFAALTGSLENHPDGFNFKV'
                                    'FGRCELAAAMKRHGLDNYRGYSLGNWVCAAKFESNFNTQATNR'
                                    'NTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDI'
                                    'TASVNCAKKIVSDGNGMNAWVAWRNRCKGTDVQAWIRGCR')
        sequence2 = Sequence('bar', '-------------------------------------------'
                                    '-------------------------------------------'
                                    '--------W------------TV--------------------'
                                    'F--C----AM---GLD-----------C--KFE-NF-------'
                                    'N-D-----G---------C-D----G--NLC-IP--------I'
                                    '--------------NG--------------D----IRGC-')
        sequence1.align_local(sequence2, id_chars=2, nonid_chars=1, gap_open_pen=-0.5, gap_ext_pen=-0.2, inplace=True)
        aligned1 = "DDLTISSLAKGETTKAAFNQMVQGHKLPAWVMKGGTYTPAQTVTLGDETYQVMSACKPHDCGSQRIAVMWSEKSN" \
                   "QMTGLFSTIDEKTSQEKLTWLNVNDALSIDGKTVLFAALTGSLENHPDGFNFKVFGRCELAAAMKRHGLDNYRGY" \
                   "SLGNWVCAAKFESNFNTQATNRNTDGSTDYGILQINSRWWCNDGRTPGSRNLCNIPCSALLSSDITASVNCAKKI" \
                   "VSDGNGMNAWVAWRNRCKGTDVQAWIRGCR"
        aligned2 = "---------------------------------------------------------------------------" \
                   "-------------------W------------TV--------------------F--C----AM---GLD-----" \
                   "------C--KFE-NF-------N-D-----G---------C-D----G--NLC-IP--------I----------" \
                   "----NG--------------D----IRGC-"
        self.assertEqual(aligned1, sequence1.seq)
        self.assertEqual(aligned2, sequence2.seq)

    def test_align_local_3(self):
        sequence1 = Sequence('foo', '------------------------------------------'
                                    '------------------------------------------'
                                    '----------W------------TV-----------------'
                                    '---F--C----AM---GLD-----------C--KFE-NF---'
                                    '----N-D-----G---------C-D----G--NLC-IP----'
                                    '----I--------------NG--------------D----IR'
                                    'GC-')
        sequence2 = Sequence('bar', '-D-------------------------------GGTYTP---'
                                    '-------------C-PHDCGS-R-------------------'
                                    '------------------------------TG--EN------'
                                    '-KV------------------------------KFESN-N-Q'
                                    'ATNR------D----Q--------------------------'
                                    '------------------------WVA--NR-----------'
                                    '---')
        sequence1.align_local(sequence2, id_chars=2, nonid_chars=1, gap_open_pen=-1.0, gap_ext_pen=-0.5, inplace=True)
        aligned1 = "---------------------------------------------------------------------------" \
                   "-------------------W------------TV--------------------F--C----AM---GLD-----" \
                   "------C--KFE-NF-------N-D-----G---------C-D----G--NLC-IP--------I----------" \
                   "----NG--------------D----IRGC-"
        aligned2 = "-D-------------------------------GGTYTP----------------C-PHDCGS-R----------" \
                   "---------------------------------------TG--EN-------KV---------------------" \
                   "---------KFESN-N-QATNR------D----Q-----------------------------------------" \
                   "---------WVA--NR--------------"
        self.assertEqual(aligned1, sequence1.seq)
        self.assertEqual(aligned2, sequence2.seq)


class TestSequenceFile(unittest.TestCase):

    def test_is_alignment_1(self):
        sequence_file = SequenceFile('test')
        sequence_file.add(Sequence('foo', 'AAAAA'))
        sequence_file.add(Sequence('bar', 'BBBBB'))
        self.assertTrue(sequence_file.is_alignment)

    def test_is_alignment_2(self):
        sequence_file = SequenceFile('test')
        sequence_file.add(Sequence('foo', 'AAAAA'))
        sequence_file.add(Sequence('bar', 'BBBB'))
        self.assertFalse(sequence_file.is_alignment)

    def test_nseqs_1(self):
        sequence_file = SequenceFile('test')
        self.assertEqual(0, sequence_file.nseqs)

    def test_nseqs_2(self):
        sequence_file = SequenceFile('test')
        sequence_file.add(Sequence('foo', 'AAAAA'))
        self.assertEqual(1, sequence_file.nseqs)

    def test_nseqs_3(self):
        sequence_file = SequenceFile('test')
        sequence_file.add(Sequence('foo', 'AAAAA'))
        sequence_file.add(Sequence('bar', 'BBBBB'))
        self.assertEqual(2, sequence_file.nseqs)

    def test_remark_1(self):
        sequence_file = SequenceFile('test')
        sequence_file.remark = 'Hello'
        self.assertEqual(['Hello'], sequence_file.remark)

    def test_remark_2(self):
        sequence_file = SequenceFile('test')
        sequence_file.remark = 'Hello'
        sequence_file.remark = 'World'
        self.assertEqual(['Hello', 'World'], sequence_file.remark)

    def test_remark_3(self):
        sequence_file = SequenceFile('test')
        sequence_file.remark = 'Hello'
        sequence_file.remark = '5'
        sequence_file.remark = 'World'
        sequence_file.remark = '!'
        self.assertEqual(['Hello', '5', 'World', '!'], sequence_file.remark)

    def test_remark_4(self):
        sequence_file = SequenceFile('test')
        self.assertEqual([], sequence_file.remark)

    def test_remark_5(self):
        sequence_file = SequenceFile('test')
        sequence_file.remark = 'hello'
        sequence = Sequence('foo', 'GSMFTPK')
        sequence.remark = 'bar'
        sequence_file.add(sequence)
        self.assertEqual(['hello'], sequence_file.remark)
        self.assertEqual(['bar'], sequence_file[0].remark)

    def test_top_sequence_1(self):
        sequence_file = SequenceFile('test')
        self.assertEqual(None, sequence_file.top_sequence)

    def test_top_sequence_2(self):
        sequence_file = SequenceFile('test')
        sequence1 = Sequence('foo', 'AAAAA')
        sequence_file.add(sequence1)
        self.assertEqual(sequence1, sequence_file.top_sequence)

    def test_top_sequence_3(self):
        sequence_file = SequenceFile('test')
        sequence1 = Sequence('foo', 'AAAAA')
        sequence2 = Sequence('bar', 'BBBBB')
        sequence_file.add(sequence1)
        sequence_file.add(sequence2)
        self.assertEqual(sequence1, sequence_file.top_sequence)

    @unittest.skipUnless(SCIPY, "SciPy not installed")
    def test_calculate_meff_1(self):
        sequence_file = SequenceFile('test')
        for s in [Sequence('foo', 'AAAAAAA'), Sequence('bar', 'AAAAAAA'),
                  Sequence('cho', 'AAAAAAA'), Sequence('baz', 'AAAAAAA')]:
            sequence_file.add(s)
        m_eff = sequence_file.calculate_meff(identity=0.7)
        self.assertTrue(isinstance(m_eff, int))
        self.assertEqual(1, m_eff)

    @unittest.skipUnless(SCIPY, "SciPy not installed")
    def test_calculate_meff_2(self):
        sequence_file = SequenceFile('test')
        for s in [Sequence('foo', 'AAAAAAA'), Sequence('bar', 'AAAAAAA'),
                  Sequence('cho', 'AAAAAAA'), Sequence('baz', 'BBBBBBB')]:
            sequence_file.add(s)
        m_eff = sequence_file.calculate_meff(identity=0.7)
        self.assertTrue(isinstance(m_eff, int))
        self.assertEqual(2, m_eff)

    @unittest.skipUnless(SCIPY, "SciPy not installed")
    def test_calculate_meff_3(self):
        sequence_file = SequenceFile('test')
        for s in [Sequence('foo', 'AAAAAAA'), Sequence('bar', 'A-AABA-'),
                  Sequence('cho', 'B-BAA--'), Sequence('baz', 'BBBBBBB')]:
            sequence_file.add(s)
        m_eff = sequence_file.calculate_meff(identity=0.7)
        self.assertTrue(isinstance(m_eff, int))
        self.assertEqual(4, m_eff)

    @unittest.skipUnless(SCIPY, "SciPy not installed")
    def test_calculate_meff_4(self):
        sequence_file = SequenceFile('test')
        for s in [Sequence('foo', 'AAAAAAA'), Sequence('bar', 'AAAABA-'),
                  Sequence('cho', 'B-BAA--'), Sequence('baz', 'BBBBBBB')]:
            sequence_file.add(s)
        m_eff = sequence_file.calculate_meff(identity=0.7)
        self.assertTrue(isinstance(m_eff, int))
        self.assertEqual(3, m_eff)

    @unittest.skipUnless(SCIPY, "SciPy not installed")
    def test_calculate_meff_5(self):
        sequence_file = SequenceFile('test')
        for s in [Sequence('foo', 'AAAAAAA'), Sequence('bar', 'AA-ABA-'),
                  Sequence('cho', 'B-BAA--'), Sequence('baz', 'BBBBBBB')]:
            sequence_file.add(s)
        m_eff = sequence_file.calculate_meff(identity=0.6)
        self.assertTrue(isinstance(m_eff, int))
        self.assertEqual(4, m_eff)
        self.assertNotEqual(3, m_eff)
        self.assertNotEqual(3, m_eff)

    @unittest.skipUnless(SCIPY, "SciPy not installed")
    def test_calculate_meff_6(self):
        sequence_file = SequenceFile('test')
        for s in [Sequence('foo', 'AAAAAAA'), Sequence('bar', 'AA-ABA-'),
                  Sequence('cho', 'AAACBAA'), Sequence('doo', 'B-BAA--'),
                  Sequence('miu', 'BBBBBBB'), Sequence('nop', 'AAAAAAB')]:
            sequence_file.add(s)
        m_eff = sequence_file.calculate_meff(identity=0.6)
        self.assertTrue(isinstance(m_eff, int))
        self.assertEqual(4, m_eff)

    def test_calculate_freq_1(self):
        sequence_file = SequenceFile('test')
        for s in [Sequence('foo', 'AAAAAAA'), Sequence('bar', 'A-AAAA-'), Sequence('cho', '--AAA--')]:
            sequence_file.add(s)
        calculated_freqs = [round(i, 6) for i in sequence_file.calculate_freq()]
        self.assertEqual([0.666667, 0.333333, 1.0, 1.0, 1.0, 0.666667, 0.333333], calculated_freqs)

    def test_calculate_freq_2(self):
        sequence_file = SequenceFile('test')
        for s in [Sequence('foo', '-------'), Sequence('bar', '-------'), Sequence('cho', '-------')]:
            sequence_file.add(s)
        calculated_freqs = [round(i, 6) for i in sequence_file.calculate_freq()]
        self.assertEqual([0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0], calculated_freqs)

    def test_calculate_freq_3(self):
        sequence_file = SequenceFile('test')
        for s in [Sequence('foo', 'AAAAAAA'), Sequence('bar', 'AAAAAAA'), Sequence('cho', 'AAAAAAA')]:
            sequence_file.add(s)
        calculated_freqs = [round(i, 6) for i in sequence_file.calculate_freq()]
        self.assertEqual([1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0], calculated_freqs)

    def test_sort_1(self):
        sequence_file = SequenceFile('test')
        for seq in [Sequence('foo', 'AAAAA'), Sequence('bar', 'BBBBB'), Sequence('doe', 'CCCCC')]:
            sequence_file.add(seq)
        sequence_file_sorted = sequence_file.sort('id', reverse=False, inplace=False)
        self.assertEqual(['bar', 'doe', 'foo'], [s.id for s in sequence_file_sorted])
        self.assertEqual(['BBBBB', 'CCCCC', 'AAAAA'], [s.seq for s in sequence_file_sorted])
        self.assertNotEqual(sequence_file, sequence_file_sorted)

    def test_sort_2(self):
        sequence_file = SequenceFile('test')
        for seq in [Sequence('foo', 'AAAAA'), Sequence('bar', 'BBBBB'), Sequence('doe', 'CCCCC')]:
            sequence_file.add(seq)
        sequence_file_sorted = sequence_file.sort('id', reverse=True, inplace=False)
        self.assertEqual(['foo', 'doe', 'bar'], [s.id for s in sequence_file_sorted])
        self.assertEqual(['AAAAA', 'CCCCC', 'BBBBB'], [s.seq for s in sequence_file_sorted])
        self.assertNotEqual(sequence_file, sequence_file_sorted)

    def test_sort_3(self):
        sequence_file = SequenceFile('test')
        for seq in [Sequence('foo', 'AAAAA'), Sequence('bar', 'BBBBB'), Sequence('doe', 'CCCCC')]:
            sequence_file.add(seq)
        sequence_file_sorted = sequence_file.sort('seq', reverse=False, inplace=True)
        self.assertEqual(['foo', 'bar', 'doe'], [s.id for s in sequence_file_sorted])
        self.assertEqual(['AAAAA', 'BBBBB', 'CCCCC'], [s.seq for s in sequence_file_sorted])
        self.assertEqual(sequence_file, sequence_file_sorted)

    def test_sort_4(self):
        sequence_file = SequenceFile('test')
        for seq in [Sequence('foo', 'AAAAA'), Sequence('bar', 'BBBBB'), Sequence('doe', 'CCCCC')]:
            sequence_file.add(seq)
        sequence_file_sorted = sequence_file.sort('seq', reverse=True, inplace=True)
        self.assertEqual(['doe', 'bar', 'foo'], [s.id for s in sequence_file_sorted])
        self.assertEqual(['CCCCC', 'BBBBB', 'AAAAA'], [s.seq for s in sequence_file_sorted])
        self.assertEqual(sequence_file, sequence_file_sorted)

    def test_trim_1(self):
        sequence_file = SequenceFile('test')
        for seq in [Sequence('foo', 'AAAAA'), Sequence('bar', 'BBBBB'), Sequence('doe', 'CCCCC')]:
            sequence_file.add(seq)
        sequence_file_trimmed = sequence_file.trim(1, 5)
        self.assertEqual(['foo', 'bar', 'doe'], [s.id for s in sequence_file_trimmed])
        self.assertEqual(['AAAAA', 'BBBBB', 'CCCCC'], [s.seq for s in sequence_file_trimmed])
        self.assertNotEqual(sequence_file, sequence_file_trimmed)

    def test_trim_2(self):
        sequence_file = SequenceFile('test')
        for seq in [Sequence('foo', 'AAAAA'), Sequence('bar', 'BBBBB'), Sequence('doe', 'CCCCC')]:
            sequence_file.add(seq)
        sequence_file_trimmed = sequence_file.trim(3, 5)
        self.assertEqual(['foo', 'bar', 'doe'], [s.id for s in sequence_file_trimmed])
        self.assertEqual(['AAA', 'BBB', 'CCC'], [s.seq for s in sequence_file_trimmed])
        self.assertNotEqual(sequence_file, sequence_file_trimmed)

    def test_trim_3(self):
        sequence_file = SequenceFile('test')
        for seq in [Sequence('foo', 'ABCDE'), Sequence('bar', 'BCDEF'), Sequence('doe', 'CDEFG')]:
            sequence_file.add(seq)
        sequence_file_trimmed = sequence_file.trim(1, 3)
        self.assertEqual(['foo', 'bar', 'doe'], [s.id for s in sequence_file_trimmed])
        self.assertEqual(['ABC', 'BCD', 'CDE'], [s.seq for s in sequence_file_trimmed])
        self.assertNotEqual(sequence_file, sequence_file_trimmed)

    def test_trim_4(self):
        sequence_file = SequenceFile('test')
        for seq in [Sequence('foo', 'ABCDE'), Sequence('bar', 'BCDEF'), Sequence('doe', 'CDEFG')]:
            sequence_file.add(seq)
        sequence_file_trimmed = sequence_file.trim(2, 3)
        self.assertEqual(['foo', 'bar', 'doe'], [s.id for s in sequence_file_trimmed])
        self.assertEqual(['BC', 'CD', 'DE'], [s.seq for s in sequence_file_trimmed])
        self.assertNotEqual(sequence_file, sequence_file_trimmed)


if __name__ == "__main__":
    unittest.main(verbosity=2)
