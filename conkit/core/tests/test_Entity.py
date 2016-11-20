"""Testing facility for conkit.core.entity"""

from builtins import range
__author__ = "Felix Simkovic"
__date__ = "12 Aug 2016"

from conkit.core.Entity import Entity

import unittest


class Test(unittest.TestCase):

    def test_contains(self):
        # ======================================================
        # Test Case 1
        entity = Entity('test')
        entity.add(Entity('foo'))
        self.assertTrue('foo' in entity)
        # ======================================================
        # Test Case 2
        entity = Entity('test')
        entity.add(Entity('foo'))
        entity.add(Entity('bar'))
        self.assertTrue('foo' in entity)
        self.assertTrue('bar' in entity)
        self.assertFalse('cho' in entity)

    def test_delitem(self):
        # ======================================================
        # Test Case 1
        entity = Entity('test')
        entity.add(Entity('foo'))
        del entity['foo']
        self.assertFalse('foo' in entity)
        # ======================================================
        # Test Case 2
        entity = Entity('test')
        entity.add(Entity('foo'))
        entity.add(Entity('bar'))
        del entity['foo']
        self.assertFalse('foo' in entity)
        self.assertTrue('bar' in entity)
        del entity['bar']
        self.assertFalse('bar' in entity)
        # ======================================================
        # Test Case 3
        entity = Entity('test')
        with self.assertRaises(KeyError):
            del entity['foo']

    def test_getitem(self):
        # ======================================================
        # Test Case 1
        entity = Entity('test')
        child_entity = Entity('foo')
        entity.add(child_entity)
        self.assertEqual(child_entity, entity['foo'])
        # ======================================================
        # Test Case 2
        entity = Entity('test')
        child_entity1 = Entity('foo')
        child_entity2 = Entity('bar')
        entity.add(child_entity1)
        entity.add(child_entity2)
        self.assertEqual(child_entity1, entity['foo'])
        self.assertEqual(child_entity2, entity['bar'])
        # ======================================================
        # Test Case 3
        entity = Entity('test')
        child_entity1 = Entity('foo')
        child_entity2 = Entity('bar')
        entity.add(child_entity1)
        entity.add(child_entity2)
        self.assertEqual(child_entity1, entity['foo'])
        self.assertEqual(child_entity2, entity['bar'])
        # ======================================================
        # Test Case 4
        entity = Entity('test')
        child_entity1 = Entity('foo')
        entity.add(child_entity1)
        self.assertEqual(child_entity1, entity[0])
        # ======================================================
        # Test Case 4
        entity = Entity('test')
        child_entity1 = Entity('foo')
        child_entity2 = Entity('bar')
        child_entity3 = Entity('cho')
        entity.add(child_entity1)
        entity.add(child_entity2)
        entity.add(child_entity3)
        self.assertEqual(child_entity1, entity[0])
        self.assertEqual(child_entity2, entity[1])
        self.assertEqual(child_entity2, entity[-2])
        self.assertEqual(child_entity3, entity[-1])
        # ======================================================
        # Test Case 5
        entity = Entity('test')
        for i in range(10):
            entity.add(Entity('foo_{0}'.format(i)))
        new_entity = entity[:5]
        self.assertEqual(type(entity), type(new_entity))
        self.assertEqual(5, len(new_entity))
        self.assertEqual(['foo_0', 'foo_1', 'foo_2', 'foo_3', 'foo_4'], [e.id for e in new_entity])
        # ======================================================
        # Test Case 6
        entity = Entity('test')
        for i in range(10):
            entity.add(Entity('foo_{0}'.format(i)))
        new_entity = entity[1::2]
        self.assertEqual(type(entity), type(new_entity))
        self.assertEqual(5, len(new_entity))
        self.assertEqual(['foo_1', 'foo_3', 'foo_5', 'foo_7', 'foo_9'], [e.id for e in new_entity])

    def test_iter(self):
        # ======================================================
        # Test Case 1
        entity = Entity('test')
        for i in range(10):
            entity.add(Entity('foo_{0}'.format(i)))
        for i, e in enumerate(entity):
            self.assertEqual('foo_{0}'.format(i), e.id)

    def test_len(self):
        # ======================================================
        # Test Case 1
        entity = Entity('test')
        for i in range(10):
            entity.add(Entity('foo_{0}'.format(i)))
        self.assertEqual(10, len(entity))
        # ======================================================
        # Test Case 2
        entity = Entity('test')
        self.assertEqual(0, len(entity))

    def test_reversed(self):
        # ======================================================
        # Test Case 1
        entity = Entity('test')
        for i in range(10):
            entity.add(Entity('foo_{0}'.format(i)))
        rev = list(reversed(list(range(10))))
        for i, e in enumerate(reversed(entity)):
            self.assertNotEqual('foo_{0}'.format(i), e.id)
            self.assertEqual('foo_{0}'.format(rev[i]), e.id)

    def test_child_list(self):
        # ======================================================
        # Test Case 1
        entity = Entity('test')
        child_entity1 = Entity('foo')
        entity.add(child_entity1)
        self.assertEqual(1, len(entity.child_list))
        self.assertEqual([child_entity1], entity.child_list)
        # ======================================================
        # Test Case 2
        entity = Entity('test')
        child_entity1 = Entity('foo')
        child_entity2 = Entity('bar')
        entity.add(child_entity1)
        entity.add(child_entity2)
        self.assertEqual(2, len(entity.child_list))
        self.assertEqual([child_entity1, child_entity2], entity.child_list)

    def test_child_dict(self):
        # ======================================================
        # Test Case 1
        entity = Entity('test')
        child_entity1 = Entity('foo')
        entity.add(child_entity1)
        self.assertDictEqual({'foo': child_entity1}, entity.child_dict)
        # ======================================================
        # Test Case 2
        entity = Entity('test')
        child_entity1 = Entity('foo')
        child_entity2 = Entity('bar')
        entity.add(child_entity1)
        entity.add(child_entity2)
        self.assertDictEqual({'foo': child_entity1, 'bar': child_entity2}, entity.child_dict)

    def test_full_id(self):
        # ======================================================
        # Test Case 1
        entity = Entity('test')
        entity.add(Entity('foo'))
        self.assertEqual(('test', 'foo'), entity[0].full_id)
        # ======================================================
        # Test Case 2
        entity = Entity('test')
        self.assertEqual(('test',), entity.full_id)
        # ======================================================
        # Test Case 3
        entity = Entity('test')
        entity.add(Entity('foo'))
        entity[0].add(Entity('bar'))
        self.assertEqual(('test', 'foo', 'bar'), entity[0][0].full_id)
        # ======================================================
        # Test Case 4
        entity = Entity('test')
        entity.add(Entity('foo'))
        entity[0].add(Entity('bar'))
        entity[0].add(Entity('cho'))
        self.assertEqual(('test', 'foo', 'bar'), entity[0][0].full_id)
        self.assertEqual(('test', 'foo', 'cho'), entity[0][1].full_id)

    def test_id(self):
        # ======================================================
        # Test Case 1
        entity = Entity('test')
        self.assertEqual('test', entity.id)
        # ======================================================
        # Test Case 2
        # self.assertRaises(TypeError, Entity, 3)
        # self.assertRaises(TypeError, Entity, 1.0)
        # ======================================================
        # Test Case 3
        entity = Entity((1, 2))
        self.assertEqual((1, 2), entity.id)
        # ======================================================
        # Test Case 4
        entity = Entity([1, 2])
        self.assertEqual((1, 2), entity.id)
        # ======================================================
        # Test Case 5
        entity = Entity((1., 2.))
        self.assertEqual((1., 2.), entity.id)

    def test_parent(self):
        # ======================================================
        # Test Case 1
        entity = Entity('test')
        child_entity = Entity('foo')
        entity.add(child_entity)
        child_child_entity1 = Entity('bar')
        child_child_entity2 = Entity('cho')
        child_entity.add(child_child_entity1)
        child_entity.add(child_child_entity2)
        self.assertEqual(entity, entity[0].parent)
        self.assertEqual(child_entity, entity[0][0].parent)
        self.assertEqual(child_entity, entity[0][1].parent)

    def test__inplace(self):
        # ======================================================
        # Test Case 1
        entity = Entity('foo')
        entity_inplace = entity._inplace(True)
        self.assertEqual(entity, entity_inplace)
        # ======================================================
        # Test Case 2
        entity = Entity('foo')
        entity_inplace = entity._inplace(False)
        self.assertNotEqual(entity, entity_inplace)

    def test__sort(self):
        # ======================================================
        # Test Case 1
        entity = Entity('test')
        entity.add(Entity('foo'))
        entity.add(Entity('bar'))
        entity._sort('id', False)
        self.assertEqual(['bar', 'foo'], [e.id for e in entity])
        # ======================================================
        # Test Case 2
        entity = Entity('test')
        entity.add(Entity('foo'))
        entity.add(Entity('bar'))
        entity._sort('id', True)
        self.assertEqual(['foo', 'bar'], [e.id for e in entity])
        # ======================================================
        # Test Case 3
        entity = Entity('test')
        entity.add(Entity('foo'))
        entity.add(Entity('bar'))
        with self.assertRaises(ValueError):
            entity._sort('test', True)

    def test_add(self):
        # ======================================================
        # Test Case 1
        entity = Entity('test')
        self.assertFalse('foo' in entity.child_dict)
        child_entity = Entity('foo')
        entity.add(child_entity)
        self.assertTrue(child_entity in entity.child_list)
        self.assertTrue('foo' in entity.child_dict)
        # ======================================================
        # Test Case 2
        entity = Entity('test')
        child_entity = Entity('foo')
        entity.add(child_entity)
        self.assertTrue(child_entity in entity.child_list)
        self.assertTrue('foo' in entity.child_dict)
        with self.assertRaises(ValueError):
            entity.add(Entity('foo'))

    def test_copy(self):
        # ======================================================
        # Test Case 1
        entity = Entity('test')
        entity.add(Entity('foo'))
        shallow = entity.copy()
        self.assertNotEqual(entity, shallow)
        self.assertTrue(entity.id, shallow.id)
        self.assertNotEqual(entity[0], shallow[0])
        self.assertEqual(entity[0].id, shallow[0].id)
        # ======================================================
        # Test Case 2
        entity = Entity('test')
        entity.add(Entity('foo'))
        shallow = entity[0].copy()
        self.assertNotEqual(entity[0], shallow)
        self.assertEqual(entity[0].id, shallow.id)
        self.assertIsNone(shallow.parent)
        self.assertEqual('foo', entity[0].id)

    def test_deepcopy(self):
        # ======================================================
        # Test Case 1
        entity = Entity('test')
        entity.add(Entity('foo'))
        deep = entity.deepcopy()
        self.assertNotEqual(entity, deep)
        self.assertTrue(entity.id, deep.id)
        self.assertNotEqual(entity[0], deep[0])
        self.assertEqual(entity[0].id, deep[0].id)
        # ======================================================
        # Test Case 2
        entity = Entity('test')
        entity.add(Entity('foo'))
        deep = entity[0].deepcopy()
        self.assertNotEqual(entity[0], deep)
        self.assertEqual(entity[0].id, deep.id)
        self.assertIsNone(deep.parent)
        self.assertEqual('foo', entity[0].id)

    def test_remove(self):
        # ======================================================
        # Test Case 1
        entity = Entity('test')
        entity.add(Entity('foo'))
        entity.remove('foo')
        self.assertFalse('foo' in entity)
        # ======================================================
        # Test Case 2
        entity = Entity('test')
        entity.add(Entity('foo'))
        entity.add(Entity('bar'))
        entity.remove('foo')
        self.assertFalse('foo' in entity)
        self.assertTrue('bar' in entity)
        entity.remove('bar')
        self.assertFalse('bar' in entity)
        # ======================================================
        # Test Case 3
        entity = Entity('test')
        with self.assertRaises(KeyError):
            entity.remove('foo')


if __name__ == "__main__":
    unittest.main(verbosity=2)
