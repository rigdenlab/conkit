"""Testing facility for conkit.core.contactmap"""

__author__ = "Felix Simkovic"
__date__ = "12 Aug 2016"

from conkit.core.ContactFile import ContactFile
from conkit.core.ContactMap import ContactMap

import unittest
        

class Test(unittest.TestCase):

    def test_author(self):
        # ======================================================
        # Test Case 1
        contact_file = ContactFile('test')
        contact_file.author = "John Doe"
        self.assertEqual("John Doe", contact_file.author)
        # ======================================================
        # Test Case 2
        contact_file = ContactFile('test')
        contact_file.author = "John Doe"
        contact_file.author = "Jane Roe"
        self.assertEqual("Jane Roe", contact_file.author)

    def test_method(self):
        # ======================================================
        # Test Case 1
        contact_file = ContactFile('test')
        contact_file.method = 'Hello'
        self.assertEqual(['Hello'], contact_file.method)
        # ======================================================
        # Test Case 2
        contact_file = ContactFile('test')
        contact_file.method = 'Hello'
        contact_file.method = 'World'
        self.assertEqual(['Hello', 'World'], contact_file.method)
        # ======================================================
        # Test Case 3
        contact_file = ContactFile('test')
        contact_file.method = 'Hello'
        contact_file.method = '5'
        contact_file.method = 'World'
        contact_file.method = '!'
        self.assertEqual(['Hello', '5', 'World', '!'], contact_file.method)
        # ======================================================
        # Test Case 4
        contact_file = ContactFile('test')
        self.assertEqual([], contact_file.method)
        # ======================================================
        # Test Case 5
        contact_file = ContactFile('test')
        contact_file.method = 'hello'
        contact_map = ContactMap('foo')
        contact_file.add(contact_map)
        self.assertEqual(['hello'], contact_file.method)

    def test_remark(self):
        # ======================================================
        # Test Case 1
        contact_file = ContactFile('test')
        contact_file.remark = 'Hello'
        self.assertEqual(['Hello'], contact_file.remark)
        # ======================================================
        # Test Case 2
        contact_file = ContactFile('test')
        contact_file.remark = 'Hello'
        contact_file.remark = 'World'
        self.assertEqual(['Hello', 'World'], contact_file.remark)
        # ======================================================
        # Test Case 3
        contact_file = ContactFile('test')
        contact_file.remark = 'Hello'
        contact_file.remark = '5'
        contact_file.remark = 'World'
        contact_file.remark = '!'
        self.assertEqual(['Hello', '5', 'World', '!'], contact_file.remark)
        # ======================================================
        # Test Case 4
        contact_file = ContactFile('test')
        self.assertEqual([], contact_file.remark)
        # ======================================================
        # Test Case 5
        contact_file = ContactFile('test')
        contact_file.remark = 'hello'
        contact_map = ContactMap('foo')
        contact_file.add(contact_map)
        self.assertEqual(['hello'], contact_file.remark)

    def test_target(self):
        # ======================================================
        # Test Case 1
        contact_file = ContactFile('test')
        contact_file.target = "John Doe"
        self.assertEqual("John Doe", contact_file.target)
        # ======================================================
        # Test Case 2
        contact_file = ContactFile('test')
        contact_file.target = "John Doe"
        contact_file.target = "Jane Roe"
        self.assertEqual("Jane Roe", contact_file.target)

    def test_top_map(self):
        # ======================================================
        # Test Case 1
        contact_file = ContactFile('test')
        self.assertEqual(None, contact_file.top_map)
        # ======================================================
        # Test Case 2
        contact_file = ContactFile('test')
        contact_map = ContactMap('foo')
        contact_file.add(contact_map)
        self.assertEqual(contact_map, contact_file.top_map)
        # ======================================================
        # Test Case 3
        contact_file = ContactFile('test')
        contact_map1 = ContactMap('foo')
        contact_map2 = ContactMap('bar')
        contact_file.add(contact_map1)
        contact_file.add(contact_map2)
        self.assertEqual(contact_map1, contact_file.top_map)
    
    def test_sort(self):
        # ======================================================
        # Test Case 1
        contact_file = ContactFile('test')
        for map in [ContactMap('foo'), ContactMap('bar'), ContactMap('doe')]:
            contact_file.add(map)
        contact_file_sorted = contact_file.sort('id', reverse=False, inplace=False)
        self.assertEqual(['bar', 'doe', 'foo'], [m.id for m in contact_file_sorted])
        self.assertNotEqual(contact_file, contact_file_sorted)
        # ======================================================
        # Test Case 2
        contact_file = ContactFile('test')
        for map in [ContactMap('foo'), ContactMap('bar'), ContactMap('doe')]:
            contact_file.add(map)
        contact_file_sorted = contact_file.sort('id', reverse=True, inplace=False)
        self.assertEqual(['foo', 'doe', 'bar'], [m.id for m in contact_file_sorted])
        self.assertNotEqual(contact_file, contact_file_sorted)
        # ======================================================
        # Test Case 3
        contact_file = ContactFile('test')
        for map in [ContactMap('foo'), ContactMap('bar'), ContactMap('doe')]:
            contact_file.add(map)
        contact_file_sorted = contact_file.sort('id', reverse=False, inplace=True)
        self.assertEqual(['bar', 'doe', 'foo'], [m.id for m in contact_file_sorted])
        self.assertEqual(contact_file, contact_file_sorted)
        # ======================================================
        # Test Case 4
        contact_file = ContactFile('test')
        for map in [ContactMap('foo'), ContactMap('bar'), ContactMap('doe')]:
            contact_file.add(map)
        contact_file_sorted = contact_file.sort('id', reverse=True, inplace=True)
        self.assertEqual(['foo', 'doe', 'bar'], [m.id for m in contact_file_sorted])
        self.assertEqual(contact_file, contact_file_sorted)


if __name__ == "__main__":
    unittest.main()
