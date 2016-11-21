"""Testing facility for conkit.io._iotools"""

__author__ = "Felix Simkovic"
__date__ = "21 Nov 2016"

from conkit.io import _iotools

import os
import sys
import unittest


class Test(unittest.TestCase):

    def test_create_tmp_f(self):
        # ===============================================
        # Test Case 1
        fname = _iotools.create_tmp_f()
        self.assertTrue(os.path.isfile(fname))
        os.unlink(fname)
        # ===============================================
        # Test Case 2
        content = 'Hello, World!'
        fname = _iotools.create_tmp_f(content=content, mode='w')
        self.assertTrue(os.path.isfile(fname))
        with open(fname, 'r') as f_in:
            written_content = f_in.read()
        self.assertEqual(content, written_content)
        os.unlink(fname)
        # ===============================================
        # Test Case 3 
        content = 'Hello, World!'
        content_bytes = content.encode('utf-8')
        fname = _iotools.create_tmp_f(content=content_bytes, mode='wb')
        self.assertTrue(os.path.isfile(fname))
        with open(fname, 'rb') as f_in:
            written_content_bytes = f_in.read()
        writte_content = written_content_bytes.decode('utf-8')
        self.assertEqual(content, written_content)
        os.unlink(fname)

    def test_is_str_like(self):
        self.assertTrue(_iotools.is_str_like('foo'))    # str
        self.assertFalse(_iotools.is_str_like(1))       # int
        self.assertFalse(_iotools.is_str_like(1.))      # float
        self.assertFalse(_iotools.is_str_like([]))      # list
        self.assertFalse(_iotools.is_str_like(()))      # tuple
        self.assertFalse(_iotools.is_str_like({}))      # dict
        self.assertFalse(_iotools.is_str_like(set()))   # set
    
    def test_open_f_handle(self):
        # ===============================================
        # Test Case 1
        fname = _iotools.create_tmp_f()
        with _iotools.open_f_handle(fname, 'append') as fhandle:
            self.assertEqual('a' + _iotools.MODE_APPENDIX, fhandle.mode)
        f_in_handle =  _iotools.open_f_handle(fname, 'append')
        with _iotools.open_f_handle(f_in_handle, 'append') as fhandle:
            self.assertEqual('a' + _iotools.MODE_APPENDIX, fhandle.mode)
        os.unlink(fname)
        # ===============================================
        # Test Case 2
        fname = _iotools.create_tmp_f()
        with _iotools.open_f_handle(fname, 'read') as fhandle:
            self.assertEqual('r' + _iotools.MODE_APPENDIX, fhandle.mode)
        f_in_handle =  _iotools.open_f_handle(fname, 'read')
        with _iotools.open_f_handle(f_in_handle, 'read') as fhandle:
            self.assertEqual('r' + _iotools.MODE_APPENDIX, fhandle.mode)
        os.unlink(fname)
        # ===============================================
        # Test Case 3
        fname = _iotools.create_tmp_f()
        with _iotools.open_f_handle(fname, 'write') as fhandle:
            self.assertEqual('w' + _iotools.MODE_APPENDIX, fhandle.mode)
        f_in_handle =  _iotools.open_f_handle(fname, 'write')
        with _iotools.open_f_handle(f_in_handle, 'write') as fhandle:
            self.assertEqual('w' + _iotools.MODE_APPENDIX, fhandle.mode)
        os.unlink(fname)
        # ===============================================
        # Test Case 4
        with self.assertRaises(TypeError):
            _iotools.open_f_handle(1, 'read')
        with self.assertRaises(TypeError):
            _iotools.open_f_handle(1.0, 'write')
        # ===============================================
        # Test Case 5 
        fname = _iotools.create_tmp_f()
        with self.assertRaises(ValueError):
            _iotools.open_f_handle(fname, 'foo')
        with self.assertRaises(ValueError):
            _iotools.open_f_handle(fname, 'bar')


if __name__ == "__main__":
    unittest.main(verbosity=2)
