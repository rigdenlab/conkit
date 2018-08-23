"""Testing facility for conkit.misc.__init__"""

__author__ = "Felix Simkovic"
__date__ = "10 Jan 2018"

import unittest

from conkit.misc import *


class Test1(unittest.TestCase):
    def test_normalize_1(self):
        self.assertListEqual([0.0, 0.5, 1.0], normalize([1, 2, 3]))

    def test_normalize_2(self):
        self.assertListEqual([0.0, 0.5, 1.0], normalize([0.0, 0.5, 1.0]))

    def test_normalize_3(self):
        self.assertListEqual([0.0, 0.5, 1.0], normalize([-3, -2, -1]))

    def test_normalize_4(self):
        self.assertListEqual([0.0, 1.0], normalize([1, 2]))

    def test_normalize_5(self):
        self.assertListEqual([-1.0, 1.0], normalize([1, 2], vmin=-1))

    def test_normalize_6(self):
        self.assertListEqual([0.0, 2.0], normalize([1, 2], vmax=2))

    def test_normalize_7(self):
        self.assertListEqual([0.0, -1.0], normalize([1, 2], vmax=-1))

    def test_normalize_8(self):
        self.assertListEqual([0.2, 0.8], normalize([1, 2], vmin=0.2, vmax=0.8))

    def test_normalize_9(self):
        self.assertListEqual([0.2, 0.5, 0.8], normalize([1, 2, 3], vmin=0.2, vmax=0.8))

    def test_normalize_10(self):
        self.assertListEqual([1.0, 1.0, 1.0], normalize([2, 2, 2]))

    def test_normalize_11(self):
        self.assertListEqual([0.8, 0.8, 0.8], normalize([2, 2, 2], vmin=0.2, vmax=0.8))


class Test2(unittest.TestCase):
    def test_deprecated_1(self):
        @deprecate('0.0.0')
        def f():
            return True

        self.assertTrue(f())

    def test_deprecated_2(self):
        @deprecate('0.0.0', msg="hello world")
        def f():
            return True

        self.assertTrue(f())

    def test_deprecated_3(self):
        @deprecate('0.0.0')
        def f(a, b):
            return a + b

        self.assertEqual(2, f(1, 1))

    def test_deprecated_4(self):
        @deprecate('0.0.0')
        class Obj(object):
            pass

        self.assertTrue(Obj())

    def test_deprecated_5(self):
        class Obj(object):
            @deprecate('0.0.0')
            def f(self, a, b):
                return a + b

        self.assertEqual(2, Obj().f(1, 1))

    def test_deprecated_6(self):
        class Obj(object):
            @staticmethod
            @deprecate('0.0.0')
            def f(a, b):
                return a + b

        self.assertEqual(2, Obj.f(1, 1))

    def test_deprecated_7(self):
        class Obj(object):
            @classmethod
            @deprecate('0.0.0')
            def f(cls, a, b):
                return a + b

        self.assertEqual(2, Obj().f(1, 1))

    def test_deprecated_8(self):
        class Obj(object):
            @property
            @deprecate('0.0.0')
            def x(self):
                return 1

        self.assertEqual(1, Obj().x)

    def test_deprecated_9(self):
        class Obj(object):
            @property
            def x(self):
                return self._x

            @x.setter
            @deprecate('0.0.0')
            def x(self, x):
                self._x = x

        o = Obj()
        o.x = 2
        self.assertEqual(2, o.x)

    def test_deprecated_10(self):
        class Obj(object):
            @deprecate('0.0.0')
            @staticmethod
            def f(a, b):
                return a + b

        with self.assertRaises(Exception):
            Obj.f(1, 1)

    def test_deprecated_11(self):
        class Obj(object):
            @deprecate('0.0.0')
            @classmethod
            def f(cls, a, b):
                return a + b

        with self.assertRaises(AttributeError):
            Obj().f(1, 1)


if __name__ == "__main__":
    unittest.main(verbosity=2)
