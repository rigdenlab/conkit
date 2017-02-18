"""Testing facility for conkit.plot._Figure"""

__author__ = "Felix Simkovic"
__date__ = "08 Feb 2017"

from conkit.plot._Figure import Figure

import unittest


class Test(unittest.TestCase):

    def test__init__1(self):
        f = Figure()
        self.assertEqual(300, f.dpi)
        self.assertEqual("png", f.format)
        self.assertEqual("conkit", f.prefix)

    def test__init__2(self):
        f = Figure(dpi=20, format="svg")
        self.assertEqual(20, f.dpi)
        self.assertEqual("svg", f.format)
        self.assertEqual("conkit", f.prefix)

        self.assertRaises(ValueError, Figure, format="ooo")

    def test_dpi_1(self):
        f = Figure()
        self.assertEqual(300, f.dpi)
        f.dpi = 20
        self.assertEqual(20, f.dpi)

    def test_filename_1(self):
        f = Figure()
        self.assertEqual("conkit.png", f.file_name)
        f.file_name = "test.eps"
        self.assertEqual("test.eps", f.file_name)
        self.assertEqual("eps", f.format)
        self.assertEqual("test", f.prefix)
        f.format = "pdf"
        self.assertEqual("test.pdf", f.file_name)
        self.assertEqual("pdf", f.format)
        self.assertEqual("test", f.prefix)
        self.assertRaises(ValueError, setattr, f, "file_name", "test.ooo")

    def test_format_1(self):
        f = Figure()
        self.assertEqual("png", f.format)
        f.format = "pdf"
        self.assertEqual("pdf", f.format)
        f.format = "ps"
        self.assertEqual("ps", f.format)
        f.format = "png"
        self.assertEqual("png", f.format)
        self.assertRaises(ValueError, setattr, f, "format", "bbb")

    def test_prefix_1(self):
        f = Figure()
        self.assertEqual("conkit", f.prefix)
        f.prefix = "test"
        self.assertEqual("test", f.prefix)
        f.prefix = "/path/to/prefix"
        self.assertEqual("/path/to/prefix", f.prefix)


if __name__ == "__main__":
    unittest.main(verbosity=2)
