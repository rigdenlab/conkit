"""Bunch of helper stuff for conkit.io tests."""

__author__ = "Felix Simkovic"
__date__ = "19 Sep 2019"

import os
import unittest

from conkit.io._iotools import create_tmp_f


class ParserTestCase(unittest.TestCase):

    def tearDown(self):
        self.doCleanups()

    def tempfile(self, *args, **kwargs):
        fname = create_tmp_f(*args, **kwargs)
        self.addCleanup(os.remove, fname)
        return fname
