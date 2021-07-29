"""Testing facility for conkit.core.DistanceFile"""


import unittest

from conkit.core.distogram import Distogram
from conkit.core.distancefile import DistanceFile


class TestDistanceFile(unittest.TestCase):

    def test_original_file_format(self):
        distance_file = DistanceFile("test")
        distance_file.original_file_format = "PDB"
        distogram = Distogram("test")
        distance_file.add(distogram)
        self.assertTrue(distogram in distance_file.child_list)
        self.assertEqual("PDB", distogram.original_file_format)
