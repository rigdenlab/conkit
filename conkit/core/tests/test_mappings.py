"""Testing facility for conkit.core.mappings"""

__author__ = "Felix Simkovic"
__date__ = "16 Mar 2018"

import unittest

from conkit.core.mappings import AminoAcidMapping, AminoAcidThreeToOne, AminoAcidOneToThree


class AminoAcidMappingTest(unittest.TestCase):
    def test_A_1(self):
        self.assertEqual(1, AminoAcidMapping["A"].value)

    def test_C_1(self):
        self.assertEqual(2, AminoAcidMapping["C"].value)

    def test_D_1(self):
        self.assertEqual(3, AminoAcidMapping["D"].value)

    def test_E_1(self):
        self.assertEqual(4, AminoAcidMapping["E"].value)

    def test_F_1(self):
        self.assertEqual(5, AminoAcidMapping["F"].value)

    def test_G_1(self):
        self.assertEqual(6, AminoAcidMapping["G"].value)

    def test_H_1(self):
        self.assertEqual(7, AminoAcidMapping["H"].value)

    def test_I_1(self):
        self.assertEqual(8, AminoAcidMapping["I"].value)

    def test_K_1(self):
        self.assertEqual(9, AminoAcidMapping["K"].value)

    def test_L_1(self):
        self.assertEqual(10, AminoAcidMapping["L"].value)

    def test_M_1(self):
        self.assertEqual(11, AminoAcidMapping["M"].value)

    def test_N_1(self):
        self.assertEqual(12, AminoAcidMapping["N"].value)

    def test_P_1(self):
        self.assertEqual(13, AminoAcidMapping["P"].value)

    def test_Q_1(self):
        self.assertEqual(14, AminoAcidMapping["Q"].value)

    def test_R_1(self):
        self.assertEqual(15, AminoAcidMapping["R"].value)

    def test_S_1(self):
        self.assertEqual(16, AminoAcidMapping["S"].value)

    def test_T_1(self):
        self.assertEqual(17, AminoAcidMapping["T"].value)

    def test_V_1(self):
        self.assertEqual(18, AminoAcidMapping["V"].value)

    def test_W_1(self):
        self.assertEqual(19, AminoAcidMapping["W"].value)

    def test_X_1(self):
        self.assertEqual(21, AminoAcidMapping["X"].value)

    def test_Y_1(self):
        self.assertEqual(20, AminoAcidMapping["Y"].value)

    def test_Z_1(self):
        self.assertEqual(21, getattr(AminoAcidMapping, "Z", AminoAcidMapping.X.value))


class AminoAcidThreeToOneTest(unittest.TestCase):
    def test_single_1(self):
        self.assertTrue("ALA" in AminoAcidThreeToOne.__members__)

    def test_single_2(self):
        self.assertFalse("ala" in AminoAcidThreeToOne.__members__)

    def test_single_3(self):
        self.assertFalse("Ala" in AminoAcidThreeToOne.__members__)

    def test_single_4(self):
        self.assertEqual("A", AminoAcidThreeToOne.ALA.value)

    def test_single_5(self):
        self.assertEqual("A", AminoAcidThreeToOne["ALA"].value)


class AminoAcidOneToThreeTest(unittest.TestCase):
    def test_single_1(self):
        self.assertTrue("A" in AminoAcidOneToThree.__members__)

    def test_single_2(self):
        self.assertFalse("a" in AminoAcidOneToThree.__members__)

    def test_single_3(self):
        self.assertEqual("ALA", AminoAcidOneToThree.A.value)

    def test_single_4(self):
        self.assertEqual("ALA", AminoAcidOneToThree["A"].value)


if __name__ == "__main__":
    unittest.main(verbosity=2)
