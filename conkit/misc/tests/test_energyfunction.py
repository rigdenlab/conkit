__author__ = "Felix Simkovic"

import unittest

from conkit.misc.energyfunction import RosettaFunctionConstructs

TEMPLATE = dict(
    atom1="CB",
    res1_seq=1,
    atom2="CB",
    res2_seq=2,
    lower_bound=0,
    upper_bound=2,
    scalar_score=0.1,
    sigmoid_cutoff=0.2,
    sigmoid_slope=0.3,
    energy_bonus=-15.0,
    raw_score=1.0,
)


class TestRosettaFunctionConstructs(unittest.TestCase):
    def test_1(self):
        output = RosettaFunctionConstructs().BOUNDED_default.format(**TEMPLATE)
        self.assertEqual(output, "AtomPair CB    1 CB    2 BOUNDED 0.000 2.000 1 0.5 #")

    def test_2(self):
        output = RosettaFunctionConstructs().BOUNDED_gremlin.format(**TEMPLATE)
        self.assertEqual(output, "AtomPair CB    1 CB    2 SCALARWEIGHTEDFUNC  0.100 BOUNDED 0 0.000 1 0.5")

    def test_3(self):
        output = RosettaFunctionConstructs().FADE.format(**TEMPLATE)
        self.assertEqual(output, "AtomPair CB    1 CB    2 FADE -10 19 10 -15.00 0")

    def test_4(self):
        output = RosettaFunctionConstructs().FADE_default.format(**TEMPLATE)
        self.assertEqual(output, "AtomPair CB    1 CB    2 FADE -10 19 10 -15.00 0")

    def test_5(self):
        output = RosettaFunctionConstructs().SIGMOID_default.format(**TEMPLATE)
        self.assertEqual(output, "AtomPair CB    1 CB    2 SIGMOID 8.00 1.00 #ContactMap: 1.000")

    def test_6(self):
        output = RosettaFunctionConstructs().SIGMOID_gremlin.format(**TEMPLATE)
        self.assertEqual(
            output,
            "AtomPair CB    1 CB    2 SCALARWEIGHTEDFUNC  0.100 SUMFUNC 2 SIGMOID  0.200  0.300 CONSTANTFUNC -0.5",
        )


if __name__ == "__main__":
    unittest.main(verbosity=2)
