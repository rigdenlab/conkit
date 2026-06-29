"""Testing facility for conkit.command_line.conkit_validate"""
import os
import tempfile
import unittest
from unittest.mock import patch

import pandas as pd

from conkit.command_line.conkit_validate import calculate_dnatco
from conkit.misc import DNATCO_CATEGORIES


# Minimal DNATCO report fixture.
#
# calculate_dnatco scans the lines for two markers and derives start/end indices:
#   "All dinucleotides" at line 1  →  start = 1 + 4 = 5
#   "Dinucleotide outliers" at line 14  →  end  = 14 - 6 = 8
#   range(5, 8) processes exactly lines 5, 6, 7 (the three data lines)
#
# The 5 blank lines between the last data line (7) and the /--- box (13) are
# required to make end = 8 rather than something smaller.
_DNATCO_FIXTURE = (
    "/------------------------------------------------------------------------------\\\n"  # 0
    "|                               All dinucleotides                              |\n"  # 1
    "|                                                                              |\n"  # 2
    "\\------------------------------------------------------------------------------/\n"  # 3
    "Chain      Step      NtC  CANA RMSD\n"                                               # 4
    "A (Cif: A) G1 G2     AA02 AAA  0.530\n"                                             # 5
    "A (Cif: A) G2 U3     AA04 AAA  0.437\n"                                             # 6
    "A (Cif: A) U3 C4     AB1S SYN  0.507\n"                                             # 7
    "\n"                                                                                   # 8
    "\n"                                                                                   # 9
    "\n"                                                                                   # 10
    "\n"                                                                                   # 11
    "\n"                                                                                   # 12
    "/------------------------------------------------------------------------------\\\n"  # 13
    "|                             Dinucleotide outliers                            |\n"  # 14
)


class TestCalculateDnatco(unittest.TestCase):

    def _run_with_fixture(self):
        """Write fixture to custom_report.txt in a temp dir and call the function.

        calculate_dnatco always reads from the hardcoded name 'custom_report.txt' in
        the current directory, so we change directory to a temporary location where we
        have placed our fixture, then restore the original directory afterwards.
        """
        with tempfile.TemporaryDirectory() as tmpdir:
            with open(os.path.join(tmpdir, 'custom_report.txt'), 'w') as fh:
                fh.write(_DNATCO_FIXTURE)

            original_dir = os.getcwd()
            try:
                os.chdir(tmpdir)
                with patch('conkit.command_line.conkit_validate.subprocess.run'):
                    result = calculate_dnatco('test.cif', 'mmcif', '/path/to/dnatco')
            finally:
                os.chdir(original_dir)

        return result

    def test_returns_dataframe(self):
        result = self._run_with_fixture()
        self.assertIsInstance(result, pd.DataFrame)

    def test_columns_are_dnatco_categories(self):
        result = self._run_with_fixture()
        self.assertEqual(DNATCO_CATEGORIES, result.columns.tolist())

    def test_index_contains_residue_numbers(self):
        # Residues 1-4 appear across the three dinucleotide steps
        result = self._run_with_fixture()
        self.assertEqual([1, 2, 3, 4], sorted(result.index.tolist()))

    def test_residue_1_aaa_score(self):
        # G1 appears only in step G1-G2 (RMSD=0.530); score = 1 / 0.530
        result = self._run_with_fixture()
        self.assertAlmostEqual(1.0 / 0.530, result.loc[1, 'AAA'], places=4)

    def test_residue_2_aaa_score_is_averaged(self):
        # G2 appears in G1-G2 (RMSD=0.530) AND G2-U3 (RMSD=0.437)
        # The function averages those RMSDs before taking 1/avg
        result = self._run_with_fixture()
        avg_rmsd = (0.530 + 0.437) / 2
        self.assertAlmostEqual(1.0 / avg_rmsd, result.loc[2, 'AAA'], places=4)

    def test_residue_3_syn_score(self):
        # U3 appears as res1 in step U3-C4 with CANA=SYN, RMSD=0.507
        result = self._run_with_fixture()
        self.assertAlmostEqual(1.0 / 0.507, result.loc[3, 'SYN'], places=4)

    def test_dnatco_tot_rmsd_accumulates(self):
        # DNATCO_TOT_RMSD for residue 2: sum of every RMSD it appears in = 0.530 + 0.437
        result = self._run_with_fixture()
        self.assertAlmostEqual(0.530 + 0.437, result.loc[2, 'DNATCO_TOT_RMSD'], places=4)

    def test_unknown_type_returns_none(self):
        # A type that is neither 'mmcif' nor 'pdb' causes an early return with None
        with patch('conkit.command_line.conkit_validate.subprocess.run'):
            result = calculate_dnatco('test.xyz', 'xyz', '/path/to/dnatco')
        self.assertIsNone(result)


if __name__ == '__main__':
    unittest.main()
