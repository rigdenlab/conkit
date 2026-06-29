"""Testing facility for conkit.misc.renumbering_tools"""
import os
import tempfile
import unittest

from Bio.PDB import PDBParser as BioPDBParser
from Bio.PDB.Chain import Chain
from Bio.PDB.Residue import Residue
from Bio.PDB.Selection import unfold_entities

from conkit.core.sequence import Sequence
from conkit.misc.renumbering_tools import (
    construct_seq_from_chain,
    get_alignment_map_dict,
    write_renumbered_version_of_chain_in_struct,
)


# Minimal 3-residue RNA PDB.  Residues are numbered 5, 6, 7 (not 1, 2, 3) so
# that renumbering is genuinely exercised — the function should map them back to
# 1, 2, 3 when aligned against the 'GAU' reference sequence.
_RNA_PDB = """\
ATOM      1  C1'   G A   5       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C1'   A A   6       5.000   0.000   0.000  1.00  0.00           C
ATOM      3  C1'   U A   7      10.000   0.000   0.000  1.00  0.00           C
"""


class TestWriteRenumberedVersionOfChainInStruct(unittest.TestCase):

    def _run(self, tmpdir):
        """Write the RNA PDB fixture into tmpdir and call the function."""
        pdb_path = os.path.join(tmpdir, 'test_rna.pdb')
        with open(pdb_path, 'w') as fh:
            fh.write(_RNA_PDB)
        seq = Sequence('test', 'GAU')
        return write_renumbered_version_of_chain_in_struct(
            pdb_path, 'pdb', seq, selected_chain='A', moltype='RNA'
        )

    def test_returns_three_tuple(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            result = self._run(tmpdir)
        self.assertIsInstance(result, tuple)
        self.assertEqual(3, len(result))

    def test_output_file_is_created(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_name, _, _ = self._run(tmpdir)
            self.assertTrue(os.path.exists(out_name))

    def test_output_filename_contains_chain_and_stem(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_name, _, _ = self._run(tmpdir)
        self.assertIn('renumbered_A_test_rna.pdb', out_name)

    def test_alignment_dict_is_identity(self):
        # 'GAU' aligns perfectly to 'GAU', so every position maps to itself (0-indexed)
        with tempfile.TemporaryDirectory() as tmpdir:
            _, alignment_dict, _ = self._run(tmpdir)
        self.assertEqual({0: 0, 1: 1, 2: 2}, alignment_dict)

    def test_residues_renumbered_from_offset(self):
        # PDB residues 5, 6, 7 should become 1, 2, 3 after alignment to 'GAU'
        with tempfile.TemporaryDirectory() as tmpdir:
            out_name, _, _ = self._run(tmpdir)
            struct = BioPDBParser(QUIET=True).get_structure('out', out_name)
            resnums = [r.get_id()[1] for r in unfold_entities(struct[0]['A'], 'R')]
        self.assertEqual([1, 2, 3], resnums)

    def test_unknown_filetype_returns_none(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            pdb_path = os.path.join(tmpdir, 'test_rna.pdb')
            with open(pdb_path, 'w') as fh:
                fh.write(_RNA_PDB)
            seq = Sequence('test', 'GAU')
            result = write_renumbered_version_of_chain_in_struct(
                pdb_path, 'xyz', seq, selected_chain='A', moltype='RNA'
            )
        self.assertIsNone(result)


class TestConstructSeqFromChain(unittest.TestCase):

    def _make_chain(self, residue_names, start_num=1):
        """Build a minimal biopython Chain with no atoms — construct_seq_from_chain
        only reads residue names and sequence IDs, so atoms are not needed."""
        chain = Chain('A')
        for i, resname in enumerate(residue_names):
            res = Residue((' ', start_num + i, ' '), resname, '')
            chain.add(res)
        return chain

    def test_standard_rna_residues_return_correct_sequence(self):
        chain = self._make_chain(['G', 'A', 'U'])
        seq, _, _ = construct_seq_from_chain(chain, alphabet='RNA')
        self.assertEqual('GAU', seq)

    def test_borders_reflect_residue_numbering(self):
        # Residues numbered 5, 6, 7 — borders should be 5 and 7
        chain = self._make_chain(['G', 'A', 'U'], start_num=5)
        _, first, last = construct_seq_from_chain(chain, alphabet='RNA')
        self.assertEqual(5, first)
        self.assertEqual(7, last)

    def test_return_borders_false_gives_string_only(self):
        chain = self._make_chain(['G', 'A', 'U'])
        result = construct_seq_from_chain(chain, return_borders=False, alphabet='RNA')
        self.assertIsInstance(result, str)
        self.assertEqual('GAU', result)

    def test_modified_nucleotide_translates_via_mod_rescodes(self):
        # PSU (pseudouridine) is in mod_nuclist and maps to 'U'
        chain = self._make_chain(['G', 'PSU', 'A'])
        seq, _, _ = construct_seq_from_chain(chain, alphabet='RNA')
        self.assertEqual('GUA', seq)

    def test_unknown_residue_name_uses_placeholder(self):
        # 'XYZ' is not in any rescodes dict, so it becomes '?'
        chain = self._make_chain(['G', 'XYZ', 'A'])
        seq, _, _ = construct_seq_from_chain(chain, alphabet='RNA')
        self.assertEqual('G?A', seq)


class TestGetAlignmentMapDict(unittest.TestCase):

    def test_returns_two_dicts_by_default(self):
        result = get_alignment_map_dict('GAU', 'GAU')
        self.assertIsInstance(result, tuple)
        self.assertEqual(2, len(result))

    def test_identical_sequences_produce_identity_map(self):
        fwd, _ = get_alignment_map_dict('GAU', 'GAU')
        self.assertEqual({0: 0, 1: 1, 2: 2}, fwd)

    def test_reverse_dict_is_inverse_of_forward(self):
        fwd, rev = get_alignment_map_dict('GAU', 'GAU')
        for moving_idx, static_idx in fwd.items():
            self.assertEqual(moving_idx, rev[static_idx])

    def test_return_score_gives_three_values(self):
        result = get_alignment_map_dict('GAU', 'GAU', return_score=True)
        self.assertEqual(3, len(result))
        _, _, score = result
        self.assertIsInstance(score, float)

    def test_return_both_directions_false_gives_single_dict(self):
        result = get_alignment_map_dict('GAU', 'GAU', return_both_directions=False)
        self.assertIsInstance(result, dict)
        self.assertEqual({0: 0, 1: 1, 2: 2}, result)


if __name__ == '__main__':
    unittest.main()
