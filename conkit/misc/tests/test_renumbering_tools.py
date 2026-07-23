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
    detect_numbering_anomalies,
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

# FASTA = 'GAU'.  Structure: G(1), A(2), A(2,'A'), U(3).
# The insertion-code residue A(2,'A') is a genuine extra structural residue with
# no FASTA slot.  It should be kept in the output with base seq_num (2) and
# original icode ('A').
_PDB_CORRECTLY_USED_ICODE = """\
ATOM      1  C1'   G A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C1'   A A   2       5.000   0.000   0.000  1.00  0.00           C
ATOM      3  C1'   A A   2A      5.000   5.000   0.000  1.00  0.00           C
ATOM      4  C1'   U A   3      10.000   0.000   0.000  1.00  0.00           C
"""

# FASTA = 'GAU'.  Structure canonical residues: G(1), U(2); plus A on insertion
# code (1,'A').  Alignment maps G→G and U→U, leaving FASTA slot 1 (A) free.
# The insertion-code A is therefore misused — it should receive a clean number.
# Expected output: G(1,' '), A(2,' '), U(3,' ') with no insertion codes.
_PDB_MISUSED_ICODE = """\
ATOM      1  C1'   G A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C1'   A A   1A      5.000   0.000   0.000  1.00  0.00           C
ATOM      3  C1'   U A   2      10.000   0.000   0.000  1.00  0.00           C
"""

# FASTA = 'GAU'.  Structure: G(1), A(2), U(3), C(4) — one extra canonical
# residue (C) with no FASTA counterpart.
# Expected: C is kept as insertion code off the last aligned FASTA position,
# i.e. (3,'A').  The three FASTA-matched residues have clean seq_ids 1, 2, 3.
_PDB_EXTRA_CANONICAL = """\
ATOM      1  C1'   G A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C1'   A A   2       5.000   0.000   0.000  1.00  0.00           C
ATOM      3  C1'   U A   3      10.000   0.000   0.000  1.00  0.00           C
ATOM      4  C1'   C A   4      15.000   0.000   0.000  1.00  0.00           C
"""


def _write_and_renumber(tmpdir, pdb_content, fasta_seq):
    """Write PDB to tmpdir and call write_renumbered_version_of_chain_in_struct."""
    pdb_path = os.path.join(tmpdir, 'test.pdb')
    with open(pdb_path, 'w') as fh:
        fh.write(pdb_content)
    seq = Sequence('test', fasta_seq)
    return write_renumbered_version_of_chain_in_struct(
        pdb_path, 'pdb', seq, selected_chain='A', moltype='RNA'
    )


def _residue_ids(struct_path):
    """Return list of (seq_id, icode) for all residues in the first model, chain A."""
    struct = BioPDBParser(QUIET=True).get_structure('out', struct_path)
    return [(r.get_id()[1], r.get_id()[2]) for r in unfold_entities(struct[0]['A'], 'R')]


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

    def test_returns_four_tuple(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            result = self._run(tmpdir)
        self.assertIsInstance(result, tuple)
        self.assertEqual(4, len(result))

    def test_output_file_is_created(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_name, _, _, _ = self._run(tmpdir)
            self.assertTrue(os.path.exists(out_name))

    def test_output_filename_contains_chain_and_stem(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_name, _, _, _ = self._run(tmpdir)
        self.assertIn('renumbered_A_test_rna.pdb', out_name)

    def test_alignment_dict_is_identity(self):
        # 'GAU' aligns perfectly to 'GAU', so every position maps to itself (0-indexed)
        with tempfile.TemporaryDirectory() as tmpdir:
            _, alignment_dict, _, _ = self._run(tmpdir)
        self.assertEqual({0: 0, 1: 1, 2: 2}, alignment_dict)

    def test_residues_renumbered_from_offset(self):
        # PDB residues 5, 6, 7 should become 1, 2, 3 after alignment to 'GAU'
        with tempfile.TemporaryDirectory() as tmpdir:
            out_name, _, _, _ = self._run(tmpdir)
            struct = BioPDBParser(QUIET=True).get_structure('out', out_name)
            resnums = [r.get_id()[1] for r in unfold_entities(struct[0]['A'], 'R')]
        self.assertEqual([1, 2, 3], resnums)

    def test_original_map_records_provenance(self):
        # PDB residues 5, 6, 7 renumbered to 1, 2, 3: original_map should map
        # each new (seq_id, ' ') back to the original (seq_id, ' ', resname).
        with tempfile.TemporaryDirectory() as tmpdir:
            _, _, _, original_map = self._run(tmpdir)
        self.assertIn((1, ' '), original_map)
        self.assertIn((2, ' '), original_map)
        self.assertIn((3, ' '), original_map)
        orig_seq_id, orig_icode, resname = original_map[(1, ' ')]
        self.assertEqual(5, orig_seq_id)
        self.assertEqual(' ', orig_icode)
        self.assertEqual('G', resname)

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

    # --- insertion-code and extra-canonical cases ---

    def test_correctly_used_icode_kept_in_output(self):
        # Structure: G(1), A(2), A(2,'A'), U(3). FASTA='GAU'.
        # A(2,'A') is a genuine extra structural residue with no FASTA slot.
        # All four residues (including the insertion-code one) should appear in output.
        with tempfile.TemporaryDirectory() as tmpdir:
            out_name, _, _, _ = _write_and_renumber(tmpdir, _PDB_CORRECTLY_USED_ICODE, 'GAU')
            ids = _residue_ids(out_name)
        self.assertEqual(4, len(ids))
        self.assertIn((1, ' '), ids)
        self.assertIn((2, ' '), ids)
        self.assertIn((3, ' '), ids)
        # Insertion-code residue kept with base seq_num and original icode.
        self.assertIn((2, 'A'), ids)

    def test_correctly_used_icode_canonical_residues_have_blank_icodes(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_name, _, _, _ = _write_and_renumber(tmpdir, _PDB_CORRECTLY_USED_ICODE, 'GAU')
            ids = _residue_ids(out_name)
        blank_icode_count = sum(1 for _, icode in ids if not icode.strip())
        self.assertEqual(3, blank_icode_count, "Three canonical residues should have blank icodes")

    def test_misused_icode_gets_clean_sequential_number(self):
        # Structure canonical: G(1), U(2); insertion code: A(1,'A'). FASTA='GAU'.
        # G aligns to G (FASTA pos 0), U aligns to U (FASTA pos 2), leaving FASTA
        # slot 1 (A) free.  The misused insertion-code A should fill that slot.
        with tempfile.TemporaryDirectory() as tmpdir:
            out_name, _, _, _ = _write_and_renumber(tmpdir, _PDB_MISUSED_ICODE, 'GAU')
            ids = _residue_ids(out_name)
        self.assertEqual(3, len(ids))
        self.assertTrue(
            all(icode.strip() == '' for _, icode in ids),
            "All residues should have clean (blank) icodes after misused-icode correction",
        )
        self.assertEqual([1, 2, 3], sorted(seq_id for seq_id, _ in ids))

    def test_extra_canonical_residue_kept_with_insertion_code(self):
        # Structure: G(1), A(2), U(3), C(4). FASTA='GAU'.
        # C(4) has no FASTA counterpart and must be retained in the output with an
        # insertion code rather than silently excluded.
        with tempfile.TemporaryDirectory() as tmpdir:
            out_name, _, _, _ = _write_and_renumber(tmpdir, _PDB_EXTRA_CANONICAL, 'GAU')
            ids = _residue_ids(out_name)
        self.assertEqual(4, len(ids), "Extra canonical residue should appear in output")
        icode_count = sum(1 for _, icode in ids if icode.strip())
        self.assertEqual(1, icode_count, "Exactly one residue should have an insertion code")
        clean_ids = sorted(seq_id for seq_id, icode in ids if not icode.strip())
        self.assertEqual([1, 2, 3], clean_ids, "FASTA-matched residues should have seq_ids 1–3")

    def test_extra_canonical_residue_emits_warning(self):
        # The function should log a WARNING so the user knows about the anomaly.
        with tempfile.TemporaryDirectory() as tmpdir:
            with self.assertLogs('conkit.misc.renumbering_tools', level='WARNING') as cm:
                _write_and_renumber(tmpdir, _PDB_EXTRA_CANONICAL, 'GAU')
        self.assertTrue(
            any('no FASTA counterpart' in msg for msg in cm.output),
            "Expected a warning about a residue with no FASTA counterpart",
        )

    def test_original_map_misused_icode_has_correct_provenance(self):
        # The misused-icode A(1,'A') gets new id (2,' '). original_map should record
        # the original insertion-code origin.
        with tempfile.TemporaryDirectory() as tmpdir:
            _, _, _, original_map = _write_and_renumber(tmpdir, _PDB_MISUSED_ICODE, 'GAU')
        # Find the entry whose orig_icode is non-blank (the misused icode residue).
        misused = [(new_id, v) for new_id, v in original_map.items() if v[1].strip()]
        self.assertEqual(1, len(misused), "Exactly one misused-icode entry expected")
        new_id, (orig_seq_id, orig_icode, resname) = misused[0]
        self.assertEqual(' ', new_id[1], "New icode should be blank after correction")
        self.assertNotEqual(' ', orig_icode, "Original icode should be non-blank")

    def test_original_map_extra_canonical_has_correct_provenance(self):
        # The extra canonical C(4,' ') gets new id (3,'A'). original_map should record
        # the original clean-integer origin.
        with tempfile.TemporaryDirectory() as tmpdir:
            _, _, _, original_map = _write_and_renumber(tmpdir, _PDB_EXTRA_CANONICAL, 'GAU')
        extra = [(new_id, v) for new_id, v in original_map.items() if new_id[1].strip()]
        self.assertEqual(1, len(extra), "Exactly one extra-canonical entry expected")
        new_id, (orig_seq_id, orig_icode, resname) = extra[0]
        self.assertNotEqual(' ', new_id[1], "New icode should be non-blank")
        self.assertEqual(' ', orig_icode, "Original icode should be blank")


class TestDetectNumberingAnomalies(unittest.TestCase):

    def test_empty_map_returns_empty_list(self):
        self.assertEqual([], detect_numbering_anomalies({}))

    def test_no_anomaly_for_plain_renumbering(self):
        # A residue that simply moved from seq_id 5 to 1 (both blank icode) is not anomalous.
        original_map = {(1, ' '): (5, ' ', 'G')}
        self.assertEqual([], detect_numbering_anomalies(original_map))

    def test_misused_icode_detected(self):
        # orig had icode 'A', new is clean → MISUSED_ICODE
        original_map = {(2, ' '): (1, 'A', 'A')}
        anomalies = detect_numbering_anomalies(original_map)
        self.assertEqual(1, len(anomalies))
        self.assertEqual('MISUSED_ICODE', anomalies[0][5])
        self.assertEqual(2, anomalies[0][0])   # new_seq_id
        self.assertEqual(' ', anomalies[0][1]) # new_icode

    def test_extra_canonical_detected(self):
        # orig had blank icode, new has icode 'A' → EXTRA_CANONICAL
        original_map = {(3, 'A'): (4, ' ', 'C')}
        anomalies = detect_numbering_anomalies(original_map)
        self.assertEqual(1, len(anomalies))
        self.assertEqual('EXTRA_CANONICAL', anomalies[0][5])
        self.assertEqual(3, anomalies[0][0])   # new_seq_id
        self.assertEqual('A', anomalies[0][1]) # new_icode

    def test_anomalies_sorted_by_new_seq_id_then_icode(self):
        original_map = {
            (5, 'A'): (6, ' ', 'G'),   # EXTRA_CANONICAL
            (2, ' '): (1, 'B', 'A'),   # MISUSED_ICODE
            (3, 'A'): (4, ' ', 'C'),   # EXTRA_CANONICAL
        }
        anomalies = detect_numbering_anomalies(original_map)
        keys = [(a[0], a[1]) for a in anomalies]
        self.assertEqual(sorted(keys), keys)


class TestConstructSeqFromChain(unittest.TestCase):

    def _make_chain(self, residue_names, start_num=1):
        """Build a minimal biopython Chain with no atoms — construct_seq_from_chain
        only reads residue names and sequence IDs, so atoms are not needed."""
        chain = Chain('A')
        for i, resname in enumerate(residue_names):
            res = Residue((' ', start_num + i, ' '), resname, '')
            chain.add(res)
        return chain

    def _make_chain_with_icodes(self, residues):
        """Build a chain from a list of (resname, seq_id, icode) tuples."""
        chain = Chain('A')
        for resname, seq_id, icode in residues:
            res = Residue((' ', seq_id, icode), resname, '')
            chain.add(res)
        return chain

    def test_standard_rna_residues_return_correct_sequence(self):
        chain = self._make_chain(['G', 'A', 'U'])
        seq, *_ = construct_seq_from_chain(chain, alphabet='RNA')
        self.assertEqual('GAU', seq)

    def test_borders_reflect_residue_numbering(self):
        # Residues numbered 5, 6, 7 — borders should be 5 and 7
        chain = self._make_chain(['G', 'A', 'U'], start_num=5)
        _, first, last, _, _ = construct_seq_from_chain(chain, alphabet='RNA')
        self.assertEqual(5, first)
        self.assertEqual(7, last)

    def test_return_borders_false_gives_seq_as_first_element(self):
        # return_borders=False returns (seq, canonical_by_chain_pos, insertion_residues).
        chain = self._make_chain(['G', 'A', 'U'])
        result = construct_seq_from_chain(chain, return_borders=False, alphabet='RNA')
        self.assertIsInstance(result, tuple)
        self.assertIsInstance(result[0], str)
        self.assertEqual('GAU', result[0])

    def test_modified_nucleotide_translates_via_mod_rescodes(self):
        # PSU (pseudouridine) is in mod_nuclist and maps to 'U'
        chain = self._make_chain(['G', 'PSU', 'A'])
        seq, *_ = construct_seq_from_chain(chain, alphabet='RNA')
        self.assertEqual('GUA', seq)

    def test_unknown_residue_name_uses_placeholder(self):
        # 'XYZ' is not in any rescodes dict, so it becomes '?'
        chain = self._make_chain(['G', 'XYZ', 'A'])
        seq, *_ = construct_seq_from_chain(chain, alphabet='RNA')
        self.assertEqual('G?A', seq)

    def test_insertion_code_residues_excluded_from_canonical_seq(self):
        # G(1,' '), A(2,' '), A(2,'A'), U(3,' ') — only the three canonical
        # residues should contribute to the sequence string.
        chain = self._make_chain_with_icodes([
            ('G', 1, ' '), ('A', 2, ' '), ('A', 2, 'A'), ('U', 3, ' ')
        ])
        seq, first, last, canonical_by_chain_pos, insertion_residues = \
            construct_seq_from_chain(chain, alphabet='RNA')
        self.assertEqual('GAU', seq)
        self.assertEqual(3, len(canonical_by_chain_pos))
        self.assertEqual(1, len(insertion_residues))

    def test_insertion_residues_list_contains_icode_residue(self):
        # Verify the insertion_residues list holds the correct residue object.
        chain = self._make_chain_with_icodes([
            ('G', 1, ' '), ('A', 1, 'A'), ('U', 2, ' ')
        ])
        _, _, _, _, insertion_residues = construct_seq_from_chain(chain, alphabet='RNA')
        self.assertEqual(1, len(insertion_residues))
        res = insertion_residues[0]
        self.assertEqual('A', res.resname)
        self.assertEqual(1, res.get_id()[1])
        self.assertEqual('A', res.get_id()[2])

    def test_canonical_by_chain_pos_does_not_include_icode_residues(self):
        # The canonical_by_chain_pos dict must only map canonical (blank-icode) residues.
        chain = self._make_chain_with_icodes([
            ('G', 1, ' '), ('A', 2, ' '), ('A', 2, 'A'), ('U', 3, ' ')
        ])
        _, _, _, canonical_by_chain_pos, _ = construct_seq_from_chain(chain, alphabet='RNA')
        for res in canonical_by_chain_pos.values():
            self.assertEqual(' ', res.get_id()[2], "Only blank-icode residues in canonical_by_chain_pos")


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
