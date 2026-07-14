"""Testing facility for conkit.plot.modelvalidation.ModelValidationFigure"""
import os
import tempfile
import unittest
from unittest.mock import MagicMock, patch

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

import conkit.io
from conkit.core.sequence import Sequence
from conkit.plot.modelvalidation import ModelValidationFigure


# Minimal 5-residue RNA PDB with C1' atoms spaced 5 Å apart in a line.
# All non-neighboring pairs are within the 20 Å cutoff used for reading.
_RNA_PDB_5 = """\
ATOM      1  C1'   G A   1       0.000   0.000   0.000  1.00  0.00           C
ATOM      2  C1'   A A   2       5.000   0.000   0.000  1.00  0.00           C
ATOM      3  C1'   U A   3      10.000   0.000   0.000  1.00  0.00           C
ATOM      4  C1'   C A   4      15.000   0.000   0.000  1.00  0.00           C
ATOM      5  C1'   G A   5      20.000   0.000   0.000  1.00  0.00           C
END
"""
_RNA_SEQ = 'GAUCG'


def _read_rna_distogram(pdb_path):
    return conkit.io.read(pdb_path, 'pdb', distance_cutoff=20.0, atom_type="C1'").top


class _FixtureBase(unittest.TestCase):
    """Shared setUpClass: write the PDB once and read the distogram."""

    @classmethod
    def setUpClass(cls):
        cls._tmpdir = tempfile.TemporaryDirectory()
        pdb_path = os.path.join(cls._tmpdir.name, 'rna.pdb')
        with open(pdb_path, 'w') as fh:
            fh.write(_RNA_PDB_5)
        cls.distogram = _read_rna_distogram(pdb_path)
        cls.sequence = Sequence('test', _RNA_SEQ)

    @classmethod
    def tearDownClass(cls):
        cls._tmpdir.cleanup()

    def tearDown(self):
        plt.close('all')

    def _make_fig(self):
        return ModelValidationFigure(self.distogram, self.distogram, self.sequence)


# ---------------------------------------------------------------------------
# Construction
# ---------------------------------------------------------------------------

class TestModelValidationFigureConstruction(_FixtureBase):

    def test_short_sequence_raises_value_error(self):
        short_seq = Sequence('short', 'GAUC')  # 4 residues — below the 5-residue minimum
        with self.assertRaises(ValueError):
            ModelValidationFigure(self.distogram, self.distogram, short_seq)

    def test_invalid_model_type_raises_type_error(self):
        with self.assertRaises(TypeError):
            ModelValidationFigure('not_a_distogram', self.distogram, self.sequence)

    def test_init_creates_data_dataframe(self):
        fig = self._make_fig()
        self.assertIsInstance(fig.data, pd.DataFrame)

    def test_init_data_has_resnum_column(self):
        fig = self._make_fig()
        self.assertIn('RESNUM', fig.data.columns)

    def test_init_data_residues_match_sequence(self):
        fig = self._make_fig()
        self.assertEqual(sorted(fig.data['RESNUM'].tolist()), [1, 2, 3, 4, 5])

    def test_init_misaligned_column_defaults_to_false(self):
        fig = self._make_fig()
        self.assertTrue((fig.data['MISALIGNED'] == False).all())

    def test_init_alignment_dict_is_empty(self):
        fig = self._make_fig()
        self.assertEqual({}, fig.alignment)


# ---------------------------------------------------------------------------
# calculate_features
# ---------------------------------------------------------------------------

class TestModelValidationFigureCalculateFeatures(_FixtureBase):

    def setUp(self):
        self.fig = self._make_fig()
        self.fig.calculate_features()

    def test_adds_wrmsd_column(self):
        self.assertIn('WRMSD', self.fig.data.columns)

    def test_adds_accuracy_column(self):
        self.assertIn('ACCURACY', self.fig.data.columns)

    def test_adds_zscore_columns(self):
        for col in ('ZSCORE_WRMSD', 'ZSCORE_ACCURACY', 'ZSCORE_FN'):
            self.assertIn(col, self.fig.data.columns)

    def test_row_count_unchanged(self):
        self.assertEqual(5, len(self.fig.data))

    def test_resnum_column_still_present(self):
        self.assertIn('RESNUM', self.fig.data.columns)


# ---------------------------------------------------------------------------
# svm — RNA paths (model-name selection)
# ---------------------------------------------------------------------------

class TestModelValidationFigureSvmRNA(_FixtureBase):

    def setUp(self):
        # Fresh figure per test: svm() merges columns so calling it twice on
        # the same object causes duplicate-column errors.
        self.fig = self._make_fig()

    def _call_svm(self, prediction_type, sec_struc_info):
        mock_clf = MagicMock()
        mock_scaler = MagicMock()
        with patch('conkit.plot.modelvalidation.load_specific_validation_model',
                   return_value=(mock_clf, mock_scaler)), \
             patch.object(self.fig, '_predict_score', return_value=0.5):
            self.fig.svm(ext_info=None, moltype='RNA',
                         prediction_type=prediction_type,
                         sec_struc_info=sec_struc_info)

    def test_dist_no_dnatco_svm_name(self):
        self._call_svm('DIST', None)
        self.assertEqual('RNA_AF3_dist_', self.fig.svm_name)

    def test_struct_no_dnatco_svm_name(self):
        self._call_svm('STRUCT', None)
        self.assertEqual('RNA_AF3_struct_', self.fig.svm_name)

    def test_dist_with_dnatco_svm_name(self):
        self._call_svm('DIST', 'DNATCO')
        self.assertEqual('RNA_DNATCO_AF3_dist_', self.fig.svm_name)

    def test_struct_with_dnatco_svm_name(self):
        self._call_svm('STRUCT', 'DNATCO')
        self.assertEqual('RNA_DNATCO_AF3_struct_', self.fig.svm_name)

    def test_svm_adds_score_column(self):
        self._call_svm('DIST', None)
        self.assertIn('SCORE', self.fig.data.columns)

    def test_ext_info_none_adds_acc_column_of_zeros(self):
        self._call_svm('DIST', None)
        self.assertIn('ACC', self.fig.data.columns)
        self.assertTrue((self.fig.data['ACC'] == 0).all())


# ---------------------------------------------------------------------------
# svm — Protein path
# ---------------------------------------------------------------------------

class TestModelValidationFigureSvmProtein(_FixtureBase):

    def setUp(self):
        self.fig = self._make_fig()

    def _call_protein_svm(self):
        mock_clf = MagicMock()
        mock_scaler = MagicMock()
        with patch('conkit.plot.modelvalidation.load_validation_model',
                   return_value=(mock_clf, mock_scaler)), \
             patch.object(self.fig, '_predict_score', return_value=0.3):
            self.fig.svm(ext_info=None, moltype='Protein')

    def test_protein_svm_name(self):
        self._call_protein_svm()
        self.assertEqual('Protein_DSSP_AF2_dist_', self.fig.svm_name)

    def test_protein_ext_info_none_adds_coil_column(self):
        self._call_protein_svm()
        self.assertIn('COIL', self.fig.data.columns)

    def test_protein_ext_info_none_adds_helix_column(self):
        self._call_protein_svm()
        self.assertIn('HELIX', self.fig.data.columns)

    def test_protein_ext_info_none_dssp_columns_are_zero(self):
        self._call_protein_svm()
        for col in ('COIL', 'HELIX', 'SHEET', 'ACC'):
            self.assertTrue((self.fig.data[col] == 0).all(),
                            msg=f'Expected {col} to be all zeros')


# ---------------------------------------------------------------------------
# svm_error_calling
# ---------------------------------------------------------------------------

class TestModelValidationFigureSvmErrorCalling(_FixtureBase):

    def setUp(self):
        self.fig = self._make_fig()
        # Inject synthetic SCORE without calling the full svm pipeline.
        self.fig.data['SCORE'] = [0.8, 0.8, 0.8, 0.1, 0.1]

    def test_adds_svm_called_error_column(self):
        self.fig.svm_error_calling()
        self.assertIn('SVM_CALLED_ERROR', self.fig.data.columns)

    def test_high_scores_above_threshold_marked_as_errors(self):
        # Three consecutive scores of 0.8 → all three should be called errors
        self.fig.svm_error_calling(min_err_size=1, score_threshold=0.5)
        called = list(self.fig.data['SVM_CALLED_ERROR'])
        self.assertEqual([True, True, True, False, False], called)

    def test_min_err_size_suppresses_short_runs(self):
        # Only two consecutive low-score residues; with min_err_size=3, the run
        # of three high scores still qualifies but if we flip the scores it won't.
        self.fig.data['SCORE'] = [0.8, 0.1, 0.8, 0.1, 0.1]
        self.fig.svm_error_calling(min_err_size=2, score_threshold=0.5)
        called = list(self.fig.data['SVM_CALLED_ERROR'])
        # Residue 0 alone (score 0.8) is a run of 1 → not called with min_err_size=2
        self.assertFalse(called[0])
        # Residues 3 and 4 (scores 0.1, 0.1) — below threshold; residue 2 (0.8) alone → not called
        self.assertFalse(called[2])


# ---------------------------------------------------------------------------
# map_align
# ---------------------------------------------------------------------------

class TestModelValidationFigureMapAlign(_FixtureBase):

    def test_no_exe_adds_misaligned_column(self):
        fig = self._make_fig()
        fig.map_align(map_align_exe=None)
        self.assertIn('MISALIGNED', fig.data.columns)

    def test_no_exe_sets_all_misaligned_false(self):
        fig = self._make_fig()
        fig.map_align(map_align_exe=None)
        self.assertTrue((fig.data['MISALIGNED'] == False).all())

    def test_with_exe_stores_exe_path(self):
        fig = self._make_fig()
        with patch.object(fig, '_get_cmap_alignment'):
            fig.map_align(map_align_exe='/path/to/map_align')
        self.assertEqual('/path/to/map_align', fig.map_align_exe)

    def test_with_exe_misaligned_reflects_alignment_keys(self):
        # Simulate _get_cmap_alignment telling us residue 1 is misaligned.
        fig = self._make_fig()

        def _fake_alignment(tempdirname=None):
            fig.alignment = {1: 2}

        with patch.object(fig, '_get_cmap_alignment', side_effect=_fake_alignment):
            fig.map_align(map_align_exe='/path/to/map_align')

        misaligned = fig.data.set_index('RESNUM')['MISALIGNED'].to_dict()
        self.assertTrue(misaligned[1])
        self.assertFalse(misaligned[2])
        self.assertFalse(misaligned[3])


# ---------------------------------------------------------------------------
# count_contacts
# ---------------------------------------------------------------------------

class TestModelValidationFigureCountContacts(_FixtureBase):

    def test_adds_contacts_column(self):
        fig = self._make_fig()
        fig.count_contacts()
        self.assertIn('CONTACTS', fig.data.columns)

    def test_contacts_are_non_negative(self):
        fig = self._make_fig()
        fig.count_contacts()
        self.assertTrue((fig.data['CONTACTS'] >= 0).all())

    def test_contacts_length_matches_residue_count(self):
        fig = self._make_fig()
        fig.count_contacts()
        self.assertEqual(5, len(fig.data['CONTACTS']))


# ---------------------------------------------------------------------------
# add_plddt
# ---------------------------------------------------------------------------

class TestModelValidationFigureAddPlddt(_FixtureBase):

    def test_adds_plddt_column(self):
        fig = self._make_fig()
        fig.add_plddt()
        self.assertIn('PLDDT', fig.data.columns)

    def test_plddt_length_matches_residue_count(self):
        fig = self._make_fig()
        fig.add_plddt()
        self.assertEqual(5, len(fig.data['PLDDT']))


# ---------------------------------------------------------------------------
# Run_combined_filter
# ---------------------------------------------------------------------------

class TestModelValidationFigureRunCombinedFilter(_FixtureBase):

    def test_missing_required_features_returns_one(self):
        # Without CONTACTS/PLDDT/Q_IN_ERROR, the filter aborts and returns 1.
        fig = self._make_fig()
        result = fig.Run_combined_filter(filter_type='CMO')
        self.assertEqual(1, result)

    def test_missing_required_features_adds_zero_filter_column(self):
        fig = self._make_fig()
        fig.Run_combined_filter(filter_type='CMO')
        self.assertIn('CMO_FILTER', fig.data.columns)
        self.assertTrue((fig.data['CMO_FILTER'] == 0).all())

    def test_unknown_filter_type_returns_one(self):
        fig = self._make_fig()
        result = fig.Run_combined_filter(filter_type='UNKNOWN')
        self.assertEqual(1, result)


# ---------------------------------------------------------------------------
# draw
# ---------------------------------------------------------------------------

class TestModelValidationFigureDraw(_FixtureBase):

    def test_draw_with_all_features_off_does_not_raise(self):
        fig = self._make_fig()
        try:
            fig.draw(RUN_SVM=False, RUN_MAP_ALIGN=False, RUN_FILTERS=False)
        except Exception as e:
            self.fail(f'draw() raised unexpectedly: {e}')

    def test_draw_sets_x_label(self):
        fig = self._make_fig()
        fig.draw(RUN_SVM=False, RUN_MAP_ALIGN=False, RUN_FILTERS=False)
        self.assertEqual('Residue Number', fig.ax.get_xlabel())

    def test_draw_sets_y_label(self):
        fig = self._make_fig()
        fig.draw(RUN_SVM=False, RUN_MAP_ALIGN=False, RUN_FILTERS=False)
        self.assertEqual('Smoothed score', fig.ax.get_ylabel())


if __name__ == '__main__':
    unittest.main()
