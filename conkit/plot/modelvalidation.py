# BSD 3-Clause License
#
# Copyright (c) 2016-21, University of Liverpool
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice, this
#   list of conditions and the following disclaimer.
#
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# * Neither the name of the copyright holder nor the names of its
#   contributors may be used to endorse or promote products derived from
#   this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
# SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
# OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
"""A module to produce a model validation plot"""

from __future__ import division
from __future__ import print_function

import os
from Bio.PDB.DSSP import DSSP
import numpy as np
import pandas as pd
import tempfile

from conkit.applications import MapAlignCommandline
from conkit.core.distance import Distance
import conkit.io
from conkit.misc import load_validation_model, SELECTED_VALIDATION_FEATURES, ALL_VALIDATION_FEATURES
from conkit.plot.figure import Figure
import conkit.plot.tools as tools

LINEKWARGS = dict(linestyle="--", linewidth=1.0, alpha=0.5, color=tools.ColorDefinitions.MISMATCH, zorder=1)
MARKERKWARGS = dict(marker='|', linestyle='None')
_MARKERKWARGS = dict(marker='s', linestyle='None')



class ModelValidationFigure(Figure):
    """A Figure object specifc for a model validation. This figure represents the proabbility that each given residue
    in the model is involved in a model error. This is donw by feeding a trained classfier the differences observed
    between the predicted distogram and the observed inter-residue contacts and distances at the PDB model.

    Attributes
    ----------
    model: :obj:`~conkit.core.distogram.Distogram`
       The PDB model that will be validated
    prediction: :obj:`~conkit.core.distogram.Distogram`
       The distogram with the residue distance predictions
    sequence: :obj:`~conkit.core.sequence.Sequence`
       The sequence of the structure
    dssp: :obj:`Bio.PDB.DSSP.DSSP`
        The DSSP output for the PDB model that will be validated
    map_align_exe: str
        The path to map_align executable [default: None]
    dist_bins: list, tuple
        A list of tuples with the boundaries of the distance bins to use in the calculation [default: CASP2 bins]
    l_factor: float
        The L/N factor used to filter the contacts before finding the False Negatives [default: 0.5]
    absent_residues: set
        The residues not observed in the model that will be validated (only if in PDB format)

    Examples
    --------
    >>> from Bio.PDB import PDBParser
    >>> from Bio.PDB.DSSP import DSSP
    >>> p = PDBParser()
    >>> structure = p.get_structure('TOXD', 'toxd/toxd.pdb')[0]
    >>> dssp = DSSP(structure, 'toxd/toxd.pdb', dssp='mkdssp', acc_array='Wilke')
    >>> import conkit
    >>> sequence = conkit.io.read('toxd/toxd.fasta', 'fasta').top
    >>> model = conkit.io.read('toxd/toxd.pdb', 'pdb').top_map
    >>> prediction = conkit.io.read('toxd/toxd.npz', 'rosettanpz').top_map
    >>> conkit.plot.ModelValidationFigure(model, prediction, sequence, dssp)

    """

    def __init__(self, model, prediction, sequence, dssp, map_align_exe=None, dist_bins=None, l_factor=0.5, **kwargs):
        """A new model validation plot

        Parameters
        ----------
        model: :obj:`~conkit.core.distogram.Distogram`
            The PDB model that will be validated
        prediction: :obj:`~conkit.core.distogram.Distogram`
            The distogram with the residue distance predictions
        sequence: :obj:`~conkit.core.sequence.Sequence`
            The sequence of the structure
        dssp: :obj:`Bio.PDB.DSSP.DSSP`
            The DSSP output for the PDB model that will be validated
        map_align_exe: str
            The path to map_align executable [default: None]
        dist_bins: list, tuple
            A list of tuples with the boundaries of the distance bins to use in the calculation [default: CASP2 bins]
        l_factor: float
            The L/N factor used to filter the contacts before finding the False Negatives [default: 0.5]

        **kwargs
           General :obj:`~conkit.plot.figure.Figure` keyword arguments

        """
        super(ModelValidationFigure, self).__init__(**kwargs)
        self._model = None
        self._prediction = None
        self._sequence = None
        self._dssp = None
        self._distance_bins = None
        self.features = None

        self.map_align_exe = map_align_exe
        self.l_factor = l_factor
        self.dist_bins = dist_bins
        self.model = model
        self.prediction = prediction
        self.sequence = sequence
        self.dssp = dssp
        self.classifier, self.scaler = load_validation_model()
        self.absent_residues = self._get_absent_residues()

        self.draw()

    def __repr__(self):
        return self.__class__.__name__

    @property
    def dist_bins(self):
        return self._dist_bins

    @dist_bins.setter
    def dist_bins(self, dist_bins):
        if dist_bins is None:
            self._dist_bins = ((0, 4), (4, 6), (6, 8), (8, 10), (10, 12), (12, 14),
                               (14, 16), (16, 18), (18, 20), (20, np.inf))
        else:
            Distance._assert_valid_bins(dist_bins)
            self._dist_bins = dist_bins

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, sequence):
        if sequence and tools._isinstance(sequence, "Sequence"):
            self._sequence = sequence
        else:
            raise TypeError("Invalid hierarchy type for sequence: %s" % sequence.__class__.__name__)

    @property
    def dssp(self):
        return self._dssp

    @dssp.setter
    def dssp(self, dssp):
        if dssp and tools._isinstance(dssp, DSSP):
            self._dssp = dssp
        else:
            raise TypeError("Invalid hierarchy type for dssp: %s" % dssp.__class__.__name__)

    @property
    def prediction(self):
        return self._prediction

    @prediction.setter
    def prediction(self, prediction):
        if prediction and tools._isinstance(prediction, "Distogram"):
            self._prediction = prediction
        else:
            raise TypeError("Invalid hierarchy type for prediction: %s" % prediction.__class__.__name__)

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model):
        if model and tools._isinstance(model, "Distogram"):
            self._model = model
        else:
            raise TypeError("Invalid hierarchy type for model: %s" % model.__class__.__name__)

    def _get_absent_residues(self):
        """Get a set of residues absent from the :attr:`~conkit.plot.ModelValidationFigure.model` and
        :attr:`~conkit.plot.ModelValidationFigure.prediction`. Only distograms originating from PDB files
        are considered."""

        absent_residues = []
        if self.model.original_file_format == "PDB":
            absent_residues += self.model.get_absent_residues(len(self.sequence))
        if self.prediction.original_file_format == "PDB":
            absent_residues += self.prediction.get_absent_residues(len(self.sequence))
        return set(absent_residues)

    def _prepare_distogram(self, distogram):
        """General operations to prepare a :obj:`~conkit.core.distogram.Distogram` instance before plotting."""
        distogram.get_unique_distances(inplace=True)
        distogram.sequence = self.sequence
        distogram.set_sequence_register()

        if distogram.original_file_format != "PDB":
            distogram.reshape_bins(self.dist_bins)

        return distogram

    def _prepare_contactmap(self, distogram):
        """General operations to prepare a :obj:`~conkit.core.contactmap.ContactMap` instance before plotting."""
        contactmap = distogram.as_contactmap()
        contactmap.sequence = self.sequence
        contactmap.set_sequence_register()
        contactmap.remove_neighbors(inplace=True)

        if distogram.original_file_format != "PDB":
            contactmap.sort("raw_score", reverse=True, inplace=True)
            contactmap.slice_map(seq_len=len(self.sequence), l_factor=self.l_factor, inplace=True)

        return contactmap

    def _parse_dssp(self):
        """Parse :obj:`Bio.PDB.DSSP.DSSP` into a :obj:`pandas.DataFrame` with secondary structure information
        about the model"""
        dssp = []
        for residue in sorted(self.dssp.keys(), key=lambda x: x[1][1]):
            resnum = residue[1][1]
            if resnum in self.absent_residues:
                dssp.append((resnum, np.nan, np.nan, np.nan, np.nan))
                continue
            acc = self.dssp[residue][3]
            if self.dssp[residue][2] in ('-', 'T', 'S'):
                ss2 = (1, 0, 0)
            elif self.dssp[residue][2] in ('H', 'G', 'I'):
                ss2 = (0, 1, 0)
            else:
                ss2 = (0, 0, 1)
            dssp.append((resnum, *ss2, acc))

        dssp = pd.DataFrame(dssp)
        dssp.columns = ['RESNUM', 'COIL', 'HELIX', 'SHEET', 'ACC']

        return dssp

    def get_cmap_alignment(self):
        """Obtain a contact map alignment between :attr:`~conkit.plot.ModelValidationFigure.model` and
        :attr:`~conkit.plot.ModelValidationFigure.prediction` and get the misaligned residues"""
        with tempfile.TemporaryDirectory() as tmpdirname:
            contact_map_a = os.path.join(tmpdirname, 'contact_map_a.mapalign')
            contact_map_b = os.path.join(tmpdirname, 'contact_map_b.mapalign')
            conkit.io.write(contact_map_a, 'mapalign', self.prediction)
            conkit.io.write(contact_map_b, 'mapalign', self.model)

            map_align_cline = MapAlignCommandline(
                cmd=self.map_align_exe,
                contact_map_a=contact_map_a,
                contact_map_b=contact_map_b)

            stdout, stderr = map_align_cline()
            alignments = tools.parse_map_align_stdout(stdout)
            misaligned_resnums = {res for misalignment in alignments for res in misalignment.residue_range}
            return misaligned_resnums

    def get_feature_df(self, predicted_dict, dssp, *metrics):
        """Create a :obj:`pandas.DataFrame` with the features of the residues in the model"""
        features = []
        for residue_features in zip(sorted(predicted_dict.keys()), *metrics):
            features.append((*residue_features,))

        features = pd.DataFrame(features)
        features.columns = ALL_VALIDATION_FEATURES
        features = features.merge(dssp, how='inner', on=['RESNUM'])

        if self.map_align_exe is not None:
            misaligned_residues = self.get_cmap_alignment()
            features['MISALIGNED'] = features.RESNUM.isin(misaligned_residues)
        else:
            features['MISALIGNED'] = False

        return features

    def _add_legend(self):
        """Adds legend to the :obj:`~conkit.plot.ModelValidationFigure`"""
        _error = self.ax.plot([], [], c=tools.ColorDefinitions.VALIDATION_ERROR, label='Error', **_MARKERKWARGS)
        _correct = self.ax.plot([], [], c=tools.ColorDefinitions.VALIDATION_CORRECT, label='Correct', **_MARKERKWARGS)
        _threshold_line = [self.ax.axvline(0, ymin=0, ymax=0, label="Threshold", **LINEKWARGS)]
        _score_plot = self.ax.plot([], [], color=tools.ColorDefinitions.VALIDATION_SCORE, label='Smoothed Score')
        plots = _score_plot + _threshold_line + _correct + _error

        if self.map_align_exe is not None:
            _misaligned = self.ax.plot([], [], c=tools.ColorDefinitions.MISALIGNED, label='Misaligned', **_MARKERKWARGS)
            _aligned = self.ax.plot([], [], c=tools.ColorDefinitions.ALIGNED, label='Aligned', **_MARKERKWARGS)
            plots += _misaligned + _aligned

        labels = [l.get_label() for l in plots]
        self.ax.legend(plots, labels, bbox_to_anchor=(0.0, 1.02, 1.0, 0.102), loc=3,
                       ncol=3, mode="expand", borderaxespad=0.0, scatterpoints=1)

    def _predict_score(self, resnum):
        """Predict whether a given residue is part of a model error or not.

        Parameters
        ----------
        resnum : int
           The residue number of the residue that will be validated

        Returns
        -------
        score: float
           The predicted probability that the residue of interest is part of a model error
        """
        residue_features = self.features.loc[self.features.RESNUM == resnum][SELECTED_VALIDATION_FEATURES]
        if (self.absent_residues and resnum in self.absent_residues) or residue_features.isnull().values.any():
            return np.nan
        scaled_features = self.scaler.transform(residue_features.values)
        return self.classifier.predict_proba(scaled_features)[0, 1]

    def draw(self):
        model_distogram = self._prepare_distogram(self.model.copy())
        prediction_distogram = self._prepare_distogram(self.prediction.copy())
        model_cmap = self._prepare_contactmap(self.model.copy())
        model_dict = model_cmap.as_dict()
        prediction_cmap = self._prepare_contactmap(self.prediction.copy())
        predicted_dict = prediction_cmap.as_dict()

        dssp = self._parse_dssp()
        cmap_metrics, cmap_metrics_smooth = tools.get_cmap_validation_metrics(model_dict, predicted_dict, self.sequence, self.absent_residues)
        rmsd, rmsd_smooth = tools.get_rmsd(prediction_distogram, model_distogram)
        zscore_metrics = tools.get_zscores(model_distogram, predicted_dict, self.absent_residues, rmsd, *cmap_metrics)
        self.features = self.get_feature_df(predicted_dict, dssp, rmsd_smooth, *cmap_metrics, *cmap_metrics_smooth, *zscore_metrics)

        misaligned_residues = set(self.features.loc[self.features.MISALIGNED, 'RESNUM'])
        residue_scores = []

        for resnum in sorted(predicted_dict.keys()):
            _score = self._predict_score(resnum)
            residue_scores.append(_score)
            if _score > 0.5:
                color = tools.ColorDefinitions.VALIDATION_ERROR
            else:
                color = tools.ColorDefinitions.VALIDATION_CORRECT
            self.ax.plot(resnum - 1, -0.01, mfc=color, c=color, **MARKERKWARGS)
            if self.map_align_exe is not None:
                if resnum in misaligned_residues:
                    color = tools.ColorDefinitions.MISALIGNED
                else:
                    color = tools.ColorDefinitions.ALIGNED
                self.ax.plot(resnum - 1, -0.05, mfc=color, c=color, **MARKERKWARGS)

        residue_scores_smooth = tools.convolution_smooth_values(np.nan_to_num(residue_scores))
        self.ax.axhline(0.5, **LINEKWARGS)
        self.ax.plot(residue_scores_smooth, color=tools.ColorDefinitions.VALIDATION_SCORE)
        self.ax.set_xlabel('Residue Number')
        self.ax.set_ylabel('Score')

        if self.legend:
            self._add_legend()

        # TODO: deprecate this in 0.10
        if self._file_name:
            self.savefig(self._file_name, dpi=self._dpi)
