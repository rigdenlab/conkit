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

from Bio.PDB.DSSP import DSSP
import numpy as np
import pandas as pd

from conkit.core.distance import Distance
from conkit.misc import load_validation_model, SELECTED_VALIDATION_FEATURES, ALL_VALIDATION_FEATURES
from conkit.plot.figure import Figure
from conkit.plot.tools import ColorDefinitions
from conkit.plot.tools import _isinstance, get_rmsd, get_cmap_validation_metrics, get_zscores


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
    distance_bins: list, tuple
        A list of tuples with the boundaries of the distance bins to use in the calculation [default: CASP2 bins]
    l_factor: float
        The L/N factor used to filter the contacts before finding the False Negatives [default: 0.5]

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

    def __init__(self, model, prediction, sequence, dssp, distance_bins=None, l_factor=0.5, **kwargs):
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
        distance_bins: list, tuple
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

        self.l_factor = l_factor
        self.distance_bins = distance_bins
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
    def distance_bins(self):
        return self._distance_bins

    @distance_bins.setter
    def distance_bins(self, distance_bins):
        if distance_bins is None:
            self._distance_bins = ((0, 4), (4, 6), (6, 8), (8, 10), (10, 12), (12, 14),
                                   (14, 16), (16, 18), (18, 20), (20, np.inf))
        else:
            Distance._assert_valid_bins(distance_bins)
            self._distance_bins = distance_bins

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, sequence):
        if sequence and _isinstance(sequence, "Sequence"):
            self._sequence = sequence
        else:
            raise TypeError("Invalid hierarchy type for sequence: %s" % sequence.__class__.__name__)

    @property
    def dssp(self):
        return self._dssp

    @dssp.setter
    def dssp(self, dssp):
        if dssp and _isinstance(dssp, DSSP):
            self._dssp = dssp
        else:
            raise TypeError("Invalid hierarchy type for dssp: %s" % dssp.__class__.__name__)

    @property
    def prediction(self):
        return self._prediction

    @prediction.setter
    def prediction(self, prediction):
        if prediction and _isinstance(prediction, "Distogram"):
            self._prediction = prediction
        else:
            raise TypeError("Invalid hierarchy type for prediction: %s" % prediction.__class__.__name__)

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model):
        if model and _isinstance(model, "Distogram"):
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
            distogram.reshape_bins(self.distance_bins)

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
        dssp = []
        for residue in sorted(self.dssp.keys(), key=lambda x: x[1][1]):
            resnum = residue[1][1]
            if resnum in self.absent_residues:
                dssp.append((resnum, np.nan, np.nan, np.nan, np.nan))
            acc = self.dssp[residue][3]
            ss2 = (1, 0, 0) if self.dssp[residue][2] in ('-', 'T', 'S') else (0, 1, 0) if self.dssp[residue][2] in ('H', 'G', 'I') else (0, 0, 1)
            dssp.append((resnum, *ss2, acc))

        dssp = pd.DataFrame(dssp)
        dssp.columns = ['RESNUM', 'COIL', 'HELIX', 'SHEET', 'ACC']

        return dssp

    @staticmethod
    def get_feature_df(predicted_dict, dssp, *metrics):
        features = []
        for residue_features in zip(sorted(predicted_dict.keys()), *metrics):
            features.append((*residue_features,))

        features = pd.DataFrame(features)
        features.columns = ALL_VALIDATION_FEATURES
        features = features.merge(dssp, how='inner', on=['RESNUM'])

        return features

    def draw(self):

        model_distogram = self._prepare_distogram(self.model.copy())
        prediction_distogram = self._prepare_distogram(self.prediction.copy())
        model_cmap = self._prepare_contactmap(self.model.copy())
        model_cmap_dict = model_cmap.as_dict()
        prediction_cmap = self._prepare_contactmap(self.prediction.copy())
        predicted_cmap_dict = prediction_cmap.as_dict()

        dssp = self._parse_dssp()
        rmsd, rmsd_smooth = get_rmsd(prediction_distogram, model_distogram)
        cmap_metrics, cmap_metrics_smooth = get_cmap_validation_metrics(model_cmap_dict, predicted_cmap_dict, self.absent_residues, len(self.sequence))
        zscore_metrics = get_zscores(model_distogram, predicted_cmap_dict, self.absent_residues, rmsd, *cmap_metrics)
        features_df = self.get_feature_df(predicted_cmap_dict, dssp, rmsd_smooth, *cmap_metrics_smooth, *zscore_metrics)

        y_pred = []
        y_pred_proba = []
        for resnum in sorted(predicted_cmap_dict.keys()):
            print(resnum)
            if self.absent_residues and resnum in self.absent_residues:
                y_pred.append(np.nan)
                y_pred_proba.append(np.nan)
                continue
            residue_features = features_df.loc[features_df.RESNUM == resnum][SELECTED_VALIDATION_FEATURES].values
            scaled_features = self.scaler.transform(residue_features)
            y_pred.append(self.classifier.predict(scaled_features)[0])
            y_pred_proba.append(self.classifier.predict_proba(scaled_features)[0, 1])

        line_kwargs = dict(linestyle="--", linewidth=1.0, alpha=0.5, color=ColorDefinitions.MISMATCH, zorder=1)
        self.ax.axhline(0.5, **line_kwargs)
        outlier_plot = [self.ax.axvline(0, ymin=0, ymax=0, label="Outlier", **line_kwargs)]

        self.ax.set_xlabel('Residue Number')
        self.ax.set_ylabel('Score')

        score_plot = self.ax.plot(y_pred_proba, color='#3299a8', label='Score')
        if self.legend:
            plots = score_plot + outlier_plot
            labels = [l.get_label() for l in plots]
            self.ax.legend(plots, labels, bbox_to_anchor=(0.0, 1.02, 1.0, 0.102), loc=3,
                           ncol=3, mode="expand", borderaxespad=0.0, scatterpoints=1)

        # TODO: deprecate this in 0.10
        if self._file_name:
            self.savefig(self._file_name, dpi=self._dpi)
