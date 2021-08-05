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
"""A module to produce a domain boundary plot"""

from __future__ import division
from __future__ import print_function

import numpy as np

from conkit.core.distance import Distance
from conkit.core.distogram import Distogram
from conkit.plot.figure import Figure
from conkit.plot.tools import ColorDefinitions
from conkit.plot.tools import _isinstance, convolution_smooth_values, find_validation_outliers, get_fn_profile


class ModelValidationFigure(Figure):
    """A Figure object specifically for a model validation illustration. This figure represents the
    RMSD profile along the sequence between the predicted distances at the distogram and the observed
    distances at the PDB model.

    Attributes
    ----------
    model: :obj:`~conkit.core.distogram.Distogram`
       The PDB model that will be validated
    prediction: :obj:`~conkit.core.distogram.Distogram`
       The distogram with the residue distance predictions
    sequence: :obj:`~conkit.core.sequence.Sequence`
       The sequence of the structure
    distance_bins: list, tuple
        A list of tuples with the boundaries of the distance bins to use in the calculation [default: CASP2 bins]
    use_weights: bool
        If set then the wRMSD is calculated using the confidence scores for predicted distances [default: True]
    l_factor: float
            The L/N factor used to filter the contacts before finding the False Negatives [default: 0.5]

    Examples
    --------
    >>> import conkit
    >>> sequence = conkit.io.read('toxd/toxd.fasta', 'fasta').top
    >>> model = conkit.io.read('toxd/toxd.pdb', 'pdb').top_map
    >>> prediction = conkit.io.read('toxd/toxd.npz', 'rosettanpz').top_map
    >>> conkit.plot.ModelValidationFigure(model, prediction, sequence)

    """

    def __init__(self, model, prediction, sequence, distance_bins=None, use_weights=False, l_factor=0.5, **kwargs):
        """A new model validation plot

        Parameters
        ----------
        model: :obj:`~conkit.core.distogram.Distogram`
            The PDB model that will be validated
        prediction: :obj:`~conkit.core.distogram.Distogram`
            The distogram with the residue distance predictions
        sequence: :obj:`~conkit.core.sequence.Sequence`
            The sequence of the structure
        distance_bins: list, tuple
            A list of tuples with the boundaries of the distance bins to use in the calculation [default: CASP2 bins]
        use_weights: bool
            If set then the wRMSD is calculated using the confidence scores for predicted distances [default: True]
        l_factor: float
            The L/N factor used to filter the contacts before finding the False Negatives [default: 0.5]

        **kwargs
           General :obj:`~conkit.plot.figure.Figure` keyword arguments

        """
        super(ModelValidationFigure, self).__init__(**kwargs)
        self._model = None
        self._prediction = None
        self._sequence = None
        self._distance_bins = None

        self.use_weights = use_weights
        self.l_factor = l_factor
        self.distance_bins = distance_bins
        self.model = model
        self.prediction = prediction
        self.sequence = sequence

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
            raise TypeError("Invalid hierarchy type: %s" % sequence.__class__.__name__)

    @property
    def prediction(self):
        return self._prediction

    @prediction.setter
    def prediction(self, prediction):
        if prediction and _isinstance(prediction, "Distogram"):
            self._prediction = prediction
        else:
            raise TypeError("Invalid hierarchy type: %s" % prediction.__class__.__name__)

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, model):
        if model and _isinstance(model, "Distogram"):
            self._model = model
        else:
            raise TypeError("Invalid hierarchy type: %s" % model.__class__.__name__)

    def _get_absent_residues(self):
        """Get a list of the residues absent from the :attr:`~conkit.plot.ModelValidationFigure.model` and
        :attr:`~conkit.plot.ModelValidationFigure.prediction`. Only distograms originating from PDB files
        are considered."""

        absent_residues = []
        if self.model.original_file_format == "PDB":
            absent_residues += self.model.get_absent_residues(len(self.sequence))
        if self.prediction.original_file_format == "PDB":
            absent_residues += self.prediction.get_absent_residues(len(self.sequence))
        return tuple(absent_residues)

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

    def draw(self):

        model_distogram = self._prepare_distogram(self.model.copy())
        prediction_distogram = self._prepare_distogram(self.prediction.copy())
        rmsd_raw = Distogram.calculate_rmsd(prediction_distogram, model_distogram, calculate_wrmsd=self.use_weights)
        rmsd_smooth = convolution_smooth_values(rmsd_raw, 10)
        rmsd_profile = np.nan_to_num(rmsd_smooth)

        model_contactmap = self._prepare_contactmap(self.model.copy())
        prediction_contactmap = self._prepare_contactmap(self.prediction.copy())
        fn_raw = get_fn_profile(model_contactmap, prediction_contactmap, ignore_residues=self._get_absent_residues())
        fn_smooth = convolution_smooth_values(fn_raw, 5)
        fn_profile = np.nan_to_num(fn_smooth)

        outliers = find_validation_outliers(rmsd_raw, rmsd_profile, fn_raw, fn_profile)
        line_kwargs = dict(linestyle="--", linewidth=1.0, alpha=0.5, color=ColorDefinitions.MISMATCH, zorder=1)
        for outlier in outliers:
            self.ax.axvline(outlier, **line_kwargs)
        outlier_plot = [self.ax.axvline(0, ymin=0, ymax=0, label="Outlier", **line_kwargs)]

        self.ax.set_xlabel('Residue Number')
        self.ax.set_ylabel('wRMSD' if self.use_weights else 'RMSD')
        rmsd_plot = self.ax.plot(rmsd_profile, color='#3299a8', label='wRMSD' if self.use_weights else 'RMSD')
        ax2 = self.ax.twinx()
        ax2.set_ylabel('FN count')
        fn_plot = ax2.plot(fn_profile, color='#325da8', label='FN count')
        if self.legend:
            plots = fn_plot + rmsd_plot + outlier_plot
            labels = [l.get_label() for l in plots]
            self.ax.legend(plots, labels, bbox_to_anchor=(0.0, 1.02, 1.0, 0.102), loc=3,
                           ncol=3, mode="expand", borderaxespad=0.0, scatterpoints=1)

        # TODO: deprecate this in 0.10
        if self._file_name:
            self.savefig(self._file_name, dpi=self._dpi)
