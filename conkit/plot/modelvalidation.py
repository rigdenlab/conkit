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
"""A module to produce a model validation plot

It uses one external program:

   map_align for contact map alignment

*** This program needs to be installed separately from https://github.com/sokrypton/map_align***
"""

from __future__ import division
from __future__ import print_function

import logging
import os
import time

from matplotlib.patches import Patch
import matplotlib.transforms as mtransforms

logger = logging.getLogger(__name__)
from Bio.PDB.DSSP import DSSP
from Bio.PDB import PDBParser, MMCIFParser
import numpy as np
import pandas as pd
import tempfile

from conkit.applications import MapAlignCommandline
from conkit.core.distance import Distance
import conkit.io
from conkit.misc import load_specific_validation_model, load_validation_model, load_filter, SELECTED_VALIDATION_FEATURES, SELECTED_VALIDATION_FEATURES_DICT, ALL_VALIDATION_FEATURES
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

    def __init__(self, model, prediction, sequence, map_align_exe=None, dist_bins=None, l_factor=0.5, **kwargs):
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
        self._distance_bins = None
        self.data = None
        self.alignment = {}
        self.filter_threshold = {}
        self.svm_name = None
        self.sorted_scores = None
        self.smooth_scores = None
        self.map_align_exe = None
        self.out_dir = None
        self._map_align_ran = False

        if len(sequence) < 5:
            raise ValueError('Cannot validate a model with less than 5 residues')

        self.l_factor = l_factor
        self.dist_bins = dist_bins
        self.model = model
        self.prediction = prediction
        self.sequence = sequence
        self.absent_residues = self._get_absent_residues()
        model_cmap = self._prepare_contactmap(self.model.copy())
        model_dict = model_cmap.as_dict()

        self.data = pd.DataFrame()

        self.data['RESNUM'] = model_dict.keys()
        self.data['MISALIGNED'] = False        
        self.data['SCORE'] = 0
        #self.data['CONTACTS'] = 0        
        #self.data['PLDDT'] = 0
        #self.data['Q_IN_ERROR'] = ''  

    def calculate_features(self,z_radius = 10, max_distance = None):

        model_distogram = self._prepare_distogram(self.model.copy())
        prediction_distogram = self._prepare_distogram(self.prediction.copy())
        model_cmap = self._prepare_contactmap(self.model.copy())
        model_dict = model_cmap.as_dict()
        prediction_cmap = self._prepare_contactmap(self.prediction.copy())
        predicted_dict = prediction_cmap.as_dict()

        cmap_metrics, cmap_metrics_smooth = tools.get_cmap_validation_metrics(model_dict, predicted_dict,
                                                                              self.sequence, self.absent_residues)
        rmsd, rmsd_smooth = tools.get_rmsd(prediction_distogram, model_distogram, max_distance = max_distance )
        zscore_metrics = tools.get_zscores(model_distogram, predicted_dict, self.absent_residues, rmsd, *cmap_metrics, population_radius = z_radius)

        self._parse_data(predicted_dict, rmsd, rmsd_smooth, *cmap_metrics, *cmap_metrics_smooth, *zscore_metrics)

        #self.draw()

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
        if self.model.original_file_format == "pdb":
            absent_residues += self.model.get_absent_residues(len(self.sequence))
        if self.prediction.original_file_format == "pdb":
            absent_residues += self.prediction.get_absent_residues(len(self.sequence))
        return set(absent_residues)

    def _prepare_distogram(self, distogram):
        """General operations to prepare a :obj:`~conkit.core.distogram.Distogram` instance before plotting."""
        distogram.get_unique_distances(inplace=True)
        distogram.sequence = self.sequence
        distogram.set_sequence_register()

        if distogram.original_file_format != "pdb":
            distogram.reshape_bins(self.dist_bins)

        return distogram

    def _prepare_contactmap(self, distogram):
        """General operations to prepare a :obj:`~conkit.core.contactmap.ContactMap` instance before plotting."""
        contactmap = distogram.as_contactmap()
        contactmap.sequence = self.sequence
        contactmap.set_sequence_register()
        contactmap.remove_neighbors(inplace=True)

        if distogram.original_file_format != "pdb":
            contactmap.sort("raw_score", reverse=True, inplace=True)
            contactmap.slice_map(seq_len=len(self.sequence), l_factor=self.l_factor, inplace=True)

        return contactmap

    def _parse_dssp(self, dssp):
        """Parse :obj:`Bio.PDB.DSSP.DSSP` into a :obj:`pandas.DataFrame` with secondary structure information
        about the model"""

        if not tools._isinstance(dssp, DSSP):
            raise TypeError("Invalid hierarchy type for dssp: %s" % dssp.__class__.__name__)

        _dssp_list = []
        for residue in sorted(dssp.keys(), key=lambda x: x[1][1]):
            resnum = residue[1][1]
            if resnum in self.absent_residues:
                _dssp_list.append((resnum, np.nan, np.nan, np.nan, np.nan))
                continue
            acc = dssp[residue][3]
            if dssp[residue][2] in ('-', 'T', 'S'):
                ss2 = (1, 0, 0)
            elif dssp[residue][2] in ('H', 'G', 'I'):
                ss2 = (0, 1, 0)
            else:
                ss2 = (0, 0, 1)
            _dssp_list.append((resnum, *ss2, acc))

        dssp = pd.DataFrame(_dssp_list)
        dssp.columns = ['RESNUM', 'COIL', 'HELIX', 'SHEET', 'ACC']
        return dssp

    def _get_cmap_alignment(self,tempdirname=None):
        """Obtain a contact map alignment between :attr:`~conkit.plot.ModelValidationFigure.model` and
        :attr:`~conkit.plot.ModelValidationFigure.prediction` and get the misaligned residues"""

        prediction_cmap = self._prepare_contactmap(self.prediction.copy())
        model_cmap = self._prepare_contactmap(self.model.copy())

        if tempdirname:
            contact_map_a = os.path.join(tempdirname, 'contact_map_a.mapalign')
            contact_map_b = os.path.join(tempdirname, 'contact_map_b.mapalign')
            conkit.io.write(contact_map_a, 'mapalign', prediction_cmap)
            conkit.io.write(contact_map_b, 'mapalign', model_cmap)
            time.sleep(2)
            map_align_cline = MapAlignCommandline(
                cmd=self.map_align_exe,
                contact_map_a=contact_map_a,
                contact_map_b=contact_map_b)

            stdout, stderr = map_align_cline()
            self.alignment = tools.parse_map_align_stdout(stdout)

        else:
            with tempfile.TemporaryDirectory() as tmpdirname:
                contact_map_a = os.path.join(tmpdirname, 'contact_map_a.mapalign')
                contact_map_b = os.path.join(tmpdirname, 'contact_map_b.mapalign')
                conkit.io.write(contact_map_a, 'mapalign', self.prediction)
                conkit.io.write(contact_map_b, 'mapalign', self.model)
                time.sleep(2)
                map_align_cline = MapAlignCommandline(
                    cmd=self.map_align_exe,
                    contact_map_a=contact_map_a,
                    contact_map_b=contact_map_b)

                stdout, stderr = map_align_cline()
                self.alignment = tools.parse_map_align_stdout(stdout)


    def _parse_data(self, predicted_dict, *metrics):
        """Create a :obj:`pandas.DataFrame` with the features of the residues in the model"""
        _features = []
        for residue_features in zip(sorted(predicted_dict.keys()), *metrics):
            _features.append((*residue_features,))

        feature_df = pd.DataFrame(_features)
        feature_df.columns = ALL_VALIDATION_FEATURES
        
        self.data = self.data.merge(feature_df, how='inner', on =['RESNUM'])


    def _add_legend(self, RUN_SVM=True, RUN_MAP_ALIGN=True, RUN_FILTERS=True,
                    n_contacts_per_res=2, plddt_threshold=65,
                    has_rf_filter=False, has_cmo_filter=False, combine_filters=False):
        """Adds legend to the :obj:`~conkit.plot.ModelValidationFigure`"""

        handles = []

        if RUN_SVM:
            handles += self.ax.plot([], [], color=tools.ColorDefinitions.SCORE, label='Smoothed Score')
            handles += [self.ax.axvline(0, ymin=0, ymax=0, label='Score Threshold', **LINEKWARGS)]
            handles += self.ax.plot([], [], c=tools.ColorDefinitions.CORRECT, label='Predicted Correct', **_MARKERKWARGS)
            handles += self.ax.plot([], [], c=tools.ColorDefinitions.ERROR, label='Predicted Error', **_MARKERKWARGS)
            if has_rf_filter:
                handles += [Patch(facecolor=tools.ColorDefinitions.ERROR, hatch='///',
                                  edgecolor='white', label='Predicted Error (filter uncertain)')]

        if RUN_MAP_ALIGN and self._map_align_ran:
            handles += self.ax.plot([], [], c=tools.ColorDefinitions.ALIGNED, label='Aligned', **_MARKERKWARGS)
            handles += self.ax.plot([], [], c=tools.ColorDefinitions.MISALIGNED, label='Misaligned', **_MARKERKWARGS)
            if has_cmo_filter:
                handles += [Patch(facecolor=tools.ColorDefinitions.MISALIGNED, hatch='///',
                                  edgecolor='white', label='Misaligned (filter uncertain)')]

        if RUN_FILTERS and not combine_filters:
            if 'CONTACTS' in self.data.columns:
                handles += self.ax.plot([], [], c=tools.ColorDefinitions.SUFFICIENT_CONTACTS,
                                        label='Sufficient contacts', **_MARKERKWARGS)
                handles += self.ax.plot([], [], c=tools.ColorDefinitions.LOW_CONTACTS,
                                        label='Low contacts <' + str(n_contacts_per_res), **_MARKERKWARGS)

            if 'PLDDT' in self.data.columns:
                color_scheme = tools.ColorDefinitions.PLDDT_COLORS
                for th in sorted(color_scheme.keys()):
                    handles += self.ax.plot([], [], c=color_scheme[th], label='pLDDT <' + str(th), **_MARKERKWARGS)

            if 'Q_IN_ERROR' in self.data.columns:
                color_scheme = tools.ColorDefinitions.Q_COLORS
                thresholds = sorted(color_scheme.keys())
                handles += self.ax.plot([], [], c=color_scheme[thresholds[2]], label='Q > 0.5', **_MARKERKWARGS)
                handles += self.ax.plot([], [], c=color_scheme[thresholds[1]], label='Q < 0.5', **_MARKERKWARGS)
                handles += self.ax.plot([], [], c=color_scheme[thresholds[0]], label='Gesamt failed to align', **_MARKERKWARGS)

        self.ax.legend(handles, [h.get_label() for h in handles],
                       bbox_to_anchor=(0.5, 1.02), loc='lower center',
                       ncol=2, borderaxespad=0.0, scatterpoints=1)

    def _predict_score(self, resnum):
        """Predict whether a given residue is part of a model error or not"""
        residue_features = self.data.loc[self.data.RESNUM == resnum][SELECTED_VALIDATION_FEATURES_DICT[self.svm_name]]

        if (self.absent_residues and resnum in self.absent_residues) or residue_features.isnull().values.any():
            return np.nan
        scaled_features = self.scaler.transform(residue_features.values)
        return self.classifier.predict_proba(scaled_features)[0, 1]


    def svm(self,ext_info,moltype='Protein',prediction_type='DIST',sec_struc_info='DSSP'):

        if moltype == 'Protein':
            self.svm_name = 'Protein_DSSP_AF2_dist_'
            self.classifier, self.scaler = load_validation_model()        
            if ext_info==None: 
                self.ext_info=pd.DataFrame()
                self.ext_info['RESNUM'] = self.data['RESNUM'].copy()
                self.ext_info['COIL'], self.ext_info['HELIX'], self.ext_info['SHEET'], self.ext_info['ACC'] = 0, 0, 0, 0
            else: 
                self.ext_info = self._parse_dssp(ext_info)


            self.data = self.data.merge(self.ext_info, how='inner', on=['RESNUM'])

        elif moltype == 'RNA':
            if sec_struc_info == 'DNATCO':
                if prediction_type == 'DIST':
                    self.svm_name = 'RNA_DNATCO_AF3_dist_'
                if prediction_type == 'STRUCT':
                    self.svm_name = 'RNA_DNATCO_AF3_struct_'
            else :
                if prediction_type == 'DIST':
                    self.svm_name = 'RNA_AF3_dist_'
                if prediction_type == 'STRUCT':
                    self.svm_name = 'RNA_AF3_struct_'

            self.classifier, self.scaler = load_specific_validation_model(self.svm_name)   # load the correct SVM for the given style of input (Distogram or structure and with or without dnatco) should  potentially make protein side compliant to this
                 
            if ext_info==None: 
                self.ext_info=pd.DataFrame()
                self.ext_info['RESNUM'] = self.data['RESNUM'].copy()
                self.ext_info['ACC'] = 0
            else: 
                self.ext_info = pd.DataFrame.from_dict(ext_info)

            self.data = self.data.merge(self.ext_info, how='inner', on=['RESNUM'])

        self.data['SCORE'] = self.data['RESNUM'].apply(lambda x: self._predict_score(x))


    def svm_error_calling(self,min_err_size=1,score_threshold=0.5):

        score_list = list(self.data['SCORE'])
        index_threshold_mask = [True if s >= score_threshold else False for s in score_list]
        error_called_indices = tools.all_consecutive_true_indices(index_threshold_mask, count=min_err_size)
        self.data['SVM_CALLED_ERROR'] = [True if i in error_called_indices else False for i in range(len(score_list))]


    def map_align(self,map_align_exe=None,temp_dir_name=None):

        self.map_align_exe = map_align_exe

        if self.map_align_exe is not None:
            try:
                self._get_cmap_alignment(tempdirname=temp_dir_name)
                self.data['MISALIGNED'] = self.data.RESNUM.isin(self.alignment.keys())
                self._map_align_ran = True
            except Exception as e:
                logger.warning("map_align failed; CMO alignment bar will be omitted from the figure. Error: %s", e)
                self.alignment = {}
        else:
            self.data['MISALIGNED'] = False

    def count_contacts(self):

        model_cmap = self._prepare_contactmap(self.model.copy())
        model_dict = model_cmap.as_dict()
        self.data['CONTACTS'] = self.data['RESNUM'].apply(lambda x: len(model_dict[int(x)]))

    def add_plddt(self,externally_supplied_plddts = {}):

        if externally_supplied_plddts == {}:
            prediction_end = len(self.prediction.plddt)
            self.data['PLDDT'] = self.data['RESNUM'].apply(lambda x: self.prediction.plddt[int(x)] if int(x) < prediction_end else 0)
        else:
            # TODO: implement using externally_supplied_plddts once the format is decided
            logger.warning("externally_supplied_plddts was provided but external pLDDT support is not yet implemented; pLDDT filter will be skipped.")

    def Run_gesamt_filter(self, experimentfile, predictionfile, gesamt_exe, moltype='Protein', experimentfiletype='pdb'):

        suggested_correspondece = self.alignment.copy()

        map_align_raw = self.data['MISALIGNED']
        svm_raw = self.data['SCORE']
        resnums_raw = self.data['RESNUM']
        absent_residues = self.absent_residues

        seen = set()
        resnums = []
        svm = []
        map_align = []

        for r, s, m in zip(resnums_raw, svm_raw, map_align_raw):
            if (not r in seen) and (not r in absent_residues):
                seen.add(r)
                resnums.append(r)
                svm.append(s)
                map_align.append(m)
                if not r in suggested_correspondece.keys():
                    suggested_correspondece[r]=r

        # identify potential errors
        flagged_experiment_regions, flagged_predcition_regions = tools.get_error_borders(svm, map_align, resnums, suggested_correspondece)

        # run gesamt for every region
        if experimentfiletype == 'pdb':
            p = PDBParser()
        if experimentfiletype == 'mmcif':
            p = MMCIFParser()
        model = p.get_structure('structure', experimentfile)[0]  ## hot fix to get chain name needed for gesamt while we are running single chain only this block needs fixing with external chain selection when multi chain handeling is introduced
        for chain in model:
            chain_experiment = chain.get_id()

        for r_exp, r_pred in zip(flagged_experiment_regions, flagged_predcition_regions):

            Q_region = tools.Gesamt_Q_score(predictionfile, r_pred, experimentfile, r_exp, gesamt_exe = gesamt_exe, chain_experiment = chain_experiment, chain_prediction = 'A', moltype=moltype)
            self.data.loc[ (self.data['RESNUM'] <= r_exp[1]) & (self.data['RESNUM'] >= r_exp[0]), 'Q_IN_ERROR'] = Q_region
        

        return 0

    def _calculate_filter_features(self):

        self.data['PLDDT_Smooth'] = tools.convolution_smooth_values(self.data['PLDDT'],window=5)
        self.data['CONTACTS_Smooth'] = tools.convolution_smooth_values(self.data['CONTACTS'],window=5)

        self.data['PLDDT_Diff'] = tools.convolution_diff_values(self.data['PLDDT'],window=3)
        self.data['CONTACTS_Diff'] = tools.convolution_diff_values(self.data['CONTACTS'],window=3)

        return 0

    def _apply_filter(self,filter_type = 'CMO'):

        if filter_type == 'CMO':
            Filter_Feature_names = ['CONTACTS','CONTACTS_Smooth','CONTACTS_Diff','PLDDT','PLDDT_Smooth','PLDDT_Diff','Q_IN_ERROR']

        elif filter_type == 'RF':
            Filter_Feature_names = ['CONTACTS','CONTACTS_Smooth','CONTACTS_Diff','PLDDT','PLDDT_Smooth','PLDDT_Diff','Q_IN_ERROR','SCORE']
        
        else:
            logger.warning("Unknown filter type %r requested; doing nothing.", filter_type)
            return 1

        _filter, _scaler = load_filter(filter_type)

        self.data[f'{filter_type}_FILTER'] = self.data['RESNUM'].apply(lambda x: self._calculate_filter_score(x,_filter,_scaler,Filter_Feature_names))


    def _calculate_filter_score(self, resnum, _filter, _scaler, _features):
        """calculate the filtering score for a resdue using the provided filter, scalar and features"""
        residue_features = self.data.loc[self.data.RESNUM == resnum][_features]

        if (self.absent_residues and resnum in self.absent_residues) or residue_features.isnull().values.any():
            return np.nan
        scaled_features = _scaler.transform(residue_features.values)

        return _filter.predict(scaled_features)[0]


    def Run_combined_filter(self, filter_type = 'CMO', filter_th=0.5):

        self.filter_threshold[filter_type] = filter_th

        if not filter_type in ['CMO','RF']:
            logger.warning("Unknown filter type %r requested; doing nothing.", filter_type)
            return 1

        if not {'CONTACTS','PLDDT','Q_IN_ERROR'}.issubset(set(self.data.columns)):
            logger.warning("Combined filtering requested but not all required features are available; setting all %s filter values to 0.", filter_type)
            self.data[f'{filter_type}_FILTER'] = 0
            return 1

        if filter_type=='RF' and 'SCORE' not in self.data.columns:
            logger.warning("Combined RF filter requested but no classifier scores found; run svm() first. Setting all RF filter values to 0.")
            self.data[f'{filter_type}_FILTER'] = 0
            return 1
        
        self._calculate_filter_features()

        self._apply_filter(filter_type)

        self.data[f'PASSED_{filter_type}_FILTER'] = self.data[f'{filter_type}_FILTER'].apply(lambda x: x >= filter_th)

        return 0

    def draw(self, RUN_SVM=True, RUN_MAP_ALIGN=True, RUN_FILTERS=True, n_contacts_per_res=2, plddt_threshold=65, svm_threshold=0.5, moltype='Protein'):

        misaligned_residues = set(self.alignment.keys())
        residues = self.data['RESNUM']
        absent = self.absent_residues
        present_residues = [r for r in sorted(residues) if r not in absent]

        # Pre-check whether filter columns exist so we can choose bar style before drawing.
        has_rf_filter = RUN_FILTERS and RUN_SVM and ('PASSED_RF_FILTER' in self.data.columns)
        has_cmo_filter = RUN_FILTERS and RUN_MAP_ALIGN and ('PASSED_CMO_FILTER' in self.data.columns)
        combine_filters = has_rf_filter or has_cmo_filter

        # Each bar row is a filled band of height _BAR_HEIGHT, separated by _BAR_GAP.
        # bar_top tracks the top edge of the next row to be drawn.
        # Step sized so a 6pt label inside the bar has comfortable clearance.
        _BAR_TOP = -0.01
        _BAR_HEIGHT = 0.034
        _BAR_GAP = 0.014
        _BAR_STEP = _BAR_HEIGHT + _BAR_GAP
        bar_top = _BAR_TOP

        classifier_label = 'RF' if moltype == 'RNA' else 'SVM'
        _blend = mtransforms.blended_transform_factory(self.ax.transAxes, self.ax.transData)

        def _color_bar(resnums, color, hatch=None):
            """Draw a single-color (optionally hatched) segment of the current bar row."""
            if resnums:
                kwargs = dict(width=1.0, color=color, linewidth=0, align='center')
                if hatch:
                    kwargs['hatch'] = hatch
                    kwargs['edgecolor'] = 'white'
                self.ax.bar(resnums, _BAR_HEIGHT, bottom=bar_top - _BAR_HEIGHT, **kwargs)

        def _label_bar(text):
            """Place a row label just inside the left edge of the axes, over the bars."""
            self.ax.text(0.005, bar_top - _BAR_HEIGHT / 2, text,
                         transform=_blend, ha='left', va='center',
                         fontsize=6, color='#333333', clip_on=False,
                         bbox=dict(boxstyle='square,pad=0.15', fc='white', alpha=0.55, ec='none'))

        if RUN_SVM:
            scores = self.data.set_index('RESNUM')['SCORE'].to_dict()
            called_errors = self.data.set_index('RESNUM')['SVM_CALLED_ERROR'].to_dict()

            # Build a full residue range so absent/missing positions become NaN,
            # which causes matplotlib to draw natural line breaks at gaps.
            full_range = list(range(min(scores), max(scores) + 1))
            score_values = np.array([
                np.nan if (r in absent or r not in scores) else scores[r]
                for r in full_range
            ])
            # Smooth using zeros at gap positions, but mask gaps back to NaN for display.
            smooth_input = np.where(np.isnan(score_values), 0.0, score_values)
            smooth_scores = tools.convolution_smooth_values(smooth_input)
            smooth_display = np.where(np.isnan(score_values), np.nan, smooth_scores)
            self.sorted_scores = smooth_input
            self.smooth_scores = smooth_display
            self.ax.plot(full_range, smooth_display, color=tools.ColorDefinitions.SCORE)

            correct_rns = [r for r in present_residues if not called_errors.get(r, False)]
            error_rns   = [r for r in present_residues if called_errors.get(r, False)]

            _color_bar(correct_rns, tools.ColorDefinitions.CORRECT)
            if has_rf_filter:
                filt = self.data.set_index('RESNUM')['PASSED_RF_FILTER'].to_dict()
                _color_bar([r for r in error_rns if filt.get(r, False)],
                           tools.ColorDefinitions.ERROR)
                _color_bar([r for r in error_rns if not filt.get(r, False)],
                           tools.ColorDefinitions.ERROR, hatch='///')
            else:
                _color_bar(error_rns, tools.ColorDefinitions.ERROR)
            _label_bar(classifier_label)
            bar_top -= _BAR_STEP

        if RUN_MAP_ALIGN and self._map_align_ran:
            aligned_rns    = [r for r in present_residues if r not in misaligned_residues]
            misaligned_rns = [r for r in present_residues if r in misaligned_residues]

            _color_bar(aligned_rns, tools.ColorDefinitions.ALIGNED)
            if has_cmo_filter:
                filt = self.data.set_index('RESNUM')['PASSED_CMO_FILTER'].to_dict()
                _color_bar([r for r in misaligned_rns if filt.get(r, False)],
                           tools.ColorDefinitions.MISALIGNED)
                _color_bar([r for r in misaligned_rns if not filt.get(r, False)],
                           tools.ColorDefinitions.MISALIGNED, hatch='///')
            else:
                _color_bar(misaligned_rns, tools.ColorDefinitions.MISALIGNED)
            _label_bar('CMO')
            bar_top -= _BAR_STEP

        # Individual filter bars shown as fallback when combined filter data is unavailable.
        if RUN_FILTERS and not combine_filters:
            if 'CONTACTS' in self.data.columns:
                n_contacts = self.data.set_index('RESNUM')['CONTACTS'].to_dict()
                _color_bar([r for r in present_residues if n_contacts.get(r, 0) < n_contacts_per_res],
                           tools.ColorDefinitions.LOW_CONTACTS)
                _color_bar([r for r in present_residues if n_contacts.get(r, 0) >= n_contacts_per_res],
                           tools.ColorDefinitions.SUFFICIENT_CONTACTS)
                _label_bar('Contacts')
                bar_top -= _BAR_STEP

            if 'PLDDT' in self.data.columns:
                plddts = self.data.set_index('RESNUM')['PLDDT'].to_dict()
                color_scheme = tools.ColorDefinitions.PLDDT_COLORS
                thresholds = sorted(color_scheme.keys(), reverse=True)
                plddt_groups = {}
                for r in present_residues:
                    color = color_scheme[thresholds[0]]
                    for th in thresholds:
                        if plddts[r] < th:
                            color = color_scheme[th]
                    plddt_groups.setdefault(color, []).append(r)
                for color, rns in plddt_groups.items():
                    _color_bar(rns, color)
                _label_bar('pLDDT')
                bar_top -= _BAR_STEP

            if 'Q_IN_ERROR' in self.data.columns:
                Qs = self.data.set_index('RESNUM')['Q_IN_ERROR'].to_dict()
                color_scheme = tools.ColorDefinitions.Q_COLORS
                thresholds = sorted(color_scheme.keys(), reverse=True)
                q_groups = {}
                for r in present_residues:
                    if Qs.get(r, '') == '':
                        continue
                    color = color_scheme[thresholds[0]]
                    for th in thresholds:
                        if Qs[r] < th:
                            color = color_scheme[th]
                    q_groups.setdefault(color, []).append(r)
                for color, rns in q_groups.items():
                    _color_bar(rns, color)
                _label_bar('Q-score')
                bar_top -= _BAR_STEP

        # Add 10 % headroom above the actual score peak so the curve doesn't hug the top edge.
        # Guard against nan/inf from an all-absent sequence (np.nanmax on all-NaN → nan).
        if RUN_SVM:
            peak = float(np.nanmax(self.smooth_scores))
            score_top = (peak * 1.10) if np.isfinite(peak) else 1.1
        else:
            score_top = None
        self.ax.set_ylim(bottom=bar_top, top=score_top)
        self.ax.axhline(svm_threshold, **LINEKWARGS)
        self.ax.set_xlabel('Residue Number')
        self.ax.set_ylabel('Smoothed score')

        if self.legend:
            self._add_legend(RUN_SVM=RUN_SVM, RUN_MAP_ALIGN=RUN_MAP_ALIGN, RUN_FILTERS=RUN_FILTERS,
                             n_contacts_per_res=n_contacts_per_res, plddt_threshold=plddt_threshold,
                             has_rf_filter=has_rf_filter, has_cmo_filter=has_cmo_filter,
                             combine_filters=combine_filters)

        # TODO: deprecate this in 0.14
        if self._file_name:
            self.savefig(self._file_name, dpi=self._dpi)
