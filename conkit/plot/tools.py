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
"""Internal utility functions"""

__author__ = "Felix Simkovic"
__date__ = "16 Feb 2017"
__version__ = "0.13"

import numpy as np
import os

from conkit.core.contact import Contact
from conkit.core.contactmap import ContactMap
from conkit.core.contactfile import ContactFile
from conkit.core.sequence import Sequence
from conkit.core.sequencefile import SequenceFile
from conkit.core.distogram import Distogram
from conkit.core.distancefile import DistanceFile
from conkit.misc import deprecate

HierarchyIndex = {
    "Contact": Contact,
    "ContactMap": ContactMap,
    "ContactFile": ContactFile,
    "Sequence": Sequence,
    "SequenceFile": SequenceFile,
    "Distogram": Distogram,
    "DistanceFile": DistanceFile
}


class ColorDefinitions(object):
    """A class storing all color definitions for the various plots
    for fast and easy handling
    """

    GENERAL = "#000000"
    MATCH = "#0F0B2C"
    MISMATCH = "#DC4869"
    STRUCTURAL = "#D8D6D6"
    L5CUTOFF = "#3F4587"
    L20CUTOFF = "#B5DD2B"
    PRECISION50 = L5CUTOFF
    FACTOR1 = L20CUTOFF
    SCORE = '#3299a8'
    ERROR = '#f54242'
    CORRECT = '#40eef7'
    ALIGNED = '#3d8beb'
    MISALIGNED = '#f7ba40'
    AA_ENCODING = {
        "A": "#882D17",
        "C": "#F3C300",
        "D": "#875692",
        "E": "#F38400",
        "F": "#A1CAF1",
        "G": "#BE0032",
        "H": "#C2B280",
        "I": "#848482",
        "K": "#008856",
        "L": "#E68FAC",
        "M": "#0067A5",
        "N": "#F99379",
        "P": "#604E97",
        "Q": "#F6A600",
        "R": "#B3446C",
        "S": "#DCD300",
        "T": "#8DB600",
        "V": "#654522",
        "W": "#E25822",
        "Y": "#2B3D26",
        "X": "#000000",
    }


def find_minima(data, order=1):
    """Find the minima in a 1-D list

    Parameters
    ----------
    data : list, tuple
       A list of values
    order : int, optional
       The order, i.e. number of points next to point to consider

    Returns
    -------
    list
       A list of indices for minima

    Warning
    -------
    For multi-dimensional problems, see :func:`~scipy.signal.argrelmin`.

    Raises
    ------
    :exc:`ValueError`
       Order needs to be >= 1!
    :exc:`ValueError`
       More than two elements required!

    """
    if order < 1:
        raise ValueError("Order needs to be >= 1!")
    data = np.asarray(data)
    nelements = data.shape[0]
    if nelements < 2:
        raise ValueError("More than two elements required!")
    results = np.zeros(nelements, dtype=np.bool_)
    for i in np.arange(1, nelements - 1):
        start = 0 if i - order < 0 else i - order
        end = nelements if i + order + 1 > nelements else i + order + 1
        results[i] = np.all(data[start:i] > data[i]) and np.all(data[i] < data[i + 1: end])
    return np.where(results)[0].tolist()


def get_adjusted_aspect(ax, aspect_ratio):
    """Adjust the aspect ratio

    Parameters
    ----------
    ax : :obj:`~matplotlib.axes.Axes`
       A :obj:`~matplotlib.axes.Axes` instance
    aspect_ratio : float
       The desired aspect ratio for :obj:`~matplotlib.axes.Axes`

    Returns
    -------
    float
       The required aspect ratio to achieve the desired one

    Warning
    -------
    This function only works for non-logarithmic axes.

    """
    default_ratio = (ax.get_xlim()[1] - ax.get_xlim()[0]) / (ax.get_ylim()[1] - ax.get_ylim()[0])
    return float(default_ratio * aspect_ratio)


@deprecate("0.11", msg="Use get_points_on_circle instead")
def points_on_circle(*args, **kwargs):
    return get_points_on_circle(*args, **kwargs)


def get_points_on_circle(radius, h=0, k=0):
    """Calculate points on a circle with even spacing

    Parameters
    ----------
    radius : int
       The radius of the circle
    h : int, optional
       The x coordinate of the origin
    k : int, optional
       The y coordinate of the origin

    Returns
    -------
    list
       The list of coordinates for each point

    """
    if radius == 0:
        return [[]]
    else:
        space = 2 * np.pi / radius
        coords = np.zeros((radius, 2))
        for i in np.arange(radius):
            coords[i] = [round(h + radius * np.cos(space * i), 6), round(k + radius * np.sin(space * i), 6)]
        return coords.tolist()


def get_radius_around_circle(p1, p2):
    """Obtain the radius around a given circle

    Parameters
    ----------
    p1 : list, tuple
       Point 1
    p2 : list, tuple
       Point 2 adjacent `p1`

    Returns
    -------
    float
       The radius for points so p1 and p2 do not intersect
    """
    dist = np.linalg.norm(np.array(p1) - np.array(p2))
    return dist / 2.0 - dist * 0.1


def _isinstance(hierarchy, hierarchy_type):
    """Confirm the data structure to be a ConKit definition"""
    if isinstance(hierarchy_type, str) and hierarchy_type in HierarchyIndex:
        return isinstance(hierarchy, HierarchyIndex[hierarchy_type])
    else:
        return isinstance(hierarchy, hierarchy_type)


def is_executable(executable):
    """Check if a given program can be executed

    Parameters
    ----------
    executable : str
       The path or name for an executable

    Returns
    -------
    str
       The absolute path to the executable

    Raises
    ------
    ValueError
        The executable cannot be accessed

    Credits
    -------
    https://stackoverflow.com/a/377028/3046533
    """
    if executable is None:
        return

    fpath, fname = os.path.split(executable)

    if fpath:
        if os.path.isfile(executable) and os.access(executable, os.X_OK):
            return executable

    else:
        for path in os.environ["PATH"].split(os.pathsep):
            exe_file = os.path.join(path, executable)
            if os.path.isfile(exe_file) and os.access(exe_file, os.X_OK):
                return exe_file

    raise ValueError('Executable {} cannot be accessed'.format(executable))


def convolution_smooth_values(x, window=5):
    """Use convolutions to smooth a list of numeric values

    Parameters
    ----------
    x : list, tuple
       A list with the numeric values to be smoothed
    window : int
       The residue window to be used to smooth values [default: 5]

    Returns
    -------
    list
       A list with the smoothed numeric values
    """
    box = np.ones(window) / window
    x_smooth = np.convolve(x, box, mode='same')
    return x_smooth


def get_rmsd(distogram_1, distogram_2, calculate_wrmsd=True):
    """Calculate the RMSD between two different distograms

    Parameters
    ----------
    distogram_1 : :obj:`~conkit.core.distogram.Distogram`
       Distogram 1
    distogram_2 : :obj:`~conkit.core.distogram.Distogram`
       Distogram 2
    calculate_wrmsd: bool
        If True then the WRMSD is calculated using the confidence scores from distogram 1

    Returns
    -------
    tuple
       Two lists with the raw/smoothed RMSD values at each residue position
    """
    rmsd_raw = Distogram.calculate_rmsd(distogram_1, distogram_2, calculate_wrmsd=calculate_wrmsd)
    rmsd_smooth = convolution_smooth_values(np.nan_to_num(rmsd_raw), 10)
    return rmsd_raw, rmsd_smooth


def get_cmap_validation_metrics(model_cmap_dict, predicted_cmap_dict, sequence, absent_residues):
    """For a given observed contact map and predicted contact map calculate a series of validation metrics at each
    residue position (Accuracy, FN, FNR, FP, FPR, Sensitivity, Specificity)

    Parameters
    ----------
    model_cmap_dict : dict
       Dictionary representation of the contact map observed in the model
    predicted_cmap_dict : dict
       Dictionary representation of the predicted contact map
    sequence: :obj:`~conkit.core.sequence.Sequence`
        The sequence of the model to be validated
    absent_residues: list, tuple, set
        The residues that are missing from the model
    Returns
    -------
    tuple
        Two lists with the raw/smoothed values of each validation metric at each residue position
    """

    def accuracy(tp, fp, tn, fn):
        if (tp + fp + tn + fn) > 0:
            return (tp + tn) / (tp + fp + tn + fn)
        return 0

    def fn(*args):
        return args[-1]

    def fn_rate(tp, fp, tn, fn):
        if (fn + tn) > 0:
            return fn / (fn + tn)
        return 0

    def fp(*args):
        return args[1]

    def fp_rate(tp, fp, tn, fn):
        if (fp + tp) > 0:
            return fp / (fp + tp)
        return 0

    def sensitivity(tp, fp, tn, fn):
        if (fn + tp) > 0:
            return tp / (fn + tp)
        return 0

    def specificity(tp, fp, tn, fn):
        if (fp + tn) > 0:
            return tn / (fp + tn)
        return 0

    nresidues = len(sequence) - len(absent_residues)
    cmap_metrics = [[] for i in range(7)]
    metrics_list = (accuracy, fn, fn_rate, fp, fp_rate, sensitivity, specificity)

    for resnum in sorted(predicted_cmap_dict.keys()):
        if absent_residues and resnum in absent_residues:
            for metric in cmap_metrics:
                metric.append(np.nan)
            continue

        predicted_contact_set = {c for c in predicted_cmap_dict[resnum] if
                                 c[0] not in absent_residues and c[1] not in absent_residues}
        model_contact_set = {c for c in model_cmap_dict[resnum] if
                             c[0] not in absent_residues and c[1] not in absent_residues}

        _fn = len(predicted_contact_set - model_contact_set)
        _tp = len(predicted_contact_set & model_contact_set)
        _fp = len(model_contact_set - predicted_contact_set)
        _tn = nresidues - _fn - _tp - _fp

        for idx, metric in enumerate(metrics_list):
            cmap_metrics[idx].append(metric(_tp, _fp, _tn, _fn))

    smooth_cmap_metrics = []
    for metric in cmap_metrics:
        smooth_cmap_metrics.append(convolution_smooth_values(np.nan_to_num(metric), 5))

    return cmap_metrics, smooth_cmap_metrics


def get_zscores(model_distogram, predicted_cmap_dict, absent_residues, *metrics):
    """Calculate the Z-Scores for a series of metrics at each residue position
    using the population of residues within 10A

    Parameters
    ----------
    model_distogram : :obj:`~conkit.core.distogram.Distogram`
       Distogram of the model that will be validated
    predicted_cmap_dict : dict
       Dictionary representation of the predicted contact map
    absent_residues: list, tuple, set
        The residues that are missing from the model
    *metrics: list
        The mertics for which the Z-Scores will be calculated

    Returns
    -------
    list
       A list of lists where each sublist contains the Z-Scores for the input metrics across all the residues. The
       sublists containing the Z-Scores are ordered in the same original order as in the input *metrics
    """
    zscore_cmap_metrics = [[] for i in metrics]

    for resnum in sorted(predicted_cmap_dict.keys()):
        if absent_residues and resnum in absent_residues:
            for zscore_metric in zscore_cmap_metrics:
                zscore_metric.append(np.nan)
            continue

        neighbour_residues = model_distogram.find_residues_within(resnum, 10)
        for cmap_metric, zscore_metric in zip(metrics, zscore_cmap_metrics):
            population_scores = [cmap_metric[resid - 1] for resid in neighbour_residues]
            observed_score = cmap_metric[resnum - 1]
            zscore_metric.append(calculate_zscore(observed_score, population_scores))

    return zscore_cmap_metrics


def calculate_zscore(observed_score, population_scores):
    """Calculate the Z-Score for a given population of values. Z-Score = (score - mean) / stdev

    Parameters
    ----------
    observed_score : float
       The observed score used to calculate the Z-Score
    population_scores : list
       A list containing the scores observed across the samples in the population

    Returns
    -------
    float
       The calculated Z-Score
    """

    if len(population_scores) < 2:
        return 0
    stdev = np.std(population_scores).astype(float)
    if stdev == 0:
        return 0
    mean = np.mean(population_scores).astype(float)
    zscore = (observed_score - mean) / stdev
    return zscore


def get_residue_ranges(numbers):
    """Given a list of integers, creates a list of ranges with the consecutive numbers found in the list.

    Parameters
    ----------
    numbers: list
       A list of integers

    Returns
    ------
    list
        A list with the ranges of consecutive numbers found in the list
    """
    nums = sorted(set(numbers))
    gaps = [[s, e] for s, e in zip(nums, nums[1:]) if s + 3 < e]
    edges = iter(nums[:1] + sum(gaps, []) + nums[-1:])
    return list(zip(edges, edges))


def parse_map_align_stdout(stdout):
    """Parse the stdout of map_align and extract the alignment of residues.

    Parameters
    ----------
    stdout : str
       Standard output created with map_align

    Returns
    ------
    dict
        A dictionary where aligned residue numbers in map_b are the keys and residue numbers in map_a values. Only
        misaligned regions are included.
    """

    alignment_dict = {}

    for line in stdout.split('\n'):
        if line and line.split()[0] == "MAX":
            line = line.rstrip().lstrip().split()
            for residue_pair in line[8:]:
                residue_pair = residue_pair.split(":")
                if residue_pair[0] != residue_pair[1]:
                    alignment_dict[int(residue_pair[1])] = int(residue_pair[0])

    return alignment_dict
