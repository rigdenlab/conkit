# BSD 3-Clause License
#
# Copyright (c) 2016-18, University of Liverpool
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

"""A collection of bandwidth estimators for Kernel Density Estimation"""

from __future__ import division
from __future__ import print_function

__author__ = "Felix Simkovic"
__date__ = "31 Mar 2018"
__version__ = "1.0.1"

import abc
import numpy as np

ABC = abc.ABCMeta('ABC', (object,), {})

SQRT_PI = np.sqrt(np.pi)
SQRT_2PI = np.sqrt(2. * np.pi)


class BandwidthBase(ABC):
    """Abstract class for bandwidth calculations"""

    @abc.abstractproperty
    def bandwidth(self):
        return 0.0

    @property
    def bw(self):
        return self.bandwidth


class AmiseBW(BandwidthBase):
    """Asymptotic Mean Integrated Squared Error (AMISE)

    This particular choice of bandwidth recovers all the important features whilst maintaining smoothness.
    It is a direct implementation of the method used by [#]_.

    .. [#] Sadowski, M.I. (2013). Prediction of protein domain boundaries from inverse covariances.

    """
    def __init__(self, data, niterations=25, eps=1e-3):
        """Instantiate a new bandwith calculator"""
        self._data = np.asarray(data)
        self._niterations = niterations
        self._eps = eps

    @property
    def bandwidth(self):
        data = self._data
        x0 = BowmanBW(data).bandwidth
        y0 = AmiseBW.optimal_bandwidth(data, x0)
        x = 0.8 * x0
        y = AmiseBW.optimal_bandwidth(data, x)
        for _ in np.arange(self._niterations):
            x -= y * (x0 - x) / (y0 - y)
            y = AmiseBW.optimal_bandwidth(data, x)
            if abs(y) < (self._eps * y0):
                break
        return x

    @staticmethod
    def curvature(data, x, w):
        """
        See Also
        --------
        gauss_curvature
        """
        import warnings
        warnings.warn("This function will be removed in a future release! Use gauss_curvature() instead")
        return AmiseBW.gauss_curvature(data, x, w)

    @staticmethod
    def gauss_curvature(data, x, w):
        M, N = data.shape
        z = ((x - data) / w)**2
        return (N * (z - 1.0) * (np.exp(-0.5 * z) / (w * SQRT_2PI)) / w**2).sum() / M

    @staticmethod
    def extended_range(min_, max_, bandwidth, ext=3):
        d = bandwidth * ext
        return min_ - d, max_ + d 

    @staticmethod
    def optimal_bandwidth(data, bandwidth):
        M, N = data.shape
        alpha = 1. / (2. * SQRT_PI)
        sigma = 1.0
        integral = AmiseBW.stiffness_integral(data, bandwidth)
        return bandwidth - ((M * integral * sigma ** 4) / alpha) ** (-1.0 / (N + 4))

    @staticmethod
    def stiffness_integral(data, bandwidth, eps=1e-4):
        min_, max_ = AmiseBW.extended_range(data.min(), data.max(), bandwidth, ext=3)
        dx = 1.0 * (max_ - min_)
        maxn = dx / np.sqrt(eps)
        if maxn > 2048:
            maxn = 2048
        yy = 0.5 * dx * (
            AmiseBW.gauss_curvature(data, min_, bandwidth)**2 
            + AmiseBW.gauss_curvature(data, max_, bandwidth)**2
        )
        n = 2
        while n <= maxn:
            dx /= 2.
            y = 0
            for i in np.arange(1, n, 2):
                y += AmiseBW.gauss_curvature(data, min_ + i * dx, bandwidth) ** 2
            yy = 0.5 * yy + y * dx
            if n > 8 and abs(y * dx - 0.5 * yy) < eps * yy:
                break
            n *= 2
        return yy


class BowmanBW(BandwidthBase):
    """Bowman & Azzalini [#]_ bandwidth calculation

    .. math::

       \\sqrt{\\frac{\\sum{X}^2}{n}-(\\frac{\\sum{X}}{n})^2}*(\\frac{(d+2)*n}{4})^\\frac{-1}{d+4}

    .. [#] Bowman, A.W. & Azzalini, A. (1997). Applied Smoothing Techniques for Data Analysis.

    """
    def __init__(self, data):
        """Instantiate a new bandwith calculator"""
        self._data = np.asarray(data)

    @property
    def bandwidth(self):
        data = self._data
        M, N = data.shape
        return np.sqrt((data ** 2).sum() / M - (data.sum() / M) ** 2) * ((((N + 2) * M) / 4.) ** (-1. / (N + 4)))


class LinearBW(BandwidthBase):
    """Linear [#]_ implementation

    .. math::

       \\frac{N_{max}}{t}

    .. [#] Sadowski, M.I. (2013). Prediction of protein domain boundaries from inverse covariances.

    """
    def __init__(self, data, threshold=15):
        self._data = np.asarray(data)
        self._threshold = threshold

    @property
    def bandwidth(self):
        return float(self._data.max() / self._threshold)


class ScottBW(BandwidthBase):
    """Scott's [#]_ implementation

    .. math::

       1.059*\\sigma*n^\\frac{-1}{d+4}

    .. [#] Scott, D.W. (1992). Multivariate Density Estimation: Theory, Practice, and Visualization.

    """
    def __init__(self, data):
        """Instantiate a new bandwith calculator"""
        self._data = np.asarray(data)

    @property
    def bandwidth(self):
        data = self._data
        M, N = data.shape
        sigma = np.minimum(np.std(data, axis=0, ddof=1), (np.percentile(data, 75) - np.percentile(data, 25)) / 1.349)[0]
        return 1.059 * sigma * M ** (-1. / (N + 4))


class SilvermanBW(BandwidthBase):
    """Silverman's [#]_ implementation

    .. math::

       0.9*\\sigma*(n*\\frac{d+2}{4})^\\frac{-1}{d+4}

    .. [#] Silverman, B.W. (1986). Density Estimation for Statistics and Data Analysis.

    """
    def __init__(self, data):
        """Instantiate a new bandwith calculator"""
        self._data = np.asarray(data)

    @property
    def bandwidth(self):
        data = self._data
        M, N = data.shape
        sigma = np.minimum(np.std(data, axis=0, ddof=1), (np.percentile(data, 75) - np.percentile(data, 25)) / 1.349)[0]
        return 0.9 * sigma * (M * (N + 2) / 4.) ** (-1. / (N + 4))


def bandwidth_factory(method):
    """Obtain the bandwidth as defined by user method"""
    if method == "amise":
        return AmiseBW
    elif method == "bowman":
        return BowmanBW
    elif method == "linear":
        return LinearBW
    elif method == "scott":
        return ScottBW
    elif method == "silverman":
        return SilvermanBW
    else:
        msg = "Undefined bandwidth method: {0}".format(method)
        raise ValueError(msg)
