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
__date__ = "02 Aug 2017"
__version__ = "1.0"

import abc
import numpy as np

ABC = abc.ABCMeta('ABC', (object,), {})


class BandwidthBase(ABC):
    """Abstract class for bandwidth calculations"""

    @abc.abstractproperty
    def bandwidth(self):
        return 0.0

    # REM: abstractproperty requires us to re-declare it in every child class
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
        y0 = AmiseBW.optimal_bandwidth_equation(data, x0)
        x = 0.8 * x0
        y = AmiseBW.optimal_bandwidth_equation(data, x)
        for _ in np.arange(self._niterations):
            x -= y * (x0 - x) / (y0 - y)
            y = AmiseBW.optimal_bandwidth_equation(data, x)
            if abs(y) < (self._eps * y0):
                break
        return x

    @staticmethod
    def curvature(p, x, w):
        z = (x - p) / w
        y = (1 * (z ** 2 - 1.0) * np.exp(-0.5 * z * z) / (w * np.sqrt(2. * np.pi)) / w ** 2).sum()
        return y / p.shape[0]

    @staticmethod
    def extended_range(mn, mx, bw, ext=3):
        return mn - ext * bw, mx + ext * bw

    @staticmethod
    def optimal_bandwidth_equation(p, default_bw):
        alpha = 1. / (2. * np.sqrt(np.pi))
        sigma = 1.0
        n = p.shape[0]
        q = AmiseBW.stiffness_integral(p, default_bw)
        return default_bw - ((n * q * sigma ** 4) / alpha) ** (-1.0 / (p.shape[1] + 4))

    @staticmethod
    def stiffness_integral(p, default_bw, eps=1e-4):
        mn, mx = AmiseBW.extended_range(p.min(), p.max(), default_bw, ext=3)
        n = 1
        dx = (mx - mn) / n
        yy = 0.5 * dx * (AmiseBW.curvature(p, mn, default_bw) ** 2 +
                         AmiseBW.curvature(p, mx, default_bw) ** 2)
        # The trapezoidal rule guarantees a relative error of better than eps
        # for some number of steps less than maxn.
        maxn = (mx - mn) / np.sqrt(eps)
        # Cap the total computation spent
        maxn = 2048 if maxn > 2048 else maxn
        n = 2
        while n <= maxn:
            dx /= 2.
            y = 0
            for i in np.arange(1, n, 2):
                y += AmiseBW.curvature(p, mn + i * dx, default_bw) ** 2
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
        sigma = np.sqrt((data ** 2).sum() / data.shape[0] - (data.sum() / data.shape[0]) ** 2)
        return sigma * ((((data.shape[1] + 2) * data.shape[0]) / 4.) ** (-1. / (data.shape[1] + 4)))


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
        sigma = np.minimum(np.std(data, axis=0, ddof=1), (np.percentile(data, 75) - np.percentile(data, 25)) / 1.349)[0]
        return 1.059 * sigma * data.shape[0] ** (-1. / (data.shape[1] + 4))


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
        sigma = np.minimum(np.std(data, axis=0, ddof=1), (np.percentile(data, 75) - np.percentile(data, 25)) / 1.349)[0]
        return 0.9 * sigma * (data.shape[0] * (data.shape[1] + 2) / 4.) ** (-1. / (data.shape[1] + 4))


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
