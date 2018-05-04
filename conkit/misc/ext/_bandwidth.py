# coding=utf-8
#
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
"""Extension to conkit.core.contactmap module"""

__author__ = "Felix Simkovic"
__date__ = "04 May 2018"
__version__ = "1.0"

import numba
import numpy as np

SQRT_PI = np.sqrt(np.pi)
SQRT_2PI = np.sqrt(2. * np.pi)


@numba.jit(nopython=True, parallel=True)
def optimize_bandwidth(A, v):
    alpha = 1. / (2. * SQRT_PI)
    sigma = 1.0
    integral = stiffness_integral(A, v, 0.0001)
    return v - ((A.shape[0] * integral * sigma**4) / alpha)**(-1.0 / (A.shape[1] + 4))


@numba.jit(nopython=True, parallel=True)
def stiffness_integral(A, v, eps):
    min_ = A.min() - v * 3
    max_ = A.max() + v * 3
    dx = 1.0 * (max_ - min_)
    maxn = dx / np.sqrt(eps)
    if maxn > 2048:
        maxn = 2048
    yy = 0.5 * dx * (
        gauss_curvature(A, min_, v)**2 + gauss_curvature(A, max_, v)**2)
    n = 2
    while n <= maxn:
        dx /= 2.
        y = 0
        for i in range(1, n, 2):
            y += gauss_curvature(A, min_ + i * dx, v)**2
        yy = 0.5 * yy + y * dx
        if n > 8 and abs(y * dx - 0.5 * yy) < eps * yy:
            break
        n *= 2
    return yy


@numba.jit(nopython=True, parallel=True)
def gauss_curvature(A, x, w):
    w_sq = w**2
    w_sqrt_2pi = w * SQRT_2PI
    curvature = 0.0
    for i in range(A.shape[0]):
        for j in range(A.shape[1]):
            z = ((x - A[i, j]) / w)**2
            curvature += (A.shape[1] * (z - 1.0) * (np.exp(-0.5 * z) / w_sqrt_2pi) / w_sq)
    return curvature / A.shape[0]
