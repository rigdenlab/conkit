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
"""Extensions to conkit.core package"""

__author__ = "Felix Simkovic"
__date__ = "03 May 2018"
__version__ = "1.0"

import numba
import numpy as np


@numba.jit(nopython=True, parallel=True)
def nb_calculate_weights(X, identity):
    N_float = float(X.shape[1])
    threshold = 1.0 - identity
    hamming = np.zeros(X.shape[0], dtype=np.float64)
    for i in numba.prange(X.shape[0]):
        for j in range(X.shape[0]):
            dist = 0
            for k in range(X.shape[1]):
                dist += X[i, k] != X[j, k]
            hamming[i] += (dist / N_float) < threshold
        hamming[i] = 1. / hamming[i]
    return hamming


@numba.jit(nopython=True, parallel=True)
def nb_filter(X, min_id, max_id):
    throwables = np.zeros(X.shape[0], dtype=np.uint8)
    for i in range(X.shape[0]):
        for j in range(i+1, X.shape[0]):
            if not throwables[j]:
                dist = 0.0
                for k in range(X.shape[1]):
                    dist += X[i, k] != X[j, k]
                ident = 1.0 - dist / X.shape[1] 
                throwables[j] = min_id > ident or ident > max_id
    return throwables


@numba.jit(nopython=True, parallel=True)
def nb_filter_gapped(X, symbol, min_prop, max_prop):
    throwables = np.zeros(X.shape[0], dtype=np.uint8)
    for i in range(X.shape[0]):
        prop = 0.0
        for k in range(X[i].shape[0]):
            prop += X[i, k] == symbol
        prop = prop / X[i].shape[0]
        throwables[i] = prop < min_prop or prop > max_prop
    return throwables
