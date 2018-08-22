#cython: boundscheck=False, cdivision=True, wraparound=False

cimport cython
import numpy as np
cimport numpy as np

from cython.parallel import prange

np.import_array()


def c_get_frequency(np.ndarray[np.int64_t, ndim=2] X, Py_ssize_t symbol, np.ndarray[np.int64_t, ndim=1] frequencies):
    cdef Py_ssize_t i, j
    for j in prange(X.shape[1], nogil=True): 
        for i in xrange(X.shape[0]):
            frequencies[j] += X[i, j] == symbol


def c_get_weights(np.ndarray[np.int64_t, ndim=2] X, double identity, np.ndarray[double, ndim=1] hamming):
    cdef Py_ssize_t i, j, k
    cdef double threshold, dist
    threshold = (1.0 - identity) * X.shape[1]
    for i in prange(X.shape[0], nogil=True):
        for j in xrange(X.shape[0]):
            dist = 0.0
            for k in xrange(X.shape[1]):
                dist += X[i, k] != X[j, k]
            hamming[i] += dist < threshold
        hamming[i] = 1.0 / hamming[i]


def c_filter(np.ndarray[np.int64_t, ndim=2] X, double min_id, double max_id, np.ndarray[np.uint8_t, ndim=1, cast=True] throwables):
    cdef Py_ssize_t i, j, k
    cdef double dist
    for i in xrange(X.shape[0]):
        for j in xrange(i + 1, X.shape[0]):
            if not throwables[j]:
                dist = 0
                for k in xrange(X.shape[1]):
                    dist += X[i, k] != X[j, k]
                ident = 1.0 - dist / X.shape[1]
                throwables[j] = (ident < min_id) or (ident > max_id)


def c_filter_symbol(np.ndarray[np.int64_t, ndim=2] X, double min_prop, double max_prop, Py_ssize_t symbol, np.ndarray[np.uint8_t, ndim=1, cast=True] throwables):
    cdef Py_ssize_t i, k
    cdef double prop
    for i in xrange(X.shape[0]):
        prop = 0
        for k in xrange(X.shape[1]):
            prop += X[i, k] == symbol
        prop /= X.shape[1]
        throwables[i] = (prop < min_prop) or (prop > max_prop)
