#cython: boundscheck=False, wraparound=False

cimport cython
import numpy as np
cimport numpy as np

from libc.stdlib cimport abs

np.import_array()


def c_singletons(np.ndarray[np.int64_t, ndim=2] X, double threshold, np.ndarray[np.uint8_t, ndim=1, cast=True] throwables):
    cdef Py_ssize_t i, j
    for i in xrange(X.shape[0]):
        for j in xrange(i + 1, X.shape[0]):
            if throwables[j]:
                continue
            if abs(int(X[j, 0] - X[i, 0])) <= threshold and abs(int(X[j, 1] - X[i, 1])) <= threshold:
                throwables[i] = True 
                throwables[j] = True 
