#cython: boundscheck=False, cdivision=True, wraparound=False

cimport cython
import numpy as np
cimport numpy as np

from libc.math cimport exp, fabs, sqrt, M_PI

np.import_array()

cdef double SQRT_PI = sqrt(M_PI)
cdef double SQRT_2PI = sqrt(2.0 * M_PI)


def c_optimize_bandwidth(np.ndarray[np.int64_t, ndim=2] A, double v):
    cdef double alpha, sigma, integral
    alpha = 1.0 / (2.0 * SQRT_PI)
    sigma = 1.0
    integral = c_get_stiffness_integral(A, v, 0.0001)
    return v - ((A.shape[0] * integral * sigma**4) / alpha)**(-1.0 / (A.shape[1] + 4))


def c_get_stiffness_integral(np.ndarray[np.int64_t, ndim=2] A, double v, double eps):
    cdef Py_ssize_t i, j, n
    cdef double min_, max_, dx, maxn, yy, y1, y2, y3, y
    min_ = A.min() - v * 3
    max_ = A.max() + v * 3
    dx = 1.0 * (max_ - min_)
    maxn = dx / sqrt(eps)
    if maxn > 2048:
        maxn = 2048
    y1 = c_get_gauss_curvature(A, min_, v)
    y2 = c_get_gauss_curvature(A, max_, v)
    yy = 0.5 * dx * (y1 * y1 + y2 * y2)
    n = 2
    
    while n <= maxn:
        dx /= 2.0
        y = 0.0
        for i in xrange(1, n, 2):
            y3 = c_get_gauss_curvature(A, min_ + i * dx, v)
            y += (y3 * y3)
        yy = 0.5 * yy + y * dx
        if n > 8 and fabs(y * dx - 0.5 * yy) < eps * yy:
            break
        n *= 2

    return yy


def c_get_gauss_curvature(np.ndarray[np.int64_t, ndim=2] A, double x, double w):
    cdef Py_ssize_t i, j
    cdef double w_sq, w_sqrt_2pi, curvature, z
    w_sq = w*w
    w_sqrt_2pi = w * SQRT_2PI
    curvature = 0.0
    for i in xrange(A.shape[0]):
        for j in xrange(A.shape[1]):
            z = (x - A[i, j]) / w
            z *= z
            curvature += (A.shape[1] * (z - 1.0) * (exp(-0.5 * z) / w_sqrt_2pi) / w_sq)
    return curvature / A.shape[0]
