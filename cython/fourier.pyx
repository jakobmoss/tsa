# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False
# distutils: extra_compile_args = -fopenmp -ffast-math -funroll-loops
# distutils: extra_link_args = -fopenmp

###############################################################################
# Modules
###############################################################################
# General
from __future__ import print_function, with_statement, division
import numpy as np

# Cython
from cython cimport boundscheck, wraparound, nonecheck, cdivision
from cython.parallel cimport prange

cdef extern from "math.h" nogil:
    double sin(double x)
    double cos(double x)

###############################################################################
# Auxiliary functions for arrays (using typed memoryviews)
#  -->  NB: Performed in-place for speed
###############################################################################
cdef double arrsum(double[::1] mv) nogil:
    """Sum the contents of the argument"""
    cdef:
        double ss = 0.0
        int i, N

    N = mv.shape[0]
    for i in range(N):
        ss += mv[i]
    return ss


cdef double arrsqsum(double[::1] mv) nogil:
    """Square the elements of the argument and sum"""
    cdef:
        double ss = 0.0
        int i, N

    N = mv.shape[0]
    for i in range(N):
        ss += mv[i] * mv[i]
    return ss


cdef void arrcos(double[::1] mvi, double[::1] mvo, int N) nogil:
    """Take cos of all elements in array"""
    cdef:
        int i

    for i in range(N):
        mvo[i] = cos(mvi[i])


cdef void arrsin(double[::1] mvi, double[::1] mvo, int N) nogil:
    """Take cos of all elements in array"""
    cdef:
        int i

    for i in range(N):
        mvo[i] = sin(mvi[i])


cdef void arrmult(double[::1] mvi1, double[::1] mvi2, double[::1] mvo,
                         int N) nogil:
    """Multiply two arrays"""
    cdef:
        int i

    for i in range(N):
        mvo[i] = mvi1[i] * mvi2[i]


cdef void arrsca(double[::1] mvi, double a, double[::1] mvo,
                 int N) nogil:
    """Multiply array by scalar"""
    cdef:
        int i

    for i in range(N):
        mvo[i] = a * mvi[i]


cdef double[::1] arrsca_safe(double[::1] mvi, double a):
    """Multiply array by scalar -- NOT IN-PLACE"""
    cdef:
        double[::1] mvo = np.zeros(mvi.shape[0])
        int i, N

    N = mvi.shape[0]
    for i in range(N):
        mvo[i] = a * mvi[i]
    return mvo


cdef void arraddsq(double[::1] mvi1, double[::1] mvi2, double[::1] mvo,
                         int N) nogil:
    """Element-by-element summation of square of two arrays"""
    cdef:
        int i

    for i in range(N):
        mvo[i] = mvi1[i]*mvi1[i] + mvi2[i]*mvi2[i]


###############################################################################
# Primary auxiliary function
###############################################################################
@cdivision(True)
def powerspec(double[::1] time, double[::1] flux,
              double low, double high, double rate):
    """
    Calculate the fourier power spectrum using a least mean square method.
    Arguments:
    - `time`: Array with the values of the time
    - `flux`: Array with the measured flux
    - `low` : The lowest test frequency
    - `high`: The highest test frequency
    - `rate`: The sampling rate (spacing between frequencies)
    """
    # Cython typing
    cdef:
        int i
        double s, c, ss, cc, sc
        double PI2 = 2 * np.pi
        double[::1] freq = np.arange(low, high, rate)
        int N = freq.shape[0]
        int M = time.shape[0]
        double[::1] ny = np.zeros(N)
        double[::1] alpha = np.zeros(N)
        double[::1] beta = np.zeros(N)
        double[::1] powers = np.zeros(N)
        double[::1] pcos = np.zeros(M)
        double[::1] psin = np.zeros(M)
        double[::1] tny = np.zeros(M)
        double[::1] v1 = np.zeros(M)
        double[::1] v2 = np.zeros(M)
        double[::1] v3 = np.zeros(M)
    
    # Convert test cyclic frequencies to angular
    arrsca(freq, PI2, ny, N)
    
    # The loop over frequencies (least mean square)
#    for i in prange(N, nogil=True):
    for i in range(N):
        # Pre-calculate sin, cos arrays
        arrsca(time, ny[i], tny, M)
        arrcos(tny, pcos, M)
        arrsin(tny, psin, M)

        # Calculate summed cos, sin terms
        arrmult(flux, pcos, v1, M)
        arrmult(flux, psin, v2, M)
        c = arrsum(v1)
        s = arrsum(v2)

        # Calculate squared terms
        cc = arrsqsum(pcos)
        ss = arrsqsum(psin)

        # Calculate cross-term
        arrmult(psin, pcos, v3, M)
        sc = arrsum(v3)

        # Calculate alpha and beta
        alpha[i] = (s*cc - c*sc) / (ss*cc - sc*sc)
        beta[i] = (c*ss - s*sc) / (ss*cc - sc*sc)

    # Calculate power
    arraddsq(alpha, beta, powers, N)
    
    # Return an array of test (cyclic) frequencies and the calculated power
    return freq, powers


###############################################################################
# Script
###############################################################################
def calc():
    cdef:
        double[::1] t, time, flux, freq, powers
        double ms, low, high, rate
    
    # Initial setup
    datdir = 'testdata/'
    outdir = 'output/'

    # Set up of power spectrum
    low = 1900.0
    high = 4100.0
    rate = 0.1

    # Load data
    infile = 'ts_7days.txt'
    t, flux = np.ascontiguousarray(np.loadtxt(datdir + infile, unpack=True))

    # Convert time to megaseconds
    ms = 1e-6
    time = arrsca_safe(t, ms)

    # Run power spectrum
    freq, powers = powerspec(time, flux, low, high, rate)

    # Return
    return freq, powers
