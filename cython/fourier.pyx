# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False
# distutils: extra_compile_args = -fopenmp -ffast-math -funroll-loops
# distutils: extra_link_args = -fopenmp

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Time Series Analysis -- Cython module
#
# Author: Jakob Rørsted Mosumgaard
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


###############################################################################
# Modules
###############################################################################
# General
from __future__ import print_function, with_statement, division
import numpy as np

# Cython
from cython cimport boundscheck, wraparound, nonecheck

# External handwritten C-module
cdef extern from "tsfourier.h" nogil:
    void fourier(double[] time, double[] flux, double[] ny, size_t N, size_t M,
                 double[] power)


###############################################################################
# Auxiliary functions for Cython arrays (using typed memoryviews)
###############################################################################
cdef void arr_sca(double[::1] mvi, double a, double[::1] mvo,
                 int N) nogil:
    """
    Multiply array by scalar.
    NOTE: Cannot interact with Python objets.
    """
    cdef:
        int i

    for i in range(N):
        mvo[i] = a * mvi[i]


cdef double[::1] arr_sca_copy(double[::1] mvi, double a):
    """
    Multiply array by scalar -- NOT IN-PLACE.
    NOTE: Can interact with Python objects.
    """
    cdef:
        double[::1] mvo = np.zeros(mvi.shape[0])
        int i, N

    N = mvi.shape[0]
    for i in range(N):
        mvo[i] = a * mvi[i]
    return mvo


###############################################################################
# Primary auxiliary function
###############################################################################
def powerspec(double[::1] time, double[::1] flux,
              double low, double high, double rate):
    """
    Calculate the power spectrum (fourier transform) using a least mean square
    method.

    This function is a wrapper for the underlying pure-C function.

    Arguments:
    - `time`: Array with the values of the time
    - `flux`: Array with the measured flux
    - `low` : The lowest test frequency
    - `high`: The highest test frequency
    - `rate`: The sampling rate (spacing between frequencies)
    """
    # Cython typing
    cdef:
        double PI2 = 2 * np.pi
        double[::1] freq = np.arange(low, high, rate)
        int N = time.shape[0]
        int M = freq.shape[0]
        double[::1] ny = np.zeros(M)
        double[::1] powers = np.zeros(M)
    
    # Convert test cyclic frequencies to angular
    arr_sca(freq, PI2, ny, M)

    # Calculate power using external C-function
    #  --> Temporary deactivate the Python memory management for performance
    with nogil:
        fourier(&time[0], &flux[0], &ny[0], N, M, &powers[0])

    # Return an array of test (cyclic) frequencies and the calculated power
    return freq, powers


###############################################################################
# Main function
###############################################################################
def calc(infile, freq_start, freq_stop, freq_rate):
    """
    Interface to calculation of the power spectrum.

    This function prepares the input for the fast Cython/C-function.

    Arguments:
    - `infile`: File to read in the format (t [seconds] , data).
    - `freq_start` : The lowest test frequency (in microHertz).
    - `freq_stop`: The highest test frequency (in microHertz).
    - `freq_rate`: The sampling rate (spacing between frequencies).
    """
    cdef:
        double[::1] t, time, flux, freq, powers
        double ms, low, high, rate
    
    # Set up of power spectrum
    low = freq_start
    high = freq_stop
    rate = freq_rate

    # Load data
    t, flux = np.ascontiguousarray(np.loadtxt(infile, unpack=True))

    # Convert time to megaseconds
    ms = 1e-6
    time = arr_sca_copy(t, ms)

    # Run power spectrum
    freq, powers = powerspec(time, flux, low, high, rate)

    # Return
    return freq, powers
