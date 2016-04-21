# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: nonecheck=False
# cython: cdivision=True
# distutils: extra_compile_args = -fopenmp -ffast-math -funroll-loops
# distutils: extra_link_args = -fopenmp

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Time Series Analysis -- Cython module
#
# *** Available functions (callable by Python) ***
# - `calc`: Calculate power spectrum
#
# *** Usage ***
# Compile the module and then just `import fourier`.
#
# Author: Jakob Rørsted Mosumgaard
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


###############################################################################
# Modules
###############################################################################
# General
from __future__ import print_function, with_statement, division
import numpy as np

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
def _powerspec(double[::1] time, double[::1] flux,
               double low, double high, double rate):
    """
    Cython wrapper for calcuation of fourier transform using fast C code.

    Note: ALL variables are declared as static Cython types for proper
    interface with C!

    Arguments:
    - `time`: Array with the values of the time
    - `flux`: Array with the measured flux
    - `low` : The lowest test frequency
    - `high`: The highest test frequency
    - `rate`: The sampling rate (spacing between frequencies)
    """
    # Generate test cyclic frequencies
    cdef:
        double[::1] freq = np.arange(low, high, rate)

    # Arrays for data storage and their length
    cdef:
        double PI2 = 2 * np.pi
        int N = time.shape[0]
        int M = freq.shape[0]
        double[::1] powers = np.zeros(M)
    
    # Calculate power using external C-function
    #  --> Temporary deactivate the Python memory management for performance
    with nogil:
        fourier(&time[0], &flux[0], &freq[0], N, M, &powers[0])

    # Return an array of test (cyclic) frequencies and the calculated power
    return freq, powers


###############################################################################
# Main function
###############################################################################
def calc(infile, freq_start, freq_stop, freq_rate):
    """
    Calculate the power spectrum using a least mean square method. Returns
    arrays with test frequencies and corresponding power.
    The actual calculation is performed by a fast Cython/C function.

    Usage:
    frequencies, powers = calc( ... )

    Arguments:
    - `infile`: File to read in the format (t [seconds] , data).
    - `freq_start` : The lowest test frequency (in microHertz).
    - `freq_stop`: The highest test frequency (in microHertz).
    - `freq_rate`: The sampling rate (spacing between frequencies).
    """
    # Make Cython-typing
    cdef:
        double[::1] time, flux, freq, powers
        double ms, low, high, rate

    # PRETTY PRINT
    print('Calculating the power spectrum of \"{0}\" ...'.format(infile))
    
    # Set up of power spectrum -- convert input arguments to static types
    low = freq_start
    high = freq_stop
    rate = freq_rate

    # Load data into C-style memory block
    time, flux = np.ascontiguousarray(np.loadtxt(infile, unpack=True))

    # Call the Cython-wrapper to fast calculation of the power spectrum
    freq, powers = _powerspec(time, flux, low, high, rate)

    # PRETTY PRINT
    print('Done!\n')
    
    # Return
    return freq, powers
