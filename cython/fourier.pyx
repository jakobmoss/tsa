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
def calc(infile, freq_start, freq_stop, freq_rate, unit='s', prep=True):
    """
    Calculate the power spectrum using a least mean square method. Returns
    arrays with test frequencies and corresponding power.
    The actual calculation is performed by a fast Cython/C function.

    Usage:
    frequencies, powers = calc( ... )

    Arguments:
    - `infile`: File to read in the format (t, data).
    - `freq_start` : The lowest test frequency (in microHertz).
    - `freq_stop`: The highest test frequency (in microHertz).
    - `freq_rate`: The sampling rate (spacing between frequencies).
    - `unit`: Unit of the time in data (allowed: 's' [default], 'day', 'ms').
    - `prep`: Subtract mean of data (default: True).
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
    t, f = np.ascontiguousarray(np.loadtxt(infile, unpack=True))

    # Convert to correct unit (seconds) and store in memoryview.
    unit = unit.lower()
    if (unit == 'day' or unit == 'days' or unit == 'd'):
        time = t * 86400.0
    elif (unit == 'ms' or unit == 'megasecond' or unit == 'megaseconds'):
        time = t * 1.0e6
    else:
        time = t

    # Subtract mean from data and store in memoryview.
    if prep:
        flux = f - np.mean(f)
    else:
        flux = f
        
    # Call the Cython-wrapper to fast calculation of the power spectrum
    freq, powers = _powerspec(time, flux, low, high, rate)

    # PRETTY PRINT
    print('Done!\n')
    
    # Return
    return freq, powers
