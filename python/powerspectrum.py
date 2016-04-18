# -*- coding: utf-8 -*-

###############################################################################
# Modules
###############################################################################
from __future__ import print_function, with_statement, division
import numpy as np


###############################################################################
# Auxiliary functions
###############################################################################
def powerspec(time, flux, low, high, rate):
    """
    Calculate the fourier power spectrum using a least mean square method.
    Arguments:
    - `time`: Array with the values of the time
    - `flux`: Array with the measured flux
    - `low` : The lowest test frequency
    - `high`: The highest test frequency
    - `rate`: The sampling rate (spacing between frequencies)
    """
    # Generate test cyclic frequencies and convert to angular
    freq = np.arange(low, high, rate)
    ny = 2 * np.pi * freq

    # Empty array to store calculated power
    powers = np.zeros(shape=freq.shape)

    # The loop over frequencies (least mean square)
    for i in range(len(ny)):
        pcos = np.cos(ny[i] * time)
        psin = np.sin(ny[i] * time)
        s = np.sum(flux * psin)
        c = np.sum(flux * pcos)
        ss = np.sum(np.square(psin))
        cc = np.sum(np.square(pcos))
        sc = np.sum(psin * pcos)
        alpha = (s*cc - c*sc) / (ss*cc - sc**2)
        beta = (c*ss - s*sc) / (ss*cc - sc**2)
        freq_power = alpha**2 + beta**2
        powers[i] = freq_power

    # Return an array of test (cyclic) frequencies and the calculated power
    return freq, powers


###############################################################################
# Script
###############################################################################
# Initial setup
datdir = 'testdata/'
outdir = 'output/'
compare = False

# Load data
infile = 'ts_7days.txt'
time, flux = np.loadtxt(datdir + infile, unpack=True)

# Convert time to megaseconds
time *= 1e-6

# Run power spectrum
freq, powers = powerspec(time, flux, 1900.0, 4100.0, 0.1)

# Compare to the true oscillations?
if compare:
    # Load module
    import matplotlib.pyplot as plt

    # Load correct oscillations
    oscfile = 'oscillations.dat'
    l, n, nu, A, delta = np.loadtxt(datdir + oscfile, unpack=True)

    # Plot
    plt.figure()
    plt.plot(freq, powers, 'r-')
    plt.plot(nu, A**2, 'g*')
    plt.title('Power spectrum')
    plt.xlabel('nu [muHz]')
    plt.ylabel('|V(t)|^2')
    plt.savefig(outdir + 'test1.pdf')
