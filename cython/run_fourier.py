# -*- coding: utf-8 -*-

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Time Series Analysis -- Python interface to Cython module.
#
# Author: Jakob Rørsted Mosumgaard
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


###############################################################################
# Modules
###############################################################################
# General
from __future__ import print_function, with_statement, division

# Cython module
import fourier


###############################################################################
# Set-up
###############################################################################
# Flags
compare = False

# Initial setup
datdir = '../testdata/'
outdir = '../output/'

# Settings
infile = 'ts_7days.txt'
unit = 's'

# Sampling of power spectrum (in microHertz)
low = 1900.0
high = 4100.0
rate = 0.1

###############################################################################
# Script
###############################################################################
# Run power spectrum
freq, powers = fourier.calc(datdir + infile, unit, low, high, rate)

# Compare to the true oscillations?
if compare:
    # Load modules
    import numpy as np
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
    plt.savefig(outdir + 'test2.pdf')
