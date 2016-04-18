# -*- coding: utf-8 -*-

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Time Series Analysis -- Python interface to Cython module.
#
# Author: Jakob Rørsted Mosumgaard
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


###############################################################################
# Modules
###############################################################################
import fourier


###############################################################################
# Script
###############################################################################
# Initial setup
datdir = '../testdata/'
outdir = '../output/'
compare = False

# Run power spectrum
freq, powers = fourier.calc()

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
