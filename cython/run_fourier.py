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
import numpy as np

# Cython module
import fourier


###############################################################################
# Set-up
###############################################################################
compare = False

# Initial setup
datdir = '../testdata/'
outdir = '../output/'

# Settings
infile = 'ts_7days.txt'
unit = 's'
prep = False

# Sampling of power spectrum (in microHertz)
low = 1900.0
high = 4100.0
rate = 0.1


###############################################################################
# Prepare for Plotting
###############################################################################
if compare:
    import matplotlib as mpl  # Pyplot is loaded after setup!

    def matplotlib_setup():
        fig_width_pt = 328
        inches_per_pt = 1.0 / 72.27
        golden_mean = (np.sqrt(5.0) - 1) / 2.0
        fig_width = fig_width_pt * inches_per_pt
        fig_height = fig_width * golden_mean
        fig_size = [fig_width, fig_height]
        mpl.rc('text', usetex=True)
        mpl.rc('figure', figsize=fig_size)
        mpl.rc('font', size=8, family='serif')
        mpl.rc('axes', labelsize=8)
        mpl.rc('legend', fontsize=8)
        mpl.rc('xtick', labelsize=8)
        mpl.rc('ytick', labelsize=8)
        mpl.rc('text.latex',
               preamble=r'\usepackage[T1]{fontenc}\usepackage{libertine}\usepackage[libertine]{newtxmath}')

    matplotlib_setup()
    import matplotlib.pyplot as plt

###############################################################################
# Script
###############################################################################
# Run power spectrum
freq, powers = fourier.calc(datdir + infile, low, high, rate, unit, prep)

# Compare to the true oscillations?
if compare:
    print('Plotting....')

    # Load correct oscillations
    oscfile = 'oscillations.dat'
    l, n, nu, A, delta = np.loadtxt(datdir + oscfile, unpack=True)

    # Plot
    plt.figure()
    plt.plot(freq, powers, 'r-', label=r'Calculated')
    plt.plot(nu, A**2, 'g.', label=r'True')
    plt.title('Power spectrum')
    plt.xlabel(r'$\nu$ [$\mu\textup{Hz}$]')
    plt.ylabel(r'Power')
    plt.legend()
    plt.savefig(outdir + 'test2.pdf', bbox_inches='tight')

    print('Done!')
