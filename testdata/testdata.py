# -*- coding: utf-8 -*-

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Time Series Analysis -- Generate artificial data for testing
#
# Author: Jakob Rørsted Mosumgaard
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###############################################################################
# Modules
###############################################################################
from __future__ import print_function, with_statement, division
import numpy as np


###############################################################################
# Functions
###############################################################################
def makets(days=14, step=60):
    """
    Generate time-series from stellar oscillation data and save to file.

    Arguments:
    - `days`: How many days of data to generate data (default: 14, corresponds
              to approximately 20000 points).
    - `step`: Sampling rate in seconds (default: 60).
    """

    # Pretty print
    print('Generating time series for {0:3g} days of data'.format(days))

    # Load oscillations from file
    l, n, nu, A, delta = np.loadtxt('oscillations.dat', unpack=True)
    omega = 2*np.pi * nu * 1e-6  # Angular frequency in Hertz

    # Define sampling
    tstart = 0
    tend = days * 86400 + 60
    t = np.arange(tstart, tend, step)

    # Make time-series data
    v = np.empty(len(t))
    for i in range(len(t)):
        v[i] = np.sum(A * np.cos(omega * t[i] + delta))

    # Display info
    print('Number of data points: {0:6g}'.format(len(t)))
    nyquist = np.divide(1, 2*np.median(np.diff(t)))
    print('Nyquist frequency: {0:.3f} muHz'.format(1e6*nyquist))

    # Save
    outfile = 'ts_' + str(days) + 'days.txt'
    np.savetxt(outfile, np.transpose([t, v]), fmt='%.15e', delimiter='\t')

    # Done!
    print('Done!\n')


###############################################################################
# Script
###############################################################################
# Define the time series to generate
days = [1, 7, 14, 30]

# Run
for length in days:
    makets(days=length)
