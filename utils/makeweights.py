# -*- coding: utf-8 -*-

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Time Series Analysis -- Generate statistical weigts from scatter
#
# Author: Jakob Rørsted Mosumgaard
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###############################################################################
# Modules
###############################################################################
from __future__ import print_function, with_statement, division
import numpy as np
import bottleneck as bn


###############################################################################
# Functions
###############################################################################
def genweight(datname, dpath, wpath):
    """
    Combine time series with statistical weights calculated from scatter

    Arguments:
    - `datname`: Identifier of data file
    - `dpath`  : Path to data file (time series).
    - `wpath`  : Path to scatter file (with same time points!)
    """

    # Pretty print
    print('Generating weights for {0} !'.format(dpath))

    # Load data and weights
    t, d = np.loadtxt(dpath, unpack=True)
    tt, sig = np.loadtxt(wpath, unpack=True)

    # Check that times are indeed the same
    tdif = t - tt
    if tdif.any() != 0:
        print('Error! Not the same time points! Quitting!')
        exit()

    # Moving variance (Hans: M = 50 - 100)
    M = 70
    movstd = bn.move_std(sig, M, min_count=1)
    movvar = np.square(movstd)

    # Remove first point
    x = 1
    t = t[x:]
    d = d[x:]
    movvar = movvar[x:]

    # Calculate weights from scatter (1 / variance)
    w = np.divide(1.0, movvar)

    # Save
    outfile = star + '_with-weights.txt'
    np.savetxt(outfile, np.transpose([t, d, w]), fmt='%.15e', delimiter='\t')

    # Done!
    print('Done!\n')


###############################################################################
# Script
###############################################################################
if __name__ == "__main__":
    # Definitions
    datdir = '../../data/'
    ext = '.txt'
    append = '-high'

    # Run for star 1
    star = 'star01'
    genweight(star, datdir + star + ext, star + append + ext)

    # Run for star 2
    star = 'star02'
    genweight(star, datdir + star + ext, star + append + ext)
