# README #

Software for performing Time Series Analysis written in C, Cython and Python.

Author: Jakob Rørsted Mosumgaard (jakob@phys.au.dk).

### Features ###

* Complete program to make a power spectrum (the fourier transform of a time series) written in C. Includes I/O and a call to the underlying fourier-function. The actual calculation is running in parallel with OpenMP.
* Stand-alone Cython-module, which is providing a Python interface to the fast C function.
* Pure Python/NumPy implementation of the algorithm for comparison.


### Installation ###

Requirements: 
* Python 3 (tested with v3.5.1)
* Cython (tested with v0.24)
* GCC (tested with v5.3.0)
* Gnuplot (tested with v5.0; not a strict requirement, only for running tests)

Setup: The whole project is controlled by several Makefile. The targets available from the root-dir are:
* `make` will build the C-executable.
* `make data` will create artificial data for testing.
* `make test` will run the abovementioned targets and make a test-run and a plot.
* `make cython` will make the stand-alone Cython module.


### File structure ###
The project contains the following directories:
* cython: Cython module and a Python wrapper.
* python: Pure Python implementation.
* source: C source code.
* test: Default test-target.
* testdata: Generation and storage of artificial data.

Furthermore, the following directories should be created:
* output: Default place for output from all programs. Automatically created by `make test`.