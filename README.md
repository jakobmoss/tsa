# Time Series Analysis #

Software for performing Time Series Analysis (targeted at astrophysical applications) written in C, Cython and Python.

Author: Jakob Rørsted Mosumgaard (jakob@phys.au.dk).


### Features ###

* Complete program to make a power spectrum (the fourier transform of a time series) written in C. The software is using OpenMP for a performance boost using multithreading.
* Stand-alone Cython-module, which is providing a Python interface to the fast C function. The interface has a very low overhead and almost as fast runtimes as the pure C.
* Pure Python/NumPy implementation of the algorithm for comparison (and for easy-to-read overview of the algorithm).
* Python script to generate artificial test data for easy verification.


### Installation ###

Requirements: 
* GCC (tested with v5.3.0)
* Python 3 (tested with v3.5.1)
* Cython (tested with v0.24)
* NumPy (tested with v1.11.0)
* Gnuplot (tested with v5.0; not a strict requirement, only for running tests)

Setup: The whole project is controlled by several Makefile. Note that the Makefiles assumes working in a virtual environment, where `python` is calling `python3`. The targets available from the root-dir are:
* `make` will build the C-executable.
* `make data` will create artificial data for testing.
* `make test` will run the abovementioned targets and make a test-run and a plot.
* `make cython` will make the stand-alone Cython module.


### Usage ###
The usage of the software is documented in the different files; it should be straigtforward for both C and Python/Cython.


### File structure ###
The project contains the following directories:
* cython: Cython module and a Python wrapper.
* python: Pure Python implementation.
* source: C source code.
* test: Default test-target.
* testdata: Generation and storage of artificial data.

Furthermore, the following directories should be created:
* output: Default place for output from all programs. Automatically created by `make test`.