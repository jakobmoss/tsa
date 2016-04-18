# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Time Series Analysis: Cython module
#
# Author: Jakob Rørsted Mosumgaard
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Default target
all: fourier.so

# Module to do fourier transform
fourier.so: setup.py fourier.pyx
	python $< build_ext --inplace


# Housekeeping
clean:
	$(RM) -r fourier.c fourier*.so build
