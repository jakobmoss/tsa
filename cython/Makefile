# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Time Series Analysis: Cython module
#
# Author: Jakob Rørsted Mosumgaard
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Default target
all: fourier.so

# Cython module
fourier.so: setup.py fourier.pyx tsfourier.c
	python $< build_ext --inplace

tsfourier.c:
	cp ../source/tsfourier.[c,h] .

# Housekeeping
clean:
	$(RM) -r fourier.c fourier*.so build tsfourier.[c,h]
