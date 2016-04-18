# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Time Series Analysis -- Master Makefile
#
# Author: Jakob Rørsted Mosumgaard
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Default target
default: c

# Generate test data
.PHONY: data
data:
	$(MAKE) -C testdata all

# Build C executable
c:
	$(MAKE) -C source all

# Housekeeping
.PHONY: clean
clean:
	$(RM) powerspec.x
	$(MAKE) -C source clean
	$(MAKE) -C testdata clean
