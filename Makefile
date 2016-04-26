# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Time Series Analysis -- Master Makefile
#
# Author: Jakob Rørsted Mosumgaard
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Settings
EXEC = powerspec.x
DAYS = 7


# Default target
default: $(EXEC)

# Generate test data
.PHONY: data
data:
	$(MAKE) -C testdata all

# Build C executable
$(EXEC): source/*.c
	$(MAKE) -C source all

# Build Cython module
.PHONY: cython
cython:
	$(MAKE) -C cython all

# Run test
.PHONY: test
test: $(EXEC) data
	@mkdir -p output
	@printf "\n" 
	time ./$(EXEC) -f 1900 4100 0.1 testdata/ts_$(DAYS)days.txt output/ctest.txt
	@printf "\n" 
	cp output/ctest.txt test/
	$(MAKE) -C test default

# Housekeeping
.PHONY: clean
clean:
	$(RM) $(EXEC)
	$(RM) output/*.txt output/*.pdf
	$(MAKE) -C source clean
	$(MAKE) -C testdata clean
	$(MAKE) -C test clean
	$(MAKE) -C cython clean
