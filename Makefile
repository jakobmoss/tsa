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
$(EXEC):
	$(MAKE) -C source all

# Run test
.PHONY: test
test: $(EXEC)
	@printf "\n" 
	time ./$(EXEC) testdata/ts_$(DAYS)days.txt output/ctest.txt
	@printf "\n" 
	cp output/ctest.txt test/
	$(MAKE) -C test default


# Housekeeping
.PHONY: clean
clean:
	$(RM) $(EXEC)
	$(RM) output/*.txt
	$(MAKE) -C source clean
	$(MAKE) -C testdata clean
	$(MAKE) -C test clean
