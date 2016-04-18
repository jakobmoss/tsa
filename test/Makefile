# Makefile for producing test-plot

default: test.pdf

test.pdf: test.plt ctest.txt
	gnuplot $<
	$(RM) $<

test.plt: Makefile
	@echo set term pdfcairo > $@
	@echo set out \'test.pdf\' >> $@
	@echo set title \'Power Spectrum\' >> $@
	@echo set xlabel \'ny [microHz]\' >> $@
	@echo set ylabel \'Power\' >> $@
	@echo plot \'ctest.txt\' with lines title \'Calculated\' \\>> $@
	@echo , \'../testdata/oscillations.dat\' using 3:\(\$$$ 4**2\) title \'True\' >> $@

clean:
	$(RM) *.pdf *.txt
