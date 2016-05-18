/*  ~~~ Time Series Analysis -- Data Filtering ~~~
 *
 * Usage:
 * filter.x [options] mode sampling inputfile outputfile
 *
 * Mode: -{band freq1 freq2 | low freq | high freq}
 *   band: Bandpass filter between freq1 and freq2 (in microHz).
 *   low : Lowpass filter with freq (in microHz) as limit.
 *   high: Highpass filter with freq (in microHz) as limit.
 * 
 * Sampling: -f {auto | low high rate}
 *   auto: Calculate power spectrum from 5 microHertz to Nyquist frequency
 *         with four times oversampling (auto is a key word, use as "-f auto").
 *   low high rate: Values for sampling in microHz (e.g. "-f 1500 4000 0.1", to
 *                  sample from 1500 to 4000 microHz in steps of 0.1 microHz).
 *
 * Options:
 *  -w: Calculate weighted power spectrum -- requires an extra column in the
 *      input file containing weight per data point.
 *  -q: Quiet-mode. No output to console.
 *  -t{sec|day|ms}: Unit of input file (seconds [default], days, megaseconds).
 *  -noprep: Do not subtract the mean of time series (for artificial data where
 *           the mean is 0).
 *  -fast: Fast-mode. Disable Nyquist calculation (and hence automatic
 *         sampling range) for lower runtime. Activates quiet-mode. Use for
 *         benchmarking the pure I/O + algorithm.
 *
 * Note:
 * Using multi-threading with OpenMP. Set number of threads used by the shell
 * variable "OMP_NUM_THREADS".
 *
 * Author: Jakob Rørsted Mosumgaard
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "fileio.h"
#include "arrlib.h"
#include "tsfourier.h"
#include "window.h"

#define PI2micro 6.28318530717958647692528676655900576839433879875e-6


int main(int argc, char *argv[])
{
    /* Important definitions */
    // Lengths
    size_t N = 0;  // Length of time series
    size_t M = 0;  // Length of sampling vector (number of frequencies)

    // Filenames
    char inname[100];
    char outname[100];

    // Sampling
    double low, high, rate;

    // Options
    int quiet = 0;
    int unit = 1;
    int prep = 1;
    int autosamp = 0;
    int fast = 0;
    int useweight = 0;
    int Nclean = 0;
    int filter = 1;  // 1 is init, 2 is bandpass, 3 is low, 4 is high

    // Filtering frequencies
    double fstart = 0;
    double fstop = 0;

    
    /* Process command line arguments and return line count of the input file */
    N = cmdarg(argc, argv, inname, outname, &quiet, &unit, &prep, &low, &high,\
               &rate, &autosamp, &fast, &useweight, NULL, NULL, &Nclean,\
               &filter, &fstart, &fstop);


    printf("f1 = %lf\nf2 = %lf\n", fstart, fstop);
    
    /* Done! */
    if ( quiet == 0 || fast ==1 ) printf("Done!\n\n");
    return 0; 
}
