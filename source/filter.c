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
#include "pass.h"

#define PI2micro 6.28318530717958647692528676655900576839433879875e-6


int main(int argc, char *argv[])
{
    /* Important definitions */
    // Lengths
    size_t N = 0;     // Length of time series

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

    // Pretty print
    if ( quiet == 0 || fast == 1){
        if ( useweight != 0 )
            printf("\nFiltering the time series \"%s\" using weights...\n",\
                   inname);
        else
            printf("\nFiltering the time series \"%s\" without weights...\n",\
                   inname);
    }

    
    /* Read data (and weights) from the input file */
    if ( quiet == 0 ) printf(" - Reading input\n");
    double* time = malloc(N * sizeof(double));
    double* flux = malloc(N * sizeof(double));
    double* weight = malloc(N * sizeof(double));
    readcols(inname, time, flux, weight, N, useweight, unit, quiet);

    // Do if fast-mode is not activated
    if ( fast == 0 ) {
        // Calculate Nyquist frequency
        double* dt = malloc(N-1 * sizeof(double));
        double nyquist;
        arr_diff(time, dt, N);
        nyquist = 1.0 / (2.0 * arr_median(dt, N-1)) * 1e6; // microHz !
        free(dt);

        // Calculate suggested sampling (4 times oversampling)
        double minsamp;
        minsamp = 1.0e6 / (4 * (time[N-1] - time[0])); // microHz !
    
        // Display info?
        if ( quiet == 0 ){
            printf(" -- INFO: Length of time series = %li\n", N);
            printf(" -- INFO: Nyquist frequency = %.2lf microHz\n", nyquist);
            printf(" -- INFO: Suggested minimum sampling = %.3lf microHz\n",\
                   minsamp);
        }

        // Apply automatic sampling?
        if ( autosamp != 0 ) {
            low = 5.0;
            high = nyquist;
            rate = minsamp;
        }
    }

    
    /* Run the desired filter */
    // Init output array
    double* filt = malloc(N * sizeof(double));

    if ( filter == 2 ) {
        if ( quiet == 0 ) {
            printf(" - Calculating bandpass filter between %.2lf and %.2lf"\
                   " microHz\n", fstart, fstop);
        }
        bandpass(time, flux, weight, N, fstart, fstop, low, high, rate, filt,\
                 useweight, quiet);
    }
    else if ( filter == 3 ) {
        if ( quiet == 0 ) {
            printf(" - Calculating lowpass filter up to %.2lf microHz\n",\
                   fstop);
        }
        lowpass(time, flux, weight, N, fstop, low, high, rate, filt,\
                useweight, quiet);
    }
    else if ( filter == 4 ) {
        if ( quiet == 0 ) {
            printf(" - Calculating highpass filter from %.2lf microHz\n",\
                   fstop);
        }
        highpass(time, flux, weight, N, fstop, low, high, rate, filt,\
                 useweight, quiet);
    }
    else {
        fprintf(stderr, "ERROR: Unknown filter chosen !");
    }

    /* Write filtered time series to file */
    if ( quiet == 0 ) printf(" - Saving to file \"%s\"\n", outname);

    // Save to file
    writecols3(outname, time, filt, weight, N, useweight, unit);

    
    /* Free data */
    free(time);
    free(flux);
    free(weight);
    free(filt);


    /* Done! */
    if ( quiet == 0 || fast ==1 ) printf("Done!\n\n");
    return 0; 
}
