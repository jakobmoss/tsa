/*  ~~~ Time Series Analysis -- Power Spectrum using Least Squares ~~~
 *
 * Usage:
 * powerspec.x [options] sampling inputfile outputfile
 *
 * Sampling: -f {auto | low high rate}
 *   auto: Calculate power spectrum from 5 microHertz to Nyquist frequency
 *         with four times oversampling (auto is a key word, use as "-f auto").
 *   low high rate: Values for sampling in microHz (e.g. "-f 1500 4000 0.1", to
 *                  sample from 1500 to 4000 microHz in steps of 0.1 microHz).
 *
 * Options:
 *  -q: Quiet-mode. No output to console.
 *  -t{sec|day|ms}: Unit of input file (seconds [default], days, megaseconds).
 *  -noprep: Do not subtract the mean of time series (for artificial data where
 *           the mean is 0).
 *  -fast: Fast-mode. Disable Nyquist calculation (and hence automatic
 *         sampling) for lower runtime. Activates quiet-mode automatically. Use
 *         for benchmarking the pure I/O + algorithm.
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

    
    /* Process command line arguments and return line count of the input file */
    N = cmdarg(argc, argv, inname, outname, &quiet, &unit, &prep, &low, &high,\
               &rate, &autosamp, &fast);
    
    // Pretty print
    if ( quiet == 0 || fast == 1)
        printf("\nCalculating the power spectrum of \"%s\" ...\n", inname);


    /* Read data from the file */
    if ( quiet == 0 ) printf(" - Reading input\n");
    double* time = malloc(N * sizeof(double));
    double* flux = malloc(N * sizeof(double));
    readcols(inname, time, flux, N, unit, quiet);

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
            printf(" -- INFO: Suggested minimum sampling = %.2lf microHz\n", minsamp);
        }

        // Apply automatic sampling?
        if ( autosamp != 0 ) {
            low = 5.0;
            high = nyquist;
            rate = minsamp;
        }
    }
    
    /* Prepare for power spectrum */
    // Get length of sampling vector
    M = arr_util_getstep(low, high, rate);

    // Fill sampling vector with cyclic frequencies
    double* freq = malloc(M * sizeof(double));
    arr_init_linspace(freq, low, rate, M);

    // Initialise arrays for data storage
    double* power = malloc(M * sizeof(double));

    // Subtract the mean to avoid "zero-frequency" problems
    if ( prep != 0 ) {
        if ( quiet == 0 ) printf(" - Subtracting the mean from time series\n");
        arr_sca_add(flux, -arr_mean(flux, N), N);
    }
    else {
        if ( quiet == 0 ) printf(" - Time series used *without* mean subtraction!\n");
    }
    
    /* Calculate power spectrum */
    if ( quiet == 0 ){
        printf(" - Calculating fourier transform\n");
        if ( autosamp != 0 ) {
            printf(" -- NB: Using automatic sampling!\n");
            printf(" -- INFO: Auto-sampling (in microHz): %.2lf to %.2lf in steps of %.2lf\n", low, high, rate);
        }
        else
            printf(" -- INFO: Sampling (in microHz): %.2lf to %.2lf in steps of %.2lf\n", low, high, rate);
        printf(" -- INFO: Number of sampling frequencies = %li\n", M);
    }
    fourier(time, flux, freq, N, M, power);

    
    /* Write data to file */
    if ( quiet == 0 ) printf(" - Saving to file \"%s\"\n", outname);
    writecols(outname, freq, power, M);

    /* Free data */
    free(time);
    free(flux);
    free(freq);
    free(power);
        
    /* Done! */
    if ( quiet == 0 || fast ==1 ) printf("Done!\n\n");
    return 0; 
}
