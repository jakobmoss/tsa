/*  ~~~ Time Series Analysis -- Frequency CLEAN ~~~
 *
 * Usage:
 * fclean.x [options] -n number sampling inputfile outputfile
 *
 * Number: Select the number of frequencies to CLEAN. E.g. -n 3
 * 
 * Sampling: -f {auto | low high rate}
 *   auto: Calculate power spectrum from 5 microHertz to Nyquist frequency
 *         with four times oversampling (auto is a key word, use as "-f auto").
 *   low high rate: Values for sampling in microHz (e.g. "-f 1500 4000 0.1", to
 *                  sample from 1500 to 4000 microHz in steps of 0.1 microHz).
 *   limit rate: ONLY IN THE CASE OF window-function-mode. Sample the window
 *               function in the range +/- limit in steps of rate.
 *
 * Options:
 *  -w: Calculate weighted power spectrum -- requires an extra column in the
 *      input file containing weight per data point.
 *  -q: Quiet-mode. No output to console.
 *  -t{sec|day|ms}: Unit of input file (seconds [default], days, megaseconds).
 *  -noprep: Do not subtract the mean of time series (for artificial data where
 *           the mean is 0).
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
    int useweight = 0;

    
    /* Process command line arguments and return line count of the input file */
    N = cmdarg(argc, argv, inname, outname, &quiet, &unit, &prep, &low, &high,\
               &rate, &autosamp, NULL, &useweight, NULL, NULL, 1);

    // Pretty print
    if ( quiet == 0 ){
        if ( useweight != 0 )
            printf("\nCLEANing frequencies from \"%s\" using weights...\n",\
                   inname);
        else
            printf("\nCLEANing frequencies from \"%s\" without weights...\n",\
                   inname);
    }

    
    /* Read data (and weights) from the input file */
    if ( quiet == 0 ) printf(" - Reading input\n");
    double* time = malloc(N * sizeof(double));
    double* flux = malloc(N * sizeof(double));
    double* weight = malloc(N * sizeof(double));
    readcols(inname, time, flux, weight, N, useweight, unit, quiet);

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
        printf(" -- INFO: Suggested minimum sampling = %.3lf microHz\n", \
               minsamp);
    }

    // Apply automatic sampling?
    if ( autosamp != 0 ) {
        low = 5.0;
        high = nyquist;
        rate = minsamp;
    }


    /* Prepare for power spectrum */
    // Get length of sampling vector
    M = arr_util_getstep(low, high, rate);

    // Fill sampling vector with cyclic frequencies
    double* freq = malloc(M * sizeof(double));
    arr_init_linspace(freq, low, rate, M);

    // Initialise arrays for data storage
    double* power = malloc(M * sizeof(double));

    
    /* Free data */
    free(time);
    free(flux);
    free(weight);
    free(freq);
    free(power);


    
    /* Done! */
    if ( quiet == 0 ) printf("Done!\n\n");
    return 0; 
}
