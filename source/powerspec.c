/*  ~~~ Time Series Analysis -- Power Spectrum using Least Squares ~~~
 *
 * Usage:
 * powerspec.x [options] sampling inputfile outputfile
 *
 * Sampling: -f {auto | low high rate | limit rate}
 *   auto: Calculate power spectrum from 5 microHertz to Nyquist frequency
 *         with four times oversampling (auto is a key word, use as "-f auto").
 *   low high rate: Values for sampling in microHz (e.g. "-f 1500 4000 0.1", to
 *                  sample from 1500 to 4000 microHz in steps of 0.1 microHz).
 *   limit rate: ONLY IN THE CASE OF window-function-mode. Sample the window
 *               function in the range +/- limit in steps of rate.
 *
 * Special options:
 *  -window f0: Activate window-function-mode. Calculate the (power spectrum of)
 *              the window function at frequency f0 (in microHertz). Uses the
 *              provided input times and weights, as well as the provided
 *              frequency sampling.
 *              NOTE: Remember to specify the frequency f0! Unit is microHertz.
 *              NOTE 2: Remember to specify a correct sampling (see above).
 *              NOTE 3: This switch will *not* produce a power spectrum of the
 *                    input data.
 *              NOTE 4: The -noprep switch is meaningless with this option,
 *                      since the input data is not used (only the times).
 *              NOTE 5: Will activate pseudo fast-mode and hence disable
 *                      automatic sampling and Nyquist calculation.

 * Options:
 *  -w: Calculate weighted power spectrum -- requires an extra column in the
 *      input file containing weight per data point.
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

    // Frequency of window function
    double winfreq = 0;

    // Options
    int quiet = 0;
    int unit = 1;
    int prep = 1;
    int autosamp = 0;
    int fast = 0;
    int useweight = 0;
    int windowmode = 0;
    int Nclean = 0;
    int filter = 0;

    
    /* Process command line arguments and return line count of the input file */
    N = cmdarg(argc, argv, inname, outname, &quiet, &unit, &prep, &low, &high,\
               &rate, &autosamp, &fast, &useweight, &windowmode, &winfreq,\
               &Nclean, &filter, NULL, NULL);
    
    // Pretty print
    if ( quiet == 0 || fast == 1 ){
        if ( windowmode == 0 && useweight != 0 )
            printf("\nCalculating the weighted power spectrum of \"%s\" ...\n",\
                   inname);
        else if ( windowmode == 0 )
            printf("\nCalculating the power spectrum of \"%s\" ...\n", inname);
        else
            printf("\nCalculating the window function of \"%s\" ...\n", inname);
    }
    

    /* Read data (and weights) from the input file */
    if ( quiet == 0 ) printf(" - Reading input\n");
    double* time = malloc(N * sizeof(double));
    double* flux = malloc(N * sizeof(double));
    double* weight = malloc(N * sizeof(double));
    readcols(inname, time, flux, weight, N, useweight, unit, quiet);
    
    // Do if fast-mode and window-mode is not activated
    if ( fast == 0 && windowmode == 0 ) {
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

    
    /* Prepare for power spectrum */
    // Calculate proper frequency range for window-function-mode
    double limit = 0;
    if ( windowmode != 0 ) {
        limit = low;
        low = winfreq - limit;
        high = winfreq + limit;
    }
    
    // Get length of sampling vector
    M = arr_util_getstep(low, high, rate);

    // Fill sampling vector with cyclic frequencies
    double* freq = malloc(M * sizeof(double));
    arr_init_linspace(freq, low, rate, M);

    // Initialise arrays for data storage
    double* power = malloc(M * sizeof(double));
    double* alpha = malloc(M * sizeof(double));
    double* beta = malloc(M * sizeof(double));


    /* Calculate power spectrum OR window function */
    if ( windowmode == 0 ) {
        // Subtract the mean to avoid "zero-frequency" problems
        if ( prep != 0 ) {
            if ( quiet == 0 ){
                printf(" - Subtracting the mean from time series\n");
            }
            arr_sca_add(flux, -arr_mean(flux, N), N);
        }
        else {
            if ( quiet == 0 )
                printf(" - Time series used *without* mean subtraction!\n");
        }
    
        // Display info
        if ( quiet == 0 ){
            printf(" - Calculating fourier transform\n");
            if ( autosamp != 0 ) {
                printf(" -- NB: Using automatic sampling!\n");
                printf(" -- INFO: Auto-sampling (in microHz): %.2lf to %.2lf"\
                       " in steps of %.4lf\n", low, high, rate);
            }
            else {
                printf(" -- INFO: Sampling (in microHz): %.2lf to %.2lf in"\
                       " steps of %.4lf\n", low, high, rate);
            }
            printf(" -- INFO: Number of sampling frequencies = %li\n", M);
        }

        // Calculate power spectrum with or without weights
        fourier(time, flux, weight, freq, N, M, power, alpha, beta, useweight);
    }
    else {
        if ( quiet == 0 ){
            printf(" - Calculating window function\n");
            printf(" -- INFO: Window frequency = %.2lf microHz\n", winfreq);
            printf(" -- INFO: Sampling in the range +/- %.2lf microHz in" \
                   " steps of %.4lf microHz\n", limit, rate);
            printf(" -- INFO: Number of sampling frequencies = %li\n", M);
        }

        // Calculate spectral window with or without weights
        windowfunction(time, freq, weight, N, M, winfreq, power, useweight);

        if ( quiet == 0 )
            printf(" - Sum of spectral window = %.4lf\n", arr_sum(power, M));

        // Move frequencies to the origin
        arr_sca_add(freq, -winfreq, M);
    }

        
    /* Write data to file */
    if ( quiet == 0 ) printf(" - Saving to file \"%s\"\n", outname);
    writecols(outname, freq, power, M);

    
    /* Free data */
    free(time);
    free(flux);
    free(weight);
    free(freq);
    free(power);
    free(alpha);
    free(beta);


    /* Done! */
    if ( quiet == 0 || fast ==1 ) printf("Done!\n\n");
    return 0; 
}
