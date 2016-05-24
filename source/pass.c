/*  ~~~ Time Series Analysis -- Auxiliary ~~~
 *
 * Routines for calculating (band, low, high)-pass filters
 * 
 * Author: Jakob Rørsted Mosumgaard
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "arrlib.h"
#include "window.h"
#include "tsfourier.h"

#define PI2micro 6.28318530717958647692528676655900576839433879875e-6



/* Bandpass filter
 *
 * Arguments:
 *  - `time`       : Array of times. In seconds!
 *  - `flux`       : Array of data.
 *  - `weight`     : Array of statistical weights.
 *  - `N`          : Length of the time series
 *  - `f1`, `f2`   : Frequency interval of filter (f1 < f2)
 *  - `low`, `high`: Frequency interval to calculate the spectrum
 *  - `rate`       : Frequency sampling
 *  - `result`     : OUTPUT -- Array containing filtered data
 *  - `useweight`  : Flag to signal whether to use weights or not (0 = no weights)
 *  - `quiet`      : Flag. 0 = verbose output. 1 = no output to console
 */
void bandpass(double time[], double flux[], double weight[], size_t N,\
              double f1, double f2, double low, double high, double rate,\
              double result[], int useweight, int quiet)
{
    // Calculate the (sum of the) window function at central frequency
    if ( quiet == 0 )
        printf(" -- TASK: Calculating window function ... \n");
    double fwin = (low + high)/2.0;
    double sumwin = windowsum(fwin, low, high, rate, time, weight, N,\
                              useweight, quiet);
    if ( quiet == 0 ) printf("      ... Done!\n");

    // Fill sampling vector with cyclic frequencies
    size_t M = arr_util_getstep(f1, f2, rate);
    double* freq = malloc(M * sizeof(double));
    arr_init_linspace(freq, f1, rate, M);
    if ( quiet == 0 )
        printf(" -- INFO: Number of sampling frequencies = %li\n", M);

    // Arrays for data storage
    double* power = malloc(M * sizeof(double));
    double* alpha = malloc(M * sizeof(double));
    double* beta = malloc(M * sizeof(double));

    // Subtract the mean to avoid "zero-frequency" problems
    double fmean = arr_mean(flux, N);
    arr_sca_add(flux, -fmean, N);

    // Calculate power spectrum and save alphas and betas
    if ( quiet == 0 ) printf(" -- TASK: Calculating power spectrum ... \n");
    fourier(time, flux, weight, freq, N, M, power, alpha, beta, useweight);
    if ( quiet == 0 ) printf("      ... Done!\n");

    // Generate new time series
    if ( quiet == 0 ) printf(" -- TASK: Calculating new time series ... \n");
    double sumfilt, ny;
    #pragma omp parallel default(shared) private(sumfilt, ny)
    {
        #pragma omp for schedule(static)
        for (size_t i = 0; i < N; ++i) {
            sumfilt = 0;
            for (size_t j = 0; j < M; ++j) {
                ny = freq[j] * PI2micro;
                sumfilt += alpha[j]*sin(ny*time[i]) + beta[j]*cos(ny*time[i]);
            }
            result[i] = sumfilt / sumwin;
        }
    }
    if ( quiet == 0 ) printf("      ... Done!\n");

    // Add the mean again -- also to the filter
    arr_sca_add(flux, fmean, N);
    arr_sca_add(result, fmean, N);

    // Done!
    free(freq);
    free(power);
    free(alpha);
    free(beta);
}


/* Lowpass filter
 *  --> Basically just a wrapper for the bandpass filter
 *
 * Arguments:
 *  - `time`       : Array of times. In seconds!
 *  - `flux`       : Array of data.
 *  - `weight`     : Array of statistical weights.
 *  - `N`          : Length of the time series
 *  - `flow`       : Frequency limit of filter
 *  - `low`, `high`: Frequency interval to calculate the spectrum
 *  - `rate`       : Frequency sampling
 *  - `result`     : OUTPUT -- Array containing filtered data
 *  - `useweight`  : Flag to signal whether to use weights or not (0 = no weights)
 *  - `quiet`      : Flag. 0 = verbose output. 1 = no output to console
 */
void lowpass(double time[], double flux[], double weight[], size_t N,\
             double flow, double low, double high, double rate,       \
             double result[], int useweight, int quiet)
{
    // Call bandpass filter from zero to lowpass frequency
    double fzero = rate; // Not defined for exactly zero!
    bandpass(time, flux, weight, N, fzero, flow, low, high, rate, result,\
             useweight, quiet);
}


/* Highpass filter
 *  --> Basically just a wrapper for the highpass filter
 *
 * Arguments:
 *  - `time`       : Array of times. In seconds!
 *  - `flux`       : Array of data.
 *  - `weight`     : Array of statistical weights.
 *  - `N`          : Length of the time series
 *  - `flow`       : Frequency limit of filter
 *  - `low`, `high`: Frequency interval to calculate the spectrum
 *  - `rate`       : Frequency sampling
 *  - `result`     : OUTPUT -- Array containing filtered data
 *  - `useweight`  : Flag to signal whether to use weights or not (0 = no weights)
 *  - `quiet`      : Flag. 0 = verbose output. 1 = no output to console
 */
void highpass(double time[], double flux[], double weight[], size_t N,\
              double fhigh, double low, double high, double rate,     \
              double result[], int useweight, int quiet)
{
    // Make temporary array
    double* temp = malloc(N * sizeof(double));

    // Run lowpass filter
    lowpass(time, flux, weight, N, fhigh, low, high, rate, temp, useweight,\
            quiet);
    
    // Calculate highpass
    for (size_t i = 0; i < N; ++i) {
        result[i] = flux[i] - temp[i];
    }

    // Done!
    free(temp);
}

