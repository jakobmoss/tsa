/*  ~~~ Time Series Analysis -- Auxiliary ~~~
 *
 * Routines for calculating the spectral window
 * 
 * Author: Jakob Rørsted Mosumgaard
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "arrlib.h"

#define PI2micro 6.28318530717958647692528676655900576839433879875e-6

void windowalpbet(double time[], double datasin[], double datacos[], size_t N,\
                  double ny, double *alphasin, double *betasin,
                  double *alphacos, double *betacos);

void windowalpbetW(double time[], double weight[], double datasin[],\
                   double datacos[], size_t N, double ny, double wsum,\
                   double *alphasin, double *betasin, double *alphacos,\
                   double *betacos);


/* Calculate the window function of a time series
 *
 * Arguments:
 *  - `time`     : Array of times. In seconds!
 *  - `freq`     : Array of cyclic frequencies to sample.
 *  - `weight`   : Array of statistical weights.
 *  - `N`        : Length of the time series
 *  - `M`        : Length of the sampling vector
 *  - `window`   : OUTPUT -- Array with power of the window
 *  - `useweight`: Flag to signal whether to use weights or not (0 = no weights)
 */
void windowfunction(double time[], double freq[], double weight[], size_t N,\
                    size_t M, double f0, double window[], int useweight)
{
    // Sample the time series using cos and sin at frequency f0
    double* datsin = malloc(N * sizeof(double));
    double* datcos = malloc(N * sizeof(double));
    double omega0 = f0 * PI2micro;
    for (size_t k = 0; k < N; ++k) {
        datsin[k] = sin(omega0 * time[k]);
        datcos[k] = cos(omega0 * time[k]);
    }
    
    // Initialise local variables
    double alphasin = 0;
    double betasin = 0;
    double alphacos = 0;
    double betacos = 0;
    double ny = 0;
    size_t i;

    // Call functions with or without weights
    if ( useweight == 0 ) {
        // Make parallel loop over all test frequencies
        #pragma omp parallel default(shared) private(alphasin, betasin, alphacos, betacos, ny)
        {
            #pragma omp for schedule(static)
            for (i = 0; i < M; ++i) {
                // Current frequency
                ny = freq[i] * PI2micro;

                // Calculate alpha and beta for cos and sin data
                windowalpbet(time, datsin, datcos, N, ny, &alphasin, &betasin, \
                             &alphacos, &betacos);
                
                // Store power
                window[i] = 0.5 * ( (alphasin*alphasin + betasin*betasin) + \
                                    (alphacos*alphacos + betacos*betacos)    );
            }
        }
    }
    else {
        // Sum of all weights
        double sumweights = arr_sum(weight, N);
        
        // Make parallel loop over all test frequencies
        #pragma omp parallel default(shared) private(alphasin, betasin, alphacos, betacos, ny)
        {
            #pragma omp for schedule(static)
            for (i = 0; i < M; ++i) {
                // Current frequency
                ny = freq[i] * PI2micro;

                // Calculate alpha and beta for cos and sin data
                windowalpbetW(time, weight, datsin, datcos, N, ny, sumweights, \
                              &alphasin, &betasin, &alphacos, &betacos);
                
                // Store power
                window[i] = 0.5 * ( (alphasin*alphasin + betasin*betasin) + \
                                    (alphacos*alphacos + betacos*betacos)    );
            }
        }
    }

    // Done
    free(datsin);
    free(datcos);
}


// Calculate alpha and beta coefficients
void windowalpbet(double time[], double datasin[], double datacos[], size_t N,\
            double ny, double *alphasin, double *betasin, double *alphacos,\
            double *betacos)
{
    // Auxiliary
    double sn, cn, D;
    
    // Sums: Individual terms
    double ssin = 0;
    double csin = 0;
    double scos = 0;
    double ccos = 0;

    // Sums: Common terms
    double cc = 0;
    double sc = 0;
    double ss;

    // Loop over the time series
    for (size_t i = 0; i < N; ++i) {
        // Pre-calculate sin, cos of point
        sn = sin(ny * time[i]);
        cn = cos(ny * time[i]);

        // Calculate sin, cos terms for both data series
        ssin += datasin[i] * sn;
        csin += datasin[i] * cn;
        scos += datacos[i] * sn;
        ccos += datacos[i] * cn;

        // Calculate common squared and cross terms
        cc += cn * cn;
        sc += sn * cn;
    }

    // Calculate ss from cc
    ss = N - cc;

    // Calculate alpha and beta for both 
    D = ss*cc - sc*sc;
    *alphasin = (ssin * cc - csin * sc)/D;
    *betasin  = (csin * ss - ssin * sc)/D;
    *alphacos = (scos * cc - ccos * sc)/D;
    *betacos  = (ccos * ss - scos * sc)/D;
}


// Calculate alpha and beta coefficients WITH WEIGHTS
void windowalpbetW(double time[], double weight[], double datasin[],\
                   double datacos[], size_t N, double ny, double wsum,\
                   double *alphasin, double *betasin, double *alphacos,\
                   double *betacos)
{
    // Auxiliary
    double sn, cn, D;
    
    // Sums: Individual terms
    double ssin = 0;
    double csin = 0;
    double scos = 0;
    double ccos = 0;

    // Sums: Common terms
    double cc = 0;
    double sc = 0;
    double ss;

    // Loop over the time series
    for (size_t i = 0; i < N; ++i) {
        // Pre-calculate sin, cos of point
        sn = sin(ny * time[i]);
        cn = cos(ny * time[i]);

        // Calculate sin, cos terms for both data series
        // NOTE: The weights are already taken into account in the data!
        ssin += weight[i] * datasin[i] * sn;
        csin += weight[i] * datasin[i] * cn;
        scos += weight[i] * datacos[i] * sn;
        ccos += weight[i] * datacos[i] * cn;

        // Calculate common squared and cross terms
        cc += weight[i] * cn * cn;
        sc += weight[i] * sn * cn;
    }

    // Calculate ss from cc
    ss = wsum - cc;

    // Calculate alpha and beta for both 
    D = ss*cc - sc*sc;
    *alphasin = (ssin * cc - csin * sc)/D;
    *betasin  = (csin * ss - ssin * sc)/D;
    *alphacos = (scos * cc - ccos * sc)/D;
    *betacos  = (ccos * ss - scos * sc)/D;
}



/* Calculate the sum of the spectral window
 *
 * Arguments:
 * - `f0`         : Desired frequency of the window function
 * - `low`, `high`: Sampling limits 
 * - `rate`       : Frequency step of sampling
 * - `time`       : Array of times (from the time series). In seconds!
 * - `weight`     : Statistical weights per data point (pass NULL if no weights).
 * - `N`          : Length of the time series
 * - `useweight`  : If != 0 weights will be used.
 * - `quiet`      : If != 0 no output will be displayed to console
 */
double windowsum(double f0, double low, double high, double rate, double time[],
                 double weight[], size_t N, int useweight, int quiet)
{
    // Init
    double result = 0;

    // Calculate length of sampling vector
    size_t M = arr_util_getstep(low, high, rate);
    if ( quiet == 0 )
        printf(" -- INFO: Number of frequencies in the window = %li\n", M);


    // Initialise arrays and generate sampling frequencies
    double* freq = malloc(M * sizeof(double));
    double* window = malloc(M * sizeof(double));
    arr_init_linspace(freq, low, rate, M);

    // Calculate spectral window with or without weights
    windowfunction(time, freq, weight, N, M, f0, window, useweight);

    // Calculate the sum
    result = arr_sum(window, M);
    
    // Done
    free(freq);
    free(window);
    return result;
}
