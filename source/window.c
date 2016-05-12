/*  ~~~ Time Series Analysis -- Auxiliary ~~~
 *
 * Routines for calculating the spectral window
 * 
 * Author: Jakob Rørsted Mosumgaard
 */

#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "arrlib.h"

#define PI2micro 6.28318530717958647692528676655900576839433879875e-6
#define PI2 6.28318530717958647692528676655900576839433879875

void windowalpbet(double time[], double datasin[], double datacos[], double ny,\
            size_t N, double *alphasin, double *betasin, double *alphacos,\
            double *betacos);


/* Calculate the window function of a time series
 *
 * Arguments:
 *  - `time` : Array of times. In seconds!
 *  - `freq` : Array of cyclic frequencies to sample.
 *  - `N`    : Length of the time series
 *  - `M`    : Length of the sampling vector
 *  - `window`: OUTPUT -- Array with power of the window
 */
void windowfunction(double time[], double freq[], size_t N, size_t M, \
                    double f0, double window[])
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

    // Make parallel loop over all test frequencies
    #pragma omp parallel default(shared) private(alphasin, betasin, alphacos, betacos, ny)
    {
        #pragma omp for schedule(static)
        for (i = 0; i < M; ++i) {
            // Current frequency
            ny = freq[i] * PI2micro;

            // Calculate alpha and beta for cos and sin data
            windowalpbet(time, datsin, datcos, ny, N, &alphasin, &betasin,\
                   &alphacos, &betacos);
                
            // Store power
            window[i] = 0.5 * ( (alphasin*alphasin + betasin*betasin) + \
                                (alphacos*alphacos + betacos*betacos)    );
        }
    }

    // Done
    free(datsin);
    free(datcos);
}


// Calculate alpha and beta coefficients
void windowalpbet(double time[], double datasin[], double datacos[], double ny,\
            size_t N, double *alphasin, double *betasin, double *alphacos,\
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
