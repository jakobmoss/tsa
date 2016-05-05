/*  ~~~ Time Series Analysis -- auxiliary ~~~
 *
 * Actual calculation of the fourier transform
 * 
 * Author: Jakob Rørsted Mosumgaard
 */

#include <stdlib.h>
#include <math.h>
#include <omp.h>

#define PI2micro 6.28318530717958647692528676655900576839433879875e-6

void alpbet(double time[], double flux[], double ny, size_t N, double *alpha,\
            double *beta);


/* Calculate the fourier transform of time series
 *
 * Arguments:
 *  - `time`: Array of times. In seconds!
 *  - `flux`: Array of data.
 *  - `freq`: Array of cyclic frequencies to sample.
 *  - `N`: Length of the time series
 *  - `M`: Length of the sampling vector
 *  - `power`: OUTPUT: Array with powers
 */
void fourier(double time[], double flux[], double freq[], size_t N, size_t M, \
             double power[])
{
    // Local variables
    double alpha = 0;
    double beta = 0;
    double ny = 0;
    size_t i;

    #pragma omp parallel default(shared) private(alpha, beta, ny)
    {
        #pragma omp for schedule(static)
        for (i = 0; i < M; ++i) {
            // Current frequency
            ny = freq[i] * PI2micro;

            // Calculate alpha and beta
            alpbet(time, flux, ny, N, &alpha, &beta);
                
            // Store power
            power[i] = alpha*alpha + beta*beta;
        }
    }
}


// Calculate alpha and beta coefficients
void alpbet(double time[], double flux[], double ny, size_t N, double *alpha,\
            double *beta)
{
    // Auxiliary
    double sn, cn, D;
    
    // Sums
    double s = 0;
    double c = 0;
    double cc = 0;
    double sc = 0;
    double ss;

    // TEMPORARY (until weights are used)
    double wsum = N;

    // Loop over the time series
    for (size_t i = 0; i < N; ++i) {
        // Pre-calculate sin, cos of point
        sn = sin(ny * time[i]);
        cn = cos(ny * time[i]);

        // Calculate sin, cos terms
        s += flux[i] * sn;
        c += flux[i] * cn;

        // Calculate squared and cross terms
        cc += cn * cn;
        sc += sn * cn;
    }

    // Calculate ss from cc
    ss = wsum - cc;

    // Calculate coefficients
    D = ss*cc - sc*sc;
    *alpha = (s * cc - c * sc)/D;
    *beta  = (c * ss - s * sc)/D;
}
