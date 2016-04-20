#include <stdio.h>
#include <math.h>
#include <stdlib.h>

#include <omp.h>

#define PI 3.141592653589793
#define PI2 6.283185307179586

#define NUM_THREADS 4

// Auxialiary: Calculate alpha and beta coefficients
void alpbet(double time[], double flux[], double nu, size_t N, double *alpha, double *beta)
{
    // Auxiliary
    double sn, cn, D;
    
    // Sums
    double s = 0;
    double c = 0;
    double cc = 0;
    double ss = 0;
    double sc = 0;

    // Loop over the time series
    for (size_t i = 0; i < N; ++i) {
        // Pre-calculate sin, cos of point
        sn = sin(nu * time[i]);
        cn = cos(nu * time[i]);

        // Calculate sin, cos terms
        s += flux[i] * sn;
        c += flux[i] * cn;

        // Calculate squared and cross terms
        ss += sn * sn;
        cc += cn * cn;
        sc += sn * cn;
        
    }

    // Calculate coefficients
    D = ss*cc - sc*sc;
    *alpha = (s * cc - c * sc)/D;
    *beta  = (c * ss - s * sc)/D;
}



// Calculate fourier transform
void fourier(double time[], double flux[], double ny[], size_t N, size_t M, \
             double power[])
{
    // Local variables
    double alpha = 0;
    double beta = 0;
    double nu = 0;
    size_t i;

#pragma omp parallel num_threads(NUM_THREADS) default(shared) private(alpha, beta, nu)
    {
#pragma omp for schedule(static)
        for (i = 0; i < M; ++i) {
            // Current frequency
            nu = ny[i];

            // Calculate alpha and beta
            alpbet(time, flux, nu, N, &alpha, &beta);
                
            // Store power
            power[i] = alpha*alpha + beta*beta;
        }
    }
}
