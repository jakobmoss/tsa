/*  ~~~ Time Series Analysis -- auxiliary ~~~
 *
 * Actual calculation of the fourier transform
 * 
 * Author: Jakob Rørsted Mosumgaard
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "arrlib.h"
#include "fmin.h"

#define PI2micro 6.28318530717958647692528676655900576839433879875e-6
#define EPS 1.0e-9

void alpbet(double time[], double flux[], size_t N, double ny, double *alpha, \
            double *beta);

void alpbetW(double time[], double flux[], double weight[], size_t N,\
             double ny, double wsum, double *alpha, double *beta);


/* Calculate the fourier transform of time series
 *
 * Arguments:
 *  - `time`     : Array of times. In seconds!
 *  - `flux`     : Array of data.
 *  - `weight`   : Array of statistical weights.
 *  - `freq`     : Array of cyclic frequencies to sample.
 *  - `N`        : Length of the time series
 *  - `M`        : Length of the sampling vector
 *  - `power`    : OUTPUT -- Array with powers
 *  - `alpha`    : OUTPUT -- Array with alphas
 *  - `beta`     : OUTPUT -- Array with betas
 *  - `useweight`: Flag to signal whether to use weights or not (0 = no weights)
 */
void fourier(double time[], double flux[], double weight[], double freq[],\
             size_t N, size_t M, double power[], double alpha[], double beta[],\
             int useweight)
{
    // Local variables
    double alp = 0;
    double bet = 0;
    double ny = 0;
    size_t i;

    // Call functions with or without weights
    if ( useweight == 0 ) {
        // Make parallel loop over all test frequencies
        #pragma omp parallel default(shared) private(alp, bet, ny)
        {
            #pragma omp for schedule(static)
            for (i = 0; i < M; ++i) {
                // Current frequency
                ny = freq[i] * PI2micro;

                // Calculate alpha and beta
                alpbet(time, flux, N, ny, &alp, &bet);
                
                // Store alpha, beta and power
                alpha[i] = alp;
                beta[i] = bet;
                power[i] = alp*alp + bet*bet;
            }
        }
    }
    else {
        // Sum of all weights
        double sumweights = arr_sum(weight, N);

        // Make parallel loop over all test frequencies
        #pragma omp parallel default(shared) private(alp, bet, ny)
        {
            #pragma omp for schedule(static)
            for (i = 0; i < M; ++i) {
                // Current frequency
                ny = freq[i] * PI2micro;

                // Calculate alpha and beta
                alpbetW(time, flux, weight, N, ny, sumweights, &alp, &bet);
                
                // Store alpha, beta and power
                alpha[i] = alp;
                beta[i] = bet;
                power[i] = alp*alp + bet*bet;
            }
        }
    }
}


// Calculate alpha and beta coefficients
void alpbet(double time[], double flux[], size_t N, double ny, double *alpha, \
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
    ss = N - cc;

    // Calculate coefficients
    D = ss*cc - sc*sc;
    *alpha = (s * cc - c * sc)/D;
    *beta  = (c * ss - s * sc)/D;
}


// Calculate alpha and beta coefficients  -- USING WEIGHTS
void alpbetW(double time[], double flux[], double weight[], size_t N,\
             double ny, double wsum, double *alpha, double *beta)
{
    // Auxiliary
    double sn, cn, D;
    
    // Sums
    double s = 0;
    double c = 0;
    double cc = 0;
    double sc = 0;
    double ss;

    // Loop over the time series
    for (size_t i = 0; i < N; ++i) {
        // Pre-calculate sin, cos of point
        sn = sin(ny * time[i]);
        cn = cos(ny * time[i]);

        // Calculate sin, cos terms
        s += weight[i] * flux[i] * sn;
        c += weight[i] * flux[i] * cn;

        // Calculate squared and cross terms
        cc += weight[i] * cn * cn;
        sc += weight[i] * sn * cn;
    }

    // Calculate ss from cc
     ss = wsum - cc;

    // Calculate coefficients
    D = ss*cc - sc*sc;
    *alpha = (s * cc - c * sc)/D;
    *beta  = (c * ss - s * sc)/D;
}


/* Calculate the fourier transform of time series and find the highest peak
 *  --> Helper routine for CLEAN
 *
 * Arguments:
 *  - `time`     : Array of times. In seconds!
 *  - `flux`     : Array of data.
 *  - `weight`   : Array of statistical weights.
 *  - `freq`     : Array of cyclic frequencies to sample.
 *  - `N`        : Length of the time series
 *  - `M`        : Length of the sampling vector
 *  - `fmax`     : OUTPUT -- Frequency of maximum power
 *  - `alpmax`   : OUTPUT -- Alpha of that frequency
 *  - `betmax`   : OUTPUT -- Beta of that frequency
 *  - `useweight`: Flag to signal whether to use weights or not (0 = no weights)  
 */
void fouriermax(double time[], double flux[], double weight[], double freq[],\
                size_t N, size_t M, double *fmax, double *alpmax,\
                double *betmax, int useweight)
{
    // Local variables
    double alpha = 0;
    double beta = 0;
    double ny = 0;
    size_t i;

    // Local variables for finding the peak
    double p = 0;
    double pmaxlocal;
    double nymaxlocal;

    // Maximum power (global)
    double pmax = 0;
    double nymax = 0;

    // Call functions with or without weights
    if ( useweight == 0 ) {
        // Function for minimisation (nested for variable access)
        double powopt(double optny)
        {
            double optalpha, optbeta, optpower;
            alpbet(time, flux, N, optny, &optalpha, &optbeta);
            optpower = optalpha*optalpha + optbeta*optbeta;
            return -optpower;
        }
        
        // Make parallel loop over all test frequencies
        #pragma omp parallel default(shared) private(alpha, beta, ny, p, pmaxlocal, nymaxlocal)
        {
            // Reset varibles
            pmaxlocal = 0;
            nymaxlocal = 0;

            // Do the loop (nowait -> each threads can move on to comparison)
            #pragma omp for schedule(static) nowait
            for (i = 0; i < M; ++i) {
                // Current frequency
                ny = freq[i] * PI2micro;

                // Calculate alpha, beta and power
                alpbet(time, flux, N, ny, &alpha, &beta);
                p = alpha*alpha + beta*beta;

                // Compare to current maximum power
                if ( p > pmaxlocal ) {
                    pmaxlocal = p;
                    nymaxlocal = ny;
                }
            }

            // Make sure we use the maximum from all the threads
            // NOTE: Double check, since the critical region is slow and should
            //       only be entered when necessary (and value can be changed
            //       by several threads, see: goo.gl/lwnzTn)!
            if ( pmaxlocal > pmax ) {
                #pragma omp critical
                {
                    if ( pmaxlocal > pmax ) {
                        pmax = pmaxlocal;
                        nymax = nymaxlocal;
                    }
                }
            }
        }

        // Search around found peak for the "true" minimum
        double df = PI2micro * (freq[1] - freq[0]);
        pmax = - fmin_golden(powopt, nymax-df, nymax+df, EPS, &nymax);

        // Store the optimised values
        alpbet(time, flux, N, nymax, alpmax, betmax);
        *fmax = nymax/PI2micro;
    }
    else {
        // Sum of all weights
        double sumweights = arr_sum(weight, N);

        // Function for minimisation (nested for variable access)
        double powopt(double optny)
        {
            double optalpha, optbeta, optpower;
            alpbetW(time, flux, weight, N, optny, sumweights, &optalpha, &optbeta);
            optpower = optalpha*optalpha + optbeta*optbeta;
            return -optpower;
        }

        // Make parallel loop over all test frequencies
        #pragma omp parallel default(shared) private(alpha, beta, ny, p, pmaxlocal, nymaxlocal)
        {
            // Reset varibles
            pmaxlocal = 0;
            nymaxlocal = 0;

            // Do the loop (nowait -> each threads can move on to comparison)
            #pragma omp for schedule(static) nowait
            for (i = 0; i < M; ++i) {
                // Current frequency
                ny = freq[i] * PI2micro;

                // Calculate alpha, beta and power
                alpbetW(time, flux, weight, N, ny, sumweights, &alpha, &beta);
                p = alpha*alpha + beta*beta;

                // Compare to current maximum power
                if ( p > pmaxlocal ) {
                    pmaxlocal = p;
                    nymaxlocal = ny;
                }
            }

            // Make sure we use the maximum from all the threads
            // NOTE: Double check, since the critical region is slow and should
            //       only be entered when necessary (and value can be changed
            //       by several threads, see: goo.gl/lwnzTn)!
            if ( pmaxlocal > pmax ) {
                #pragma omp critical
                {
                    if ( pmaxlocal > pmax ) {
                        pmax = pmaxlocal;
                        nymax = nymaxlocal;
                    }
                }
            }
        }

        // Search around found peak for the "true" minimum
        double df = PI2micro * (freq[1] - freq[0]);
        pmax = - fmin_golden(powopt, nymax-df, nymax+df, EPS, &nymax);

        // Store the optimised values
        alpbetW(time, flux, weight, N, nymax, sumweights, alpmax, betmax);
        *fmax = nymax/PI2micro;
    }
    // Done!
}
