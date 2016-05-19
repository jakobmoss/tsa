/*  ~~~ Time Series Analysis -- Auxiliary ~~~
 *
 * Routines for calculating (band, low, high)-pass filters
 * 
 * Author: Jakob Rørsted Mosumgaard
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "arrlib.h"
#include "window.h"

#define PI2micro 6.28318530717958647692528676655900576839433879875e-6



/* Bandpass filter */
void bandpass(double time[], double flux[], double weight[], size_t N,\
              double f1, double f2, double low, double high, double rate,\
              double result[], int useweight, int quiet)
{
    // Calculate the (sum of the) window function at central frequency
    if ( quiet == 0 ) printf(" -- TASK: Calculating window function ... \n");
    double fwin = (f1 + f2)/2.0;
    double sumwin = windowsum(fwin, low, high, rate, time, weight, N, useweight);
    if ( quiet == 0 ) printf("      ... Done!\n");
    
    // Fill sampling vector with cyclic frequencies
    size_t M = arr_util_getstep(f1, f2, rate);
    double* freq = malloc(M * sizeof(double));
    arr_init_linspace(freq, f1, f2, M);
    if ( quiet == 0 )
        printf(" -- INFO: Number of sampling frequencies = %li\n", M);

    // Arrays for data storage
    double* power = malloc(M * sizeof(double));
    double* alpha = malloc(M * sizeof(double));
    double* beta = malloc(M * sizeof(double));

    // Calculate power spectrum and save alphas and betaas
    if ( quiet == 0 ) printf(" -- TASK: Calculating power spectrum ... \n");
    if ( quiet == 0 ) printf("      ... Done!\n");

    
    // Done!
    free(freq);
    free(power);
    free(alpha);
    free(beta);
    
}
