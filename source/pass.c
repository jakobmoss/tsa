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
              double result[], size_t Nres, int useweight, int quiet)
{
    // Calculate the (sum of the) window function

    // Fill sampling vector with cyclic frequencies
    size_t M = arr_util_getstep(low, high, rate);
    double* freq = malloc(M * sizeof(double));
    arr_init_linspace(freq, low, rate, M);
    if ( quiet == 0 )
        printf(" -- INFO: Number of sampling frequencies = %li\n", M);


    // Done!
    free(freq);
}
