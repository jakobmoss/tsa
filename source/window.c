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
                    double window[])
{
    ;
    
}

