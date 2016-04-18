#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "fileio.h"
#include "arrlib.h"
#include "tsfourier.h"

#define PI 3.141592653589793
#define PI2 6.283185307179586


int main(int argc, char *argv[])
{
    /* Important definitions */
    // Lengths
    size_t N = 0;  // Length of time series
    size_t M = 0;  // Length of sampling vector (number of frequencies)

    // Sampling
    double low = 1900.0;
    double high = 4100.0;
    double rate = 0.1;

    
    /* Get filename from command line and count the number of lines */
    if (argc != 3) {
        printf("Usage: %s  path/to/data  path/to/output\n", argv[0]);
    }
    else {
        N = countlines(argc, argv);
    }

    
    /* Read data from the file */
    double* t = malloc(N * sizeof(double));
    double* time = malloc(N * sizeof(double));
    double* flux = malloc(N * sizeof(double));
    readcols(argv[1], t, flux, N);

    // Convert to mega seconds
    double s_to_ms = 1e-6;
    arr_scale(t, s_to_ms, time, N);

    
    /* Prepare for power spectrum */
    // Get length of sampling vector
    M = arr_util_getstep(low, high, rate);

    // Fill sampling vector (in both cyclic and angular frequency)
    double* freq = malloc(M * sizeof(double)); // Cyclic
    arr_init_linspace(freq, low, rate, M);
    double* ny = malloc(M * sizeof(double)); // Angular
    arr_scale(freq, PI2, ny, M);

    // Initialise arrays for data storage
    double* alpha = malloc(M * sizeof(double));
    double* beta = malloc(M * sizeof(double));
    double* power = malloc(M * sizeof(double));

    
    /* Calculate power spectrum */
    fourier(time, flux, ny, N, M, power);

    /* Write data to file */
    writecols(argv[2], freq, power, M);

    /* Free data */
    free(t);
    free(time);
    free(flux);
    free(freq);
    free(ny);
    free(alpha);
    free(beta);
    free(power);
        
    
    return 0; 
}
