#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

#include "fileio.h"
#include "arrlib.h"
#include "tsfourier.h"


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
        fprintf(stderr, "usage: %s  path/to/data  path/to/output\n", argv[0]);
    }
    else {
        N = countlines(argc, argv);
    }

    // Pretty print
    printf("\nCalculating the power spectrum of \"%s\" ...\n", argv[1]);

    
    /* Read data from the file */
    double* time = malloc(N * sizeof(double));
    double* flux = malloc(N * sizeof(double));
    readcols(argv[1], time, flux, N);

    
    /* Prepare for power spectrum */
    // Get length of sampling vector
    M = arr_util_getstep(low, high, rate);

    // Fill sampling vector with cyclic frequencies
    double* freq = malloc(M * sizeof(double));
    arr_init_linspace(freq, low, rate, M);

    // Initialise arrays for data storage
    double* power = malloc(M * sizeof(double));

    
    /* Calculate power spectrum */
    fourier(time, flux, freq, N, M, power);

    
    /* Write data to file */
    writecols(argv[2], freq, power, M);

    /* Free data */
    free(time);
    free(flux);
    free(freq);
    free(power);
        
    /* Done! */
    printf("Done!\n\n");
    return 0; 
}
