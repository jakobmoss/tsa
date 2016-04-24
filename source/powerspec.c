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

    // Filenames
    char inname[100];
    char outname[100];

    
    /* Process command line arguments and return line count of the input file */
    N = cmdarg(argc, argv, inname, outname);
    
    // Pretty print
    printf("\nCalculating the power spectrum of \"%s\" ...\n", argv[1]);

    // TESTING
    printf("\nInput: \"%s\" \nOutput: \"%s\"\n\n", inname, outname);

    
    /* Read data from the file */
    double* time = malloc(N * sizeof(double));
    double* flux = malloc(N * sizeof(double));
    readcols(inname, time, flux, N);

    
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
    writecols(outname, freq, power, M);

    /* Free data */
    free(time);
    free(flux);
    free(freq);
    free(power);
        
    /* Done! */
    printf("Done!\n\n");
    return 0; 
}
