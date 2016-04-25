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

    // Options
    int quiet = 0;
    int unit = 1;
    int prep = 1;

    
    /* Process command line arguments and return line count of the input file */
    N = cmdarg(argc, argv, inname, outname, &quiet, &unit, &prep);
    
    // Pretty print
    if ( quiet == 0 )
        printf("\nCalculating the power spectrum of \"%s\" ...\n", inname);


    /* Read data from the file */
    if ( quiet == 0 ) printf(" - Reading input\n");
    double* time = malloc(N * sizeof(double));
    double* flux = malloc(N * sizeof(double));
    readcols(inname, time, flux, N, unit);

    
    /* Prepare for power spectrum */
    // Get length of sampling vector
    M = arr_util_getstep(low, high, rate);

    // Fill sampling vector with cyclic frequencies
    double* freq = malloc(M * sizeof(double));
    arr_init_linspace(freq, low, rate, M);

    // Initialise arrays for data storage
    double* power = malloc(M * sizeof(double));

    // Subtract the mean to avoid "zero-frequency" problems
    if ( prep != 0 ) {
        if ( quiet == 0 ) printf(" - Subtracting the mean from time series\n");
        arr_sca_add(flux, -arr_mean(flux, N), N);
    }
    else {
        if ( quiet == 0 ) printf(" - Time series used *without* mean subtraction!\n");
    }
    
    /* Calculate power spectrum */
    if ( quiet == 0 ) printf(" - Calculating fourier transform\n");
    fourier(time, flux, freq, N, M, power);

    
    /* Write data to file */
    if ( quiet == 0 ) printf(" - Saving to file\n");
    writecols(outname, freq, power, M);

    /* Free data */
    free(time);
    free(flux);
    free(freq);
    free(power);
        
    /* Done! */
    if ( quiet == 0 ) printf("Done!\n\n");
    return 0; 
}
