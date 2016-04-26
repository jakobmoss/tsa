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

    // Filenames
    char inname[100];
    char outname[100];

    // Sampling
    double low, high, rate;

    // Options
    int quiet = 0;
    int unit = 1;
    int prep = 1;

    
    /* Process command line arguments and return line count of the input file */
    N = cmdarg(argc, argv, inname, outname, &quiet, &unit, &prep, &low, &high,\
               &rate);
    
    // Pretty print
    if ( quiet == 0 )
        printf("\nCalculating the power spectrum of \"%s\" ...\n", inname);


    /* Read data from the file */
    if ( quiet == 0 ) printf(" - Reading input\n");
    double* time = malloc(N * sizeof(double));
    double* flux = malloc(N * sizeof(double));
    readcols(inname, time, flux, N, unit);

    // Calculate and print Nyquist frequency
    if ( quiet == 0 ) {
        double* dt = malloc(N-1 * sizeof(double));
        double nyquist;
        arr_diff(time, dt, N);
        nyquist = 1.0 / (2.0 * arr_median(dt, N-1));
        printf(" -- INFO: Nyquist frequency = %.2lf microHz\n", 1e6*nyquist);
        free(dt);
    }
    
    
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
    if ( quiet == 0 ){
        printf(" - Calculating fourier transform\n");
        printf(" -- INFO: Sampling (in microHz): %.2lf to %.2lf in steps of %.2lf\n", low, high, rate);
    }
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
