/*  ~~~ Time Series Analysis -- Frequency CLEAN ~~~
 *
 * Usage:
 * fclean.x [options] -n number sampling inputfile outputfile
 *
 * Number: Select the number of frequencies to CLEAN. E.g. -n 3
 * 
 * Sampling: -f {auto | low high factor}
 *   auto: Calculate power spectrum from 5 microHertz to the Nyquist frequency
 *         with four times oversampling (auto is a key word, use as "-f auto").
 *   low high factor: Values for sampling in microHz (e.g. "-f 1500 4000 2", to
 *             sample from 1500 to 4000 microHz with 2 times oversampling).
 *
 * Options:
 *  -w: Calculate weighted power spectrum -- requires an extra column in the
 *      input file containing weight per data point.
 *  -q: Quiet-mode. No output to console.
 *  -t{sec|day|ms}: Unit of input file (seconds [default], days, megaseconds).
 *  -noprep: Do not subtract the mean of time series (for artificial data where
 *           the mean is 0).
 *  -fast: Fast-mode. Disable Nyquist calculation (and hence automatic
 *         sampling range) for lower runtime. Activates quiet-mode. Use for
 *         benchmarking the pure I/O + algorithm.
 *
 * Note:
 * Using multi-threading with OpenMP. Set number of threads used by the shell
 * variable "OMP_NUM_THREADS".
 *
 * Author: Jakob Rørsted Mosumgaard
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <omp.h>

#include "fileio.h"
#include "arrlib.h"
#include "tsfourier.h"

#define PI2micro 6.28318530717958647692528676655900576839433879875e-6


int main(int argc, char *argv[])
{
    /* Important definitions */
    // Lengths
    size_t N = 0;  // Length of time series
    size_t M = 0;  // Length of sampling vector (number of frequencies)

    // Filenames
    char inname[100];
    char outname[100];
    char logname[100];

    // Sampling
    double low, high, rate;

    // Options
    int quiet = 0;
    int unit = 1;
    int prep = 1;
    int autosamp = 0;
    int fast = 0;
    int useweight = 0;
    int Nclean = 1;
    int filter = 0;

    
    /* Process command line arguments and return line count of the input file */
    N = cmdarg(argc, argv, inname, outname, &quiet, &unit, &prep, &low, &high,\
               &rate, &autosamp, &fast, &useweight, NULL, NULL, &Nclean, &filter);

    // Pretty print
    if ( quiet == 0 || fast == 1){
        if ( useweight != 0 )
            printf("\nCLEANing the time series \"%s\" using weights...\n",\
                   inname);
        else
            printf("\nCLEANing the time series \"%s\" without weights...\n",\
                   inname);
    }

    
    /* Read data (and weights) from the input file */
    if ( quiet == 0 ) printf(" - Reading input\n");
    double* time = malloc(N * sizeof(double));
    double* flux = malloc(N * sizeof(double));
    double* weight = malloc(N * sizeof(double));
    readcols(inname, time, flux, weight, N, useweight, unit, quiet);

    // Do if fast-mode is not activated
    if ( fast == 0 ) {
        // Calculate Nyquist frequency
        double* dt = malloc(N-1 * sizeof(double));
        double nyquist;
        arr_diff(time, dt, N);
        nyquist = 1.0 / (2.0 * arr_median(dt, N-1)) * 1e6; // microHz !
        free(dt);

        // Calculate sampling (N times oversampling: N is stored in 'rate')
        int oversamp = (int) rate;
        double minsamp = 1.0e6 / (oversamp * (time[N-1] - time[0])); // microHz !
    
        // Display info?
        if ( quiet == 0 ){
            printf(" -- INFO: Length of time series = %li\n", N);
            printf(" -- INFO: Nyquist frequency = %.2lf microHz\n", nyquist);
            printf(" -- INFO: Using %i times oversampling = %.3lf microHz\n",\
                   oversamp, minsamp);
        }

        // Apply sampling in a provided range or the full range
        if ( autosamp == 0 ) {
            rate = minsamp;
        }
        else { 
            low = 5.0;
            high = nyquist;
            rate = minsamp;
        }
    }
    else {
        // Only set the suggested sampling
        int oversamp = (int) rate;
        double minsamp = 1.0e6 / (oversamp * (time[N-1] - time[0])); // microHz !
        rate = minsamp;
    }

    
    /* Prepare for power spectrum */
    // Get length of sampling vector
    M = arr_util_getstep(low, high, rate);

    // Fill sampling vector with cyclic frequencies
    double* freq = malloc(M * sizeof(double));
    arr_init_linspace(freq, low, rate, M);
    if ( quiet == 0 )
        printf(" -- INFO: Number of sampling frequencies = %li\n", M);

    // Subtract the mean to avoid "zero-frequency" problems
    double fmean = 0;
    if ( prep != 0 ) {
        if ( quiet == 0 ) printf(" - Subtracting the mean from time series\n");
        fmean = arr_mean(flux, N);
        arr_sca_add(flux, -fmean, N);
    }
    else {
        if ( quiet == 0 )
            printf(" - Time series used *without* mean subtraction!\n");
    }


    /* Prepare file for writing the CLEAN-output */
    // Create log-file
    strcpy(logname, outname);
    strcat(logname, ".cleanlog");
    FILE* logfile = fopen(logname, "w");

    // Write header
    fprintf(logfile, "# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"\
            "~~~~~~~~~~\n");
    fprintf(logfile, "# Log of CLEAN on \"%s\"\n", inname);
    fprintf(logfile, "# Interval: [%.2lf, %.2lf] microHz\n", low, high);
    fprintf(logfile, "# Finding %i frequencies\n", Nclean);
    fprintf(logfile, "# \n");
    fprintf(logfile, "# %8s %11s %11s %12s %12s\n", "Number", "Frequency",\
            "Power", "Alpha", "Beta");
    fprintf(logfile, "# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"\
            "~~~~~~~~~~\n");
    

    
    /* Find and clean peaks */
    // Init variables
    double fmax, alpmax, betmax, powmax;
    
    // Display info
    if ( quiet == 0 ) {
        printf(" - CLEANing %i frequencies in the range %.1lf to %.1lf"\
               " microHz\n", Nclean, low, high);
        printf("\n %9s %11s %11s\n", "Number", "Frequency", "Power");
    }

    // Enter CLEAN-loop
    for (int i = 0; i < Nclean; ++i) {
        // Display progress
        if ( quiet == 0 ) printf(" %6i", i+1);

        // Reset variables
        fmax = 0;
        alpmax = 0;
        betmax = 0;
        powmax = 0;

        // Call with or without weights
        if ( useweight == 0 )
            fouriermax(time, flux, NULL, freq, N, M, &fmax, &alpmax, &betmax, 0);
        else
            fouriermax(time, flux, weight, freq, N, M, &fmax, &alpmax, &betmax, 1);

        // Calculate the power and write to log
        powmax = alpmax*alpmax + betmax*betmax;
        fprintf(logfile, " %6i %15.6lf %12.6lf %12.6lf %12.6lf\n", i+1, fmax,\
                powmax, alpmax, betmax);
        if ( quiet == 0) printf(" %15.6lf %12.6lf \n", fmax, powmax);

        // Remove frequency from time series
        for (int j = 0; j < N; ++j) {
            flux[j] = flux[j] - alpmax * sin( PI2micro*fmax * time[j] ) - \
                                betmax * cos( PI2micro*fmax * time[j] );
        }
    }

    // Final touch
    fclose(logfile);
    if ( quiet == 0 ) printf("\n");

    
    /* Write CLEANed time series to file */
    if ( quiet == 0 ) printf(" - Saving to file \"%s\"\n", outname);

    // Add the mean to the time series data again
    if ( prep != 0 ) arr_sca_add(flux, fmean, N);

    // Save to file
    writecols3(outname, time, flux, weight, N, useweight, unit);

    
    /* Free data */
    free(time);
    free(flux);
    free(weight);
    free(freq);


    
    /* Done! */
    if ( quiet == 0 || fast ==1 ) printf("Done!\n\n");
    return 0; 
}
