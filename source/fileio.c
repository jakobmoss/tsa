/*  ~~~ Time Series Analysis -- Auxiliary ~~~
 *
 * Routines for input/output
 * 
 * Author: Jakob Rørsted Mosumgaard
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* Check command-line argument and count lines in given file */
int cmdarg(int argc, char *argv[], char inname[], char outname[], int *quiet,\
           int *unit, int *prep, double *low, double *high, double *rate,\
           int *autosamp, int *fast, int *useweight, int *windowmode,\
           double *winfreq)
{
    // Internal
    int isamp = 0;
    int ifast = 0;
    
    // Quit if wrong number of arguments is given!
    if (argc < 5) {
        fprintf(stderr, "usage: %s  [-window] [-w] [-q] [-t{sec|day|ms}]"\
                " [-noprep] [-fast] -f {auto | low high rate} input_file"\
                " output_file\n", argv[0]);
        exit(1);
    }

    // Loop through given arguments (skipping program name)
    for (int i = 1; i < argc; ++i) {
        // Optional arguments
        // Quiet-mode
        if ( strcmp(argv[i], "-q" ) == 0 ) {
            *quiet = 1;
        }
        // Units
        else if ( strcmp(argv[i], "-tsec" ) == 0 ) {
            *unit = 1;
        }
        else if ( strcmp(argv[i], "-tday" ) == 0 ) {
            *unit = 2;
        }
        else if ( strcmp(argv[i], "-tms" ) == 0 ) {
            *unit = 3;
        }
        // Modify data
        else if ( strcmp(argv[i], "-noprep" ) == 0 ) {
            *prep = 0;
        }
        // Fast-mode
        else if ( strcmp(argv[i], "-fast" ) == 0 ) {
            *fast = 1;
            ifast = 1;
        }
        // Weights
        else if ( strcmp(argv[i], "-w" ) == 0 ) {
            *useweight = 1;
        }
        // Windowfunction-mode
        else if ( strcmp(argv[i], "-window" ) == 0 ) {
            *windowmode = 1;

            // Read frequency
            i++;
            *winfreq = atof(argv[i]);
        }
        // Sampling
        else if ( strcmp(argv[i], "-f" ) == 0 ) {
            // Go to first option
            i++;
            
            // Check for automatic sampling
            if ( strcmp(argv[i], "auto") == 0) {
                isamp = 1;
                *autosamp = 1;
            }
            // If manual, check that enough arguments is left
            else if ( i + 4 <= argc - 1) {
                isamp = 2;
                
                // Read the values and increment i
                *low = atof(argv[i]);
                i++;
                *high = atof(argv[i]);
                i++;
                *rate = atof(argv[i]);
            }
            else {
                break;
            }
        }
        // Non-optional arguments (filenames)
        else {
            // Read input file
            strcpy(inname, argv[i]);

            // Increment i and read output file
            i++;
            strcpy(outname, argv[i]);
        }
    }

    // Exit if no (or wrong) sampling provided
    if ( isamp == 0 ) {
        fprintf(stderr, "No or wrong sampling provided! Quitting!\n");
        exit(1);
    }

    // Override options if fast-mode is activated
    if ( ifast == 1 ) {
        printf(" * Fast-mode activated. Going (almost) quiet * \n");
        *quiet = 1;
        if ( isamp == 1 ) {
            fprintf(stderr, "Cannot autosample in fast mode! Quitting!\n");
            exit(1);
        }
    }
    

    // Read file and quit if it cannot be opened
    size_t ch, number_of_lines = 0;
    FILE* tmpfile = fopen(inname, "r");
    if ( tmpfile == 0 ) {
        fprintf(stderr, "Could not open file:  %s \n", inname);
        exit(1);
    }
    
    // Go through file and handle files not ending in newline
    do {
        ch = fgetc(tmpfile);
        if(ch == '\n')
            number_of_lines++;
    } while (ch != EOF);
    if(ch != '\n' && number_of_lines != 0)
        number_of_lines++;

    // Close file and save number of actual lines
    fclose(tmpfile);
    return number_of_lines - 1;
}


/* Read file with two or three columns of data */
void readcols(char *fname, double x[], double y[], double z[], size_t N,\
              int three, int unit, int quiet)
{
    // Read the file depending on the number of columns
    FILE* infile = fopen(fname, "r");
    if ( three == 0) {
        for (size_t i = 0; i < N+1; ++i) {
            if ( fscanf(infile ,"%lf%lf", &x[i], &y[i] ) != 2) break;
        }
    }
    else {
        if ( quiet == 0 ) printf(" -- INFO: Using weights\n");
        for (size_t i = 0; i < N+1; ++i) {
            if ( fscanf(infile ,"%lf%lf%lf", &x[i], &y[i], &z[i] ) != 3) break;
        }
    }
    fclose(infile);

    // Convert unit into seconds
    if ( unit != 1 ) {
        double scaling = 1;
        
        // Days
        if ( unit == 2 ) {
            scaling = 86400.0;
            if ( quiet == 0 ) printf(" -- INFO: Unit is %s\n", "days");
        }
        // Mega seconds
        else if ( unit == 3 ) {
            scaling = 1e6;
            if ( quiet == 0 ) printf(" -- INFO: Unit is %s\n", "megaseconds");
        }
        // Failure -- keep scaling of 1
        else {
            fprintf(stderr,"Error: Wrong unit. Assuming seconds.\n");
        }

        // Apply to the time vector
        for (size_t j = 0; j < N+1; ++j) {
            x[j] *= scaling;
        }
    }
    else {
        if ( quiet == 0 ) printf(" -- INFO: Unit is %s\n", "seconds");
    }
}


/* Write file with two columns of data */
void writecols(char *fname, double x[], double y[], size_t N)
{
    FILE* outfile = fopen(fname, "w");

    // Check if file is available
    if (outfile != NULL) {
        size_t i;
        for ( i = 0; i < N; ++i ) {
            fprintf(outfile, "%15.9e %18.9e\n", x[i], y[i]);
        }
        fclose(outfile);
    }
}
