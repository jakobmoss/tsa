#include <stdio.h>
#include <stdlib.h>
#include <string.h>


/* Check command-line argument and count lines in given file */
int cmdarg(int argc, char *argv[], char inname[], char outname[], int *quiet,\
           int *unit, int *prep, double *low, double *high, double *rate)
{
    // Init
    int samp = 0;
    
    // Quit if wrong number of arguments is given!
    if (argc < 3) {
        fprintf(stderr, "usage: %s  [-q] [-t{sec|day|ms}] [-noprep] [-f low high rate] input_file  output_file\n", argv[0]);
        exit(1);
    }

    // Loop through given arguments (skipping program name)
    for (int i = 1; i < argc; ++i) {
        // Optional arguments
        // Quiet-mode
        if ( strcmp(argv[i], "-q") == 0 ) {
            *quiet = 1;
        }
        // Units
        else if ( strcmp(argv[i], "-tsec") == 0 ) {
            *unit = 1;
        }
        else if ( strcmp(argv[i], "-tday") == 0 ) {
            *unit = 2;
        }
        // Modify data
        else if ( strcmp(argv[i], "-noprep") == 0 ) {
            *prep = 0;
        }
        // Sampling
        else if ( strcmp(argv[i], "-f") == 0 ) {
            samp = 1;

            // Check that enough arguments is left
            if ( i + 4 <= argc - 1) {
                // Increment i and read the values
                i++;
                *low = atof(argv[i]);
                i++;
                *high = atof(argv[i]);
                i++;
                *rate = atof(argv[i]);
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

    // AT THE MOMENT EXIT IF NO SAMPLING IS PROVIDED!
    if ( samp == 0 ) {
        fprintf(stderr, "No sampling provided! Unsupported option!\n");
        exit(1);
    }

    // Read file and quit if it cannot be opened
    size_t ch, number_of_lines = 0;
    FILE* tmpfile = fopen(inname, "r");
    if ( tmpfile == 0 ) {
        fprintf(stderr,"Could not open file:  %s \n", inname);
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


/* Read file with two columns of data */
void readcols(char *fname, double x[], double y[], size_t N, int unit)
{
    // Read the file
    FILE* infile = fopen(fname, "r");
    for (size_t i = 0; i < N+1; ++i) {
        if ( fscanf(infile ,"%lf%lf", &x[i], &y[i] ) != 2) {
            break;
        }
    }
    fclose(infile);

    // Convert unit into seconds
    if ( unit != 1 ) {
        double scaling = 1;
        
        // Days
        if ( unit == 2 ) {
            scaling = 86400.0;
        }
        // Mega seconds
        else if ( unit == 3 ) {
            scaling = 1e6;
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
