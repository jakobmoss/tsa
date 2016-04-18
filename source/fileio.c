#include <stdio.h>


/* Check command-line argument and count lines in given file */
int countlines(int argc, char *argv[])
{
    size_t ch, number_of_lines = 0;
    FILE* tmpfile = fopen(argv[1], "r");

    // Check if file can be opened
    if ( tmpfile == 0 ) {
        printf( "Could not open file:  %s \n", argv[0]);
    }
    // If yes: read through file and handle files not ending in newline
    else {
        do {
            ch = fgetc(tmpfile);
            if(ch == '\n')
                number_of_lines++;
        } while (ch != EOF);
        if(ch != '\n' && number_of_lines != 0)
            number_of_lines++;
    }
    // Close file and save number of actual lines
    fclose(tmpfile);
    return number_of_lines - 1;
}


/* Read file with two columns of data */
void readcols(char *fname, double x[], double y[], size_t N)
{
    FILE* infile = fopen(fname, "r");

    // Check if file exists and loop through
    if (infile != NULL) {
        size_t i;
        for (i = 0; i < N+1; ++i) {
            if ( fscanf(infile ,"%lf%lf", &x[i], &y[i] ) != 2) {
               break;
            }
         }
        fclose(infile);
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