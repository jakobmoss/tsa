int cmdarg(int argc, char *argv[], char inname[], char outname[], int *quiet,\
           int *unit, int *prep, double *low, double *high, double *rate);

void readcols(char *fname, double x[], double y[], size_t N, int unit);

void writecols(char *fname, double x[], double y[], size_t N);
