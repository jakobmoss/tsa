int cmdarg(int argc, char *argv[], char inname[], char outname[], int *quiet,\
           int *unit, int *prep, double *low, double *high, double *rate,\
           int *autosamp, int *fast, int *useweight);

void readcols(char *fname, double x[], double y[], double z[], size_t N,\
              int three, int unit, int quiet);

void writecols(char *fname, double x[], double y[], size_t N);
