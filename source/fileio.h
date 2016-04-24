int cmdarg(int argc, char *argv[], char inname[], char outname[], int *quiet,\
           int *unit);

void readcols(char *fname, double x[], double y[], size_t N);

void writecols(char *fname, double x[], double y[], size_t N);
