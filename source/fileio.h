int cmdarg(int argc, char *argv[], char inname[], char outname[], int *quiet,\
           int *unit, int *prep);

void readcols(char *fname, double x[], double y[], size_t N, int unit);

void writecols(char *fname, double x[], double y[], size_t N);
