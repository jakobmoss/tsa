void fourier(double time[], double flux[], double weight[], double freq[],\
             size_t N, size_t M, double power[], double alpha[], double beta[],\
             int useweight);

void fouriermax(double time[], double flux[], double weight[], double freq[],\
                size_t N, size_t M, double *fmax, double *alpmax,\
                double *betmax, int useweight);
