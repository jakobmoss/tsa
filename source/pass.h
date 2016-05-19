void bandpass(double time[], double flux[], double weight[], size_t N,\
              double f1, double f2, double low, double high, double rate,\
              double result[], int useweight, int quiet);

void lowpass(double time[], double flux[], double weight[], size_t N,\
             double flow, double low, double high, double rate,       \
             double result[], int useweight, int quiet);

void highpass(double time[], double flux[], double weight[], size_t N,\
              double fhigh, double low, double high, double rate,     \
              double result[], int useweight, int quiet);
