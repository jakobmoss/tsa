void windowfunction(double time[], double freq[], double weight[], size_t N,
                     size_t M, double f0, double window[], int useweight);

double windowsum(double f0, double low, double high, double rate, double time[],
                 double weight[], size_t N, int useweight, int quiet);


