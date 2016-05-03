void arr_init_crange(double x[], size_t N);

void arr_init_frange(double x[], size_t N);

void arr_init_linspace(double x[], double a, double rate, size_t N);


double arr_mean(double x[], size_t N);

double arr_median(double x[], size_t N);


void arr_sca_add(double x[], double a, size_t N);

void arr_diff(double x[], double y[], size_t N);


size_t arr_util_getstep(double a, double b, double rate);
