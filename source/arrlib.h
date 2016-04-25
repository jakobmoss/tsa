void arr_init_zero(double x[], size_t N);

void arr_init_crange(double x[], size_t N);

void arr_init_frange(double x[], size_t N);

void arr_init_linspace(double x[], double a, double rate, size_t N);


double arr_sum(double x[], size_t N);

double arr_sumsq(double x[], size_t N);

double arr_mean(double x[], size_t N);


void arr_cos(double x[], double y[], size_t N);

void arr_sin(double x[], double y[], size_t N);

void arr_scale(double x[], double a, double y[], size_t N);

void arr_sca_add(double x[], double a, double y[], size_t N);


void arr_mult(double x[], double y[], double z[], size_t N);

void arr_add(double x[], double y[], double z[], size_t N);

void arr_addsq(double x[], double y[], double z[], size_t N);

size_t arr_util_getstep(double a, double b, double rate);

