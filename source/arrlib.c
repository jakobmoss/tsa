#include <stdio.h>
#include <math.h>

/* ~~~~~ Initialisations ~~~~~ */

// Fill an array with zeros
void arr_init_zero(double x[], size_t N)
{
    for (int i = 0; i < N; ++i) {
        x[i] = 0;
    }
}

// Fill an array with numbers: 0, ..., N-1 ; (C-style range)
void arr_init_crange(double x[], size_t N)
{
    for (int i = 0; i < N; ++i) {
        x[i] = i;
    }
}

// Fill an array with numbers: 1, ..., N ; (FORTRAN-style range)
void arr_init_frange(double x[], size_t N)
{
    for (int i = 0; i < N; ++i) {
        x[i] = i+1;
    }
}

// Fill array x with N points of increment rate
//  - NOTE: Use utility function to calculate N
void arr_init_linspace(double x[], double a, double rate, size_t N) 
{
    for (int i = 0; i < N; ++i) {
        x[i] = a + i*rate;
    }
}



/* ~~~~~ Return value for single array ~~~~~ */

// Sum the elements of an array
double arr_sum(double x[], size_t N)
{
    double sum = 0.0;
    for (int i = 0; i < N; ++i) {
        sum += x[i];
    }
    return sum;
}

// Sum the squared elements of an array
double arr_sumsq(double x[], size_t N)
{
    double sum = 0.0;
    for (int i = 0; i < N; ++i) {
        sum += x[i]*x[i];
    }
    return sum;
}



/* ~~~~~ Array operations on single array ~~~~~ */

// Take the cos of elements x and store in y.
void arr_cos(double x[], double y[], size_t N)
{
    for (int i = 0; i < N; ++i) {
        y[i] = cos(x[i]);
    }
}

// Take the sin of elements x and store in y.
void arr_sin(double x[], double y[], size_t N)
{
    for (int i = 0; i < N; ++i) {
        y[i] = sin(x[i]);
    }
}

// Multiply array x by scalar a and store in y
void arr_scale(double x[], double a, double y[], size_t N)
{
    for (int i = 0; i < N; ++i) {
        y[i] = x[i] * a;
    }
}



/* ~~~~~ Array operations on two arrays ~~~~~ */

// Multiply arrays x, y and store in z
void arr_mult(double x[], double y[], double z[], size_t N)
{
    for (int i = 0; i < N; ++i) {
        z[i] = x[i] * y[i];
    }
}


// Add arrays x, y and store in z
void arr_add(double x[], double y[], double z[], size_t N)
{
    for (int i = 0; i < N; ++i) {
        z[i] = x[i] + y[i];
    }
}


// Add the square of the arrays x, y and store in z
void arr_addsq(double x[], double y[], double z[], size_t N)
{
    for (int i = 0; i < N; ++i) {
        z[i] = x[i]*x[i] + y[i]*y[i];
    }
}



/* ~~~~~ Mixed utilities ~~~~~ */

// FOR LINSPACE: Calculate number of steps required
size_t arr_util_getstep(double a, double b, double rate)
{
    size_t steps = 0;
    double val = a;
    do {
        val += rate;
        steps++;
    } while (val < b);
    return steps-1;
}
