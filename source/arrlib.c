#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


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

// Average (mean) of elements in array
double arr_mean(double x[], size_t N)
{
    double sum = 0.0;
    double mean;
    
    for (int i = 0; i < N; ++i) {
        sum += x[i];
    }
    mean = (double) sum / N;

    return mean;
}

// Median of all elements in array
double arr_median(double x[], size_t N)
{
    double median, temp;

    // Copy of array to preserve original
    double* y = malloc(N * sizeof(double));
    memcpy(y, x, N * sizeof(double));

    // Sort array in ascending order
    for (size_t i = 0; i < N-1; ++i) {
        for (size_t j = i+1; j < N; ++j) {
            // Swap elements
            if ( y[j] < y[i] ) {
                temp = y[i];
                y[i] = y[j];
                y[j] = temp;
            }
        }
    }

    // Median depends on even/uneven number of elements
    //  - Even: Mean of the two elements in the middle
    //  - Odd : Element in the middle
    if ( N % 2 == 0 ) {
        median = (y[N/2] + y[N/2 - 1]) / 2.0;
    }
    else {
        median = y[N/2];
    }
    
    // Done
    free(y);
    return median;
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

// Add scalar a to array x -- IN-PLACE
void arr_sca_add(double x[], double a, size_t N)
{
    for (int i = 0; i < N; ++i) {
        x[i] += a;
    }
}

// Find the difference between elements in x and store in y.
// - NB: y should have the length N-1
void arr_diff(double x[], double y[], size_t N)
{
    for (int i = 0; i < N-1; ++i) {
        y[i] = x[i+1] - x[i];
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
