#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

void quicksort(double x[], size_t first, size_t last);


/* ~~~~~ Initialisations ~~~~~ */

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

// Fill array x with N points of increment rate starting from a
//  - NOTE: Uses utility function to calculate N
void arr_init_linspace(double x[], double a, double rate, size_t N) 
{
    for (int i = 0; i < N; ++i) {
        x[i] = a + i*rate;
    }
}



/* ~~~~~ Return value for single array ~~~~~ */

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

// Median of array x
double arr_median(double x[], size_t N)
{
    double median;

    // Copy of array to preserve original
    double* y = malloc(N * sizeof(double));
    memcpy(y, x, N * sizeof(double));

    // Sort array in ascending order using quicksort
    quicksort(y, 0, N-1);

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

// FOR MEDIAN: Sort array using (recursive) 'quicksort' algorithm
void quicksort(double x[], size_t first, size_t last)
{
    double temp;
    size_t pivot, i, j;
    
    if ( first < last ){
        pivot = first;
        i = first;
        j = last;

        while( i < j ){
            while( x[i] <= x[pivot] && i < last )
                i++;
            while( x[j] > x[pivot] )
                j--;
            if( i < j ){
                temp = x[i];
                x[i] = x[j];
                x[j] = temp;
            }
        }

        temp = x[pivot];
        x[pivot] = x[j];
        x[j] = temp;
        quicksort(x, first, j-1);
        quicksort(x, j+1, last);
    }
}


