/*  ~~~ Time Series Analysis -- Auxiliary ~~~
 *
 * Routines for minimisation of scalar functions
 * 
 * Author: Jakob Rørsted Mosumgaard
 */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

// Maximum number of iterations
#define MAXITER 100

// Golden ratio and (1 - golden ratio)
#define GOLD 0.6180339887498948482046
#define IGOLD 0.3819660112501051517954


/* Find minimum using a golden section search
 *
 * Arguments:
 * - `f`: Function to find the minimum of.
 * - `a`, `b`: Interval containing the minimum (such that a < b).
 * - `eps`: Desired tolerance
 * - `xmin`: OUTPUT -- location of minimum
 *
 * The return value is function is f(xmin).
 * */
double fmin_golden(double (*f)(double), double a, double b, double eps,\
                   double *xmin)
{
    // Initialisations
    double fx;
    int success = 0;
    
    // Initial interval
    double x1 = GOLD*a + IGOLD*b;
    double x2 = IGOLD*a + GOLD*b;
    double fx1 = (*f)(x1);
    double fx2 = (*f)(x2);

    // Main loop
    for (int i = 0; i < MAXITER; ++i) {
        // Case 1: Reduce interval to [a, x2]
        if ( fx1 < fx2 ) {
            b = x2;
            x2 = x1;
            fx2 = fx1;
            x1 = GOLD*a + IGOLD*b;
            fx1 = (*f)(x1);
        }
        // Case 2: Reduce interval to [x1, b]
        else {
            a = x1;
            x1 = x2;
            fx1 = fx2;
            x2 = IGOLD*a + GOLD*b;
            fx2 = (*f)(x2);
        }

        // Check for convergence
        if ( fabs(b-a) < eps ) {
            success = 1;
            break;
        }
    }

    // Signal if no minimum is found within desired number of iterations
    if ( success == 0 ) {
        fprintf(stderr, "fmin: Accuracy not reached in %i interations!",\
                MAXITER);
        exit(1);
    }
    else {
        // Save midpoint and return function value
        *xmin = a + (b-a)/2;
        fx = (*f)(*xmin);
        return fx;
    }
}

