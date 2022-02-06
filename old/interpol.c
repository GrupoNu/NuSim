#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>

/* size of arrays grow in multiples of CHUNK */
#define CHUNK   (100)

int main() {
    int i = 0, N = 0;

    double *x = (double *) malloc(CHUNK*sizeof(double));
    double *y = (double *) malloc(CHUNK*sizeof(double));

    /* getting input */
    while (scanf("%lf%lf", &x[N], &y[N]) != EOF) {
        i++; N++;
        if (i > CHUNK-1) {
            x = (double *) realloc(x, (N+CHUNK)*sizeof(double));
            y = (double *) realloc(y, (N+CHUNK)*sizeof(double));
            i = 0;
        }
    }

    /* resizing the arrays to correct size */
    x = (double *) realloc(x, N * sizeof(double));
    y = (double *) realloc(y, N * sizeof(double));

    double Ndouble = (double) N;
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline_steffen = gsl_spline_alloc(gsl_interp_steffen, N);
    gsl_spline_init(spline_steffen, x, y, N);

    double xi, yi_steffen;
    for (i = 0; i < N; i++) {
        xi = x[0] + (i / Ndouble) * (x[N-1] - x[0]);
        yi_steffen = gsl_spline_eval(spline_steffen, xi, acc);
        printf("%.6e %.6e\n", xi, yi_steffen);
    }

    gsl_spline_free(spline_steffen);
    gsl_interp_accel_free(acc);
    free(x);
    free(y);

    return 0;
}
