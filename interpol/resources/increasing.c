#include <stdio.h>
#include <stdlib.h>

/* size of arrays grow in multiples of CHUNK */
#define CHUNK   (100)

int main() {
    int i = 0, N = 0;
    double u = -100000;

    double *x = (double *) malloc(CHUNK*sizeof(double));
    double *y = (double *) malloc(CHUNK*sizeof(double));

    /* getting input */
    while (scanf("%lf%lf", &x[N], &y[N]) != EOF) {
        if (x[N] <= u) {
            printf("parou em N = %d\n", N);
            break;
        }
        u = x[N];
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

    free(x);
    free(y);
    printf("N = %d\n", N);

    return 0;
}
