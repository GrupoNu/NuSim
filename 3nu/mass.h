/* standard libraries */
#include <stdio.h>                  /* get input and print output */
#include <stdlib.h>                 /* free, alloc, NULL pointer */
#include <math.h>                   /* basic math functions: exp */

/* GSL */
#include <gsl/gsl_complex.h>        /* definition of complex numbers */
#include <gsl/gsl_complex_math.h>   /* complex math operations */
#include <gsl/gsl_vector.h>         /* vector definitions */
#include <gsl/gsl_matrix.h>         /* matrix definitions */
#include <gsl/gsl_blas.h>           /* basic linear algebra operations */
#include <gsl/gsl_errno.h>          /* error handling: GSL_SUCCESS */
#include <gsl/gsl_odeiv2.h>         /* solve ODEs */
#include <gsl/gsl_spline.h>         /* interpolation of real data */

/* function macros */
#define MAT(M, i, j, z)     gsl_matrix_complex_set((M), (i), (j), (z))
#define REAL(x)             gsl_complex_rect((x), 0.0)
#define POLAR(r, t)         gsl_complex_polar((r), (t))
#define H0re(i, j)          gsl_matrix_get((H0_re), (i), (j))
#define H0im(i, j)          gsl_matrix_get((H0_im), (i), (j))
#define DEG_TO_RAD(ang)     ((M_PI / 180.0) * ang)
#define phi_R(i)            GSL_REAL(gsl_vector_complex_get((phi), (i)))
#define phi_I(i)            GSL_IMAG(gsl_vector_complex_get((phi), (i)))

/* constant macros */
#define GERAC   (3)
#define NUM_IT  ((int) ((T_FINAL - T_INIC) / PASSO))    /* number of iterations */
#define DIM     (2*GERAC)       /* dimension of the problem ( 2 n ) */

/* structures */
typedef struct {
    gsl_interp_accel *acc;
    gsl_spline *spline;
    gsl_matrix *H0_re;
    gsl_matrix *H0_im;
} par;  /* parameters */

/* functions */
int func(double t, const double Psi[], double f[], void *params); /* ODE step function */
double D(double t, void *params);   /* sqrt(2) * G_F * N_e, where N_e is the electron density */
void action(const gsl_matrix_complex *A, const gsl_matrix_complex *B, gsl_matrix_complex *C);   /* C = A B A^dagger */
int readalloc(double **r, double **logNe, int chunk);  /* read and alloc data, logNe is "ln(N_e)" */
void printMatrix(gsl_matrix **M);
