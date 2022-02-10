/* standard libraries */
#include <stdio.h>                  /* get input and print output */
#include <stdlib.h>                 /* free, malloc, NULL pointer */
#include <math.h>                   /* basic math functions: exp */

/* GSL */
#include <gsl/gsl_complex.h>        /* definition of complex numbers */
#include <gsl/gsl_complex_math.h>   /* complex math operations */
#include <gsl/gsl_vector.h>         /* vector definitions */
#include <gsl/gsl_matrix.h>         /* matrix definitions */
#include <gsl/gsl_blas.h>           /* basic linear algebra operations */
#include <gsl/gsl_eigen.h>          /* solve eigensystems */
#include <gsl/gsl_errno.h>          /* error handling: GSL_SUCCESS */
#include <gsl/gsl_odeiv2.h>         /* solve ODEs */
#include <gsl/gsl_spline.h>         /* interpolation of real data */

/* function macros */
#define MAT(M, i, j, z)     gsl_matrix_complex_set((M), (i), (j), (z))
#define REAL(x)             gsl_complex_rect((x), 0.0)
#define POLAR(r, t)         gsl_complex_polar((r), (t))
#define ADD(a, b)           gsl_complex_add_real((a), (b))
#define H0re(i, j)          gsl_matrix_get((H0_re), (i), (j))
#define H0im(i, j)          gsl_matrix_get((H0_im), (i), (j))
#define DEG_TO_RAD(ang)     ((M_PI / 180.0) * ang)
#define phi_R(j)            GSL_REAL(gsl_vector_complex_get((phi), (j)))
#define phi_I(j)            GSL_IMAG(gsl_vector_complex_get((phi), (j)))
#define psi_R(j)            GSL_REAL(gsl_vector_complex_get((psi), (j)))
#define psi_I(j)            GSL_IMAG(gsl_vector_complex_get((psi), (j)))

/* constant macros */
#define GERAC   (3)     /* numero de geracoes = 3 para fazer o calculo. NAO mudar */
#define NUM_IT  ((int) ((T_FINAL - T_INIC) / PASSO))    /* number of iterations */
#define DIM     (2*GERAC)       /* dimension of the problem ( 2 n ) */
#define MASS    (1)             /* for mass neutrinos */
#define PART    (0)             /* for interaction neutrinos */

/* structures */
typedef struct {
    gsl_interp_accel *acc;
    gsl_spline *spline;
    gsl_matrix_complex *h0;
    gsl_matrix *H0_re, *H0_im;
} par;  /* parameters */

typedef struct {
    gsl_matrix_complex *H, *eigvec;
    gsl_vector *eigval;
    gsl_eigen_hermv_workspace *workspace;
} eigen_problem;    /* to solve eigensystems */

typedef struct {
    gsl_vector_complex *psi, *phi;
} complex_vecs;     /* to change between mass and interaction basis */

/* functions */
par *genAllocParams(double *x, double *logNe, int N);  /* generate and alloc params */
void freeParams(par *params);   /* frees params allocated by genAllocParams */
int func(double t, const double Psi[], double f[], void *params); /* ODE step function */
double D(double t, void *params);   /* gsl_function: sqrt(2) * G_F * N_e, where N_e is the electron density */
void solveEigensys(double t, par *params, eigen_problem *eig_prob);  /* solves eigen_problem */
eigen_problem *allocEigensys();  /* generate and alloc eig_prob */
void freeEigensys(eigen_problem *eig_prob); /* frees what allocEigensys allocated */
complex_vecs *allocVecs();                  /* alloc complex vectors */
void freeVecs(complex_vecs *cplx_vecs);     /* frees memory allocated by freeVecs */
void action(const gsl_matrix_complex *A, const gsl_matrix_complex *B, gsl_matrix_complex *C);   /* C = A B A^dagger */
void multMatrixVec(CBLAS_TRANSPOSE_t TransA, const gsl_matrix_complex *A, gsl_vector_complex *x, gsl_vector_complex *y);    /* y = op(A) x */
int readalloc(double **r, double **logNe, int chunk);  /* read and alloc data, logNe is "ln(N_e)" */
