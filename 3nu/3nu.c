/* headers */
#include <stdio.h>                  /* get input and print output */
#include <stdlib.h>                 /* free, alloc, NULL pointer */
#include <math.h>                   /* basic math functions: exp */
#include <gsl/gsl_complex_math.h>   /* complex numbers */
#include <gsl/gsl_matrix.h>         /* matrix definitions */
#include <gsl/gsl_blas.h>           /* basic linear algebra operations */
#include <gsl/gsl_errno.h>          /* error handling: GSL_SUCCESS */
#include <gsl/gsl_odeiv2.h>         /* solve ODEs */
#include <gsl/gsl_spline.h>         /* interpolation of real data */

/* EDO parameters */
#define T_INIC  (0.02)
#define T_FINAL (1.0)
#define PASSO   (1e-2)
#define EPS_ABS (1e-2)
#define EPS_REL (1e-3)
#define NUM_IT  ((int) ((T_FINAL - T_INIC) / PASSO))    /* number of iterations */
#define GERAC   (3)             /* number of neutrinos ( n = 3 ) */
#define DIM     (2*GERAC)       /* dimension of the problem ( 2 n ) */

/* initial condition */
#define RE1     (1.0)
#define RE2     (0.0)
#define RE3     (0.0)
#define IM1     (0.0)
#define IM2     (0.0)
#define IM3     (0.0)

/* CONSTANTS */
#define G_F     (3.0)

/* MACROS */
#define MAT(M, i, j, z)     gsl_matrix_complex_set((M), (i), (j), (z))
#define CARTE(x, y)         gsl_complex_rect((x), (y))
#define POLAR(r, t)         gsl_complex_polar((r), (t))
#define C_ONE               (CARTE(1.0, 0.0))
#define C_ZERO              (CARTE(0.0, 0.0))
#define H0re(i, j)          gsl_matrix_get((H0_re), (i), (j))
#define H0im(i, j)          gsl_matrix_get((H0_im), (i), (j))

/* parameters structure */
typedef struct {
    gsl_interp_accel *acc;
    gsl_spline *spline;
    gsl_matrix *H0_re;
    gsl_matrix *H0_im;
} par;

/* DEFINITIONS */
double D(double t, void *params);   /* sqrt(2) * G_F * N_e, where N_e is the electron density */
void action(const gsl_matrix_complex *A, const gsl_matrix_complex *B, gsl_matrix_complex *C);
int func(double t, const double Psi[], double f[], void *params); /* ODE step function */
int readalloc(double **r, double **logNe);  /* read and alloc data, logNe is "ln(N_e)" */

int main() {
    /* indexes for matrices */
    int i, j;

    /* getting data and resizing properly */
    double *x, *logNe;
    int N = readalloc(&x, &logNe);

    /* initializing gsl structures for interpolation */
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_steffen, N);
    gsl_spline_init(spline, x, logNe, N);

    /* defining the mixing matrix */
    gsl_complex U11, U12, U13,
                U21, U22, U23,          /* here we put the actual data from NuFit */
                U31, U32, U33;
    gsl_matrix_complex *U = gsl_matrix_complex_alloc(GERAC, GERAC);
    MAT(U, 0, 0, U11);          MAT(U, 0, 1, U12);          MAT(U, 0, 2, U13);
    MAT(U, 1, 0, U21);          MAT(U, 1, 1, U22);          MAT(U, 1, 2, U23);
    MAT(U, 2, 0, U31);          MAT(U, 2, 1, U32);          MAT(U, 2, 2, U33);

    /* diagonal matrix for the masses */
    double d1, d2, d3;  /* the diagonal elements */
    gsl_matrix_complex *M = gsl_matrix_complex_alloc(GERAC, GERAC);
    gsl_matrix_complex_set_zero(M);
    MAT(M, 1, 1, CARTE(d1, 0));
                                MAT(M, 2, 2, CARTE(d2, 0));
                                                            MAT(M, 3, 3, CARTE(d3, 0));

    /* constant matrix is such that H = h0 + diag(D(t), 0, 0) */
    gsl_matrix_complex *h0 = gsl_matrix_complex_alloc(GERAC, GERAC);
    action(U, M, h0);    /* h0 = U . M . U^dagger */

    /* obtaining the matrices H0_re, H0_im such that h0 = H0_re + sqrt(-1) . H0_im */
    gsl_matrix *H0_re = gsl_matrix_alloc(GERAC, GERAC), *H0_im = gsl_matrix_alloc(GERAC, GERAC);
    for (i = 0; i < GERAC; i++) {
        for (j = 0; j < GERAC; j++) {
            gsl_matrix_set(H0_re, i, j, GSL_REAL(gsl_matrix_complex_get(h0, i, j)));
            gsl_matrix_set(H0_im, i, j, GSL_IMAG(gsl_matrix_complex_get(h0, i, j)));
        }
    }

    /* parameters data structure */
    par params = { acc, spline, H0_re, H0_im };

    /* initializing system */
    gsl_odeiv2_system ode_sys = {func, NULL, DIM, &params};
    gsl_odeiv2_driver *driver =
        gsl_odeiv2_driver_alloc_y_new(&ode_sys, gsl_odeiv2_step_rkf45, PASSO, EPS_ABS, EPS_REL);    /* Runge-Kutta-Fehlberg 45 method */

    /* initial conditon */
    double t = T_INIC;
    double Psi[DIM] = { RE1, RE2, RE3, IM1, IM2, IM3 };

                            /*************************/
                            /*    ODE CALCULATION    */
                            /*************************/
    double ti;
    for (i = 1; i <= NUM_IT; i++) {
        ti = i * PASSO + T_INIC;
        int status = gsl_odeiv2_driver_apply(driver, &t, ti, Psi);

        if (status != GSL_SUCCESS) {
            printf ("Error, return status = %d\n", status);
            break;
        }
        printf("%d %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n", i, t, Psi[0], Psi[1], Psi[2], Psi[3], Psi[4], Psi[5]);
    }

                            /*************************/
                            /* FREEING THE RESOURCES */
                            /*************************/
    /* interpolation structures */
    gsl_odeiv2_driver_free(driver); gsl_spline_free(spline); gsl_interp_accel_free(acc);
    /* data */
    free(x); free(logNe);
    /* complex matrices */
    gsl_matrix_complex_free(U); gsl_matrix_complex_free(M); gsl_matrix_complex_free(h0);
    /* real matrices */
    gsl_matrix_free(H0_re); gsl_matrix_free(H0_im);

    return 0;
}

/* evaluates C = A B A^dagger */
void action(const gsl_matrix_complex *A, const gsl_matrix_complex *B, gsl_matrix_complex *C) {
    gsl_matrix_complex *A_mult_B = gsl_matrix_complex_alloc(A->size1, B->size2);
    /* zgemm does: C = alpha op(A) op(B) + beta C */
                                             /* alpha     beta */
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, C_ONE, A, B, C_ZERO, A_mult_B);      /* A_mult_B = A . B */
    gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, C_ONE, A_mult_B, A, C_ZERO, C);    /* C = A . B . A^dagger */

    gsl_matrix_complex_free(A_mult_B);  /* freeing the memory of the auxiliary variable */
}

/* calculates D = sqrt(2) . G_F . N_e, where N_e = exp(log N_e) is interpolated */
double D(double t, void *params) {
    par *param = (par *) params;
    return M_SQRT2 * G_F * exp(gsl_spline_eval(param->spline, t, param->acc));
}

/* dPsi/dt = f(Psi, t) */
int func(double t, const double Psi[], double f[], void *params) {
    par *param = (par *) params;
    gsl_matrix *H0_re = param->H0_re;
    gsl_matrix *H0_im = param->H0_im;

                                /* notice: the DENSITY TERM always corresponds to R11 = H0re(0, 0) */
/*             J11                 J12                 J13                   R11                 R12                 R13                DENSITY TERM  */
f[0] =   H0im(0, 0)*Psi[0] + H0im(0, 1)*Psi[1] + H0im(0, 2)*Psi[2]  +  H0re(0, 0)*Psi[3] + H0re(0, 1)*Psi[4] + H0re(0, 2)*Psi[5]  -  D(t, param)*Psi[3];
/*             J21                 J22                 J23                   R21                 R22                 R23       */
f[1] =   H0im(1, 0)*Psi[0] + H0im(1, 1)*Psi[1] + H0im(1, 2)*Psi[2]  +  H0re(1, 0)*Psi[3] + H0re(1, 1)*Psi[4] + H0re(1, 2)*Psi[5];
/*             J31                 J32                 J33                   R31                 R32                 R33       */
f[2] =   H0im(2, 0)*Psi[0] + H0im(2, 1)*Psi[1] + H0im(2, 2)*Psi[2]  +  H0re(2, 0)*Psi[3] + H0re(2, 1)*Psi[4] + H0re(2, 2)*Psi[5];

/*            -R11                -R12                -R13                   J11                 J12                 J13                DENSITY TERM  */
f[3] = - H0re(0, 0)*Psi[0] - H0re(0, 1)*Psi[1] - H0re(0, 2)*Psi[2]  +  H0im(0, 0)*Psi[3] + H0im(0, 1)*Psi[4] + H0im(0, 2)*Psi[5]  -  D(t, param)*Psi[0];
/*            -R21                -R22                -R23                   J21                 J22                 J23       */
f[4] = - H0re(1, 0)*Psi[0] - H0re(1, 1)*Psi[1] - H0re(1, 2)*Psi[2]  +  H0im(1, 0)*Psi[3] + H0im(1, 1)*Psi[4] + H0im(1, 2)*Psi[5];
/*            -R31                -R32                -R33                   J31                 J32                 J33       */
f[5] = - H0re(2, 0)*Psi[0] - H0re(2, 1)*Psi[1] - H0re(2, 2)*Psi[2]  +  H0im(2, 0)*Psi[3] + H0im(2, 1)*Psi[4] + H0im(2, 2)*Psi[5];


    return GSL_SUCCESS;
}

/* reads data from standard input, storages them in x_ptr and y_ptr and returns N */
int readalloc(double **x_ptr, double **y_ptr) {
    int i = 0, N = 0;
    int chunk = 100;
    double *x = (double *) malloc(chunk * sizeof(double));
    double *y = (double *) malloc(chunk * sizeof(double));

    /* getting input */
    while (scanf("%lf%lf", &x[N], &y[N]) != EOF) {
        i++; N++;
        if (i > chunk-1) {
            x = (double *) realloc(x, (N+chunk)*sizeof(double));
            y = (double *) realloc(y, (N+chunk)*sizeof(double));
            i = 0;
        }
    }

    /* resizing the arrays to correct size */
    *x_ptr = (double *) realloc(x, N * sizeof(double));
    *y_ptr = (double *) realloc(y, N * sizeof(double));
    return N;
}
