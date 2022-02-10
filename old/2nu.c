#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>      /* Error handling */
#include <gsl/gsl_odeiv2.h>     /* Solve ODEs */
#include <gsl/gsl_spline.h>     /* Interpolation of electron density */

/* EDO parameters */
#define T_INIC  (0.0)
#define T_FINAL (1.0)
#define PASSO   (1e-2)
#define EPS_ABS (1e-2)
#define EPS_REL (1e-3)
#define NUM_IT  ((int) ((T_FINAL - T_INIC) / PASSO))    /* number of iterations */
#define GERAC   (2)             /* number of generations, 2 or 3 */
#define DIM     (2*GERAC)       /* dimension of the problem */
#define PARNUM  (2*GERAC*GERAC) /* number of hamiltonian entries */

/* initial conditions */
#define RE1     (1.0)
#define RE2     (0.0)
#define IM1     (0.0)
#define IM2     (0.0)

/* CONSTANTS */
/* from: http:www.nu-fit.org/ */
#define THETA   (M_PI / 6)
#define M12     (1.0)
#define M22     (2.0)
#define DM2     (M22 - M12)
#define ENERG   (1.0)
#define G_F     (3.0)
#define KNE     (M_SQRT2*G_F)   /* constant that multiplies Ne */

/* parameters structure */
typedef struct {
    gsl_interp_accel *acc;
    gsl_spline *spline;
    int num;
    double *real;
} par;

/* DEFINITIONS */
double Ne(double t, void *params);  /* electron density */
int func(double t, const double Psi[], double f[], void *params); /* ODE step function */
int readalloc(double **r, double **logNe);  /* read and alloc data, logNe is "ln(N_e)" */

int main()
{
//d                           /*************************/
//d                           /***     DEBUGGING     ***/
//d                           /*************************/
//d
//d   /* initializing parameters */
//d   double sin2th = pow(sin(THETA), 2), cos2th = pow(cos(THETA), 2);
//d   double real[PARNUM];    /* real matrix */
//d   /* Ar = */ real[0] = (M12*cos2th + M22*sin2th)/(2*ENERG);  /* Br = */ real[1] = (DM2*sin(2*THETA))/(4*ENERG);
//d   /* Cr = */ real[2] =                      real[1]       ;  /* Dr = */ real[3] = (M12*sin2th + M22*cos2th)/(2*ENERG);
//d   /* Ai = */ real[4] =                        0.0         ;  /* Bi = */ real[5] =           0.0;
//d   /* Ci = */ real[6] =                        0.0         ;  /* Di = */ real[7] =           0.0;
//d   printf("\nmatrix\n\n%.3lf %.3lf\n%.3lf %.3lf\n", real[0], real[1], real[2], real[3]);

    double *x, *logNe;
    /* getting data and resizing appropriately */
    int N = readalloc(&x, &logNe);
    /* initializing gsl structures */
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_steffen, N);
    gsl_spline_init(spline, x, logNe, N);

    /* initializing parameters */
    double sin2th = pow(sin(THETA), 2), cos2th = pow(cos(THETA), 2);
    double real[PARNUM];    /* real matrix */
    /* Ar = */ real[0] = (M12*cos2th + M22*sin2th)/(2*ENERG);  /* Br = */ real[1] = (DM2*sin(2*THETA))/(4*ENERG);
    /* Cr = */ real[2] =                      real[1]       ;  /* Dr = */ real[3] = (M12*sin2th + M22*cos2th)/(2*ENERG);
    /* Ai = */ real[4] =                        0.0         ;  /* Bi = */ real[5] =           0.0;
    /* Ci = */ real[6] =                        0.0         ;  /* Di = */ real[7] =           0.0;
    par params = { acc, spline, PARNUM, real };

    /* initializing system */
    gsl_odeiv2_system ode_sys = {func, NULL, DIM, &params};
    gsl_odeiv2_driver *driver =
        gsl_odeiv2_driver_alloc_y_new(&ode_sys, gsl_odeiv2_step_rkf45, PASSO, EPS_ABS, EPS_REL);    /* RKF45 method */
    /* initial conditon */
    double t = T_INIC;
    double Psi[DIM] = { RE1, RE2, IM1, IM2 };

    int i; double ti;
    for (i = 0; i <= NUM_IT; i++) {
        ti = i * PASSO + T_INIC;
        int status = gsl_odeiv2_driver_apply(driver, &t, ti, Psi);

        if (status != GSL_SUCCESS) {
            printf ("Error, return value = %d\n", status);
            break;
        }

        printf("%d  %.5e  %.5e  %.5e  %.5e  %.5e\n", i, t, Psi[0], Psi[1], Psi[2], Psi[3]);
    }

    gsl_odeiv2_driver_free(driver);
    gsl_spline_free(spline);
    gsl_interp_accel_free(acc);
    free(x); free(logNe);
    return 0;
}

/* calculates Ne = exp(logNe), where logNe was interpolated */
double Ne(double t, void *params) {
    par *param = (par *) params;
    return exp(gsl_spline_eval(param->spline, t, param->acc));
}

/* dPsi/dt = f(Psi, t) */
int func(double t, const double Psi[], double f[], void *params) {
    par *param = (par *) params;
    double *real = param->real;
          /*  Ar                                        Br */
    f[0] =  (real[0] + KNE*Ne(t, param)) * Psi[2]   +  real[1] * Psi[3];

          /*  Cr                                        Dr */
    f[1] =  real[2] * Psi[2]                        +  real[3] * Psi[3];

          /* - Ar                                       - Br */
    f[2] = - (real[0] + KNE*Ne(t, param)) * Psi[0]  -  real[1] * Psi[1];

          /* - Cr                                       - Dr */
    f[3] = - real[2] * Psi[0]                       -  real[3] * Psi[1];

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
