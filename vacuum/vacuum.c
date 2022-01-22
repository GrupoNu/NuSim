#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>      /* Error handling */
#include <gsl/gsl_matrix.h>     /* Matrix for the Jacobian */
#include <gsl/gsl_odeiv2.h>     /* Solve ODE's */

/* PARAMETERS */
#define T_INIC  (0.0)
#define T_FINAL (50.0)
#define PASSO   (1e-2)
#define EPS_ABS (1e-4)
#define EPS_REL (1e-5)
#define NUM_IT  ((int) (T_FINAL / PASSO))   /* numero de iteracoes */
#define DIM     (4)

/* INITIAL CONDITIONS */
#define RE1     (1.0)
#define RE2     (0.0)
#define IM1     (0.0)
#define IM2     (0.0)

/* CONSTANTS */
#define THETA   (M_PI / 6)
#define M12     (1.0)
#define M22     (2.0)
#define DM2     (M22 - M12)
#define ENERG   (1.0)
#define G_F     (3.0)
#define N_E     (0.0)

/* DEFINITIONS */
int func(double t, const double y[], double f[], void *params); /* ODE step function */
int jac(double t, const double y[], double *dfdy, double dfdt[], void *params); /* jacobian */

int main()
{
    /* initializing parameters */
    double sin2th = pow(sin(THETA), 2), cos2th = pow(cos(THETA), 2);
    double par[8];
    /* Ar = */ par[0] = (M12*cos2th + M22*sin2th)/(2*ENERG) + sqrt(2)*G_F*N_E;  /* Br = */ par[1] = (DM2*sin(2*THETA))/(4*ENERG);
    /* Cr = */ par[2] =                       par[1]                         ;  /* Dr = */ par[3] = (M12*sin2th + M22*cos2th)/(2*ENERG);
    /* Ai = */ par[4] =                        0.0                           ;  /* Bi = */ par[5] =           0.0;
    /* Ci = */ par[6] =                        0.0                           ;  /* Di = */ par[7] =           0.0;

    /* initializing system */
    gsl_odeiv2_system ode_sys = {func, jac, DIM, par};
    /* driver with RFK45 */
    gsl_odeiv2_driver *driver =
        gsl_odeiv2_driver_alloc_y_new(&ode_sys, gsl_odeiv2_step_rkf45, PASSO, EPS_ABS, EPS_REL);
    /* initial condition */
    double t = T_INIC;
    double y[DIM] = { RE1, RE2, IM1, IM2 };

    int i; double ti;
    for (i = 1; i <= NUM_IT; i++) {
        ti = i * PASSO + T_INIC;
        int status = gsl_odeiv2_driver_apply(driver, &t, ti, y);

        if (status != GSL_SUCCESS) {
            printf ("Error, return value = %d\n", status);
            break;
        }

        printf("%d %.5e %.5e %.5e %.5e %.5e\n", i, t, y[0], y[1], y[2], y[3]);
    }

    gsl_odeiv2_driver_free(driver);
    return 0;
}

/* dy/dt = f(y, t) */
int func(double t, const double y[], double f[], void *params) {
    double *par = (double *) params;
            /* Ai               Bi                  Ar              Br */
    f[0] =   par[4] * y[0]  +  par[5] * y[1]  +  par[0] * y[2]  +  par[1] * y[3];

            /* Ci               Di                  Cr              Dr */
    f[1] =   par[6] * y[0]  +  par[7] * y[1]  +  par[2] * y[2]  +  par[3] * y[3];

            /* Ar               Br                  Ai              Bi */
    f[2] = - par[0] * y[0]  -  par[1] * y[1]  +  par[4] * y[2]  +  par[5] * y[3];

            /* Cr               Dr                  Ci              Di */
    f[3] = - par[2] * y[0]  -  par[3] * y[1]  +  par[6] * y[2]  +  par[7] * y[3];

    return GSL_SUCCESS;
}

/* jacobian: J_{ij} = df_i(t, y(t)) / dy_j */
int jac(double t, const double y[], double *dfdy, double dfdt[], void *params) {
    double *par = (double *) params;
    gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, DIM, DIM);
    gsl_matrix *mat = &dfdy_mat.matrix;
                            /* Ai                                   Bi                              Ar                                  Br */
    gsl_matrix_set(mat, 0, 0, par[4]); gsl_matrix_set(mat, 0, 1, par[5]); gsl_matrix_set(mat, 0, 2, par[0]); gsl_matrix_set(mat, 0, 3, par[1]);
                            /* Ci                                   Di                              Cr                                  Dr */
    gsl_matrix_set(mat, 1, 0, par[6]); gsl_matrix_set(mat, 1, 1, par[7]); gsl_matrix_set(mat, 1, 2, par[2]); gsl_matrix_set(mat, 1, 3, par[3]);
                            /* Ar                                   Br                              Ai                                  Bi */
    gsl_matrix_set(mat, 2, 0, par[0]); gsl_matrix_set(mat, 2, 1, par[1]); gsl_matrix_set(mat, 2, 2, par[4]); gsl_matrix_set(mat, 2, 3, par[5]);
                            /* Cr                                   Dr                              Ci                                  Di */
    gsl_matrix_set(mat, 3, 0, par[2]); gsl_matrix_set(mat, 3, 1, par[3]); gsl_matrix_set(mat, 3, 2, par[6]); gsl_matrix_set(mat, 3, 3, par[7]);
    dfdt[0] = 0.0; dfdt[1] = 0.0; dfdt[2] = 0.0; dfdt[3] = 0.0;
    return GSL_SUCCESS;
}
