#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_errno.h>      /* Error handling */
#include <gsl/gsl_odeiv2.h>     /* Solve ODEs */

/* PARAMETERS */
#define T_INIC  (0.0)
#define T_FINAL (10.0)
#define PASSO   (1e-2)
#define EPS_ABS (1e-2)
#define EPS_REL (1e-3)
#define NUM_IT  ((int) (T_FINAL / PASSO))   /* number of iterations */
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
#define N_E     (M_E)

/* DEFINITIONS */
int func(double t, const double y[], double f[], void *params); /* ODE step function */

int main()
{
    /* initializing parameters */
    double sin2th = pow(sin(THETA), 2), cos2th = pow(cos(THETA), 2);
    double par[8];
    /* Ar = */ par[0] = (M12*cos2th + M22*sin2th)/(2*ENERG) + sqrt(2)*G_F*N_E;  /* Br = */ par[1] = (DM2*sin(2*THETA))/(4*ENERG);
    /* Cr = */ par[2] =                       par[1]                         ;  /* Dr = */ par[3] = (M12*sin2th + M22*cos2th)/(2*ENERG);
    /* Ai = */ par[4] =                        0.0                           ;  /* Bi = */ par[5] =           0.0;
    /* Ci = */ par[6] =                        0.0                           ;  /* Di = */ par[7] =           0.0;

    /* initializing system */       /* jac */
    gsl_odeiv2_system ode_sys = {func, NULL, DIM, par};
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
