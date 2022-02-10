#include <gsl/gsl_matrix.h>     /* Matrix for the Jacobian */
#include <gsl/gsl_deriv.h>      /* Numerical differentiation of N_e */

int jac(double t, const double Psi[], double *dfdy, double dfdt[], void *params); /* jacobian of the system */

/* jacobiano: J_{ij} = df_i(t, Psi(t)) / dPsi_j e J_t = df_i/dt */
int jac(double t, const double Psi[], double *dfdPsi, double dfdt[], void *params) {
    par *param = (par *) params;
    double *real = param->real;
    gsl_matrix_view dfdPsi_mat = gsl_matrix_view_array(dfdPsi, DIM, DIM);
    gsl_matrix *mat = &dfdPsi_mat.matrix;
                       /*    Ai                               Bi                                      Ar                                        Br */
    gsl_matrix_set(mat,0,0,real[4]); gsl_matrix_set(mat,0,1,real[5]); gsl_matrix_set(mat,0,2,real[0]+KNE*Ne(t, param)); gsl_matrix_set(mat,0,3,real[1]);
                       /*    Ci                               Di                                      Cr                                        Dr */
    gsl_matrix_set(mat,1,0,real[6]); gsl_matrix_set(mat,1,1,real[7]);         gsl_matrix_set(mat,1,2,real[2]);                  gsl_matrix_set(mat,1,3,real[3]);
                       /*  - Ar                             - Br                                      Ai                                        Bi */
    gsl_matrix_set(mat,2,0,-real[0]-KNE*Ne(t, param)); gsl_matrix_set(mat,2,1,-real[1]); gsl_matrix_set(mat,2,2,real[4]); gsl_matrix_set(mat,2,3,real[5]);
                       /*  - Cr                             - Dr                                      Ci                                        Di */
    gsl_matrix_set(mat,3,0,-real[2]); gsl_matrix_set(mat,3,1,-real[3]);         gsl_matrix_set(mat,3,2,real[6]);          gsl_matrix_set(mat,3,3,real[7]);

    gsl_function gslNe;
    gslNe.function = Ne;
    gslNe.params = param;
    double derivNe, err;
    gsl_deriv_central(&gslNe, t, 1e-8, &derivNe, &err);

    dfdt[0] = KNE * derivNe * Psi[2];
    dfdt[1] = 0.0;
    dfdt[2] = - KNE * derivNe * Psi[0];
    dfdt[3] = 0.0;
    return GSL_SUCCESS;
}
