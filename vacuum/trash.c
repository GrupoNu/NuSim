#include <gsl/gsl_matrix.h>     /* Matrix for the Jacobian */

int jac(double t, const double y[], double *dfdy, double dfdt[], void *params); /* jacobian */

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
