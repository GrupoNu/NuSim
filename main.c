#include "main.h"
#include "param.h"

int main() {

                            /*************************/
                            /***   READING DATA    ***/
                            /*************************/
    double *x, *logNe;
    int N = readalloc(&x, &logNe, 4800);    /* we have 4726 points of data */


                            /*************************/
                            /***   INITIALIZING    ***/
                            /*************************/

    /* structures for interpolation */
    gsl_interp_accel *acc = gsl_interp_accel_alloc();
    gsl_spline *spline = gsl_spline_alloc(gsl_interp_steffen, N);
    gsl_spline_init(spline, x, logNe, N);

    /* matrizes:  H = H0 + diag( D(t), 0, 0 ) */
    /* H0 = H0_re + i H0_im */
    gsl_matrix *H0_re, *H0_im;
    gsl_matrix_complex *h0;
    genMatrix_alloc(&H0_re, &H0_im, &h0);

    /* parameters data structure */
    par params = { acc, spline, H0_re, H0_im };

    /* sistema de EDO */            /* jac : jacobiano nao eh necessario para esse algoritmo  */
                                    /*  |  */
                                    /*  v  */
    gsl_odeiv2_system ode_sys = {func, NULL, DIM, &params};
    gsl_odeiv2_driver *driver =
        gsl_odeiv2_driver_alloc_y_new(&ode_sys, gsl_odeiv2_step_rkf45, PASSO, EPS_ABS, EPS_REL);    /* Runge-Kutta-Fehlberg 45 method */


    /* initial condition */
    double t = T_INIC;
    int i; double ti;

                            /*************************/
                            /***      MASS NU      ***/
                            /*************************/

    /* mass neutrinos */
    if (NU_TYPE == MASS)
    {
        /* structure for eigensystem */
        gsl_matrix_complex *H = gsl_matrix_complex_alloc(GERAC, GERAC);
        gsl_matrix_complex *eigvec = gsl_matrix_complex_alloc(GERAC, GERAC);
        gsl_vector *eigval = gsl_vector_alloc(GERAC);   /* eigval eh inutil, so ta aqui pq eh preciso dele para o eigvec */
        gsl_eigen_hermv_workspace *workspc = gsl_eigen_hermv_alloc(GERAC);
        eigen_problem eig_prob = { &params, h0, H, eigvec, eigval, workspc };

        gsl_vector_complex *psi_complex = gsl_vector_complex_alloc(GERAC);  /* psi_complex = (psi_e, psi_mu, psi_tau) */
        gsl_vector_complex *phi = gsl_vector_complex_alloc(GERAC);  /* phi is the complex vector (psi_1m, psi_2m, psi_3m) */
        gsl_vector_complex_set(phi, 0, gsl_complex_rect(RE1, IM1));
        gsl_vector_complex_set(phi, 1, gsl_complex_rect(RE2, IM2));
        gsl_vector_complex_set(phi, 2, gsl_complex_rect(RE3, IM3));

        /* solve the eigensystem a first time to get eigvec matrix */
        solveEigensys(t, &eig_prob);

        /* change basis to get initial condition in interaction basis */
        mult_MatrixVec(CblasConjTrans, eigvec, phi, psi_complex);   /* psi_complex = eigvec^dagger . phi */
        double Psi[DIM] = { psi_R(0), psi_R(1), psi_R(2), psi_I(0), psi_I(1), psi_I(2) };

        for (i = 0; i <= NUM_IT; i++) {
            ti = i * PASSO + T_INIC;
            int status = gsl_odeiv2_driver_apply(driver, &t, ti, Psi);

            if (status != GSL_SUCCESS) {
                printf ("Error, return status = %d\n", status);
                break;
            }

            /* getting eigvec matrix */
            solveEigensys(t, &eig_prob);

            /* pegando o psi_complex a partir do Psi */
            for (int j = 0; j < GERAC; j++)
                gsl_vector_complex_set(psi_complex, j, gsl_complex_rect(Psi[j], Psi[j+GERAC]));
            mult_MatrixVec(CblasNoTrans, eigvec, psi_complex, phi);  /* phi = eigvec . psi_complex */

            printf("%d    %.5e    %.5e      %.5e      %.5e      %.5e      %.5e      %.5e\n"  ,
                     i,     t,   phi_R(0), phi_R(1), phi_R(2), phi_I(0), phi_I(1), phi_I(2) );
        }

        /* freeing eigensystem resources */
        gsl_vector_complex_free(phi); gsl_vector_complex_free(psi_complex);
        gsl_eigen_hermv_free(workspc);
        gsl_vector_free(eigval);
        gsl_matrix_complex_free(eigvec);
        gsl_matrix_complex_free(H);
    }

                            /*************************/
                            /***      PART NU      ***/
                            /*************************/

    else
    {
        double Psi[DIM] = { RE1, RE2, RE3, IM1, IM2, IM3 };
        for (i = 0; i <= NUM_IT; i++) {
            ti = i * PASSO + T_INIC;
            int status = gsl_odeiv2_driver_apply(driver, &t, ti, Psi);

            if (status != GSL_SUCCESS) {
                printf ("Error, return status = %d\n", status);
                break;
            }
            printf("%d    %.5e    %.5e    %.5e    %.5e    %.5e    %.5e    %.5e\n",
                     i,     t,   Psi[0], Psi[1], Psi[2], Psi[3], Psi[4], Psi[5] );
        }
    }


                            /*************************/
                            /* FREEING THE RESOURCES */
                            /*************************/

    /* matrices */
    gsl_matrix_free(H0_re); gsl_matrix_free(H0_im);
    gsl_matrix_complex_free(h0);

    /* interpol */
    gsl_odeiv2_driver_free(driver); gsl_spline_free(spline); gsl_interp_accel_free(acc);

    /* data */
    free(x); free(logNe);

    return 0;
}

void genMatrix_alloc(gsl_matrix **H0r_ptr, gsl_matrix **H0i_ptr, gsl_matrix_complex **h0_ptr) {
    /* defining the mixing matrix */
    double th12 = DEG_TO_RAD(THETA12),
           th23, th13, d_CP;

    /* checking NUM_NU */
    if (NUM_NU == 3) {
        th23 = DEG_TO_RAD(THETA23), th13 = DEG_TO_RAD(THETA13), d_CP = DEG_TO_RAD(DELTACP);
    }
    else if (NUM_NU == 2) {
        th23 = th13 = d_CP = 0.0;
    }
    else {
        printf("Erro: NUM_NU nao eh 2 nem 3. NUM_NU = %d\n", NUM_NU);
        exit(EXIT_FAILURE);
    }

    /* calculating mixing matrix */
    double s12 = sin(th12), c12 = cos(th12),
           s23 = sin(th23), c23 = cos(th23),
           s13 = sin(th13), c13 = cos(th13);

    gsl_complex
    U11=REAL(c12 * c13),                        U12=REAL(s12 * c13),                        U13=POLAR(s13,-d_CP),
    U21=ADD(POLAR(-c12*s23*s13,d_CP),-s12*c23), U22=ADD(POLAR(-s12*s23*s13,d_CP),c12*c23),  U23=REAL(s23*c13),
    U31=ADD(POLAR(-c12*c23*s13,d_CP),s12*s23),  U32=ADD(POLAR(-s12*c23*s13,d_CP),-c12*s23), U33=REAL(c23*c13);

    gsl_matrix_complex *U = gsl_matrix_complex_alloc(GERAC, GERAC);
    MAT(U, 0, 0, U11);  MAT(U, 0, 1, U12);  MAT(U, 0, 2, U13);
    MAT(U, 1, 0, U21);  MAT(U, 1, 1, U22);  MAT(U, 1, 2, U23);
    MAT(U, 2, 0, U31);  MAT(U, 2, 1, U32);  MAT(U, 2, 2, U33);


    /* diagonal matrix for the masses */
    gsl_matrix_complex *M = gsl_matrix_complex_alloc(GERAC, GERAC);
    gsl_matrix_complex_set_zero(M);
    double m1 = 0.0, m2 = DM2_21, m3 = DM2_31;
    MAT(M, 0, 0, REAL(m1/2.0));
                              MAT(M, 1, 1, REAL(m2/2.0));
                                                           MAT(M, 2, 2, REAL(m3/2.0));

    *h0_ptr = gsl_matrix_complex_alloc(GERAC, GERAC);
    action(U, M, *h0_ptr);   /* h0 = U . M . U^dagger */

    /* h0 = H0_re + i H0_im */
    *H0r_ptr = gsl_matrix_alloc(GERAC, GERAC),
    *H0i_ptr = gsl_matrix_alloc(GERAC, GERAC);

    /* colocando as partes reais e imaginarias */
    for (int i = 0; i < GERAC; i++) {
        for (int j = 0; j < GERAC; j++) {
            gsl_matrix_set(*H0r_ptr, i, j, GSL_REAL(gsl_matrix_complex_get(*h0_ptr, i, j)));
            gsl_matrix_set(*H0i_ptr, i, j, GSL_IMAG(gsl_matrix_complex_get(*h0_ptr, i, j)));
        }
    }

    /* liberando a memoria das matrizes auxiliares */
    gsl_matrix_complex_free(M);
    gsl_matrix_complex_free(U);
}

/* dPsi/dt = f(Psi, t) */
int func(double t, const double Psi[], double f[], void *params) {
    par *param = (par *) params;
    gsl_matrix *H0_re = param->H0_re;
    gsl_matrix *H0_im = param->H0_im;
                                /* notice: the DENSITY TERM always corresponds to R11 = H0re(0, 0) */
/*             J11                 J12                 J13                   R11                 R12                 R13                DENSITY TERM  */
f[0] =   H0im(0, 0)*Psi[0] + H0im(0, 1)*Psi[1] + H0im(0, 2)*Psi[2]  +  H0re(0, 0)*Psi[3] + H0re(0, 1)*Psi[4] + H0re(0, 2)*Psi[5]  +  D(t, param)*Psi[3];
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

/* calculates D = sqrt(2) . G_F . N_e, where N_e = exp(log N_e) is interpolated */
double D(double t, void *params) {
    par *param = (par *) params;
    return M_SQRT2 * G_F * exp(gsl_spline_eval(param->spline, t, param->acc));
}

/* solves eig_prob depending on t */
void solveEigensys(double t, eigen_problem *eig_prob) {
    gsl_matrix *H0_re = eig_prob->params->H0_re;
    /* the matrix H is destroyed when calculating its eigenvectors */
    gsl_matrix_complex_memcpy(eig_prob->H, eig_prob->h0);
    gsl_matrix_complex_set(eig_prob->H, 0, 0, REAL(H0re(0, 0) + D(t, eig_prob->params)));
    gsl_eigen_hermv(eig_prob->H, eig_prob->eigval, eig_prob->eigvec, eig_prob->workspc);
    gsl_eigen_hermv_sort(eig_prob->eigval, eig_prob->eigvec, GSL_EIGEN_SORT_VAL_ASC);
}

/* evaluates u = op(A) v, where op(A) = A, A^t, A^dagger for TransA = CblasNoTrans, CblasTrans, CblasConjTrans */
void mult_MatrixVec(CBLAS_TRANSPOSE_t TransA, const gsl_matrix_complex *A, gsl_vector_complex *x, gsl_vector_complex *y) {
    /* zgemv does: y = alpha op(A).x + beta y */
                                /*    alpha                    beta    */
    gsl_blas_zgemv(TransA, GSL_COMPLEX_ONE, A, x, GSL_COMPLEX_ZERO, y);
}

/* evaluates C = A B A^dagger */
void action(const gsl_matrix_complex *A, const gsl_matrix_complex *B, gsl_matrix_complex *C) {
    gsl_matrix_complex *A_mult_B = gsl_matrix_complex_alloc(A->size1, B->size2);
    /* zgemm does: C = alpha op(A) op(B) + beta C */
                                                /*    alpha                          beta  */
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,   GSL_COMPLEX_ONE, A, B,        GSL_COMPLEX_ZERO, A_mult_B); /* A_mult_B = A . B */
    gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, A_mult_B, A, GSL_COMPLEX_ZERO, C);        /* C = A.B.A^dagger */

    gsl_matrix_complex_free(A_mult_B);  /* freeing the memory of the auxiliary variable */
}

/* reads data from stdin, storages them in x_ptr and y_ptr and returns number of data */
int readalloc(double **x_ptr, double **y_ptr, int chunk) {
    int i = 0, N = 0;
    double *x = (double *) malloc(chunk * sizeof(double));  /*  chunk determines size of  */
    double *y = (double *) malloc(chunk * sizeof(double));  /*  allocation for each step  */

    /* getting input */
    while (scanf("%lf%lf", &x[N], &y[N]) != EOF) {
        i++; N++;
        if (i > chunk-1) {
            x = (double *) realloc(x, (N+chunk) * sizeof(double));
            y = (double *) realloc(y, (N+chunk) * sizeof(double));
            i = 0;
        }
    }
    /* resizing the arrays to correct size */
    *x_ptr = (double *) realloc(x, N * sizeof(double));
    *y_ptr = (double *) realloc(y, N * sizeof(double));
    return N;
}

void printMatrix(gsl_matrix **M) {
    for (int i = 0; i < GERAC; i++)
        for (int j = 0; j < GERAC; j++)
            printf("%.3lf%c", gsl_matrix_get(*M, i, j), (j == GERAC-1) ? '\n' : ' ');
}
