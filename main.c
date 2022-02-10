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

    /* generating all parameters and allocating memory */
    par *params = genAllocParams(x, logNe, N);

    /* sistema de EDO */            /* jac */
    gsl_odeiv2_system ode_sys = {func, NULL, DIM, params};
    gsl_odeiv2_driver *driver =
        gsl_odeiv2_driver_alloc_y_new(&ode_sys, METHOD, PASSO, EPS_ABS, EPS_REL);

    /* initial condition */
    int i, status;
    double ti, t = T_INIC;


                            /*************************/
                            /***    MASS  BASIS    ***/
                            /*************************/

    if (NU_TYPE == MASS) {
        /* allocating memory for eigen operations */
        eigen_problem *eig_prob = allocEigensys();
        complex_vecs *cplx_vecs = allocVecs();
        gsl_vector_complex *psi = cplx_vecs->psi, *phi = cplx_vecs->phi;    /* for convenient macros */

        /* solve the eigensystem a first time to get eigvec matrix */
        solveEigensys(t, params, eig_prob);

        /* change basis to get initial condition in interaction basis */
        multMatrixVec(CblasConjTrans, eig_prob->eigvec, phi, psi);   /* psi = eigvec^dagger . phi */
        double Psi[DIM] = { psi_R(0), psi_R(1), psi_R(2), psi_I(0), psi_I(1), psi_I(2) };   /* initial condition in mass basis */

        for (i = 0; i <= NUM_IT; i++) {
            ti = i * PASSO + T_INIC;
            if ((status = gsl_odeiv2_driver_apply(driver, &t, ti, Psi)) != GSL_SUCCESS) {
                printf ("Error, return value = %d\n", status);
                break;
            }

            solveEigensys(t, params, eig_prob); /* getting eigvec matrix */
            for (int j = 0; j < GERAC; j++)     /* updating the complex psi */
                gsl_vector_complex_set(psi, j, gsl_complex_rect(Psi[j], Psi[j+GERAC]));
            multMatrixVec(CblasNoTrans, eig_prob->eigvec, psi, phi);   /* phi = eigvec . psi */

            printf("%d    %.5e    %.5e      %.5e      %.5e      %.5e      %.5e      %.5e\n" ,
                     i,     t,   phi_R(0), phi_R(1), phi_R(2), phi_I(0), phi_I(1), phi_I(2));
        }

        /* free eigen resources */
        freeVecs(cplx_vecs);
        freeEigensys(eig_prob);
    }

                            /*************************/
                            /*** INTERACTION BASIS ***/
                            /*************************/

    else {
        double Psi[DIM] = { RE1, RE2, RE3, IM1, IM2, IM3 }; /* initial condition in interaction basis */

        for (i = 0; i <= NUM_IT; i++) {
            ti = i * PASSO + T_INIC;
            if ((status = gsl_odeiv2_driver_apply(driver, &t, ti, Psi)) != GSL_SUCCESS) {
                printf ("Error, return value = %d\n", status);
                break;
            }

            printf("%d    %.5e    %.5e      %.5e      %.5e      %.5e      %.5e      %.5e\n" ,
                     i,     t,   Psi[0],   Psi[1],   Psi[2],   Psi[3],   Psi[4],   Psi[5]  );
        }
    }


                            /*************************/
                            /* FREEING THE RESOURCES */
                            /*************************/

    /* params */
    freeParams(params);
    /* EDO */
    gsl_odeiv2_driver_free(driver);
    /* data */
    free(x); free(logNe);

    if (status != GSL_SUCCESS)
        return EXIT_FAILURE;
    else
        return EXIT_SUCCESS;
}

/* gera os parametros que precisamos e retorna o pointer dele */
par *genAllocParams(double *x, double *logNe, int N) {
    par *params = malloc(sizeof(par));
    params->acc = gsl_interp_accel_alloc();
    params->spline = gsl_spline_alloc(gsl_interp_steffen, N);
    gsl_spline_init(params->spline, x, logNe, N);

    /* defining the parameters */
    double th12 = DEG_TO_RAD(THETA12),
           th23, th13, d_CP;

    /* checking number of neutrinos and setting parameters accordingly */
    if (NUM_NU == 2)
        th23 = th13 = d_CP = 0.0;
    else
        th23 = DEG_TO_RAD(THETA23), th13 = DEG_TO_RAD(THETA13), d_CP = DEG_TO_RAD(DELTACP);

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

    /* diagonal matrix for the masses/(2*ENERG) */
    gsl_matrix_complex *M = gsl_matrix_complex_alloc(GERAC, GERAC);
    gsl_matrix_complex_set_zero(M);
    double m1 = -DM2_21, m2 = 0.0, m3 = DM2_32; /* de acordo com o paper da Gonzalez */
    MAT(M, 0, 0, REAL(m1/(2.0*ENERG)));
                                        MAT(M, 1, 1, REAL(m2/(2.0*ENERG)));
                                                                            MAT(M, 2, 2, REAL(m3/(2.0*ENERG)));

    params->h0 = gsl_matrix_complex_alloc(GERAC, GERAC);
    action(U, M, params->h0);   /* h0 = U . M . U^dagger */
    /* h0 = H0_re + i H0_im */
    params->H0_re = gsl_matrix_alloc(GERAC, GERAC); params->H0_im = gsl_matrix_alloc(GERAC, GERAC);

    /* colocando as partes reais e imaginarias */
    for (int i = 0; i < GERAC; i++) {
        for (int j = 0; j < GERAC; j++) {
            gsl_matrix_set(params->H0_re, i, j, GSL_REAL(gsl_matrix_complex_get(params->h0, i, j)));
            gsl_matrix_set(params->H0_im, i, j, GSL_IMAG(gsl_matrix_complex_get(params->h0, i, j)));
        }
    }

    /* liberando a memoria das matrizes auxiliares */
    gsl_matrix_complex_free(M);
    gsl_matrix_complex_free(U);

    return params;
}

void freeParams(par *params) {
    /* matrices */
    gsl_matrix_free(params->H0_re); gsl_matrix_free(params->H0_im);
    gsl_matrix_complex_free(params->h0);
    /* interpol */
    gsl_spline_free(params->spline); gsl_interp_accel_free(params->acc);
    free(params);
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
                                /*******************************************************************/
    return GSL_SUCCESS;
}

/* calculates D = sqrt(2) . G_F . N_e, where N_e = exp(log N_e) is interpolated */
double D(double t, void *params) {  /* void pointer because this defines a gsl_function */
    par *param = (par *) params;    /* for numeric differentiation, we use gsl_function */
    return M_SQRT2 * G_F * exp(gsl_spline_eval(param->spline, t, param->acc));
}

/* solves eig_prob depending on t */
void solveEigensys(double t, par *params, eigen_problem *eig_prob) {
    gsl_matrix *H0_re = params->H0_re;
    /* the matrix H is destroyed when calculating its eigenvectors */
    gsl_matrix_complex_memcpy(eig_prob->H, params->h0);
    gsl_matrix_complex_set(eig_prob->H, 0, 0, REAL(H0re(0, 0) + D(t, params)));
    gsl_eigen_hermv(eig_prob->H, eig_prob->eigval, eig_prob->eigvec, eig_prob->workspace);
    gsl_eigen_hermv_sort(eig_prob->eigval, eig_prob->eigvec, GSL_EIGEN_SORT_VAL_ASC);
}

/* defines and eig_prob based on params */
eigen_problem *allocEigensys() {
    eigen_problem *eig_prob = malloc(sizeof(eigen_problem));
    eig_prob->H = gsl_matrix_complex_alloc(GERAC, GERAC);
    eig_prob->eigvec = gsl_matrix_complex_alloc(GERAC, GERAC);
    eig_prob->eigval = gsl_vector_alloc(GERAC);   /* eigval eh inutil, mas o GSL precisa dele para o eigvec */
    eig_prob->workspace = gsl_eigen_hermv_alloc(GERAC);
    return eig_prob;
}

/* frees everything that allocEigensys allocated. It does not free params and h0 */
void freeEigensys(eigen_problem *eig_prob) {
    gsl_eigen_hermv_free(eig_prob->workspace);
    gsl_vector_free(eig_prob->eigval);
    gsl_matrix_complex_free(eig_prob->eigvec);
    gsl_matrix_complex_free(eig_prob->H);
    free(eig_prob);
}

/* alloc the two complex vectors: psi = interaction basis, phi = mass basis */
complex_vecs *allocVecs() {
    complex_vecs *cplx_vecs = malloc(sizeof(complex_vecs));
    cplx_vecs->psi = gsl_vector_complex_alloc(GERAC);  /* psi_complex = (psi_e, psi_mu, psi_tau) */
    cplx_vecs->phi = gsl_vector_complex_alloc(GERAC);  /* phi is the complex vector (psi_1m, psi_2m, psi_3m) */
    gsl_vector_complex_set(cplx_vecs->phi, 0, gsl_complex_rect(RE1, IM1));
    gsl_vector_complex_set(cplx_vecs->phi, 1, gsl_complex_rect(RE2, IM2));
    gsl_vector_complex_set(cplx_vecs->phi, 2, gsl_complex_rect(RE3, IM3));
    return cplx_vecs;
}

/* frees everything allocated by allocVecs */
void freeVecs(complex_vecs *cplx_vecs) {
    gsl_vector_complex_free(cplx_vecs->phi);
    gsl_vector_complex_free(cplx_vecs->psi);
    free(cplx_vecs);
}

/* evaluates C = A B A^dagger */
void action(const gsl_matrix_complex *A, const gsl_matrix_complex *B, gsl_matrix_complex *C) {
    gsl_matrix_complex *A_mult_B = gsl_matrix_complex_alloc(A->size1, B->size2);
    /* zgemm does: C = alpha op(A) op(B) + beta C */
                                                  /*  alpha                          beta  */
    gsl_blas_zgemm(CblasNoTrans, CblasNoTrans,   GSL_COMPLEX_ONE, A, B,        GSL_COMPLEX_ZERO, A_mult_B); /* A_mult_B = A . B */
    gsl_blas_zgemm(CblasNoTrans, CblasConjTrans, GSL_COMPLEX_ONE, A_mult_B, A, GSL_COMPLEX_ZERO, C);        /* C = A.B.A^dagger */

    gsl_matrix_complex_free(A_mult_B);  /* freeing the memory of the auxiliary variable */
}

/* evaluates y = op(A) x, where op(A) = A, A^t, A^dagger for TransA = CblasNoTrans, CblasTrans, CblasConjTrans */
void multMatrixVec(CBLAS_TRANSPOSE_t TransA, const gsl_matrix_complex *A, gsl_vector_complex *x, gsl_vector_complex *y) {
                            /*  alpha                   beta  */
    gsl_blas_zgemv(TransA, GSL_COMPLEX_ONE, A, x, GSL_COMPLEX_ZERO, y);
}

/* reads data from standard input, storages them in x_ptr and y_ptr and returns N */
int readalloc(double **x_ptr, double **y_ptr, int chunk) {
    int i = 0, N = 0;
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
