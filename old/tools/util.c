void printMatrix(gsl_matrix *M);   /* print a real matrix, only for debugging purposes */
void printMatrixComplex(gsl_matrix_complex *M); /* print complex matrix, for debugging */

void printMatrix(gsl_matrix *M) {
    for (int i = 0; i < GERAC; i++)
        for (int j = 0; j < GERAC; j++)
            printf("%.3e%c", gsl_matrix_get(M, i, j), (j == GERAC-1) ? '\n' : ' ');
}

void printMatrixComplex(gsl_matrix_complex *M) {
    for (int i = 0; i < GERAC; i++)
        for (int j = 0; j < GERAC; j++)
            printf("%.3e + %.3e i %c", GSL_REAL(gsl_matrix_complex_get(M, i, j)),
                                       GSL_IMAG(gsl_matrix_complex_get(M, i, j)), (j == GERAC-1) ? '\n' : ' ');
}
