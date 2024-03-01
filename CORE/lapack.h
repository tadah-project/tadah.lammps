#ifndef LAPACK_H
#define LAPACK_H

#include "maths.h"
extern "C" {
    /** the Cholesky factorization of a real symmetric
     *  positive definite matrix a.*/
    void dpotrf_(char *uplo, int *n, double *a, int *lda, int *info);

    /** Compute inverse using the Cholesky factorisation from dpotrf */
    void dpotri_(char *uplo, int *n, double *a, int *lda, int *info);

    /** Cholesky decomposition for packed matrices*/
    void dpptrf_(char *uplo, int *n, double *ap, int *info);

    /** LU decomoposition of a general matrix */
    void dgetrf_(int *m, int *n, double *a, int *lda, int *ipiv, int *info);

    /** Compute inverse of the matrix from its  LU factorisation */
    void dgetri_(int *n, double *a, int *lda, int *ipiv,
            double *work, int *lwork, int *info);

    void dgesdd_(char *jobz, int *m, int *n, double *a, int *lda,
            double *s, double *u, int *ldu, double *vt, int *ldvt,
            double *work, int *lwork, int *iwork, int *info);

    void dgemv_(const char *trans, const int *m, const int *n,
            const double *alpha, const double *a, const int *lda, const double *x,
            const int *incx, const double *beta, double *y, const int *incy);

    void dgelsd_ (int *m, int* n, int* nrhs, double *a, int* lda,
            double *b, int* ldb, double *s, double* rcond, int* rank,
            double *work, int* lwork, int* iwork, int* info);

    void dgemm_(char *transa, char *transb, int *m, int *n,
            int *k, double *alpha, double *a, int *lda, double *b, int *ldb,
            double *beta, double *c, int *ldc);

    /** Balance a general real matrix a */
    void dgebal_(char *job, int *n, double *a, int *lda, int *ilo, int *ihi,
            double *scale, int *info);

    void dposv_(char *uplo, int *n, int *nrhs, double *a, int *lda, double *b,
            int *ldb, int *info);

    void dpocon_(char *uplo, int *n, double *a, int *lda, double *anorm,
            double *rcond, double *work, int *iwork, int *info);

    double dlange_(char *norm, int *m, int *n, double *a, int *lda, double *work);

}


// Some usefull operations on vectors and matrices
aed_type2 T_dgemv(Matrix &mat, aed_type2 &vec, char trans='N');
aed_type2 T_MDMT_diag(const Matrix &M, const Matrix &D);
double norm_X(Matrix &A, char norm='1');
aed_type2 solve_posv(Matrix &A, aed_type2 /*&*/T, char uplo='U');
double condition_number(Matrix &A, char uplo='U');



#endif
