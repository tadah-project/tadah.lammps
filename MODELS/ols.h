#ifndef OLS_H
#define OLS_H


#include "../CORE/core_types.h"

class OLS {
    public:
        template <typename M, typename V>
            static void solve(M &A, V &B, V &weights) {
                // A will be destroyed on exit
                // B will contain solution
                int m = A.rows();  //IN
                int n = A.cols();  //IN
                int nrhs = 1; // number of colums in B (targets T) //IN
                double *a = A.ptr(); //IN
                int lda = m;
                double *b = B.ptr();    // IN/OUT
                int ldb = std::max(m,n);
                double *s = new double[std::min(m,n)];   // size: rows
                double rcond = 1e-8;   // 1 use machine precision
                int rank;
                double *work;
                int lwork=-1; // QUERRY
                //int liwork = std::max(
                //        0,int(std::log2(std::min(m,n)/(25+1))) + 1)
                //    + 11*std::min(m,n);
                int smlsiz = 25;
                int nlvl = std::max( 0, int( std::log2( std::min(m,n)/(smlsiz+1) ) )+1 );
                int *iwork = new int[3*std::min(m,n)*nlvl + 11*std::min(m,n)];
                int info;
                double wkopt;
                dgelsd_(&m,&n,&nrhs,a,&lda,b,&ldb,s,&rcond,&rank,
                        &wkopt, &lwork, iwork, &info);
                lwork = (int)wkopt;
                work = new double[lwork];
                dgelsd_(&m,&n,&nrhs,a,&lda,b,&ldb,s,&rcond,&rank,
                        work, &lwork, iwork, &info);

                weights = B.copy(A.cols());

                delete[] s;
                delete[] work;
                delete[] iwork;
            }
};
#endif
