#ifndef RIDGE_REGRESSION_H
#define RIDGE_REGRESSION_H

#include "../CORE/core_types.h"
#include "../CORE/maths.h"
#include "../CORE/lapack.h"
#include "svd.h"

class RidgeRegression {
    public:
        template <typename V, typename W>
            static void solve(const SVD &svd, V /*&*/T, W &weights,
                    const double lambda) {

                double *u = svd.getU();
                double *s = svd.getS();
                double *vt = svd.getVT();

                double *d = new double[svd.sizeS()];

                for (size_t i=0; i<svd.sizeS(); ++i) {
                    double sv = s[i];
                    if (sv>std::numeric_limits<double>::min())
                        d[i] = (sv)/(sv*sv+lambda);
                    else
                        d[i] = 0.0;
                }

                int m = svd.shapeA().first;
                int n = svd.shapeA().second;
                int ldu = svd.shapeU().first;
                int ucol = svd.shapeU().second;
                int lda = m;
                int ldvt = svd.shapeVT().first;

                // SOLVE: V D U^T T

                //STEP1: U^T T
                char trans = 'T';
                double alpha_ = 1.0;
                double *x = T.ptr();
                double *y = new double[lda];
                int incx = 1;
                double beta_ = 0.0;
                dgemv_(&trans, &ldu, &ucol, &alpha_, u, &lda, x, &incx,
                        &beta_, y, &incx);

                // STEP2 : D (U^T T)
                for (int r=0; r<ldvt; ++r) {
                    y[r]*=d[r];
                }

                // STEP3: V (D U^T T)
                trans='T';
                dgemv_(&trans, &ldvt, &n, &alpha_, vt, &ldvt, y, &incx,
                        &beta_, x, &incx);

                weights = t_type(x, (size_t)ldvt, true);

                delete[] d;
                delete[] y;
            }
};
#endif
