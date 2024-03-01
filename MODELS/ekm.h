#ifndef EKM_H
#define EKM_H

#include "../CORE/core_types.h"
#include "../CORE/lapack.h"

template <typename Kern>
class EKM {
    private:
        Kern kernel;
        // TODO add to config? Experiment with value
        double eps = std::numeric_limits<double>::epsilon()*0;
        std::vector<size_t> indices;
    public:
        Matrix EKM_mat;
        template <typename T>
            EKM(T &t):
                kernel(t)
    {}
        void configure(Matrix &basis) {

            // check norm of basis vectors and keep indices
            for (size_t i=0; i<basis.cols(); i++) {
                double norm = kernel(basis.col(i),basis.col(i));
                if (norm > eps) indices.push_back(i);
            }

            // Copy basis vectors with non zero norm
            // to new basis b
            // basis vectors are stored as columns in basis
            Matrix b(basis.rows(), indices.size());
            for (const auto &i:indices) {
                for (size_t j=0; j<basis.rows(); ++j) {
                    b(j,i) = basis(j,i);
                }
            }
            // Update old basis with new one
            basis=b;

            Matrix K(b.cols(), b.cols());
            for (size_t i=0; i<K.cols(); i++) {
                for (size_t j=0; j<K.cols(); j++) {
                    K(i,j) = kernel(b.col(i),b.col(j));
                }
            }

            // Matrix balancing
            double *a = K.ptr();
            int n = K.cols();
            char job='B';
            int ilo; int ihi; int info;
            double *scale = new double[n];
            dgebal_(&job, &n, a, &n, &ilo, &ihi, scale, &info);
            delete [] scale;

            // Find max value in K diagonal
            // and update eps
            double max_val=0;
            for (size_t i=0; i<K.cols(); i++) {
                if (K(i,i)>max_val) max_val=K(i,i);
            }
            eps *= max_val;

            SVD svd(K);

            double *s = svd.getS();
            double *vt = svd.getVT();
            t_type S(s,svd.sizeS());
            Matrix VT(vt,svd.shapeVT().first, svd.shapeVT().second);

            // First count number of non zero singular values, then
            // filter out corresponding right singular vectors from
            // the final EKM_mat
            size_t nonzero_sv=0;
            for (size_t i=0; i<S.size(); ++i) if (s[i]>eps) nonzero_sv++;

            EKM_mat = Matrix(nonzero_sv, svd.shapeVT().second);

            for (size_t i=0; i<S.size(); ++i) {
                if (s[i]>eps) {
                    for (size_t j=0; j<EKM_mat.cols(); j++) {
                        EKM_mat(i,j) = VT(i,j); // /s[i];
                    }
                }
            }

        }
        void inverse(Matrix &A) {
            double *a = &A.data()[0];
            int n = A.rows();   // A is square
            int *ipiv = new int[n];
            int lwork = n*n;
            double *work = new double[lwork];
            int info;

            dgetrf_(&n,&n,a,&n,ipiv,&info);
            dgetri_(&n,a,&n,ipiv,work,&lwork,&info);

            delete[] ipiv;
            delete[] work;
        }
        void cholesky(Matrix &A, char uplo) {
            // matrix must be either U or L
            int n = A.cols(); // A is square
            double *a = &A.data()[0];
            int info;
            dpotrf_(&uplo, &n, a, &n, &info);
        }

        void project(phi_type &Phi) {
            // Phi'=Phi*EKM_mat
            char transa = 'N';
            char transb = 'N';
            int m = Phi.rows();
            int n = Phi.cols();
            double *a = &Phi.data()[0];
            double *b = &EKM_mat.data()[0];
            double *c = new double[m*n];
            double alpha = 1.0;
            double beta = 0.0;
            dgemm_(&transa, &transb, &m, &n, &n, &alpha, a, &m, b, &n,
                    &beta, c, &m);

            Phi = Matrix(c,m,n);

            delete[] c;

        }
        aed_type2 project(aed_type2 &aed, Matrix &basis) {
            aed_type2 temp(basis.cols());
            for (size_t b=0; b<basis.cols(); b++) {
                temp(b) = kernel(aed,basis[b]);
            }
            return  EKM_mat*temp;
        }
};
#endif
