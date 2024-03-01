#ifndef SVD_H
#define SVD_H

#include "../CORE/core_types.h"
#include "../CORE/maths.h"
#include "../CORE/lapack.h"

class SVD {
    /** Compute Compact SVD
     *
     * M = USV^T
     */
    private:
        char jobz = 'S'; // compute the first min(M,N) columns of U
        int m;
        int n;
        int lda;
        int ldu;
        int ldvt;
        int ucol;
        int lwork;
        int info;

        // arrays
        double *a = nullptr;
        double *s = nullptr;
        double *u = nullptr;
        double *vt = nullptr;
        double *work = nullptr;
        int *iwork = nullptr;

        bool work_alloc=false;

    public:
        ~SVD()
        {
            delete[] s;
            delete[] u;
            delete[] vt;

            if (work_alloc)
                delete[] iwork;

            if (work_alloc && lwork>0)
                delete[] work;
        }
        SVD(const SVD &other)
        {
            std::cout << "Copy-constructor" << std::endl;
            jobz = other.jobz;
            m = other.m;
            n = other.n;
            lda = other.lda;
            ldu = other.ldu;
            ldvt = other.ldvt;
            ucol = other.ucol;
            lwork = other.lwork;
            info = other.info;

            // allocate all
            s = new double[sizeS()];
            u = new double[sizeU()];
            vt = new double[sizeVT()];

            // We don't need this arrays
            a = nullptr;
            work = nullptr;
            iwork = nullptr;

            work_alloc=false;

            std::cout << sizeS() << std::endl;;
            std::cout << sizeU() << std::endl;;
            std::cout << sizeVT() << std::endl;;

            // hard copy
            std::memcpy(s, other.s, sizeS()*sizeof(double));
            std::memcpy(u, other.u, sizeU()*sizeof(double));
            std::memcpy(vt, other.vt, sizeVT()*sizeof(double));
        }

        //SVD &operator=(const SVD &other)
        //{
        //    //if (this == &other)
        //    //    return *this;
        //    std::cout << "Copy-assignment" << std::endl;
        //}

        SVD(Matrix &Phi):
            m(Phi.rows()),
            n(Phi.cols()),
            lda(Phi.rows()),
            ldu(Phi.rows()),
            ldvt(std::min(m,n)),
            ucol(std::min(m,n)),
            a(&Phi.data()[0]),
            s(new double[sizeS()]),
            u(new double[sizeU()]),
            vt(new double[sizeVT()]),
            iwork(new int[8*std::min(m,n)]),
            work_alloc(true)
        {
            lwork=-1; // QUERRY
            double wkopt;
            dgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
                    &wkopt, &lwork, iwork, &info);
            lwork = (int)wkopt;
            work = new double[lwork];
            dgesdd_(&jobz, &m, &n, a, &lda, s, u, &ldu, vt, &ldvt,
                    work, &lwork, iwork, &info);
            // Filter out v.small singularvalues from S and its inverse
            for (size_t i=0; i<sizeS(); ++i) {
                if (s[i]<std::numeric_limits<double>::min()) {
                    s[i] = 0.0;
                }
            }
        }

        int get_info() const{
            return info;
        }
        double *getA() const{
            return a;
        }
        double *getU() const{
            return u;
        }
        double *getVT() const{
            return vt;
        }
        double *getS() const{
            return s;
        }
        size_t sizeA() const {
            return m*n;
        }
        size_t sizeU() const {
            return ldu*ucol;
        }
        size_t sizeVT() const {
            return ldvt*n;
        }
        size_t sizeS() const {
            return std::min(m,n);
        }
        std::pair<size_t,size_t> shapeA() const{
            return std::pair<size_t, size_t>(m,n);
        }
        std::pair<size_t,size_t> shapeU() const{
            return std::pair<size_t, size_t>(ldu,ucol);
        }
        std::pair<size_t,size_t> shapeVT() const{
            return std::pair<size_t, size_t>(ldvt,n);
        }
        std::pair<size_t,size_t> shapeS() const{
            return std::pair<size_t, size_t>(std::min(m,n),1);
        }

        template <typename T>
            T predict(const T &w) const {
                /* Return U S V^T w */

                double *t = new double[shapeA().first];
                double *x = new double[shapeA().second];

                T res = predict(w,t,x);

                delete[] t;
                delete[] x;
                return res;
            }
        template <typename T>
            T predict(const T &w, double *t, double *x) const {
                /* Return U S V^T w */

                const double *w_ptr = &w.data()[0];
                // 1. V^T w
                char trans='N';
                int m = static_cast<int>(shapeA().first);
                int n = static_cast<int>(shapeA().second);
                double alpha_ = 1.0;
                double beta_ = 0.0;
                int incx=1;

                dgemv_(&trans, &n, &n, &alpha_, vt, &n, w_ptr, &incx, &beta_, x, &incx);

                // 2. S (V^T w)
                for (int i=0; i<n; ++i) {
                    x[i]*=s[i];
                }

                // 3. U (S V^T w)
                dgemv_(&trans, &m, &n, &alpha_, u, &m, x, &incx, &beta_, t, &incx);
                return T(t,m);
            }

};

#endif
