#ifndef M_KRR_TRAIN_H
#define M_KRR_TRAIN_H

#include "linear_regressor.h"
#include "../CORE/config/config.h"
#include "functions/function_base.h"
#include "functions/kernels/kern_base.h"
#include "m_train.h"
#include "m_krr_core.h"
#include "ekm.h"
#include "../CORE/maths.h"

#include <iostream>
#include <limits>

template
<class Kern=Function_Base&>
class M_KRR_Train:
    public M_KRR_Core<Kern>,    // MUST NOT BE VIRUAL!!
    public M_Train
{

    public:
        EKM<Kern> ekm;
        using M_KRR_Core<Kern>::trained;
        using M_KRR_Core<Kern>::config;
        using M_KRR_Core<Kern>::kernel;
        using M_KRR_Core<Kern>::weights;
        using M_KRR_Core<Kern>::Sigma;

        M_KRR_Train(Config &c):
            M_KRR_Core<Kern>(c),
            ekm(c)
    {}

        M_KRR_Train(Kern &kernel, Config &c):
            M_KRR_Core<Kern>(kernel, c),
            ekm(kernel)
    {}

        void train(phi_type &Phi, t_type &T) {
            // implementation using EKM
            if(trained) {
                throw std::runtime_error("This object is already trained!");
            }
            if (kernel.get_label()!="Kern_Linear") {
                ekm.project(Phi);
            }
            LinearRegressor::train(config,Phi, T,weights,Sigma);

            if (kernel.get_label()!="Kern_Linear") {
                //kernalize weights
                weights = ekm.EKM_mat*weights;
            }
            trained=true;
        }

        void train2(Matrix &B, t_type &T) {
            // The standad KRR implementation where covariance
            // matrix K is computed.
            // The training is only done on energies.
            // The added regularisation is adhoc
            //double lambda=9e-11;

            double lambda = config.template get<double>("LAMBDA");

            // Imporve ill conditioned problems
            if(lambda<0)
                throw std::runtime_error("This KRR implementation cannot use LAMBDA -1 (Evidence Approximation)");
            std::cout << "Computing Kernel Matrix" << std::endl;
            Matrix K(B.cols(), B.cols()); // define covariance matrix
            std::cout << K.rows() << "  " << K.cols() << std::endl;
            for (size_t i=0; i<K.cols(); i++) {
                for (size_t j=0; j<K.cols(); j++) {
                    K(i,j) = kernel(B.col(i),B.col(j));
                    if (i==j) K(i,j)+=lambda;
                }
            }
            std::cout << "Solving..." << std::flush;
            weights = solve_posv(K,T);
            std::cout << "Done" << std::endl;

            //bool = config.template get<bool>("MODEL")[2];
            std::cout << "Matrix condition number: " << condition_number(K) << std::endl;

            trained=true;
        }
};
#endif
