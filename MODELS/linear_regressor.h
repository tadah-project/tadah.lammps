#ifndef LINEAR_REGRESSOR_H
#define LINEAR_REGRESSOR_H

#include "../CORE/config/config.h"
#include "../CORE/core_types.h"
#include "ea.h"
#include "svd.h"
#include "ols.h"
#include "ridge_regression.h"

#include <vector>

class LinearRegressor {
            // SVD - most accurate but the slowest
            // w = V*D*S(*)*U^T*T
            // where S(*) = S.asDiagonal().inverse()
            // and D is a filter matrix with diagonal entires: S(i)^2/(S(i)^2+lambda)
    public:
        static void train(Config &config, phi_type &Phi, t_type &T,
                t_type &weights, Matrix &Sigma) {

            int verbose = config.get<int>("VERBOSE");
            double lambda=config.get<double>("LAMBDA");

            if (lambda == 0) {
                OLS::solve(Phi,T, weights);
            }
            else {
                double alpha = config.get<double>("ALPHA");
                double beta = config.get<double>("BETA");

                SVD svd(Phi);

                if (lambda<0) {
                    EA ea(config, svd, T);
                    ea.run(alpha,beta);
                    lambda=alpha/beta;

                    Sigma = ea.get_sigma(alpha, beta);
                    config_add_sigma(config,Sigma);

                    config.remove("ALPHA");
                    config.remove("BETA");
                    config.add("ALPHA",alpha);
                    config.add("BETA",beta);
                }

                if (verbose) std::cout << std::endl << "REG LAMBDA: " << lambda << std::endl;
                RidgeRegression::solve(svd,T,weights,lambda);
            }
        }
        static void config_add_sigma(Config &config, Matrix &Sigma) {
            t_type Sigma_as_vector(Sigma.data(), Sigma.cols()*Sigma.rows());
            config.add("SIGMA",Sigma.rows());
            config.add<t_type>("SIGMA",Sigma_as_vector);
        }
        static void read_sigma(Config &config, Matrix &Sigma) {
            using V=std::vector<double>;
            // get N = number of rows or cols, Sigma is NxN
            size_t N;
            try {
                N = config.get<size_t>("SIGMA");
            }
            catch (const std::runtime_error& error) {
                throw std::runtime_error("Sigma matrix is not computed. \n\
                        Hint: It is only computed for LAMBDA -1.\n");
            }
            // +1 as first number is N - a number of rows and cols
            V Sigma_as_vector(N*N+1);
            // read sigma including N, it is sotred col by col
            config.get<V>("SIGMA",Sigma_as_vector);
            Sigma_as_vector.erase(Sigma_as_vector.begin());
            Sigma.resize(N,N);
            for (size_t c=0; c<N; ++c)
                for (size_t r=0; r<N; ++r)
                    Sigma(r,c) = Sigma_as_vector.at(N*c+r);

        }
};
#endif
