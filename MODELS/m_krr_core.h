#ifndef M_KRR_Core_H
#define M_KRR_Core_H

#include "../CORE/config/config.h"
#include "functions/function_base.h"
#include "functions/kernels/kern_base.h"
#include "m_core.h"

#include <iostream>

template
<class K=Function_Base&>
class M_KRR_Core:
    public virtual M_Core
{
    public:
        Matrix Sigma;
        Config &config;
        K kernel;

        M_KRR_Core(Config &c):
            config(c),
            kernel(c)
    {
        static_assert(std::is_base_of<Kern_Base, K>::value,
                "\nThe provided Kernel K is not of a Kernel type.\n\
                For KRR use kernels defined in kernel.h.\n\
                e.g. Use Kern_Linear instead of BF_Linear\n.");

        init();
    }

        M_KRR_Core(K &kernel, Config &c):
            config(c),
            kernel(kernel)
    {
        if (dynamic_cast<Kern_Base*>(&kernel) == nullptr)
            throw std::invalid_argument("Provided object is not of Kern_Base type");

        init();
    }
        double predict(const aed_type2 &v) const{
            return kernel.epredict(weights,v);
        }

        t_type get_weights_uncertainty() const{
            double lambda=config.template get<double>("LAMBDA");
            if(lambda >= 0) throw std::runtime_error(
                    "Sigma matrix is only computed for LAMBDA < 0");
            return Sigma.diagonal();
        }


    private:
        void init() {
            if (config.exist("WEIGHTS")) {
                weights.resize(config.size("WEIGHTS"));
                config.get<t_type>("WEIGHTS",weights);
                trained=true;
            }

            verbose=(config.get<int>("VERBOSE"));
        }
};

//template
//<class K=Function_Base&>
//class M_KRR_Predict:
//    public M_KRR_Core<K>,    // MUST NOT BE VIRUAL!!
//    public virtual M_Predict
//{
//
//    public:
//        using M_KRR_Core<K>::config;
//        using M_KRR_Core<K>::kernel;
//        using M_KRR_Core<K>::weights;
//
//        M_KRR_Predict(Config &c):
//            M_KRR_Core<K>(c)
//    {}
//
//        M_KRR_Predict(K &kernel, Config &c):
//            M_KRR_Core<K>(kernel, c)
//    {}
//
//        double epredict(const aed_type2 &aed) {
//            return kernel.epredict(weights,aed);
//        };
//
//        double fpredict(const fd_type &fdij, const aed_type2 &aedi, const size_t k) {
//            return kernel.fpredict(weights,fdij,aedi,k);
//        }
//
//        force_type fpredict(const fd_type &fdij, const aed_type2 &aedi) {
//            return kernel.fpredict(weights,fdij,aedi).array();
//        }
//
//};
#endif
