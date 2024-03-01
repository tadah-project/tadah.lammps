#ifndef M_BLR_Core_H
#define M_BLR_Core_H

#include "m_core.h"
#include "functions/function_base.h"
#include "functions/basis_functions/bf_base.h"
#include "../CORE/config/config.h"

#include <iostream>

template
<class BF=Function_Base&>
class M_BLR_Core:
    public virtual M_Core
{
    public:
        Matrix Sigma;
        Config &config;
        BF bf;
        M_BLR_Core(Config &c):
            config(c),
            bf(c)
    {
        static_assert(std::is_base_of<BF_Base, BF>::value,
                "\nThe provided Basis function is not of a BasisFunction type.\n\
                For BLR use basis functions defined in basis_functions.h.\n\
                If you want to use this kernel functions use M_KRR class instead.\n\
                e.g. Use BF_Linear instead of Kern_Linear\n.");

        init();
    }

        M_BLR_Core(BF &bf, Config &c):
            config(c),
            bf(bf)
    {
        if (dynamic_cast<BF_Base*>(&bf) == nullptr)
            throw std::invalid_argument("Provided object is not of BF_Base type");

        init();
    }
        double predict(const aed_type2 &v) const{
            return bf.epredict(get_weights(),v);
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
#endif
