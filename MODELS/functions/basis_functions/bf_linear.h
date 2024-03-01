#ifndef BF_LINEAR_H
#define BF_LINEAR_H
#include "bf_base.h"

/** Linear basis function */
struct BF_Linear: public virtual BF_Base
{
    BF_Linear();
    BF_Linear(const Config &c);
    std::string get_label();
    std::string label = "BF_Linear";
    double epredict(const t_type &weights, const aed_type2& aed);
    double fpredict(const t_type &weights, const fd_type &fdij,
            const aed_type2& aedi, const size_t k);
    force_type fpredict(const t_type &weights, const fd_type &fdij,
            const aed_type2& aedi);
};
#endif
