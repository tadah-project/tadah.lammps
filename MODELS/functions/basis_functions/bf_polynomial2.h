#ifndef BF_POLYNOMIAL2_H
#define BF_POLYNOMIAL2_H
#include "bf_base.h"

/** Polynomial order 2 basis function with c=0 */
struct BF_Polynomial2: public virtual BF_Base
{
    BF_Polynomial2();
    BF_Polynomial2(const Config &c);
    std::string label = "BF_Polynomial2";
    std::string get_label();
    double epredict(const t_type &weights, const aed_type2& aed);
    double fpredict(const t_type &weights, const fd_type &fdij,
            const aed_type2& aedi, const size_t k);
    force_type fpredict(const t_type &weights, const fd_type &fdij,
            const aed_type2& aedi);
};
#endif
