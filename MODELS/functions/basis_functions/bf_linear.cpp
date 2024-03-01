#include "bf_linear.h"



BF_Linear::BF_Linear() {}
BF_Linear::BF_Linear(const Config &c): Function_Base(c) {}
std::string BF_Linear::get_label() {
    return label;
}
double BF_Linear::epredict(const t_type &weights, const aed_type2& aed)
{
    return weights*aed;
}
double BF_Linear::fpredict(const t_type &weights, const fd_type &fdij,
        const aed_type2& , const size_t k)
{
    return -1*(fdij(k) * weights);
}
force_type BF_Linear::fpredict(const t_type &weights, const fd_type &fdij,
        const aed_type2& )
{
    return -1*(fdij * weights);
}
