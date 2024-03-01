#include "kern_linear.h"


Kern_Linear::Kern_Linear() {}
Kern_Linear::Kern_Linear (const Config &c):
    Function_Base(c)
{}
double Kern_Linear::operator() (const aed_type2& b, const aed_type2& af)
{
    return b*af;;
}
aed_type2 Kern_Linear::derivative(const aed_type2& b, const aed_type2&)
{
    return b;
}
double Kern_Linear::prime(const aed_type2& b, const aed_type2& , const aed_type2& ff)
{
    // return derivative(b,af).transpose() * ff;
    return b * ff;
}
std::string Kern_Linear::get_label()
{
    return label;
}
double Kern_Linear::epredict(const t_type &weights, const aed_type2& aed)
{
    return weights*aed;
}
double Kern_Linear::fpredict(const t_type &weights, const fd_type &fdij,
        const aed_type2& , const size_t k)
{
    return -(fdij(k) * weights);
}
force_type Kern_Linear::fpredict(const t_type &weights, const fd_type &fdij,
        const aed_type2& )
{
    return -(fdij * weights);
}
