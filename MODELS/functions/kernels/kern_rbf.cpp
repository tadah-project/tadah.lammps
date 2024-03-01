#include "kern_rbf.h"


Kern_RBF::Kern_RBF () {}
Kern_RBF::Kern_RBF (const Config &c):
    Function_Base(c),
    sigma(c.get<double>("MPARAMS")),
    gamma(1.0/(2.0*sigma*sigma))
{
    if (get_verbose()) std::cout << std::endl << label << " | sigma: " << sigma
        << " | gamma: " << gamma << std::endl;

    read_basis_from_config(c,basis);
}
double Kern_RBF::operator() (const aed_type2& b, const aed_type2& af)
{
    const aed_type2 v = af-b;
    return exp(-gamma*(v*v));
}
aed_type2 Kern_RBF::derivative(const aed_type2& b, const aed_type2& af)
{
    return  -2.0*gamma*(af-b)*(*this)(af,b);
}
double Kern_RBF::prime(const aed_type2& b, const aed_type2& af, const aed_type2& ff)
{
    return derivative(b, af)*ff;
}
std::string Kern_RBF::get_label()
{
    return label;
}
