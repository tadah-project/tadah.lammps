#include "kern_polynomial.h"


Kern_Polynomial::Kern_Polynomial () {}
Kern_Polynomial::Kern_Polynomial (const Config &c):
    Function_Base(c),
    d(c.get<int>("MPARAMS",0)),
    gamma(c.get<double>("MPARAMS",1)),
    c(c.get<double>("MPARAMS",2))
{
    if (get_verbose()) std::cout << std::endl << label << " | degree: " << d
        <<" | gamma: " << gamma << " | c: " << c << std::endl;
    read_basis_from_config(c,basis);
}
double Kern_Polynomial::operator() (const aed_type2& b, const aed_type2& af)
{
    return std::pow(gamma*(b*af)+c,d);
}
aed_type2 Kern_Polynomial::derivative(const aed_type2& b, const aed_type2& af)
{
    return  d*gamma*std::pow(gamma*(b*af)+c,d-1)*b;
}
double Kern_Polynomial::prime(const aed_type2& b, const aed_type2& af, const aed_type2& ff)
{
    return derivative(b, af)*ff;
}
std::string Kern_Polynomial::get_label()
{
    return label;
}
