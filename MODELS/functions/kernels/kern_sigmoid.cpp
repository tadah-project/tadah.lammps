#include "kern_sigmoid.h"


Kern_Sigmoid::Kern_Sigmoid () {}
Kern_Sigmoid::Kern_Sigmoid (const Config &c):
    Function_Base(c),
    gamma(c.get<double>("MPARAMS",0)),
    c(c.get<double>("MPARAMS",1))
{
    if (get_verbose()) std::cout << std::endl << label << " | gamma: " << gamma
        << " | c: " << c << std::endl;

    read_basis_from_config(c,basis);
}
double Kern_Sigmoid::operator() (const aed_type2& b, const aed_type2& af)
{
    return std::tanh(gamma*(b*af)+c);
}
aed_type2 Kern_Sigmoid::derivative(const aed_type2& b, const aed_type2& af)
{
    return gamma*(1.0-std::pow((*this)(af,b),2))*b;
}
double Kern_Sigmoid::prime(const aed_type2& b, const aed_type2& af, const aed_type2& ff)
{
    return derivative(b, af)*ff;
}
std::string Kern_Sigmoid::get_label()
{
    return label;
}
