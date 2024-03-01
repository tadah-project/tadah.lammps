#include "kern_quadratic.h"


Kern_Quadratic::Kern_Quadratic () {}
Kern_Quadratic::Kern_Quadratic (const Config &c):
    Function_Base(c)
{
    if (get_verbose()) std::cout << std::endl << label << std::endl;

    read_basis_from_config(c,basis);
}
double Kern_Quadratic::operator() (const aed_type2& b, const aed_type2& af)
{
    double x = b*af;
    return  x*x;
}
aed_type2 Kern_Quadratic::derivative(const aed_type2& b, const aed_type2& af)
{
    double temp = 2.0 * b*af;
    return b*temp;;
}
double Kern_Quadratic::prime(const aed_type2& b, const aed_type2& af, const aed_type2& ff)
{
    //return dot(derivative(b, af), ff);
    // it is faster this way
    double scale = 2.0*b*af;
    return scale*ff*b;
}
std::string Kern_Quadratic::get_label()
{
    return label;
}
