#include "cutoffs.h"
#include <cmath>
#include <stdexcept>


template<> Registry<Cut_Base,double>::Map Registry<Cut_Base,double>::registry{};
Registry<Cut_Base,double>::Register<Cut_Cos> Cut_Cos_1( "Cut_Cos" );
Registry<Cut_Base,double>::Register<Cut_Tanh> Cut_Tanh_1( "Cut_Tanh" );
Registry<Cut_Base,double>::Register<Cut_Poly2> Cut_Poly_1( "Cut_Poly2" );
Registry<Cut_Base,double>::Register<Cut_Dummy> Cut_Dummy_1( "Cut_Dummy" );


Cut_Base::~Cut_Base() {}
void Cut_Base::test_rcut(const double r)
{
    if (r<0.0)
        throw std::runtime_error("Cutoff distance cannot be negative.");
}

Cut_Cos::Cut_Cos() {}
Cut_Cos::Cut_Cos(double rcut)
{
    set_rcut(rcut);
    test_rcut(rcut);
}

std::string Cut_Cos::label() {
    return lab;
}

void Cut_Cos::set_rcut(const double r) {
    test_rcut(r);
    rcut=r;
    rcut_sq=r*r;
    rcut_inv= r<=0 ? 0.0 : 1.0/r;
}
double Cut_Cos::get_rcut() {
    return rcut;
}

double Cut_Cos::get_rcut_sq() {
    return rcut_sq;
}

double Cut_Cos::calc(double r) {
    if (r>=rcut) return 0.0;
    return 0.5 * (std::cos(M_PI*r*rcut_inv) + 1.0);
}
double Cut_Cos::calc_prime(double r) {
    if (r>=rcut) return 0.0;
    double t = M_PI*rcut_inv;
    return -0.5 * t*sin(t*r);

}

Cut_Dummy::Cut_Dummy() {}
Cut_Dummy::Cut_Dummy(double rcut)
{
    set_rcut(rcut),
    test_rcut(rcut);
}

std::string Cut_Dummy::label() {
    return lab;
}

void Cut_Dummy::set_rcut(const double r) {
    test_rcut(r);
    rcut=r;
    rcut_sq=r*r;
    rcut_inv= r<=0 ? 0.0 : 1.0/r;
}
double Cut_Dummy::get_rcut() {
    return rcut;
}

double Cut_Dummy::get_rcut_sq() {
    return rcut_sq;
}

double Cut_Dummy::calc(double r) {
    if (r>=rcut) return 0.0;
    return 1.0;
}
double Cut_Dummy::calc_prime(double) {
    return 0.0;
}

Cut_Tanh::Cut_Tanh() {}
Cut_Tanh::Cut_Tanh(double rcut)
{
    set_rcut(rcut);
    test_rcut(rcut);
}

std::string Cut_Tanh::label() {
    return lab;
}

void Cut_Tanh::set_rcut(const double r) {
    test_rcut(r);
    rcut=r;
    rcut_sq=r*r;
    rcut_inv= r<=0 ? 0.0 : 1.0/r;
}
double Cut_Tanh::get_rcut() {
    return rcut;
}

double Cut_Tanh::get_rcut_sq() {
    return rcut_sq;
}

double Cut_Tanh::calc(double r) {
    if (r>=rcut) return 0.0;
    double t = std::tanh(1.0 - r*rcut_inv);
    return t*t*t;
}
double Cut_Tanh::calc_prime(double r) {
    if (r>=rcut) return 0.0;
    double a = (1.0-r*rcut_inv);
    double t = std::tanh(a);
    double s = 1.0/(std::cosh(a));
    return -3.0*t*t*rcut_inv*s*s;
}

Cut_Poly2::Cut_Poly2() {}
Cut_Poly2::Cut_Poly2(double rcut)
{
    if (rcut<=1.0)
        throw std::runtime_error(
                "Cut_Poly2 cutoff distance cannot be smaller than 1.0.");
    set_rcut(rcut);
    test_rcut(rcut);
}

std::string Cut_Poly2::label() {
    return lab;
}

void Cut_Poly2::set_rcut(const double r) {
    test_rcut(r);
    rcut=r;
    rcut_sq=r*r;
    rcut_inv= r<=0 ? 0.0 : 1.0/r;
    rcut_inner=rcut-1.0;
}
double Cut_Poly2::get_rcut() {
    return rcut;
}

double Cut_Poly2::get_rcut_sq() {
    return rcut_sq;
}

double Cut_Poly2::calc(double r) {
    if (r>=rcut) return 0.0;
    else if (r<= rcut_inner) return 1.0;
    double rs=r-rcut_inner;
    return rs*rs*rs*(rs*(15.0-6.0*rs)-10.0)+1;
}
double Cut_Poly2::calc_prime(double r) {
    if (r>=rcut) return 0.0;
    else if (r<= rcut_inner) return 0.0;
    double rs=r-rcut_inner;
    return -30.0*(rs-1.0)*(rs-1.0)*rs*rs;
}
