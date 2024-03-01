#include "kern_base.h"

Kern_Base::~Kern_Base() {}
void Kern_Base::set_basis(const Matrix b)
{
    basis = b;
}
double Kern_Base::epredict(const t_type &kweights, const aed_type2 &aed)
{
    double energy=0;
    for (size_t b=0; b<basis.cols(); ++b) {
        energy += kweights[b]*(*this)(basis.col(b),aed);
    }
    return energy;
}
double Kern_Base::fpredict(const t_type &kweights, const fd_type &fdij,
        const aed_type2 &aedi, const size_t k)
{
    double res=0.0;
    for (size_t b=0; b<basis.cols(); ++b) {
        res -= kweights[b]*(*this).prime(basis.col(b),aedi,fdij(k));
    }
    return res;
}
force_type Kern_Base::fpredict(const t_type &kweights, const fd_type &fdij,
        const aed_type2 &aedi)
{
    force_type v;
    v.set_zero();
    for (size_t b=0; b<basis.cols(); ++b) {
        for (size_t k=0; k<3; ++k) {
            v(k) -= kweights[b]*(*this).prime(basis.col(b),aedi,fdij(k));
        }
    }
    return v;
}
void Kern_Base::read_basis_from_config(const Config &config, Matrix &b) {
    if (!(config.exist("SBASIS") && config.exist("BASIS"))) return;
    size_t sbasis = config.get<size_t>("SBASIS");
    size_t vlen = config("BASIS").size()/sbasis;
    std::vector<double> temp(sbasis*vlen);
    config.get<std::vector<double>>("BASIS",temp);
    b.resize(vlen,sbasis);  // aeds are in column-major order
    size_t ii=0;
    for (size_t i=0;i<sbasis; ++i) {
        for (size_t j=0;j<vlen; ++j) {
            b(j,i) = temp[ii++];
        }
    }
}
const Matrix &Kern_Base::get_basis() const {
    return basis;
}
Matrix Kern_Base::get_basis() {
    return basis;
}
