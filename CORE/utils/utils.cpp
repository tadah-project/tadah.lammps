#include "utils.h"
#include <stdexcept>

v_type logspace(double start, double stop, int num, double base) {
    if (start==0) {
        throw std::runtime_error("log series cannot start at 0\n");
    }
    double realStart = pow(base, log(start));
    double realBase = pow(base, (log(stop)-log(start))/(num-1));

    v_type retval;
    retval.reserve(num);
    std::generate_n(std::back_inserter(retval), num, Logspace<double>(realStart,realBase));
    return retval;
}

size_t number_of_digits (size_t i)
{
    return i > 0 ? (int) log10 ((double) i) + 1 : 1;
}
size_t max_number(size_t n) {
    size_t val = 0;
    if (n == 0) return 0;
    for (size_t i=0; i<n; ++i) val+= 9*std::pow(10,i);
    return val;
}

// convert vector to string
std::string vec_to_string(std::vector<double> vec,
        std::string delim=" ") {
    std::stringstream ss;
    for (auto &v:vec) ss << std::setprecision(15)
        << v << delim;
    return ss.str().substr(0,ss.str().size()-delim.size()); // remove last delim
}

//std::vector<double> calc_rho(Structure &st) {
//    // Computes density according to Ta2 by Ravelo
//    double p=5.9913;
//    double q=8.0;
//    double rho0=0.074870;
//    double rc=5.5819;
//    double rc_sq=rc*rc;
//    double rcp=pow(rc,p);
//    double a0=3.304;
//    double r0=sqrt(3)*a0/2;
//    double r0p=pow(r0,p);
//    // assume full NN list
//    std::vector<double> rho(st.natoms(),0);
//    for (size_t i=0; i<st.natoms(); ++i) {
//        const Atom &a1 = st.atoms[i];
//        for (size_t jj=0; jj<st.nn_size(i); ++jj) {
//            const Vec3d &a2pos = st.nn_pos(i,jj);
//            Vec3d delij = a1.position - a2pos;
//            double rij_sq = delij.transpose() * delij;
//            if (rij_sq > rc_sq) continue;
//            double r=sqrt(rij_sq);
//            rho[i]+=rho0*pow((rcp-pow(r,p))/(rcp-r0p),q);
//        }
//    }
//    return rho;
//}
std::string to_upper(const std::string s) {
    std::string str = s;
    std::transform(str.begin(), str.end(),str.begin(), ::toupper);
    return str;
}
