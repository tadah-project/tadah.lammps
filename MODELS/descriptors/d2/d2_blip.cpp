#include "d2_blip.h"
#include "../d_basis_functions.h"
#include <stdexcept>


D2_Blip::D2_Blip(Config &c):
    verbose(c.get<int>("VERBOSE"))
{

    if (!c.get<bool>("INIT2B")) return;

    get_grid(c,"CGRID2B",mius);
    get_grid(c,"SGRID2B",etas);

    if (verbose) {
        std::cout << std::endl;
        std::cout << "SGRID2B: " << etas.size() << std::endl;
        for (auto e:etas) std::cout << e << "  ";

        std::cout << std::endl;
        std::cout << "CGRID2B: " << mius.size() << std::endl;
        for (auto m:mius) std::cout << m << "  ";
        std::cout << std::endl;
    }

    if (mius.size()!=etas.size()) {
        throw std::runtime_error("SGRID2B and CGRID2B arrays differ in size.\n");
    }

    s=mius.size();
}

void D2_Blip::calc_aed(
        const double rij,
        const double ,
        const double fc_ij,
        aed_type2 &aed)
{

    size_t i=fidx;
    for (size_t c=0; c<mius.size(); c++) {
        double t = B(etas[c]*(rij-mius[c]),fc_ij);
        aed(i++) += t;
    }

}
void D2_Blip::calc_dXijdri(
        const double rij,
        const double ,
        const double fc_ij,
        const double fcp_ij,
        fd_type &fd_ij)
{
    size_t i=fidx;
    for (size_t c=0; c<mius.size(); c++) {
        fd_ij(i++,0) = dB(etas[c]*(rij-mius[c]),etas[c],fc_ij,fcp_ij);
    }
}
void D2_Blip::calc_all(
        const double rij,
        const double ,
        const double fc_ij,
        const double fcp_ij,
        aed_type2 &aed,
        fd_type &fd_ij)
{

    size_t i=fidx;
    for (size_t c=0; c<mius.size(); c++) {
        aed(i) += B(etas[c]*(rij-mius[c]),fc_ij);
        fd_ij(i,0) = dB(etas[c]*(rij-mius[c]),etas[c],fc_ij,fcp_ij);
        ++i;
    }
}

size_t D2_Blip::size() {
    return s;
}
std::string D2_Blip::label() {
    return lab;
}
