#include "d2_bp.h"
#include "../d_basis_functions.h"


D2_BP::D2_BP(Config &c):
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

void D2_BP::calc_aed(
        const double rij,
        const double ,
        const double fc_ij,
        aed_type2 &aed)
{

    size_t i=fidx;
    for (size_t c=0; c<mius.size(); c++) {         
        aed(i++) += G(rij,etas[c],mius[c],fc_ij);
    }

}
void D2_BP::calc_dXijdri(
        const double rij,
        const double ,
        const double fc_ij,
        const double fcp_ij,
        fd_type &fd_ij)
{
    size_t i=fidx;
    for (size_t c=0; c<mius.size(); c++) {         
        fd_ij(i++,0) = dG(rij,etas[c],mius[c],fc_ij,fcp_ij);
    }
}
void D2_BP::calc_all(
        const double rij,
        const double ,
        const double fc_ij,
        const double fcp_ij,
        aed_type2 &aed,
        fd_type &fd_ij)
{

    size_t i=fidx;
    for (size_t c=0; c<mius.size(); c++) {         
        aed(i) += G(rij,etas[c],mius[c],fc_ij);
        fd_ij(i,0) = dG(rij,etas[c],mius[c],fc_ij,fcp_ij);
        ++i;
    }
}
size_t D2_BP::size() {
    return s;
}
std::string D2_BP::label() {
    return lab;
}
