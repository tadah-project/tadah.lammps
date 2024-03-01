#include "d2_mie.h"


D2_MIE::D2_MIE(Config &c):
    verbose(c.get<int>("VERBOSE"))
{
    if (!c.get<bool>("INIT2B")) {
        s=0;
        return;
    }
    get_grid(c,"SGRID2B",etas);
    if (etas.size()!=2) {
        throw std::runtime_error("Number of elements in SGRID2B is != 2\n\
                Mie descriptor requires SGRID2B with two positive elements.\n");
    }
    n=etas[1];
    m=etas[0];
    if (n<0 || m<0) {
        throw std::runtime_error("Both Mie exponents must by positive\n");
    }
}

void D2_MIE::calc_aed(
        const double rij,
        const double,
        const double fc_ij,
        aed_type2 &aed)
{
    double rn_inv = 1.0/pow(rij,n);
    double rm_inv = 1.0/pow(rij,m);

    aed(fidx) -= rn_inv*fc_ij;
    aed(fidx+1) += rm_inv*fc_ij;
}
void D2_MIE::calc_dXijdri(
        const double rij,
        const double,
        const double fc_ij,
        const double fcp_ij,
        fd_type &fd_ij)
{

    fd_ij(fidx,0)   = n/pow(rij,n+1)*fc_ij - fcp_ij/pow(rij,n);
    fd_ij(fidx+1,0) = -m/pow(rij,m+1)*fc_ij + fcp_ij/pow(rij,m);

}
void D2_MIE::calc_all(
        const double rij,
        const double,
        const double fc_ij,
        const double fcp_ij,
        aed_type2 &aed,
        fd_type &fd_ij)
{

    double rn_inv = 1.0/pow(rij,n);
    double rm_inv = 1.0/pow(rij,m);

    aed(fidx) -= rn_inv*fc_ij;
    aed(fidx+1) += rm_inv*fc_ij;

    fd_ij(fidx,0) = n*rn_inv/rij*fc_ij - fcp_ij*rn_inv;
    fd_ij(fidx+1,0) = -m*rm_inv/rij*fc_ij + fcp_ij*rm_inv;

}
size_t D2_MIE::size() {
    return s;
}
std::string D2_MIE::label() {
    return lab;
}
