#include "d2_lj.h"


D2_LJ::D2_LJ(Config &c):
    verbose(c.get<int>("VERBOSE"))
{
    if (!c.get<bool>("INIT2B")) {
        s=0;
    }
}

void D2_LJ::calc_aed(
        const double,
        const double rij_sq,
        const double fc_ij,
        aed_type2 &aed)
{
    double r2_inv = 1.0/rij_sq;
    double r6_inv = r2_inv*r2_inv*r2_inv;

    aed(fidx) -= r6_inv*fc_ij;
    aed(fidx+1) += r6_inv*r6_inv*fc_ij;
}
void D2_LJ::calc_dXijdri(
        const double rij,
        const double rij_sq,
        const double fc_ij,
        const double fcp_ij,
        fd_type &fd_ij)
{

    double r2_inv = 1.0/rij_sq;
    double r6_inv = r2_inv*r2_inv*r2_inv;

    fd_ij(fidx,0) = 6.0*r6_inv*r2_inv*rij*fc_ij - fcp_ij*r6_inv;
    fd_ij(fidx+1,0) = -12.0*r6_inv*r6_inv*r2_inv*rij*fc_ij + fcp_ij*r6_inv*r6_inv;
}
void D2_LJ::calc_all(
        const double rij,
        const double rij_sq,
        const double fc_ij,
        const double fcp_ij,
        aed_type2 &aed,
        fd_type &fd_ij)
{
    double r2_inv = 1.0/rij_sq;
    double r6_inv = r2_inv*r2_inv*r2_inv;

    aed(fidx) -= r6_inv*fc_ij;
    aed(fidx+1) += r6_inv*r6_inv*fc_ij;

    fd_ij(fidx,0) = 6.0*r6_inv*r2_inv*rij*fc_ij - fcp_ij*r6_inv;
    fd_ij(fidx+1,0) = -12.0*r6_inv*r6_inv*r2_inv*rij*fc_ij + fcp_ij*r6_inv*r6_inv;
}
size_t D2_LJ::size() {
    return s;
}
std::string D2_LJ::label() {
    return lab;
}
