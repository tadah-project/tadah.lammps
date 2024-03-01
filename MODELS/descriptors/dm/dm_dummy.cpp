#include "dm_dummy.h"


DM_Dummy::DM_Dummy() {}
DM_Dummy::DM_Dummy(Config&) {}
void DM_Dummy::calc_aed(
        rho_type&,
        aed_type2 & ) {}
int DM_Dummy::calc_dXijdri_dXjidri(
        const double,
        const double,
        const Vec3d &,
        const double,
        const double,
        rho_type&,
        rho_type&,
        fd_type &,
        const double,
        const double) { return 0; }
int DM_Dummy::calc_dXijdri(
        const double,
        const double,
        const Vec3d &,
        const double,
        const double,
        rho_type&,
        fd_type &) { return 0; }
size_t DM_Dummy::size() { return s;}
std::string DM_Dummy::label() {return lab;}
void DM_Dummy::init_rhoi(rho_type&) {}
void DM_Dummy::calc_rho(
        const double,
        const double,
        const double,
        const Vec3d &,
        rho_type&) {}
size_t DM_Dummy::rhoi_size() { return 0; }
size_t DM_Dummy::rhoip_size() { return 0; }

