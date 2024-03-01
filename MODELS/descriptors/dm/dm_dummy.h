#ifndef DM_DUMMY_H
#define DM_DUMMY_H
#include "dm_base.h"

/**
 * Dummy many-body descriptor.
 *
 * Use it to satisfy DescriptorsCalc requirements in case
 * when many-body descriptor is not required.
 *
 */
class DM_Dummy: public DM_Base {
    private:
        size_t s=0;
        std::string lab="DM_Dummy";
    public:
        DM_Dummy();
        DM_Dummy(Config&);
        void calc_aed(
                rho_type& rho,
                aed_type2 &aed);
        int calc_dXijdri_dXjidri(
                const double rij,
                const double rij_sq,
                const Vec3d &vec_ij,
                const double fc_ij,
                const double fcp_ij,
                rho_type& rhoi,
                rho_type& rhoj,
                fd_type &fd_ij,
                const double wi,
                const double wj);
        int calc_dXijdri(
                const double rij,
                const double rij_sq,
                const Vec3d &vec_ij,
                const double fc_ij,
                const double fcp_ij,
                rho_type& rhoi,
                fd_type &fd_ij);
        size_t size();
        std::string label();
        void init_rhoi(rho_type& rhoi);
        void calc_rho(
                const double rij,
                const double rij_sq,
                const double fc_ij,
                const Vec3d &vec_ij,
                rho_type& rho);
        size_t rhoi_size();
        size_t rhoip_size();
};
#endif
