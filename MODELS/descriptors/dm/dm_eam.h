#ifndef DM_EAM_H
#define DM_EAM_H
#include "dm_base.h"

/** \brief many-body part for the Embedded Atom Method descriptor.
 *
 * \f[
 * V_i = F\Bigg(\sum_{j \neq i} \rho(r_{ij}) \Bigg)
 * \f]
 *
 * This descriptor will load tabulated values for the density \f$ \rho \f$
 * and embedded energy \f$ F \f$ from the provided \ref SETFL file.
 *
 * This descriptor is usually used together with the two-body descriptor \ref D2_EAM
 * although this is not required and user can mix it with any other descriptors
 * or use it on its own.
 *
 * This descriptor will enforce cutoff distance as specified in a \ref SETFL file.
 * Set \ref RCUTMB to the same value to suppress the warning message.
 *
 * Required Config keys:
 * \ref INITMB \ref SETFL
 */
class DM_EAM: public DM_Base {
    private:
        int verbose;
        struct eam_file {
            std::string file_path;
            std::vector<double> frho;
            std::vector<double> rhor;
            std::vector<double> z2r;
            int nrho=0;
            double drho=0;
            int nr;
            double dr;
            double rdr;
            double rdrho;
            double rcut;
            int atomic_number;
            double atomic_mass;
            double lattice_param;
            std::string lattice;
            double rhomax;

        };
        eam_file ef;

        std::vector<std::vector<double>> frho_spline;
        std::vector<std::vector<double>> rhor_spline;
        std::vector<std::vector<double>> z2r_spline;

        size_t s=1;
        std::string lab="DM_EAM";
        void read_setfl();
        void gen_splines(int &n, double &delta, std::vector<double> &f, std::vector<std::vector<double>> &spline);

    public:
        DM_EAM(Config&);
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
