#ifndef DM_EAD_H
#define DM_EAD_H
#include "dm_base.h"

/**  Embedded Atom Descriptor
 *
 * \f[
 * V_i^{L,\eta,r_s} = \sum_{l_x,l_y,l_z}^{l_x+l_y+l_z=L} \frac{L!}{l_x!l_y!l_z!}
 * \Big( \rho_i^{\eta,r_s,l_x,l_y,l_z} \Big)^2
 * \f]
 *
 * where density \f$ \rho \f$ is calculated using Gaussian Type Orbitals:
 * \f[
 * \rho_i^{\eta,r_s,l_x,l_y,l_z} = \sum_{j \neq i} x_{ij}^{l_x}y_{ij}^{l_y}z_{ij}^{l_z}
 * \exp{\Big(-\eta(r_{ij}-r_s)^2\Big)}f_c(r_{ij})
 * \f]
 *
 * \ref CGRIDMB parameters control position \f$ r_s \f$ of the gaussian basis function.
 *
 * \ref SGRIDMB parameters control width \f$ \eta \f$ of the gaussian basis function.
 *
 * \ref AGRIDMB parameter specify maximum value for the angular momentum \f$L_{max}\f$.
 *
 * e.g. \f$L_{max}=2\f$ will calculate descriptors with \f$ L=0,1,2 \f$ (s,p,d orbitals).
 *
 * More information about this descriptor:
 *
 * <div class="csl-entry">Zhang, Y., Hu, C.,Jiang, B.
 * (2019). Embedded atom neural network potentials: efficient and accurate
 * machine learning with a physically inspired representation.
 * <i>Journal of Physical Chemistry Letters</i>, <i>10</i>(17),
 * 4962–4967. https://doi.org/10.1021/acs.jpclett.9b02037</div>
 *
 * Required Config keys:
 * \ref INITMB \ref CGRIDMB \ref SGRIDMB \ref AGRIDMB
 */
class DM_EAD: public DM_Base {
    private:
        int Lmax=0;
        size_t rhoisize=0;
        size_t s=0;
        std::string lab="DM_EAD";
        v_type sgrid;
        v_type cgrid;
        std::vector<std::vector<std::vector<int>>> orbitals;
        std::vector<std::vector<double>> fac;
        size_t gen_atomic_orbitals(int L);
        double fact(int n);
        double my_pow(double,int);
        int verbose;
        double rc;

    public:
        DM_EAD(Config&);
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
