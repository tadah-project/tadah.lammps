#ifndef D2_LJ_H
#define D2_LJ_H
#include "d2_base.h"
/** \anchor D2_LJ \brief Standard Lennard - Jones descriptor
 *
 * \f[
 * V_i = \sum_{j \neq i} 4 \epsilon \Bigg(\Big(\frac{\sigma}{r_{ij}}\Big)^{12} - \Big(\frac{\sigma}{r_{ij}}\Big)^6\Bigg) f_c(r_{ij})
 * \f]
 *
 * or equivalently:
 *
 * \f[
 * V_i = \sum_{j \neq i} \frac{C_{12}}{r_{ij}^{12}} - \frac{C_6}{r_{ij}^6} f_c(r_{ij})
 * \f]
 *
 * Note that machined learned coefficients \f$C_6\f$ and \f$C_{12}\f$
 * corresponds to \f$\sigma\f$ and \f$\epsilon\f$ through the following relation:
 *
 * \f[ \sigma = \Big(\frac{C_{12}}{C_6}\Big)^{1/6} \f]
 * \f[ \epsilon = \frac{1}{4} \frac{C_6^2}{C_{12}} w(Z) \f]
 * where \f$w(Z)\f$ is a species depended weight factor (default is an atomic number).
 *
 * The machine learned \f$\sigma\f$ and \f$\epsilon\f$ only make sense
 * (say to compare with the literature ones) when \ref BIAS false
 * and \ref NORM false and system in monatomic.
 * It is ok thought to set them to true it's just that
 * numerical values will be different.
 *
 * Required Config Key: \ref INIT2B
 */
class D2_LJ : public D2_Base {
    private:
        size_t s=2;
        std::string lab="D2_LJ";
        int verbose;
    public:
        D2_LJ(Config &config);
        void calc_aed(
                const double,
                const double rij_sq,
                const double fc_ij,
                aed_type2 &aed);
        void calc_dXijdri(
                const double rij,
                const double rij_sq,
                const double fc_ij,
                const double fcp_ij,
                fd_type &fd_ij);
        void calc_all(
                const double rij,
                const double rij_sq,
                const double fc_ij,
                const double fcp_ij,
                aed_type2 &aed,
                fd_type &fd_ij);
        size_t size();
        std::string label();
};
#endif
