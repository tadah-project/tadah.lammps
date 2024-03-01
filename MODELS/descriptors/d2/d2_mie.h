#ifndef D2_MIE_H
#define D2_MIE_H
#include "d2_base.h"
/** \brief Mie descriptor
 *
 * \f[
 * V_i = \sum_{j \neq i} C \epsilon \Bigg(\Big(\frac{\sigma}{r_{ij}}\Big)^{n} - \Big(\frac{\sigma}{r_{ij}}\Big)^m\Bigg)
 * \f]
 *
 * where
 * \f[
 * C=\frac{n}{n-m}\Big( \frac{n}{m} \Big)^{\frac{m}{n-m}}
 * \f]
 *
 *
 * Any cutoff can be used
 *
 * Required Config Key: \ref INIT2B \ref SGRID2B
 *
 * e.g. SGID2B 12 6
 *
 * will result in Lennard-Jones type descriptor
 */
class D2_MIE : public D2_Base {
    private:
        size_t s=2;
        std::string lab="D2_MIE";
        double m,n;
        v_type etas;
        int verbose;
    public:
        D2_MIE(Config &config);
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
