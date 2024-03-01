#ifndef D2_BLIP_H
#define D2_BLIP_H
#include "d2_base.h"

/** \brief Blip two-body descriptor.
 *
 * \f[
 * V_i^{\eta,r_s} =\sum_{j \neq i} B(\eta(r_{ij}-r_s))f_c(r_{ij})
 * \f]
 *
 * where \f$ f_c \f$ is a cutoff function and \f$ B \f$ is a blip
 * basis function centered at \f$r_s\f$ of width \f$4/\eta\f$.
 *
 * \ref CGRID2B parameters control position \f$ r_s \f$ of blip centres.
 *
 * \ref SGRID2B parameters control width \f$ \eta \f$ of blips.
 *
 * Blip basis function is built out of 3rd degree polynomials
 * in the four intervals [-2,-1], [-1,0], [0,1], [1,2] and is defined as:
 *
 * \f[
	\begin{equation}
	  B(r) =
		\begin{cases}
		  1-\frac{3}{2}r^2+\frac{3}{4}|r|^3 & \text{if} \qquad 0<|r|<1\\
		  \frac{1}{4}(2-|r|)^3 & \text{if} \qquad 1<|r|<2\\
		  0 & \text{if} \qquad |r|>2
		\end{cases}
	\end{equation}
 * \f]
 *
 * More details about the blip basis functions can be found in the following paper:
 *
 * <div class="csl-entry">Hernández, E., Gillan, M., Goringe, C.
 * (1997). Basis functions for linear-scaling first-principles calculations.
 * <i>Physical Review B - Condensed Matter and Materials Physics</i>,
 * <i>55</i>(20), 13485–13493. https://doi.org/10.1103/PhysRevB.55.13485</div>
 *
 * Required keys:
 * \ref INIT2B \ref CGRID2B \ref SGRID2B
 */
class D2_Blip : public D2_Base {
    private:
        size_t s=0;
        std::string lab="D2_Blip";
        v_type etas;
        v_type mius;
        double get_blip(double);
        double get_dblip(double, double);
        v_type gen_blip_grid(double, double);
        int verbose;

    public:
        D2_Blip(Config &config);
        void calc_aed(
                const double rij,
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
