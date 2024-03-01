#ifndef D2_EAM_H
#define D2_EAM_H
#include "d2_base.h"

/** \brief Pair-wise part for the Embedded Atom Method descriptor.
 *
 * \f[
 * V_i = \frac{1}{2} \sum_{j \neq i} \psi(r_{ij})
 * \f]
 *
 * This descriptor will load tabulated values for the two-body potential \f$ \phi \f$
 * from the provided \ref SETFL file.
 *
 * This descriptor is usually used together with the many-body descriptor \ref DM_EAM
 * although this is not required and user can mix it with any other descriptors
 * or use it on its own.
 *
 * This descriptor will enforce cutoff distance as specified in a \ref SETFL file.
 * Set \ref RCUT2B to the same value to suppress the warning message.
 *
 * Required Config keys:
 * \ref INIT2B \ref SETFL
 */
class D2_EAM : public D2_Base {
    private:
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
        };
		eam_file ef;

        std::vector<std::vector<double>> frho_spline;
        std::vector<std::vector<double>> rhor_spline;
        std::vector<std::vector<double>> z2r_spline;

        int verbose;

        size_t s=1;
        std::string lab="D2_EAM";
		void read_setfl();
		void gen_splines(int &n, double &delta, std::vector<double> &f, std::vector<std::vector<double>> &spline);

    public:
        D2_EAM(Config &config);
        void calc_aed(
                double rij,
                const double rij_sq,
                const double,
                aed_type2 &aed);
        void calc_dXijdri(
                double rij,
                const double rij_sq,
                const double,
                const double,
                fd_type &fd_ij);
        void calc_all(
                double rij,
                const double rij_sq,
                const double,
                const double,
                aed_type2 &aed,
                fd_type &fd_ij);
        size_t size();
        std::string label();
};
#endif
