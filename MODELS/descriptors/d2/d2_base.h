#ifndef D2_BASE_H
#define D2_BASE_H

#include "../d_base.h"

/** \brief Base class for all two-body type descriptors.
 *
 * All two-body descriptors **must** inherit this class.
 */
class D2_Base: public D_Base {
    public:
        size_t fidx=0;    // first index in aed and fd arrays for this descriptor
        virtual ~D2_Base() {};

        /** \brief Calculate \ref AED
         *
         * Calculate Atomic Energy Descriptor for the atom local environment.
         */
        virtual void calc_aed(
                const double rij,
                const double rij_sq,
                const double fc_ij,
                aed_type2& aed)=0;

        /** \brief Calculate \ref FD
         *
         * Calculate Force Descriptor for the atom local environment.
         */
        virtual void calc_dXijdri(
                const double rij,
                const double rij_sq,
                const double fc_ij,
                const double fcp_ij,
                fd_type &fd_ij)=0;

        /** \brief Calculate approximation to \ref FD
         *
         * Calculate Force Descriptor for the atom local environment
         * using central difference approximation.
         *
         * h -> something small e.g. 1e-8
         * fcijp -> cutoff(rij+h,...)
         * fcijn -> cutoff(rij-h,...)
         * aedp -> calc_aed(rij+h,...)
         * aedn -> calc_aed(rij-h,...)
         * fd_ij - descriptor derivative
         */
        virtual void calc_fd_approx(
        const double rij,
        const double fcijp,
        const double fcijn,
        aed_type2& aedp,
        aed_type2& aedn,
        fd_type &fd_ij,
        const double h);

        /** \brief Calculate \ref AED + \ref FD */
        virtual void calc_all(
                const double rij,
                const double rij_sq,
                const double fc_ij,
                const double fcp_ij,
                aed_type2& aed,
                fd_type &fd_ij)=0;
        virtual std::string label()=0;
};
//template<> inline Registry<D2_Base,Config&>::Map Registry<D2_Base,Config&>::registry{};
#endif
