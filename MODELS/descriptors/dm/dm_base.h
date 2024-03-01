#ifndef DM_BASE_H
#define DM_BASE_H

#include <string>
#include "../d_base.h"

/** \brief Base class for all many-body type descriptors.
 *
 * All many-body descriptors **must** inherit this class.
 */
class DM_Base: public D_Base {
    public:
        size_t fidx=0;    // first index in aed and fd arrays for this descriptor
        virtual ~DM_Base() {};

        /** \brief Calculate \ref AED
         *
         * Calculate Atomic Energy Descriptor for the atom local environment.
         */
        virtual void calc_aed(
                rho_type& rho,
                aed_type2 &aed)=0;

        /** \brief Calculate \ref FD
         *
         * Calculate Force Descriptor between to atoms.
         *
         * This method works for half NN lists
         * and linear models.
         *
         * \f[
         * \mathbf{fd}_{ij} \mathrel{+}= 
         * \frac{\partial \mathbf{X}_{ij}}{\partial \mathbf{r}_i} +
         * \frac{\partial \mathbf{X}_{ji}}{\partial \mathbf{r}_i}
         * \f]
         */
        virtual int calc_dXijdri_dXjidri(
                const double rij,
                const double rij_sq,
                const Vec3d &vec_ij,
                const double fc_ij,
                const double fcp_ij,
                rho_type& rhoi,
                rho_type& rhoj,
                fd_type &fd_ij,
                const double wi,
                const double wj)=0;

        /** \brief Calculate \ref FD
         *
         * Calculate Force Descriptor between to atoms.
         *
         * This method works for full NN lists
         * and all models.
         *
         * \f[
         * \mathbf{fd}_{ij} \mathrel{+}=
         * \frac{\partial \mathbf{X}_{ij}}{\partial \mathbf{r}_i}
         * \f]
         */
        virtual int calc_dXijdri(
                const double rij,
                const double rij_sq,
                const Vec3d &vec_ij,
                const double fc_ij,
                const double fcp_ij,
                rho_type& rhoi,
                fd_type &fd_ij)=0;

        /** \brief Resize arrays for a density and F' */
        virtual void init_rhoi(rho_type& rhoi)=0;

        /** \brief Calculate density */
        virtual void calc_rho(
                const double rij,
                const double rij_sq,
                const double fc_ij,
                const Vec3d &vec_ij,
                rho_type& rho)=0;

        /** \brief Return size of the density array */
        virtual size_t rhoi_size()=0;

        /** \brief Return size of the derivative of the embedding energy array */
        virtual size_t rhoip_size()=0;
        virtual size_t size()=0;
        virtual std::string label()=0;
};
//template<> inline Registry<DM_Base,Config&>::Map Registry<DM_Base,Config&>::registry{};
#endif
