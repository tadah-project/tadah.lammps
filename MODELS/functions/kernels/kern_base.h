#ifndef KERN_BASE_H
#define KERN_BASE_H

#include "../../../CORE/core_types.h"
#include "../function_base.h"

#include <iostream>

/** \brief Abstract class to be used as a base for all kernels.
 *
 *  - b = basis vector
 *  - af = atomic energy descriptor
 *  - ff = force descriptor
 *  - all derivatives are defined wrt to the second argument
 */
class Kern_Base: public virtual Function_Base {
    public:
        Matrix basis;
        virtual ~Kern_Base();

        /**
         * Calculate kernel value given two vectors
         *
         * - b = basis vector
         * - af = atomic energy descriptor
         */
        virtual double operator() (const aed_type2 &b, const aed_type2 &af)=0;

        /**
         * Calculate the kernel derivative wrt to the second argument
         *
         * - b = basis vector
         * - af = atomic energy descriptor
         */
        virtual aed_type2 derivative(const aed_type2 &b, const aed_type2 &af)=0;

        /**
         * Calculate inner product of the kernel derivative
         * wrt to the second argument (b wrt af) with the force descriptor
         *
         * - b = basis vector
         * - af = atomic energy descriptor
         * - ff = force descriptor (TODO i-th dir component of it)
         *   TODO consider calculating all 3-directions at once
         */
        virtual double prime(const aed_type2 &b, const aed_type2 &af,
                const aed_type2 &ff)=0;

        /** \brief Set basis for calculations */
        virtual void set_basis(const Matrix b);

        /** \brief Return basis stored by this object */
        virtual const Matrix &get_basis() const;
        /** \brief Return basis stored by this object */
        virtual Matrix get_basis();

        virtual double epredict(const t_type &kweights, const aed_type2 &aed);
        virtual double fpredict(const t_type &kweights, const fd_type &fdij,
                const aed_type2 &aedi, const size_t k);
        virtual force_type fpredict(const t_type &kweights, const fd_type &fdij,
                const aed_type2 &aedi);

        virtual void read_basis_from_config(const Config &config,
                Matrix &b);
};
#endif
