#ifndef M_PREDICT_H
#define M_PREDICT_H

#include "../CORE/core_types.h"

class M_Predict {
    public:

        virtual ~M_Predict() {}

        /** \brief Predict local energy of an atom or bond energy.
         *
         *  The result depends on how aed is computed.
         *
         *  If it is computed between a pair of atoms than the result is a bond energy.
         *
         *  If aed contains sum over all nearest neighbours than the result is
         *  a local atomic energy \f$ E_i \f$.
         * */
        virtual double epredict(const aed_type2 &aed)=0;

        /** \brief Predict force between a pair of atoms in a k-direction. */
        virtual double fpredict(const fd_type &fdij, const aed_type2 &aedi, size_t k)=0;

        /** \brief Predict force between a pair of atoms. */
        virtual force_type fpredict(const fd_type &fdij,
                const aed_type2 &aedi)=0;
};
#endif
