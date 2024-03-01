#ifndef BASIS_FUNCTIONS_H
#define BASIS_FUNCTIONS_H

#include "../../../CORE/core_types.h"
#include "../function_base.h"

#include <iostream>

/** \brief Base class to be used by all basis functions */
struct BF_Base: public virtual Function_Base {
    virtual ~BF_Base();
    virtual double epredict(const t_type &, const aed_type2& )=0;
    virtual double fpredict(const t_type &, const fd_type &,
            const aed_type2& , const size_t )=0;
    virtual force_type fpredict(const t_type &, const fd_type &,
            const aed_type2& )=0;
};
#endif
