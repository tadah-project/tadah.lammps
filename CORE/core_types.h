#ifndef CORE_TYPES_H
#define CORE_TYPES_H

#include "maths.h"

#include <vector>

typedef Vec3d force_type;
typedef aed_type2 t_type;
typedef aed_type2 rho_type;
typedef aeds_type2 rhos_type;
typedef Matrix phi_type;
typedef Matrix3d stress_type;

typedef std::vector<double> v_type;

////////////////////////////////////////////////////////////////////////////////
// Define dummy mat class to satisfy Sigma in core model classes. //////////////
////////////////////////////////////////////////////////////////////////////////
#ifdef M_MD_BASE_H
//#define TEMP_MAT_SIGMA
class mat {
    public:
        t_type diagonal() {
            return t_type();
        }
        t_type diagonal() const {
            return t_type();
        }

};
#endif
/////////////////////////////////////////

#endif
