#include "d2_base.h"
//template struct Registry<D2_Base, Config &>;
void D2_Base::calc_fd_approx(
        const double ,
        const double ,
        const double ,
        aed_type2& aedp,
        aed_type2& aedn,
        fd_type &fd_ij,
        const double h) {

    for (size_t i=0; i<(*this).size(); ++i)
        fd_ij(i,0) = (aedp(i)-aedn(i))/(2*h);
}
