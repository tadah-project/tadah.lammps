#ifndef D2_DUMMY_H
#define D2_DUMMY_H
#include "d2_base.h"

/**
 * Dummy two-body descriptor.
 *
 * Use it to satisfy DescriptorsCalc requirements in case
 * when two-body descriptor is not required.
 *
 */
class D2_Dummy : public D2_Base {
    private:
        size_t s=0;
        std::string lab="D2_Dummy";
        int verbose;
    public:
        D2_Dummy();
        D2_Dummy(Config &);
        void calc_aed(
                const double,
                const double,
                const double,
                aed_type2 &aed);
        void calc_dXijdri(
                const double,
                const double,
                const double,
                const double,
                fd_type &fd_ij);
        void calc_all(
                const double,
                const double,
                const double,
                const double,
                aed_type2 &aed,
                fd_type &fd_ij);
        size_t size();
        std::string label();
};
#endif
