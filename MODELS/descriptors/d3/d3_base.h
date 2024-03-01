#ifndef D3_BASE_H
#define D3_BASE_H

#include "../d_base.h"

class D3_Base: public D_Base {
    public:
        virtual ~D3_Base() {};

        virtual void calc_aed(
                const size_t fidx,
                const double rij,
                const double rik,
                const double fc_ij,
                const double fc_ik,
                aed_type2& aed)=0;
        virtual void calc_fd(
                const size_t fidx,
                const double rij,
                const double rik,
                const double fc_ij,
                const double fc_ik,
                const double fcp_ij,
                const double fcp_ik,
                fd_type &fd_ij)=0;
        virtual void calc_all(
                const size_t fidx,
                const double rij,
                const double rik,
                const double fc_ij,
                const double fc_ik,
                const double fcp_ij,
                const double fcp_ik,
                aed_type2& aed,
                fd_type &fd_ij)=0;
        virtual size_t size()=0;
        virtual std::string label()=0;
};
//template<> inline Registry<D3_Base,Config&>::Map Registry<D3_Base,Config&>::registry{};
#endif
