#include "d2_dummy.h"


D2_Dummy::D2_Dummy() {}
D2_Dummy::D2_Dummy(Config &c):
    verbose(c.get<int>("VERBOSE"))
{}

void D2_Dummy::calc_aed(
        const double,
        const double,
        const double,
        aed_type2 & ) {}
void D2_Dummy::calc_dXijdri(
        const double,
        const double,
        const double,
        const double,
        fd_type &) {}
void D2_Dummy::calc_all(
        const double,
        const double,
        const double,
        const double,
        aed_type2 & ,
        fd_type &) {}
size_t D2_Dummy::size() { return s; }
std::string D2_Dummy::label() { return lab; }
